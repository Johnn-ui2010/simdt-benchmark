/*
 * Copyright (C) 2011-2021 Intel Corporation. All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *   * Redistributions of source code must retain the above copyright
 *     notice, this list of conditions and the following disclaimer.
 *   * Redistributions in binary form must reproduce the above copyright
 *     notice, this list of conditions and the following disclaimer in
 *     the documentation and/or other materials provided with the
 *     distribution.
 *   * Neither the name of Intel Corporation nor the names of its
 *     contributors may be used to endorse or promote products derived
 *     from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#include "Enclave.h"
#include "Enclave_t.h" /* print_string */
#include <stdarg.h>
#include <stdio.h> /* vsnprintf */
#include <string.h>
//#include <user_types.h>
#include "../../simd/xmmintrin.h"
#include "../../simd/avxintrin.h"
#include "../../simd/avx2intrin.h"

/* 
 * printf: 
 *   Invokes OCALL to display the enclave buffer to the terminal.
 */
int printf(const char* fmt, ...)
{
    char buf[BUFSIZ] = { '\0' };
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print_string(buf);
    return (int)strnlen(buf, BUFSIZ - 1) + 1;
}


void ecall_scalar_op(struct table *t, uint32_t OP) {
    switch(OP){
        case ADD:
            for (uint32_t i = 0; i < t->size; i ++) {
                t->row[i] += t->row[i];
            }
            break;

        case AND:
            for (uint32_t i = 0; i < t->size; i ++) {
                t->row[i] &= NUM1;
            }
            break;

        case SHIFT_RIGHT:
            for (uint32_t i = 0; i < t->size; i ++) {
                t->row[i] = t->row[i] >> SHIFT;
            }
            break;
    }
}

void ecall_unroll_op(struct table *t, uint32_t OP) {
    switch(OP){
        case ADD:{
            for (uint32_t i = 0; i < t->size; i +=8) {
                t->row[i] += t->row[i];
                t->row[i+1] += t->row[i+1];
                t->row[i+2] += t->row[i+2];
                t->row[i+3] += t->row[i+3];
                t->row[i+4] += t->row[i+4];
                t->row[i+5] += t->row[i+5];
                t->row[i+6] += t->row[i+6];
                t->row[i+7] += t->row[i+7];
            }
            break;
        }
        case AND:{
            for (uint32_t i = 0; i < t->size; i +=8) {
                t->row[i] &= NUM1;
                t->row[i+1] &= NUM1;
                t->row[i+2] &= NUM1;
                t->row[i+3] &= NUM1;
                t->row[i+4] &= NUM1;
                t->row[i+5] &= NUM1;
                t->row[i+6] &= NUM1;
                t->row[i+7] &= NUM1;
            }
            break;
        }
        case SHIFT_RIGHT:{
            for (uint32_t i = 0; i < t->size; i +=8) {
                t->row[i] = t->row[i] >> SHIFT;
                t->row[i+1] = t->row[i+1] >> SHIFT;
                t->row[i+2] = t->row[i+2] >> SHIFT;
                t->row[i+3] = t->row[i+3] >> SHIFT;
                t->row[i+4] = t->row[i+4] >> SHIFT;
                t->row[i+5] = t->row[i+5] >> SHIFT;
                t->row[i+6] = t->row[i+6] >> SHIFT;
                t->row[i+7] = t->row[i+7] >> SHIFT;
            }
            break;
        }
    }
    
}
void ecall_simd_sse_op(struct table *t, uint32_t OP) {
    __m128i a_vec;

    switch(OP){
        case ADD:{
            for (uint32_t i = 0; i < t->size; i+=4) {
                a_vec = _mm_loadu_si128((__m128i const *) (t->row +i));
                a_vec = _mm_add_epi32(a_vec, a_vec);
                _mm_storeu_ps((float *) (t->row +i), _mm_castsi128_ps (a_vec) ); 
            }
            break;
        }
        case AND:{
            __m128i vec1 = _mm_set1_epi32(NUM1);
            for (uint32_t i = 0; i < t->size; i+=4) {
                a_vec = _mm_loadu_si128((__m128i const *) (t->row +i));
                a_vec = _mm_and_si128(a_vec, vec1);
                _mm_storeu_ps((float *) (t->row +i), _mm_castsi128_ps (a_vec) ); 
            }
            break;
        }
        case SHIFT_RIGHT:{
            for (uint32_t i = 0; i < t->size; i+=4) {
                a_vec = _mm_loadu_si128((__m128i const *) (t->row +i));
                a_vec = _mm_srli_epi32(a_vec, SHIFT);
                _mm_storeu_ps((float *) (t->row +i), _mm_castsi128_ps (a_vec) );  
            }
            break;
        }
    }
}

void ecall_simd_avx_op(struct table *t, uint32_t OP) {
    __m256i a_vec;

    switch(OP){
        case ADD:{
            for (uint32_t i = 0; i < t->size; i+=8) {
                a_vec = _mm256_loadu_si256((__m256i const *) (t->row +i));
                a_vec = _mm256_add_epi32(a_vec, a_vec);
                _mm256_storeu_ps((float *) (t->row +i), _mm256_castsi256_ps (a_vec) ); 
            }
            break;
        }
        case AND:{
            __m256i vec1 = _mm256_set1_epi32(NUM1);
            for (uint32_t i = 0; i < t->size; i+=8) {
                a_vec = _mm256_loadu_si256((__m256i const *) (t->row +i));
                a_vec = _mm256_and_si256(a_vec, vec1);
                _mm256_storeu_ps((float *) (t->row +i), _mm256_castsi256_ps (a_vec) ); 
            }
            break;
        }
        case SHIFT_RIGHT:{
            for (uint32_t i = 0; i < t->size; i+=8) {
                a_vec = _mm256_loadu_si256((__m256i const *) (t->row +i));
                a_vec = _mm256_srli_epi32(a_vec, SHIFT);
                _mm256_storeu_ps((float *) (t->row +i), _mm256_castsi256_ps (a_vec) ); 
            }
            break;
        }
    }
}

void ecall_simd_ld_st_op(struct table *t, uint32_t OP) {
    __m256i a_vec;
    for (uint32_t i = 0; i < t->size; i+=8) {
                a_vec = _mm256_loadu_si256((__m256i const *) (t->row +i));
                _mm256_storeu_ps((float *) (t->row +i), _mm256_castsi256_ps (a_vec) ); 
    }
}