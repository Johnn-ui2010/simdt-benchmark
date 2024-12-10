
#include "string.h"
#include <stdio.h>
#include <string.h>

//#include "user_types.h"
//#include <immintrin.h> // AVX 1 and 2 (256bit) and AVX512 (includes all SSE, too)

#include "../../simd/xmmintrin.h"
#include "../../simd/avxintrin.h"
#include "../../simd/avx2intrin.h"

void scalar_op(struct table *t, uint32_t OP);
void unroll_op(struct table *t, uint32_t OP);
void simd_sse_op(struct table *t, uint32_t OP);
void simd_avx_op(struct table *t, uint32_t OP);

void print256_num_uint32(__m256i var)
{  
    uint32_t val[8];
    memcpy(val, &var, sizeof(val));
    printf("Integer Numerical: %i %i %i %i %i %i %i %i \n", 
           val[0], val[1], val[2], val[3], val[4], val[5], 
           val[6], val[7]);
}

void scalar_op(struct table *t, uint32_t OP) {
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

void unroll_op(struct table *t, uint32_t OP) {
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

void simd_sse_op(struct table *t, uint32_t OP) {
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

void simd_avx_op(struct table *t, uint32_t OP) {
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

void simd_ld_st_op(struct table *t, uint32_t OP) {
    __m256i a_vec;
    for (uint32_t i = 0; i < t->size; i+=8) {
                a_vec = _mm256_loadu_si256((__m256i const *) (t->row +i));
                _mm256_storeu_ps((float *) (t->row +i), _mm256_castsi256_ps (a_vec) ); 
    }
}