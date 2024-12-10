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

#include <stdio.h>
#include <string.h>
#include <assert.h>

# include <unistd.h>
# include <pwd.h>
#include <cmath>

# define MAX_PATH FILENAME_MAX

#include <iomanip>
#include <iostream>
#include <fstream>
//#include <sys/time.h>       
#include <time.h>
#include <ctime>

using namespace std;
//#include "user_types.h"
//#include "user_helper.hpp" //Did not work

#include "sgx_urts.h"
#include "App.h"
#include "Enclave_u.h"

//#include "immintrin.h"  


/* Global EID shared by multiple threads */
sgx_enclave_id_t global_eid = 0;

typedef struct _sgx_errlist_t {
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] = {
    {
        SGX_ERROR_UNEXPECTED,
        "Unexpected error occurred.",
        NULL
    },
    {
        SGX_ERROR_INVALID_PARAMETER,
        "Invalid parameter.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_MEMORY,
        "Out of memory.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_LOST,
        "Power transition occurred.",
        "Please refer to the sample \"PowerTransition\" for details."
    },
    {
        SGX_ERROR_INVALID_ENCLAVE,
        "Invalid enclave image.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ENCLAVE_ID,
        "Invalid enclave identification.",
        NULL
    },
    {
        SGX_ERROR_INVALID_SIGNATURE,
        "Invalid enclave signature.",
        NULL
    },
    {
        SGX_ERROR_OUT_OF_EPC,
        "Out of EPC memory.",
        NULL
    },
    {
        SGX_ERROR_NO_DEVICE,
        "Invalid SGX device.",
        "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."
    },
    {
        SGX_ERROR_MEMORY_MAP_CONFLICT,
        "Memory map conflicted.",
        NULL
    },
    {
        SGX_ERROR_INVALID_METADATA,
        "Invalid enclave metadata.",
        NULL
    },
    {
        SGX_ERROR_DEVICE_BUSY,
        "SGX device was busy.",
        NULL
    },
    {
        SGX_ERROR_INVALID_VERSION,
        "Enclave version was invalid.",
        NULL
    },
    {
        SGX_ERROR_INVALID_ATTRIBUTE,
        "Enclave was not authorized.",
        NULL
    },
    {
        SGX_ERROR_ENCLAVE_FILE_ACCESS,
        "Can't open enclave file.",
        NULL
    },
};

/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist/sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++) {
        if(ret == sgx_errlist[idx].err) {
            if(NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }
    
    if (idx == ttl)
    	printf("Error code is 0x%X. Please refer to the \"Intel SGX SDK Developer Reference\" for more details.\n", ret);
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    
    /* Call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */
    ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL, &global_eid, NULL);
    if (ret != SGX_SUCCESS) {
        print_error_message(ret);
        return -1;
    }

    return 0;
}

/* OCall functions */
void ocall_print_string(const char *str)
{
    /* Proxy/Bridge will check the length and null-terminate 
     * the input string to prevent buffer overflow. 
     */
    printf("%s", str);
}

const uint32_t mb_of_data = 262144; //For the data size, number of rows. Based on TEEBench. 
const char *operations[3] = {"ADD", "AND", "SHIFT_RIGHT"};
const char *types[4] = {"scalar", "unroll", "simd_sse", "simd_avx"};

//Similar values like from the TEEBench.
uint32_t data_sizes[4] = {int(0.2*mb_of_data),
               int(6.4 * mb_of_data), 
               16 * mb_of_data,       
               100 * mb_of_data};

/* Application entry */
int SGX_CDECL main(int argc, char *argv[])
{   
    //Start total time:
    const std::clock_t c_start2 = std::clock();
    (void)(argc);
    (void)(argv);

    if(argc != 2){
        std::cout << "Don't forget the 1.argument. (1 means 'append', else normal write)." << std::endl;
        return 1;
    }

    /* Initialize the enclave */
    if(initialize_enclave() < 0){
        printf("Enter a character before exit ...\n");
        getchar();
        return -1; 
    }

    // Create/open the file
    std::ofstream File("scripts/data/simd_file2.csv", ios::app);
    if(!File){
        std::cout << "Opening the file failed." << std::endl;
        return 1;  
    }

    //Check, whether write or append is needed. 
    //(1 means "append", else normal write)
    if(atoi(argv[1])!=1){
        std::cout << "I'm here!!" << std::endl;
        File.close();
        File.open("scripts/data/simd_file2.csv", ios::out);
        // Write into the file
        File << "operation,type,mode,size,time\n";

        if(!File) {
            std::cout << "Opening the file failed." << std::endl;
            return 1;  
        }
    }

    //Define the data set.
    struct table t1;
    
    for(auto s: data_sizes){
        //Resize and init the data set.
        t1.row = NULL;
        t1.row = (uint32_t*) calloc(s, sizeof(uint32_t)) ;
        t1.size = s;
    
        for(uint32_t op = 0; op<3; op++){
            for(int t = 0; t < 4; t++){
                double sum_all_times = 0;
                for(int r = 0; r < RUNS; r++){
                    // Init the data set with values.
                    for(uint32_t a = 0; a < t1.size; a++){
                        t1.row[a] = a;
                    }
                    if(t==3 && op==0 && r==0){
                        const std::clock_t c_start0 = std::clock();
                        ecall_simd_ld_st_op(global_eid, &t1,  op);
                        const std::clock_t c_end0 = std::clock();
                        double time0 = 1000.0 * (c_end0 - c_start0) / CLOCKS_PER_SEC;
                        std::cout << "SIMD_load_store" <<  "," << types[t] << ",sgx,"  << std::setprecision(3) << (float) s/mb_of_data << "," << time0 <<  std::endl;
                        File << "SIMD_load_store" <<  "," << types[t] << ",sgx,"  << std::setprecision(3) << (float) s/mb_of_data << "," << time0 <<  std::endl;
                    }
                    const std::clock_t c_start1 = std::clock();
                    
                        switch(t){
                            case 0:
                                //ocall_scalar_op(&t1,  op);
                                ecall_scalar_op(global_eid, &t1, op);
                                break;

                            case 1:
                                //ocall_unroll_op(&t1, op);
                                ecall_unroll_op(global_eid, &t1, op);
                                break;

                            case 2:
                                //ocall_simd_op(&t1, op);
                                ecall_simd_sse_op(global_eid, &t1, op);
                                break;
                            
                            case 3:
                                ecall_simd_avx_op(global_eid, &t1, op);
                                break;
                        }
                    const std::clock_t c_end1 = std::clock();
                    double time1 = 1000.0 * (c_end1 - c_start1) / CLOCKS_PER_SEC;
                    sum_all_times +=time1;
                }
                double avg_time = sum_all_times / RUNS;

                std::cout << operations[op] <<  "," << types[t] << ",sgx,"  << std::setprecision(3) << (float) s/mb_of_data << "," << avg_time <<  std::endl;
                File << operations[op] <<  "," << types[t] << ",sgx," << std::setprecision(3) << (float) s/mb_of_data << "," << avg_time <<  std::endl;
            }
        }
    }
     
    /* Destroy the enclave */
    sgx_destroy_enclave(global_eid);

    const std::clock_t c_end2 = std::clock();
    double total_time = 1000.0 * (c_end2 - c_start2) / CLOCKS_PER_SEC;
    std::cout << "total time in native,"  << std::setprecision(3) << "," << total_time <<  std::endl;
    //File << "total time in native,"  << std::setprecision(3) << "," << total_time <<  std::endl;

    // Close the file
    File.close();
    printf("\n");
    printf("Info: SampleEnclave successfully returned.\n");
    return 0;
}

