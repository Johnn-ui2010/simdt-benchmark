#include <stdint.h>
#include "btree.h"
#include <pthread.h>
#include "data-types.h"

#ifdef NATIVE_COMPILATION

#include "Logger.h"
#include "native_ocalls.h"
#include <immintrin.h> //AVX 1 and 2 (256bit) and AVX512 (includes all SSE, too)

#else
#include "Enclave_t.h"
#include "Enclave.h"

#include "../simd/mmintrin.h"
#include "../simd/xmmintrin.h"
#include "../simd/emmintrin.h"
#include "../simd/pmmintrin.h"
#include "../simd/tmmintrin.h"
#include "../simd/avxintrin.h"
#include "../simd/avx2intrin.h"
#include "../simd/avx512fintrin.h" 
#endif

typedef struct arg_nl_t {
    tuple_t * relR;
    tuple_t * relS;

    uint64_t numR;
    uint64_t numS;

    uint64_t result;
    int32_t my_tid;
} arg_nl_t;

static void
print_timing(uint64_t total, uint64_t numtuples, int64_t result)
{
    double cyclestuple = (double) total / (double) numtuples;
    logger(DBG, "Total input tuples : %lu", numtuples);
    logger(DBG, "Result tuples : %lu", result);
    logger(DBG, "Phase Join (cycles) : %lu", total);
    logger(DBG, "Cycles-per-tuple : %.4lf", cyclestuple);
}

void * nlj_thread(void * param) {
    arg_nl_t *args = (arg_nl_t*) param;
    uint64_t results = 0;

    for (int32_t i = 0; i < args->numR; i++)
    {
        for (int32_t j = 0; j < args->numS; j++)
        {
            if (args->relR[i].key == args->relS[j].key)
            {
                results++;
                //logger(INFO, "res: %d, key: %d\n", results, args->relR[i].key);
                //logger(INFO, "res: %d\n", results);
            }
        }
    }

    args->result = results;
    return nullptr;
}

result_t* NL (struct table_t* relR, struct table_t* relS, int nthreads) {
    (void) (nthreads);

    int64_t result = 0;
    pthread_t tid[nthreads];
    arg_nl_t args[nthreads];
    uint64_t numperthr[2];
#ifndef NO_TIMING
    uint64_t timer1;
    ocall_startTimer(&timer1);
#endif
#ifdef PCM_COUNT
    ocall_set_system_counter_state("Start join phase");
#endif

    numperthr[0] = relR->num_tuples / nthreads;
    numperthr[1] = relS->num_tuples / nthreads;

    for (int i = 0; i < nthreads; i++) {
        args[i].my_tid = i;
        args[i].relR = relR->tuples + i * numperthr[0];
        args[i].relS = relS->tuples;
        args[i].numR = (i == (nthreads-1)) ?
                       (relR->num_tuples - i * numperthr[0]) : numperthr[0];
        args[i].numS = relS->num_tuples;
        args[i].result = 0;

        int rv = pthread_create(&tid[i], nullptr, nlj_thread, (void*)&args[i]);
        if (rv){
            logger(ERROR, "return code from pthread_create() is %d\n", rv);
            ocall_exit(-1);
        }
    }

    for (int i = 0; i < nthreads; i++) {
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }

#ifdef PCM_COUNT
    ocall_get_system_counter_state("Join", 0);
#endif

#ifndef NO_TIMING
    ocall_stopTimer(&timer1);
    print_timing(timer1, relR->num_tuples + relS->num_tuples, result);
#endif

    result_t * joinresult;
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = result;
    joinresult->nthreads = nthreads;
    return joinresult;
}

/* Vectorization version type: duplicate-outer nested loop join (Comparing to Zhou paper*) 
* Jingren Zhou, Kenneth A. Ross, Implementing Database Operations Using SIMD Instructions, SIGMOD '02
https://doi.org/10.1145/564691.564709 , ACM SIGMOD ’2002 June 4-6
*/

/* Adding SIMD AVX2 version of nested loop join. Load 8 keys into each vector.
* Problem: Too many calculations, so this SIMD version is still slower than the scalar version.
* But this code with int instructions is about 40% or more faster than with float. 
* The SSE code and AVX2 with float values are also much slower. So not included.
*/
void * NL_keys_thread(void * param) {
    arg_nl_t *args = (arg_nl_t*) param;
    uint64_t results = 0;
    
    //Temporary result vector of results.
    __m256i res_vector = _mm256_setzero_si256();
    //for ( int j=0; j < args->numS; j+=1){ 
    //    logger(INFO,"%d" ,args->relS[j].key);
    //}

    for (int32_t i = 0; i < args->numR; i++)
    {   
        __m256i  args_r_tmp = _mm256_set1_epi32(args->relR[i].key);
        uint64_t j = 0;
        for ( j=0; j < args->numS-8; j+=8)
        {   
            __m256i args_s_tmp =  _mm256_set_epi32(args->relS[j].key, args->relS[j+1].key, args->relS[j+2].key, args->relS[j+3].key,
                    args->relS[j+4].key, args->relS[j+5].key, args->relS[j+6].key, args->relS[j+7].key );

            __m256i eq = _mm256_cmpeq_epi32( args_r_tmp, args_s_tmp);
            res_vector = _mm256_add_epi32(eq, res_vector);
        }

        //This code part avoid seg fault if the table size isn't a multiply of eight.
        //This code part can be slower.
        //logger(INFO,"num_tuples: %d, i: %d\n", args->numS, j);
        if( j+8 != args->numS ){
            for(uint64_t m=j; m< args->numS; m++){
                if (args->relR[i].key == args->relS[m].key){
                    results++;
                }
            }
        }
        else{
            __m256i args_s_tmp =  _mm256_set_epi32(args->relS[j].key, args->relS[j+1].key, args->relS[j+2].key, args->relS[j+3].key,
                    args->relS[j+4].key, args->relS[j+5].key, args->relS[j+6].key, args->relS[j+7].key );
            __m256i eq = _mm256_cmpeq_epi32( args_r_tmp, args_s_tmp);
            res_vector = _mm256_add_epi32(eq, res_vector);
        }
        
    }
    results += -( _mm256_extract_epi32(res_vector, 0) + _mm256_extract_epi32(res_vector, 1)
    + _mm256_extract_epi32(res_vector, 2) + _mm256_extract_epi32(res_vector, 3)
    + _mm256_extract_epi32(res_vector, 4) + _mm256_extract_epi32(res_vector, 5) 
    + _mm256_extract_epi32(res_vector, 6) + _mm256_extract_epi32(res_vector, 7) );
    
    args->result = results;
    return nullptr;
}

/* Vectorization version type: duplicate-outer nested loop join (Comparing to Zhou paper*) 
* Jingren Zhou, Kenneth A. Ross, Implementing Database Operations Using SIMD Instructions, SIGMOD '02
https://doi.org/10.1145/564691.564709 , ACM SIGMOD ’2002 June 4-6
*/

/* Adding SIMD AVX2 version of nested loop join. Load 8 keys into each vector.
* Problem: Too many calculations, so this SIMD version is still slower than the scalar version.
* The SSE code and AVX2 with float values are also much slower. So not included. */
result_t* NL_keys (struct table_t* relR, struct table_t* relS, int nthreads) {
    (void) (nthreads);

    int64_t result = 0;
    pthread_t tid[nthreads];
    arg_nl_t args[nthreads];
    uint64_t numperthr[2];
#ifndef NO_TIMING
    uint64_t timer1;
    ocall_startTimer(&timer1);
#endif
#ifdef PCM_COUNT
    ocall_set_system_counter_state("Start join phase");
#endif

    numperthr[0] = relR->num_tuples / nthreads;
    numperthr[1] = relS->num_tuples / nthreads;

    for (int i = 0; i < nthreads; i++) {
        args[i].my_tid = i;
        args[i].relR = relR->tuples + i * numperthr[0];
        args[i].relS = relS->tuples;
        args[i].numR = (i == (nthreads-1)) ?
                       (relR->num_tuples - i * numperthr[0]) : numperthr[0];
        args[i].numS = relS->num_tuples;
        args[i].result = 0;

        int rv = pthread_create(&tid[i], nullptr, NL_keys_thread, (void*)&args[i]);
        if (rv){
            logger(ERROR, "return code from pthread_create() is %d\n", rv);
            ocall_exit(-1);
        }
    }

    for (int i = 0; i < nthreads; i++) {
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }

#ifdef PCM_COUNT
    ocall_get_system_counter_state("Join", 0);
#endif

#ifndef NO_TIMING
    ocall_stopTimer(&timer1);
    print_timing(timer1, relR->num_tuples + relS->num_tuples, result);
#endif

    result_t * joinresult;
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = result;
    joinresult->nthreads = nthreads;
    return joinresult;
}

void * NL_tuples_thread(void * param) {
    arg_nl_t *args = (arg_nl_t*) param;
    uint64_t results = 0;

    //Temporary result vector of results.
    __m256i res_vector = _mm256_setzero_si256();;

    /* Reading into undefined area, if the number of S tuples isn't a multiply of 4.
    for (int32_t i = 0; i < args->numR; i++)
    {   
        __m256i  args_r_tmp = _mm256_set1_epi32(args->relR[i].key);
        uint64_t j = 0;
        
        for ( j=0; j < args->numS; j += 4)
        {   
            __m256i args_s_tmp = _mm256_loadu_si256((__m256i const *) (args->relS + j) );
            __m256i eq = _mm256_cmpeq_epi32( args_r_tmp, args_s_tmp);
            res_vector = _mm256_add_epi32(eq, res_vector);
        }
    }*/

    
    //Use the code below instead. If dealing with 0 keys.     
    for (int32_t i = 0; i < args->numR; i++)
    {   
        __m256i  args_r_tmp = _mm256_set1_epi32(args->relR[i].key);
        uint64_t j = 0;
        
        for ( j=0; j < args->numS-4; j += 4)
        {   
            __m256i args_s_tmp = _mm256_loadu_si256((__m256i const *) (args->relS + j) );
            __m256i eq = _mm256_cmpeq_epi32( args_r_tmp, args_s_tmp);
            res_vector = _mm256_add_epi32(eq, res_vector);
        }
        
        //This code part deals with key_id = 0;
        //This code part can be slower.
        //logger(INFO,"num_tuples: %d, i: %d\n", args->numS, j);
        if( j+4 != args->numS ){
            //logger(INFO,"scalar loop: num_tuples: %d, j: %d\n", args->numS, j);
            for(uint64_t m=j; m< args->numS; m++){
                //logger(INFO,"num_tuples: %d, i: %d\n", args->numS, j);
                if (args->relR[i].key == args->relS[m].key){
                    results++;
                }
            }
        }
        else{
            __m256i args_s_tmp = _mm256_loadu_si256((__m256i const *) (args->relS + j) );
            __m256i eq = _mm256_cmpeq_epi32( args_r_tmp, args_s_tmp);
            res_vector = _mm256_add_epi32(eq, res_vector);
        }
    }
    
   
    //Add the extracted number of res_vector to results. results is of type uint64_t.
    //inline long long _mm256_extract_epi64(__m256i __X, int __N) -> N is the index.
    results += -( _mm256_extract_epi32(res_vector, 0) + _mm256_extract_epi32(res_vector, 2)
    + _mm256_extract_epi32(res_vector, 4) + _mm256_extract_epi32(res_vector, 6) );
    
    args->result = results;
    return nullptr;
}

/*  Adding SIMD AVX2 version of nested loop join. This time load 4 whole tuples into each vector.
*   Faster than scalar version, by higher number of tuples at least equal.
*   
*   Vectorization idea from Cagri Balkesen, Gustavo Alonso, Jens Teubner, and M. Tamer Özsu. 2013. 
*   Multi-core, main-memory joins: sort vs. hash revisited. Proc. VLDB Endow. 7, 1 (September 2013), 85–96. 
*   https://doi.org/10.14778/2732219.2732227  */
result_t* NL_tuples (struct table_t* relR, struct table_t* relS, int nthreads) {
    (void) (nthreads);

    int64_t result = 0;
    pthread_t tid[nthreads];
    arg_nl_t args[nthreads];
    uint64_t numperthr[2];
#ifndef NO_TIMING
    uint64_t timer1;
    ocall_startTimer(&timer1);
#endif
#ifdef PCM_COUNT
    ocall_set_system_counter_state("Start join phase");
#endif

    numperthr[0] = relR->num_tuples / nthreads;
    numperthr[1] = relS->num_tuples / nthreads;

    for (int i = 0; i < nthreads; i++) {
        args[i].my_tid = i;
        args[i].relR = relR->tuples + i * numperthr[0];
        args[i].relS = relS->tuples;
        args[i].numR = (i == (nthreads-1)) ?
                       (relR->num_tuples - i * numperthr[0]) : numperthr[0];
        args[i].numS = relS->num_tuples;
        args[i].result = 0;

        int rv = pthread_create(&tid[i], nullptr, NL_tuples_thread, (void*)&args[i]);
        if (rv){
            logger(ERROR, "return code from pthread_create() is %d\n", rv);
            ocall_exit(-1);
        }
    }

    for (int i = 0; i < nthreads; i++) {
        pthread_join(tid[i], NULL);
        result += args[i].result;
    }

#ifdef PCM_COUNT
    ocall_get_system_counter_state("Join", 0);
#endif

#ifndef NO_TIMING
    ocall_stopTimer(&timer1);
    print_timing(timer1, relR->num_tuples + relS->num_tuples, result);
#endif

    result_t * joinresult;
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = result;
    joinresult->nthreads = nthreads;
    return joinresult;
}



struct arg_inl_t {
    uint32_t my_tid;
    tuple_t * relR;
//    tuple_t * relS;

    uint32_t numR;
    uint32_t totalR;

    stx::btree<type_key, type_value> * indexS;

    uint32_t matches = 0;

};

void print_timing(uint64_t total_cycles, uint64_t numtuples, uint64_t join_matches,
                  uint64_t start, uint64_t end)
{
    double cyclestuple = (double) total_cycles / (double) numtuples;
    uint64_t time_usec = end - start;
    double throughput = numtuples / (double) time_usec;

    logger(ENCLAVE, "Total input tuples     : %lu", numtuples);
    logger(ENCLAVE, "Result tuples          : %lu", join_matches);
    logger(ENCLAVE, "Phase Join (cycles)    : %lu", total_cycles);
    logger(ENCLAVE, "Cycles-per-tuple       : %.4lf", cyclestuple);
    logger(ENCLAVE, "Total Runtime (us)     : %lu ", time_usec);
    logger(ENCLAVE, "Throughput (M rec/sec) : %.2lf", throughput);
}

void * inl_thread(void * param)
{
    uint32_t i, matches = 0;
    arg_inl_t * args = (arg_inl_t*) param;
//    uint32_t my_tid = args->my_tid;

    stx::btree<type_key, type_value> * index = args->indexS;

    // for each R scan S-index
    #pragma omp simd
    for (i = 0; i < args->numR; i++) {
        row_t r = args->relR[i];
        size_t count = index->count(r.key);
        //printf ("Here A %d\n", count);
        if (count) {
            auto it = index->find(r.key);
            for (size_t j = 0; j < count; j++) {
                matches++;
                //printf ("Here %d\n", it.data());
                it++;
                
            }
        }
    }
//    logger(INFO, "Thread-%d matches: %u", my_tid, matches);
    args->matches = matches;
    return nullptr;
}

result_t* INL (struct table_t* relR, struct table_t* relS, int nthreads) {
    uint64_t i, matches = 0;
    int rv;
    stx::btree<type_key, type_value> index;

    pthread_t tid[nthreads];
    arg_inl_t args[nthreads];
    uint64_t numperthr[2];

    uint64_t timer, start, end;

    numperthr[0] = relR->num_tuples / nthreads;
    numperthr[1] = relS->num_tuples / nthreads;

    // build index on S
    for (i = 0; i < relS->num_tuples; i++) {
        index.insert(std::make_pair(relS->tuples[i].key, relS->tuples[i].payload));
    }

    logger(DBG, "Index complete. Size: %zu", index.size());

    ocall_startTimer(&timer);
    ocall_get_system_micros(&start);
#ifdef PCM_COUNT
    ocall_set_system_counter_state("Start join phase");
#endif
    for (i = 0; i < nthreads; i++) {
        args[i].relR = relR->tuples + i * numperthr[0];

        args[i].numR = (i == (nthreads-1)) ?
                       (relR->num_tuples - i * numperthr[0]) : numperthr[0];
        args[i].totalR = relR->num_tuples;

        args[i].my_tid = i;
        args[i].indexS = &index;

        rv = pthread_create(&tid[i], nullptr, inl_thread, (void*)&args[i]);

        if (rv){
            logger(ERROR, "return code from pthread_create() is %d\n", rv);
            ocall_exit(-1);
        }
    }

    for (i = 0; i < nthreads; i++) {
        pthread_join(tid[i], nullptr);
        matches += args[i].matches;
    }
#ifdef PCM_COUNT
    ocall_get_system_counter_state("Join", 0);
#endif
    ocall_get_system_micros(&end);
    ocall_stopTimer(&timer);
    print_timing(timer, relR->num_tuples + relS->num_tuples, matches, start, end);

    result_t * joinresult;
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}
