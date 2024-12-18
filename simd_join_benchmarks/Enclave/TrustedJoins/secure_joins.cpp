#include "Enclave_t.h"
#include "data-types.h"
#include "nested_loop_join.h"
#include "no_partitioning_join.h"
#include "rhobli.h"
#include "radix_join.h"
#include "CHTJoinWrapper.hpp"
#include "radix_sortmerge_join.h"
#include <sgx_tseal.h>
#include <mway/sortmergejoin_multiway.h>

#include "parallel_sortmerge_join.h"
#include "util.h"
#include "rho_atomic/radix_join_atomic.h"

extern char aad_mac_text[256];
//extern result_t* rhobli_join (struct table_t *relR, struct table_t *relS, int nthreads);

void print_relation(relation_t *rel, uint32_t num, uint32_t offset)
{
    logger(DBG, "****************** Relation sample ******************");
    for (uint32_t i = offset; i < rel->num_tuples && i < num + offset; i++)
    {
        logger(DBG, "%u -> %u", rel->tuples[i].key, rel->tuples[i].payload);
    }
    logger(DBG, "******************************************************");
}

result_t* CHT(struct table_t * relR, struct table_t * relS, int nthreads)
{
    join_result_t  join_result = CHTJ<7>(relR, relS, nthreads);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTfsimd(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    //join_result_t  join_result = CHTJ<7>(relR, relS, nthreads);
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTFSIMD);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTopt(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    //join_result_t  join_result = CHTJ<7>(relR, relS, nthreads);
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTOPT);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTopt2(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    //join_result_t  join_result = CHTJ<7>(relR, relS, nthreads);
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTOPT2);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTopt3(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    //join_result_t  join_result = CHTJ<7>(relR, relS, nthreads);
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTOPT3);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPa1(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPA1);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPa2(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPA2);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPa3(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPA3);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

//Not finished. TO DO in future
result_t* CHTPa3sl(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPA3SL);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPa4(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPA4);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTb1(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTB1);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTb1s(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTB1S);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTb2(struct table_t * relR, struct table_t * relS, int nthreads)
{
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTB2);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTb2s(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTB2S);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPr1(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPR1);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

result_t* CHTPr2(struct table_t * relR, struct table_t * relS, int nthreads)
{   
    join_result_t  join_result = CHTJ_simd<7>(relR, relS, nthreads, CHTPR2);
    result_t * joinresult;
    
    joinresult = (result_t *) malloc(sizeof(result_t));
    joinresult->totalresults = join_result.matches;
    joinresult->nthreads = nthreads;
    return joinresult;
}

static struct algorithm_t sgx_algorithms[] = {
    // Here is the list of used implementetions.
    {"NL", NL},
    {"NL_keys", NL_keys}, 
    {"NL_tuples", NL_tuples}, 

    {"RHT", RHT},
    {"PRHO", PRHO},
        
    {"CHT", CHT},
    //Specific vectorized implementations of CHT.
    {"CHTfsimd", CHTfsimd},
    {"CHTopt", CHTopt},
    {"CHTopt2", CHTopt2},
    {"CHTopt3", CHTopt3},
    {"CHTPa1", CHTPa1},
    {"CHTPa2", CHTPa2},
    {"CHTPa3", CHTPa3},
    {"CHTPa3sl", CHTPa3sl},
    {"CHTPa4", CHTPa4},
    {"CHTb1", CHTb1},
    {"CHTb1s", CHTb1s},
    {"CHTb2", CHTb2},
    {"CHTb2s", CHTb2s},
    {"CHTPr1", CHTPr1},
    {"CHTPr2", CHTPr2},

    //Not used in the experiment.
    {"PHT", PHT},
    {"NPO_st", NPO_st},
    {"INL", INL},
    {"RJ", RJ},
    {"RHO", RHO},
    {"PSM", PSM},
    {"RSM", RSM},
    {"MWAY", MWAY},
};

uint8_t* unseal_rel(const uint8_t *sealed_rel, size_t size)
{
    uint32_t mac_text_len = sgx_get_add_mac_txt_len((const sgx_sealed_data_t *)sealed_rel);
    uint32_t decrypt_data_len = sgx_get_encrypt_txt_len((const sgx_sealed_data_t *) sealed_rel);
    if (mac_text_len == UINT32_MAX || decrypt_data_len == UINT32_MAX)
    {
//        return SGX_ERROR_UNEXPECTED;
        return nullptr;
    }
    if (mac_text_len > size || decrypt_data_len > size)
    {
        // SGX_ERROR_INVALID_PARAMETER
        return 0;
    }
    uint8_t *de_mac_text = (uint8_t *) malloc(mac_text_len);
    if (de_mac_text == nullptr)
    {
        //SGX_ERROR_OUT_OF_MEMORY
        return nullptr;
    }
    uint8_t *decrypt_data = (uint8_t *) malloc(decrypt_data_len);
    if (decrypt_data == nullptr)
    {
        //SGX_ERROR_OUT_OF_MEMORY
        free(de_mac_text);
        return nullptr;
    }
    sgx_status_t ret = sgx_unseal_data((const sgx_sealed_data_t *)sealed_rel,
                                       de_mac_text,
                                       &mac_text_len,
                                       decrypt_data,
                                       &decrypt_data_len);
    if (ret != SGX_SUCCESS)
    {
        free(de_mac_text);
        free(decrypt_data);
        return nullptr;
    }

    if (memcmp(de_mac_text, aad_mac_text, strlen(aad_mac_text)))
    {
        //SGX_ERROR_UNEXPECTED
        return nullptr;
    }
    free(de_mac_text);
    return decrypt_data;
}

uint8_t* sealed_buf;

uint32_t seal_relation(relation_t * rel, uint32_t seal_chunk_size)
{
    uint32_t output_size = sizeof(relation_t) + rel->num_tuples * (sizeof(row_t));
    logger(DBG, "Size of unsealed data = %.2lf MB", B_TO_MB(output_size));
    if (seal_chunk_size != 0) {
        logger(DBG, "Size of seal chunk = %d kB", seal_chunk_size);
        uint32_t chunk_size, offset = 0, temp_sealed_size;
        uint32_t * tmp;
        uint32_t available_space = seal_chunk_size * 1024 - sizeof(sgx_sealed_data_t) - (uint32_t) strlen(aad_mac_text); // reverse sgx_calc_sealed_data_size
        uint32_t chunks = (output_size - 1) / available_space + 1; // get the number of chunks rounded up
        logger(DBG, "Available space in seal_chunk_size: %d, total chunks: %d", available_space, chunks);
        uint8_t * temp_sealed_buf = (uint8_t*) malloc(seal_chunk_size * 1024);
        for (uint32_t i = 0; i < chunks; i ++) {
            chunk_size = (i == (chunks - 1)) ?
                                  (output_size - i * available_space) : available_space;
            tmp = (uint32_t *) rel;
            tmp += offset;
            temp_sealed_size = sgx_calc_sealed_data_size((uint32_t)strlen(aad_mac_text), chunk_size);
//            logger(DBG, "Chunk %d: size=%d, offset=%d, sealed_size=%d", i+1, chunk_size, offset, temp_sealed_size);
            sgx_status_t err = sgx_seal_data((uint32_t) strlen(aad_mac_text),
                                             (const uint8_t *) aad_mac_text,
                                             chunk_size,
                                             (uint8_t*) tmp,
                                             temp_sealed_size,
                                             (sgx_sealed_data_t *) temp_sealed_buf);
            if (err == SGX_SUCCESS) {
                sealed_buf = temp_sealed_buf;
                offset += chunk_size;
            } else {
                logger(ERROR, "[Chunk %d/%d] sealing error: %d", i+1, chunks, err);
                ocall_exit(1);
            }
        }
        logger(DBG, "Sealing relation in chunks successful");
        return seal_chunk_size;
    } else {
        uint32_t sealed_data_size = sgx_calc_sealed_data_size((uint32_t)strlen(aad_mac_text), output_size);
        logger(DBG, "Size of seal chunk = %d kB", sealed_data_size);
        uint8_t* temp_sealed_buf = (uint8_t *) malloc(sealed_data_size);
        sgx_status_t err = sgx_seal_data((uint32_t)strlen(aad_mac_text),
                                         (const uint8_t *) aad_mac_text,
                                         output_size,
                                         (uint8_t*) rel,
                                         sealed_data_size,
                                         (sgx_sealed_data_t *) temp_sealed_buf);
        if (err == SGX_SUCCESS)
        {
            logger(DBG, "Sealing relation successful");
            sealed_buf = temp_sealed_buf;
            return sealed_data_size;
        }
    }
    return 0;
}

sgx_status_t ecall_get_sealed_data(uint8_t* sealed_blob, uint32_t data_size)
{
    if (sealed_buf == nullptr)
    {
        printf("Nothing to return...");
        return SGX_ERROR_UNEXPECTED;
    }
    memcpy(sealed_blob, sealed_buf, data_size);
    free(sealed_buf);
    return SGX_SUCCESS;
}

result_t* ecall_join(struct table_t * relR, struct table_t * relS, char *algorithm_name, int nthreads)
{
    int i =0, found = 0;
    algorithm_t *algorithm = nullptr;
    while(sgx_algorithms[i].join)
    {
        if (strcmp(algorithm_name, sgx_algorithms[i].name) == 0)
        {
            found = 1;
            algorithm = &sgx_algorithms[i];
            break;
        }
        i++;
    }
    if (found == 0)
    {
        printf("Algorithm not found: %s", algorithm_name);
        ocall_exit(EXIT_FAILURE);
    }
    struct rusage_reduced_t usage;
    usage.ru_utime_sec = 0;
    usage.ru_utime_usec = 0;
    usage.ru_stime_sec = 0;
    usage.ru_stime_usec = 0;
    usage.ru_minflt = 0;
    usage.ru_majflt = 0;
    usage.ru_nvcsw = 0;
    usage.ru_nivcsw = 0;
    ocall_getrusage(&usage, 0);
    result_t *res = algorithm->join(relR, relS, nthreads);
    ocall_getrusage(&usage, 1);

    return res;
}

relation_t *to_relation(result_t *result) {
#ifndef JOIN_MATERIALIZE
    logger(WARN, "JOIN_MATERIALIZE not defined. to_relation might fail.");
#endif
    relation_t * output = (relation_t*) malloc(sizeof(relation_t));
    malloc_check(output);
    output->tuples = (tuple_t*) malloc(sizeof(tuple_t)*result->totalresults);
    malloc_check(output->tuples);
    uint64_t items = 0;
    for (int i = 0; i < result->nthreads; i++)
    {
        output_list_t *list = result->resultlist[i].results;
        while (list != nullptr)
        {
            memcpy(output->tuples + items, list, sizeof(type_key)+sizeof(type_value));
            items++;
            list = list->next;
        }
    }
    output->num_tuples = result->totalresults;
    output->sorted = 0;
    output->ratio_holes = 0;
    logger(DBG, "to_relation check if %lu == %lu", result->totalresults, items);
    return output;
}

uint32_t ecall_join_sealed_tables(const uint8_t *sealed_r,
                                  size_t size_r,
                                  const uint8_t *sealed_s,
                                  size_t size_s,
                                  char *algorithm,
                                  int nthreads,
                                  uint32_t seal_chunk_size)
{
    uint64_t seal_timer = 0, unseal_timer = 0, join_timer = 0;
    relation_t * output = nullptr;
    uint32_t sealed_data_size = 0;
    ocall_startTimer(&unseal_timer);
    struct table_t * relR = (struct table_t *) unseal_rel(sealed_r, size_r);
    printf("Unseal R successful");
    struct table_t * relS = (struct table_t *) unseal_rel(sealed_s, size_s);
    printf("Unseal S successful");
    ocall_stopTimer(&unseal_timer);
    ocall_startTimer(&join_timer);
    result_t* result = ecall_join(relR, relS, algorithm, nthreads);
    if (strcmp(algorithm, "RHO_seal_buffer") != 0) {
        output = to_relation(result);
    }
    ocall_stopTimer(&join_timer);

    if (strcmp(algorithm, "RHO_seal_buffer") != 0) {
        ocall_startTimer(&seal_timer);
        uint64_t start, end;
        ocall_get_system_micros(&start);
        sealed_data_size = seal_relation(output, seal_chunk_size);
        ocall_get_system_micros(&end);
        ocall_stopTimer(&seal_timer);
        logger(INFO, "seal_micros = %lu", (end-start));
    }
    logger(INFO, "uns_timer = %lu (%.2lf%%)", unseal_timer, (double) unseal_timer*100/(unseal_timer + seal_timer));
    logger(INFO, "s_timer   = %lu (%.2lf%%)", seal_timer, (double) seal_timer*100/(unseal_timer + seal_timer));
    logger(INFO, "seal_timer = %lu", (unseal_timer + seal_timer));
    logger(INFO, "join_timer = %lu", join_timer);


    free(relR);
    free(relS);
    //TODO: free result
    return sealed_data_size;
}

uint32_t ecall_three_way_join_sealed_tables(const uint8_t *sealed_r,
                                            size_t size_r,
                                            const uint8_t *sealed_s,
                                            size_t size_s,
                                            const uint8_t *sealed_t,
                                            size_t size_tt,
                                            char *algorithm,
                                            int nthreads,
                                            uint32_t seal_chunk_size) {
    uint64_t seal_timer, unseal_timer, join1_timer, join2_timer;
    ocall_startTimer(&unseal_timer);
    relation_t * relR = (relation_t *) unseal_rel(sealed_r, size_r);
    printf("Unseal R successful");
    relation_t * relS = (relation_t *) unseal_rel(sealed_s, size_s);
    printf("Unseal S successful");
    relation_t * relT = (relation_t *) unseal_rel(sealed_t, size_tt);
    printf("Unseal T successful");
    ocall_stopTimer(&unseal_timer);

    ocall_startTimer(&join1_timer);
    result_t* result1 = ecall_join(relR, relS, "RHO", nthreads);
    relation_t * output1 = to_relation(result1);
    ocall_stopTimer(&join1_timer);

    ocall_startTimer(&join2_timer);
    result_t* result2 = ecall_join(relT, output1, "RHO", nthreads);
    relation_t * output2 = to_relation(result2);
    ocall_stopTimer(&join2_timer);

    ocall_startTimer(&seal_timer);
    uint64_t start, end;
    ocall_get_system_micros(&start);
    uint32_t sealed_data_size = seal_relation(output2, seal_chunk_size);
    ocall_get_system_micros(&end);
//    free(sealed_buf);
//    seal_relation(relS);
    ocall_stopTimer(&seal_timer);
    logger(INFO, "uns_timer = %lu (%.2lf%%)", unseal_timer, (double) unseal_timer*100/(unseal_timer + seal_timer));
    logger(INFO, "s_timer   = %lu (%.2lf%%)", seal_timer, (double) seal_timer*100/(unseal_timer + seal_timer));
    logger(INFO, "seal_timer = %lu", (unseal_timer + seal_timer));
    logger(INFO, "join1_timer = %lu", join1_timer);
    logger(INFO, "join2_timer = %lu", join2_timer);
    logger(INFO, "seal_micros = %lu", (end-start));
    logger(INFO, "total sealing share = %.2lf%%",
           (double) (unseal_timer+seal_timer)*100/(unseal_timer+seal_timer+join1_timer+join2_timer));

    free(relR);
    free(relS);
    free(relT);
    //TODO: free result
    return sealed_data_size;
}
