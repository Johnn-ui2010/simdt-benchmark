#include <vector>
#include <atomic>
#include <iostream>

#include "data-types.h"
#include "CHT.hpp"
#include "Barrier.h"

#include <cstring> //memset

#ifdef NATIVE_COMPILATION
#include "Logger.h"
#include "native_ocalls.h"
#include <immintrin.h> //AVX 1 and 2 (256bit) and AVX512 (includes all SSE, too)

#else
#include "Enclave_t.h"
#include "Enclave.h"

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

//From CHTJoin.hpp
class CHTPartitionQueue{
private:
    std::atomic<int> counter;
    int nPart;
public:
    CHTPartitionQueue(int _nPart):counter(0),nPart(_nPart) {
    }

    int getPartition() {
        return counter++;
    }

    bool empty() {
        return counter.load() == nPart;
    }
};

const unsigned int TUPLES_PER_CACHELINE= 64 / sizeof(tuple_t);

typedef struct {
    tuple_t tuples[TUPLES_PER_CACHELINE - 1];
    uint64_t target;
} cht_cacheline_t;


#define STREAM_UNIT 4

class CHTJoin {
private:
    int nthreads;
    int npart;
    relation_t *relR;
    relation_t *relS;
    tuple_t *partBuffer;
    uint64_t **hist;
    uint64_t **dst;
    CHTPartitionQueue *partitions;
    CHT ht;
    const intkey_t MASK;
    const intkey_t SHIFT;
    uint64_t *matches;
    uint64_t *checksums;
    uint64_t time_usec;
    uint64_t partitionTime;
    /*
         * timer1 - entire algorithm
         * timer2 - partition phase
         * timer3 - build phase
         * timer4 - probe phase
         * */
//    uint64_t start, end, timer1, timer2, timer3, timer4;
    struct timers_t timers;

    PThreadLockCVBarrier *barrier;

    /*
    Standard radix partition
    */
    void radix_partition(const tuple_t *input,
                         tuple_t *output,
                         uint64_t *histogram,
                         size_t numTuples) {

        __attribute__((aligned(64))) cht_cacheline_t buffers[npart];
        //logger(INFO,"TUPLES_PER_CACHELINE: %d, sizeof(tuple_t): %d \n", TUPLES_PER_CACHELINE, sizeof(tuple_t));

        // we have to make sure that the size of each partitioning in terms of elements is a multiple of 4 (since 4 tuples fit into 32 bytes = 256 bits).
        for(int i = 0; i < npart; ++i){
            buffers[i].target = histogram[i];
        }

        __attribute__((aligned(64))) uint64_t bucket_num = 0;
        for(uint64_t j = 0; j < numTuples; ++j){
            bucket_num = (input[j].key >> SHIFT) & MASK;
            int slot=buffers[bucket_num].target & (TUPLES_PER_CACHELINE - 1); // Some bits are set?

            if(slot == TUPLES_PER_CACHELINE - 1) {

                uint64_t targetBkp=buffers[bucket_num].target- (TUPLES_PER_CACHELINE-1);
                //logger(INFO,"targetBkp: %d, bucket_num: %d, TUPLES_PER_CACHELINE: %d\n", targetBkp, bucket_num, TUPLES_PER_CACHELINE);
                buffers[bucket_num].tuples[slot]=input[j];
                for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
//                    _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                    memcpy(output + targetBkp, buffers[bucket_num].tuples + b, 32);
                    targetBkp += STREAM_UNIT;
                }
                buffers[bucket_num].target=targetBkp;

            } else {
                buffers[bucket_num].tuples[slot] = input[j];
                buffers[bucket_num].target++;
            }
        }

        barrier->Arrive();

        for (int i = npart - 1; i >= 0; i--) {
            uint32_t slot  = (uint32_t)buffers[i].target;
            uint32_t sz    = (slot) & (TUPLES_PER_CACHELINE - 1);
            slot          -= sz;
            uint32_t startPos = (slot < histogram[i]) ? ((uint32_t)histogram[i] - slot) : 0;
            for(uint32_t j = startPos; j < sz; j++) {
                output[slot+j]  = buffers[i].tuples[j];
            }
        }
    }

    //Using standard radix_partition with histogram.
    //The function creates npart=128 partitions.
    void partition(int threadID, relation_t chunkR)
    {
        const tuple_t * tupleR = chunkR.tuples;
        uint64_t * my_hist     = hist[threadID];
        uint64_t * my_dst      = dst[threadID];
        uint64_t sum           = 0;

        /* Compute the hash index for each key. .
        *  Count the number of key elements for each hash index.
        *  "Store" the results into a histogram aka 2D array.
        */

        //logger(INFO, "SHIFT= %d, MASK=%d \n", SHIFT,MASK);
        for (size_t i = 0; i < chunkR.num_tuples; ++i) {
            intkey_t hk = (hashKey(tupleR[i].key) >> SHIFT) & MASK;

            my_hist[hk]++;
            //if(i<5){
            //        logger(INFO,"Before: i: %d, hash_key: %d, my_hist[hash_keys[j]: %d ", i,  tupleR[i].key,my_hist[hk]);
            //}
        }


        for (int i = 0; i < npart; ++i) {
            sum += my_hist[i];
            my_hist[i] = sum;
            //logger(INFO, "my_hist[i]: %d, no. i=%d\n",my_hist[i],i);
        }
        
        barrier->Arrive();

        //my_dst[] is assigning for all the threads and partitions.
        for (int i = 0; i < threadID; ++i) {
            for (int j = 0; j < npart; ++j) {
                my_dst[j] += hist[i][j];
                
            }   
        }
        for (int i = threadID; i < nthreads; ++i) {
            for (int j = 1; j < npart; ++j) {
                my_dst[j] += hist[i][j-1];
            }
        }

        //for(int a = 0; a < npart; a++){
        //    logger(INFO,"%ld ", my_dst[a]);
        //}

        radix_partition(chunkR.tuples, partBuffer, my_dst, chunkR.num_tuples);
    }

    void build(int threadID, int part)
    {
        (void) (threadID);
        tuple_t * tuples= partBuffer + dst[0][part];
        const uint64_t num_tuples = part == npart - 1 ? relR->num_tuples - dst[0][part] : dst[0][part + 1] - dst[0][part];
        //logger(INFO,"num_tuples = %ld, dst[0][part] = %ld, part = %ld \n", num_tuples, dst[0][part], part );
        
        //Set bit in the bitmap for each key.
        for (uint64_t i = 0; i < num_tuples; ++i) {
            ht.setBit(tuples[i].key);
		}

        ht.computePartPopCount(part, (uint32_t)dst[0][part]);

        for (uint64_t i = 0; i < num_tuples; ++i)
            //Insert into the non-empty array.
            ht.setTuple(tuples[i]);
    }

    void probe(int threadID, relation_t chunkS)
    {
        uint64_t match = 0;
        uint64_t checksum = 0;
 		tuple_t *tupleS = chunkS.tuples;
		const size_t numS = chunkS.num_tuples;
		const size_t batchStartUpperBound = numS - PROBE_BATCH_SIZE;
        
        //logger(INFO,"Before probe: match: %d, checksum: %d, threadID: %d\n", match, checksum, threadID);

        for (size_t i = 0; i <= batchStartUpperBound; i += PROBE_BATCH_SIZE)
		{
			ht.batch_probe(tupleS + i, match, checksum);
		}

		const size_t leftOver = numS % PROBE_BATCH_SIZE;
		tupleS += numS - leftOver;

        //logger(INFO,"After batch probes: match: %d, checksum: %d, threadID: %d, leftover %d\n", match, checksum, threadID,leftOver);
        //logger(INFO,"match: %d, first_hit_probes = %d, second_hit_probes = %d \n \n",match, ht.first_hit_probes, ht.second_hit_probes);

        //This is the scalar part for leftovers.
        for (size_t i = 0; i < leftOver; ++i)
        {
            tuple_t *foundTuple = ht.probe(tupleS[i].key);
            if (foundTuple)
            {
                match++;
                //checksum += foundTuple->payload + tupleS[i].payload;
            }
        }
        matches[threadID] = match;
        checksums[threadID] = checksum;
        //logger(INFO,"After probes: match: %d, checksum: %d, threadID: %d\n", match, checksum, threadID);
    }

public:

    CHTJoin(int _nthreads, int _npart, relation_t *_relR, relation_t *_relS, tuple_t * _partBuffer) : nthreads(_nthreads),
                                                     npart(_npart), relR(_relR) , relS(_relS),partBuffer(_partBuffer),
                                                     partitions(new CHTPartitionQueue(npart)),
                                                     ht(relR->num_tuples, nthreads, npart),
                                                     MASK(npart-1),SHIFT(__builtin_ctz(Utils::nextPowerOfTwo(relR->num_tuples)) - __builtin_ctz(npart)),
                                                     matches(new uint64_t[nthreads]), checksums(new uint64_t[nthreads])
    {
        barrier = new PThreadLockCVBarrier(nthreads);
        hist = (uint64_t **) malloc(nthreads * sizeof(uint64_t *));
        dst = (uint64_t **) malloc(nthreads * sizeof(uint64_t *));
    }


    void join(int threadID) {
        //Init the CHT.
        ht.init(threadID);

        hist[threadID] = (uint64_t*) calloc(npart, sizeof(uint64_t));
        dst[threadID] = (uint64_t*) calloc(npart, sizeof(uint64_t));

        barrier->Arrive();

        if (threadID == 0) {
            ocall_get_system_micros(&timers.start);
            ocall_startTimer(&timers.total);
            ocall_startTimer(&timers.timer1);
#ifdef PCM_COUNT
            ocall_set_system_counter_state("Start Partition");
#endif
        }

        // partition phase
        uint64_t numRPerThread = relR->num_tuples / nthreads;
        relation_t myChunkR;
        myChunkR.tuples = relR->tuples + numRPerThread * threadID;
        myChunkR.num_tuples = threadID == nthreads - 1 ? (uint32_t)(relR->num_tuples - numRPerThread * threadID) : (uint32_t)numRPerThread;

        partition(threadID, myChunkR);

        barrier->Arrive();

        if (threadID == 0) {
            ocall_get_system_micros(&timers.end);
            ocall_stopTimer(&timers.timer1);
            ocall_startTimer(&timers.timer2);
            partitionTime=timers.end-timers.start;
		}

        // build phase
        int partID;
        while ((partID = partitions->getPartition()) < npart)
        {
            build(threadID, partID);
        }

        barrier->Arrive();

        if (threadID == 0)
        {
#ifdef PCM_COUNT
            ocall_get_system_counter_state("Partition", 0);
            ocall_set_system_counter_state("Join");
#endif
            ocall_stopTimer(&timers.timer2);
            ocall_startTimer(&timers.timer3);
        }
      //  std::cerr << "build finished" << std::endl;

        // probe phase
        uint64_t numSperThread = relS->num_tuples / nthreads;
        relation_t myChunkS;
        myChunkS.tuples = relS->tuples + numSperThread * threadID;
        myChunkS.num_tuples = threadID == nthreads - 1 ? (uint32_t)(relS->num_tuples - numSperThread * threadID) : (uint32_t)numSperThread;
        probe(threadID, myChunkS);

        barrier->Arrive();

        if (threadID == 0) {
#ifdef PCM_COUNT
            ocall_get_system_counter_state("Join", 0);
#endif
	//		std::cerr << "probe finished" << std::endl;
	        ocall_stopTimer(&timers.timer3);
            ocall_get_system_micros(&timers.end);
            ocall_stopTimer(&timers.total);
            time_usec = timers.end - timers.start;

        }
    }

    join_result_t get_join_result()
    {
        uint64_t m = 0;
        uint64_t c = 0;
        for (int i = 0; i < nthreads; ++i) {
            m += matches[i];
            c += checksums[i];
        }
        join_result_t result;
        result.time_usec = time_usec;
        result.part_usec=partitionTime;
        result.join_usec=time_usec-partitionTime;
        result.matches = m;
        result.checksum = c;
        return result;
    }

    struct timers_t get_timers()
    {
        return timers;
    }

    ~CHTJoin()
    {
		for (int i = 0; i < nthreads; ++i) {
			free(dst[i]);
			free(hist[i]);
		}
		free(dst);
		free(hist);
        delete[] checksums;
        delete[] matches;
        delete 	 barrier;
        delete   partitions;
    }
};

//CHANGE the name later to CHTJoin_simd_histogram1
class CHTJoin_simd {
private:
    int nthreads;
    int npart;
    relation_t *relR;
    relation_t *relS;
    tuple_t *partBuffer;
    uint64_t **hist;
    uint64_t **dst;
    CHTPartitionQueue *partitions;
    CHT ht;
    const intkey_t MASK;
    const intkey_t SHIFT;
    uint64_t *matches;
    uint64_t *checksums;
    uint64_t time_usec;
    uint64_t partitionTime;
    /*
         * timer1 - entire algorithm
         * timer2 - partition phase
         * timer3 - build phase
         * timer4 - probe phase
         * */
//    uint64_t start, end, timer1, timer2, timer3, timer4;
    struct timers_t timers;

    //Choose a SIMD implementation type for CHT.
    int impl_type;

    PThreadLockCVBarrier *barrier;

    /*
    Standard radix partition
    */
    void radix_partition(const tuple_t *input,
                         tuple_t *output,
                         uint64_t *histogram,
                         size_t numTuples) {

        __attribute__((aligned(64))) cht_cacheline_t buffers[npart];
        //logger(INFO,"TUPLES_PER_CACHELINE: %d, sizeof(tuple_t): %d \n", TUPLES_PER_CACHELINE, sizeof(tuple_t));

        // we have to make sure that the size of each partitioning in terms of elements is a multiple of 4 (since 4 tuples fit into 32 bytes = 256 bits).
        for(int i = 0; i < npart; ++i){
            buffers[i].target = histogram[i];
        }
        __attribute__((aligned(64))) uint64_t bucket_num = 0;
        
        if(impl_type == CHTFSIMD || impl_type == CHTPA3 || impl_type == CHTOPT ){
            //logger(INFO,"Welcome in CHTPA3 loop.");
            __m256i MASK_VEC = _mm256_set1_epi32(MASK);
            __m256i TPC_VEC = _mm256_set1_epi64x(TUPLES_PER_CACHELINE - 1);
            uint64_t j = 0;

            for(j = 0; j < numTuples -4; j+=4){
                // Equivalent to //bucket_num = (input[j].key >> SHIFT) & MASK; // Can be vectorized.                
                //
                __m256i bn_vec = _mm256_loadu_si256((__m256i const *)(input + j));
                bn_vec = _mm256_srli_epi32( bn_vec, SHIFT);
                bn_vec = _mm256_and_si256( bn_vec, MASK_VEC);

                uint64_t bucket_numbers[4];
                _mm256_storeu_pd(( double *) (bucket_numbers), _mm256_castsi256_pd (bn_vec) ); 
               
                /*if(j<64){ //Seems to be correct.
                    logger(INFO,"bucket_numbers[0] = %ld, bucket_numbers[1] = %ld, bucket_numbers[2] = %ld, bucket_numbers[3] = %ld \n",
                            bucket_numbers[0], bucket_numbers[1], bucket_numbers[2], bucket_numbers[3]);
                    
                }*/

                for(uint64_t m = 0; m < 4; m++){
                    bucket_num = (input[j+m].key >> SHIFT) & MASK;
                    int slot=buffers[bucket_numbers[m]].target & (TUPLES_PER_CACHELINE - 1);

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        uint64_t targetBkp=buffers[bucket_numbers[m]].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_numbers[m]].tuples[slot]=input[j+m];

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            //Replaced with _mm256_stream_si256 with memcpy (Both suggested original solutions).
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_numbers[m]].target=targetBkp;

                    } 
                    else {
                        buffers[bucket_numbers[m]].tuples[slot] = input[j+m];
                        buffers[bucket_numbers[m]].target++;
                    }
                }
            }
            //
            if(  j+4 != numTuples ){
                for(uint64_t m=j; m< numTuples; m++){
                    bucket_num = (input[m].key >> SHIFT) & MASK; // Can be vectorized.
                    int slot=buffers[bucket_num].target & (TUPLES_PER_CACHELINE - 1); // Some bits are simd_set? // Can be vectorized.

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        uint64_t targetBkp=buffers[bucket_num].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_num].tuples[slot]=input[m]; //Can be vectorized.

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_num].target=targetBkp;

                    } 
                    else {
                        buffers[bucket_num].tuples[slot] = input[m];
                        buffers[bucket_num].target++;
                    }
                }
            }
            else {
                // Equivalent to //bucket_num = (input[j].key >> SHIFT) & MASK; // Can be vectorized.                
                __m256i bn_vec = _mm256_loadu_si256((__m256i const *)(input + j));

                bn_vec = _mm256_srli_epi32( bn_vec, SHIFT);
                bn_vec = _mm256_and_si256( bn_vec, MASK_VEC);

                uint64_t bucket_numbers[4];
                _mm256_storeu_pd(( double *) (bucket_numbers), _mm256_castsi256_pd (bn_vec) ); 
  
                for(uint64_t m = 0; m < 4; m++){
                    bucket_num = (input[j+m].key >> SHIFT) & MASK;
                    int slot=buffers[bucket_numbers[m]].target & (TUPLES_PER_CACHELINE - 1);

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        uint64_t targetBkp=buffers[bucket_numbers[m]].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_numbers[m]].tuples[slot]=input[j+m]; //Can be vectorized.

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_numbers[m]].target=targetBkp;
                    } 
                    else {
                        buffers[bucket_numbers[m]].tuples[slot] = input[j+m];
                        buffers[bucket_numbers[m]].target++;
                    }
                }
            }        
        }
        // [TO DO] Currently no difference to CHTPA3. How to vectorize slot?
        else if(impl_type == CHTPA3SL){
            //logger(INFO,"Welcome in CHTPA3 loop.\n");
            __m256i MASK_VEC = _mm256_set1_epi32(MASK);
            __m256i TPC_VEC = _mm256_set1_epi64x(TUPLES_PER_CACHELINE - 1);
            uint64_t j = 0;

            for(j = 0; j < numTuples -4; j+=4){
                // Equivalent to //bucket_num = (input[j].key >> SHIFT) & MASK; // Can be vectorized.                
                //
                __m256i bn_vec = _mm256_loadu_si256((__m256i const *)(input + j));
                bn_vec = _mm256_srli_epi32( bn_vec, SHIFT);
                bn_vec = _mm256_and_si256( bn_vec, MASK_VEC);

                uint64_t bucket_numbers[4];
                _mm256_storeu_pd(( double *) (bucket_numbers), _mm256_castsi256_pd (bn_vec) ); 
               

                // [TO DO]
                //__m256i slot_vec = _mm256_set_epi64x(buffers[bucket_numbers[0]].target, buffers[bucket_numbers[1]].target,
                // buffers[bucket_numbers[2]].target, buffers[bucket_numbers[3]].target);
                //slot_vec = _mm256_and_si256( slot_vec, TPC_VEC);
                //uint64_t slot_numbers[4];
                //_mm256_storeu_pd(( double *) (slot_numbers), _mm256_castsi256_pd (slot_vec) ); 

                //How to vectorize this part?

                for(uint64_t m = 0; m < 4; m++){
                    bucket_num = (input[j+m].key >> SHIFT) & MASK;
                    int slot=buffers[bucket_numbers[m]].target & (TUPLES_PER_CACHELINE - 1);
                    //Vectorizing AND by doing two times simd_set() and have to deal with collisions. Is it worth?

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        //Vectorizing subtraction by doing two times simd_set(), is it worth?
                        uint64_t targetBkp=buffers[bucket_numbers[m]].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_numbers[m]].tuples[slot]=input[j+m];

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            //Replaced with _mm256_stream_si256 with memcpy (Both suggested original solutions).
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_numbers[m]].target=targetBkp;

                    } 
                    else {
                        buffers[bucket_numbers[m]].tuples[slot] = input[j+m];
                        buffers[bucket_numbers[m]].target++;
                    }
                }
            }
            //
            if(  j+4 != numTuples ){
                for(uint64_t m=j; m< numTuples; m++){
                    bucket_num = (input[m].key >> SHIFT) & MASK; // Can be vectorized.
                    int slot=buffers[bucket_num].target & (TUPLES_PER_CACHELINE - 1); // Some bits are simd_set? // Can be vectorized.

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        uint64_t targetBkp=buffers[bucket_num].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_num].tuples[slot]=input[m]; //Can be vectorized.

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_num].target=targetBkp;

                    } 
                    else {
                        buffers[bucket_num].tuples[slot] = input[m];
                        buffers[bucket_num].target++;
                    }
                }
            }
            else {
                // Equivalent to //bucket_num = (input[j].key >> SHIFT) & MASK; // Can be vectorized.                
                __m256i bn_vec = _mm256_loadu_si256((__m256i const *)(input + j));

                bn_vec = _mm256_srli_epi32( bn_vec, SHIFT);
                bn_vec = _mm256_and_si256( bn_vec, MASK_VEC);

                uint64_t bucket_numbers[4];
                _mm256_storeu_pd(( double *) (bucket_numbers), _mm256_castsi256_pd (bn_vec) ); 
  
                for(uint64_t m = 0; m < 4; m++){
                    bucket_num = (input[j+m].key >> SHIFT) & MASK;
                    int slot=buffers[bucket_numbers[m]].target & (TUPLES_PER_CACHELINE - 1);

                    if(slot == TUPLES_PER_CACHELINE - 1) {
                        uint64_t targetBkp=buffers[bucket_numbers[m]].target- (TUPLES_PER_CACHELINE-1);
                        buffers[bucket_numbers[m]].tuples[slot]=input[j+m]; //Can be vectorized.

                        //This loop has only 2 iterations.
                        for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                            _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                            targetBkp += STREAM_UNIT;
                        }
                        buffers[bucket_numbers[m]].target=targetBkp;
                    } 
                    else {
                        buffers[bucket_numbers[m]].tuples[slot] = input[j+m];
                        buffers[bucket_numbers[m]].target++;
                    }
                }
            }        
        }
        
        else{
            for(uint64_t j = 0; j < numTuples; ++j){
                bucket_num = (input[j].key >> SHIFT) & MASK;
                int slot=buffers[bucket_num].target & (TUPLES_PER_CACHELINE - 1); // Some bits are set?

                if(slot == TUPLES_PER_CACHELINE - 1) {

                    uint64_t targetBkp=buffers[bucket_num].target- (TUPLES_PER_CACHELINE-1);
                    //logger(INFO,"targetBkp: %d, bucket_num: %d, TUPLES_PER_CACHELINE: %d\n", targetBkp, bucket_num, TUPLES_PER_CACHELINE);
                    buffers[bucket_num].tuples[slot]=input[j];
                    for(uint32_t b = 0; b < TUPLES_PER_CACHELINE; b += STREAM_UNIT) {
                       _mm256_stream_si256(reinterpret_cast<__m256i*>(output + targetBkp), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[bucket_num].tuples + b))));
                        //memcpy(output + targetBkp, buffers[bucket_num].tuples + b, 32);
                        targetBkp += STREAM_UNIT;
                    }
                    buffers[bucket_num].target=targetBkp;

                } else {
                    buffers[bucket_num].tuples[slot] = input[j];
                    buffers[bucket_num].target++;
                }
            }
        }
        barrier->Arrive();

           //logger(INFO, "In radix_part(), impl_type = %d\n",impl_type);

        if(impl_type == CHTFSIMD || impl_type == CHTPA4){ //
            //logger(INFO, "In radix_part() CHTPa5, impl_type = %d\n",impl_type);
            __m256i TR_VEC = _mm256_set1_epi32(TUPLES_PER_CACHELINE - 1);
            for (int i = npart - 8; i >= 0; i-=8) {
                //Equivalent to:
                //uint32_t slot  = (uint32_t)buffers[i].target;
                __m256i slot_vec = _mm256_setr_epi32(buffers[i].target, buffers[i+1].target, 
                                    buffers[i+2].target, buffers[i+3].target,
                                    buffers[i+4].target,buffers[i+5].target,
                                    buffers[i+6].target, buffers[i+7].target) ;
                //Equivalent to:
                //uint32_t sz    = (slot) & (TUPLES_PER_CACHELINE - 1); // Can have values between 0 and 7!
                __m256i sz_vec = _mm256_and_si256(slot_vec, TR_VEC );
                
                //Equivalent to: //slot -= sz;
                slot_vec = _mm256_sub_epi32(slot_vec,sz_vec);
                uint32_t slot_numbers[8];
                _mm256_storeu_ps(( float *) (slot_numbers),  _mm256_castsi256_ps(slot_vec) );
                uint32_t sz_numbers[8];
                _mm256_storeu_ps(( float *) (sz_numbers),  _mm256_castsi256_ps(sz_vec) );

                for(int m = 0; m < 8; m++){
                    uint32_t startPos = (slot_numbers[m] < histogram[i+m]) ? ((uint32_t)histogram[i+m] - slot_numbers[m]) : 0;
                    //logger(INFO, "startPos: %d,\n",startPos);

                    if(sz_numbers[m] - startPos >= 4){
                        //logger(INFO, "I'm using _mm256_stream_si256");
                        _mm256_stream_si256(reinterpret_cast<__m256i*>(output+ slot_numbers[m]+startPos), _mm256_load_si256((reinterpret_cast<__m256i*>(buffers[i+m].tuples+startPos) ) ));
                        for(uint32_t j = 4; j < sz_numbers[m]; j++) {
                            //logger(INFO, "j = %d, ", j);
                            output[slot_numbers[m]+j]  = buffers[i+m].tuples[j];
                        }
                    }
                    else{
                        for(uint32_t j = startPos; j < sz_numbers[m]; j++) {
                            output[slot_numbers[m]+j]  = buffers[i+m].tuples[j];
                        }
                    }
                }
            }
        }    
        else{
            for (int i = npart - 1; i >= 0; i--) {
                uint32_t slot  = (uint32_t)buffers[i].target;
                uint32_t sz    = (slot) & (TUPLES_PER_CACHELINE - 1);
                slot          -= sz;
                uint32_t startPos = (slot < histogram[i]) ? ((uint32_t)histogram[i] - slot) : 0;
                for(uint32_t j = startPos; j < sz; j++) {
                    output[slot+j]  = buffers[i].tuples[j];
                }
            }   
        }
    }  //Using standard radix_partition with histogram.
    //The function creates npart=128 partitions.
    void partition(int threadID, relation_t chunkR)
    {
        const tuple_t * tupleR = chunkR.tuples;
        uint64_t * my_hist     = hist[threadID];
        uint64_t * my_dst      = dst[threadID];
        uint64_t sum           = 0;

        /* Compute the hash index for each key. .
        *  Count the number of key elements for each hash index.
        *  "Store" the results into a histogram aka 2D array.
        */

        //Vectorization: Very big hoping here.
        //logger(INFO,"here1\n");
        size_t i = 0;
        if(impl_type == CHTFSIMD || impl_type == CHTPA1 || impl_type == CHTOPT || impl_type == CHTOPT2){
        
            //Very tricky to vectorize this loop. Try later. In overall speedup not big, +- 3% .

            //Problem worse with table_sizes not of multiply of four.
            __m256i mask_vec = _mm256_set1_epi32(MASK);
            //logger(INFO,"here2\n");

            size_t i = 0;
            for (i = 0; i < chunkR.num_tuples-4; i+=4) { 
                //intkey_t hk = (hashKey(tupleR[i].key) >> SHIFT) & MASK; // Rewrite this Hash as SIMD instructions like SHIFT and AND.
                //my_hist[hk]++; //Difficult to vectorize. 
                //Especially this line is very tricky. How to store back to different places in memory? In AVX512 we have scatter.

                __m256i keyvals = _mm256_loadu_si256((__m256i const *)(tupleR + i));

                keyvals = _mm256_srli_epi32( keyvals, SHIFT);
                keyvals = _mm256_and_si256( keyvals, mask_vec);
                uint64_t hash_keys[4];
                
                _mm256_storeu_pd(( double *) (hash_keys), _mm256_castsi256_pd (keyvals) ); 
                //Difficult to vectorize. TO DO try later.
                for(uint64_t j=0; j<4; j++){
                    my_hist[hash_keys[j]]++;
                }
            }
            //logger(INFO,"num_tuples: %d, rest: %d\n", chunkR.num_tuples,rest);
            if(  i+4 != chunkR.num_tuples){
                for(uint64_t m=i; m< chunkR.num_tuples; m++){
                    int32_t hk = ((tupleR[m].key) >> SHIFT) & MASK; 
                    my_hist[hk]++;
                }
            }
            else{
                __m256i keyvals = _mm256_loadu_si256((__m256i const *)(tupleR + i));
                //Shift right "individual lanes" by SHIFT number.
                keyvals = _mm256_srli_epi32( keyvals, SHIFT);
                keyvals = _mm256_and_si256( keyvals, mask_vec);
                uint64_t hash_keys[4];
                
                _mm256_storeu_pd(( double *) (hash_keys), _mm256_castsi256_pd (keyvals) ); 
                //Difficult to vectorize.
                for(uint64_t j=0; j<4; j++){
                    my_hist[hash_keys[j]]++;
                }
            }
        }
        else{
            for (i = 0; i < chunkR.num_tuples; i++) { 
                intkey_t hk = (hashKey(tupleR[i].key) >> SHIFT) & MASK;
                my_hist[hk]++;
            }
        }

       //Not to vectorize, because of the dependency with the element before.
        for (int i = 0; i < npart; ++i) {
            sum += my_hist[i];
            my_hist[i] = sum;
        }

        barrier->Arrive();
        //my_dst[] is assigning for all the threads and partitions.

        for (int i = 0; i < threadID; ++i) {
            if (impl_type == CHTFSIMD || impl_type == CHTPA2){
                // Not much faster. In two same tuple size -> scalar faster by about 3%, if not equal simd by 3% faster.
                for (int j = 0; j < npart; j+=4) {
                    __m256i hist_int = _mm256_loadu_si256 ((__m256i const *) (my_dst+j) );
                    __m256i hist_int2 = _mm256_loadu_si256 ((__m256i const *) (hist[i]+j) );
                    hist_int = _mm256_add_epi64(hist_int,hist_int2);
                    _mm256_storeu_pd(( double *) (my_dst+j), _mm256_castsi256_pd (hist_int) );
                }  
            }
            else{
                for (int j = 0; j < npart; ++j) {
                    my_dst[j] += hist[i][j];
                }
            }
        }

        for (int i = threadID; i < nthreads; ++i) {
            if (impl_type == CHTFSIMD || impl_type == CHTPA2){
                for (int j = 1; j < npart-3; j+=4) {
                    __m256i hist_int = _mm256_loadu_si256 ((__m256i const *) (my_dst+j) );
                    __m256i hist_int2 = _mm256_loadu_si256 ((__m256i const *) (hist[i]+j-1) );
                    hist_int = _mm256_add_epi64(hist_int,hist_int2);
                    _mm256_storeu_pd(( double *) (my_dst+j), _mm256_castsi256_pd (hist_int) );
                }
                for (int j = npart-3; j < npart; j++) {
                    my_dst[j] += hist[i][j-1];
                }
            }
            else{
                for (int j = 1; j < npart; ++j) {
                    my_dst[j] += hist[i][j-1];
                }
            }
        }

        //for(int a = 0; a < npart; a++){
        //logger(INFO,"%ld ", my_dst[a]);
        //}
        radix_partition(chunkR.tuples, partBuffer, my_dst, chunkR.num_tuples);
    }

    void build(int threadID, int part)
    {
        (void) (threadID);
        tuple_t * tuples= partBuffer + dst[0][part];
        const uint64_t num_tuples = part == npart - 1 ? relR->num_tuples - dst[0][part] : dst[0][part + 1] - dst[0][part];

        if(num_tuples ==0){
            return;
        }
        if(impl_type == CHTB1){

            if(num_tuples < 4){
                for (uint64_t i = 0; i < num_tuples; i++){
                    ht.setBit(tuples[i].key);
                }
            }
            else{
                uint64_t i = 0;
                const __m256i BM_VEC = _mm256_set1_epi32(ht.get_bitMapSize() - 1); 
                //const __m256i PS_VEC = _mm256_set1_epi32(ht.get_partitionSize() - 1); 
                const __m256i BP_VEC = _mm256_set1_epi64x(ht.get_bitsPerBucket()-1); 
                const __m256i ZERO_VEC = _mm256_setzero_si256(); // Move to CHTJoin.hpp

                for ( i = 0; i < num_tuples-4; i+=4){
                    //Insert into the non-empty array.
                    ht.setBits_simd(tuples + i, BM_VEC, BP_VEC, ZERO_VEC);
                }
                //logger(INFO, "last vector: i = %d\n",i);
                if(i + 4 == num_tuples){
                    ht.setBits_simd(tuples + i, BM_VEC, BP_VEC, ZERO_VEC);
                }
                else{
                    for (uint64_t m = i; m < num_tuples; m++){
                        //Insert into the non-empty array.
                        ht.setBit(tuples[m].key);
                    }
                }
                }
            
        }
        else if(impl_type == CHTFSIMD || impl_type == CHTB1S || impl_type == CHTOPT || impl_type == CHTOPT2 || impl_type == CHTOPT3){

            if(num_tuples < 8){
                for (uint64_t i = 0; i < num_tuples; i++){
                    ht.setBit(tuples[i].key);
                }
            }
            else{
                uint64_t i = 0;
                const __m256i BM_VEC = _mm256_set1_epi32(ht.get_bitMapSize() - 1); 
                //const __m256i PS_VEC = _mm256_set1_epi32(ht.get_partitionSize() - 1); 
                const __m256i BP_VEC = _mm256_set1_epi32(ht.get_bitsPerBucket()-1); 
                const __m256i ZERO_VEC = _mm256_setzero_si256(); // Move to CHTJoin.hpp

                for ( i = 0; i < num_tuples-8; i+=8){
                    //Insert into the non-empty array.
                    //Try to vectorize here if possible.

                    ht.setBits_simd2(tuples + i, BM_VEC, BP_VEC, ZERO_VEC);
                    //ht.setBit(tuples[i].key);
                    //ht.setBit(tuples[i].key);
                }
                //logger(INFO, "last vector: i = %d\n",i);
                if(i + 8 == num_tuples){
                    //logger(INFO, "vector chosen:\n");
                    ht.setBits_simd2(tuples + i, BM_VEC, BP_VEC, ZERO_VEC);
                }
                else{
                    //logger(INFO, "scalar chosen:\n");
                    for (uint64_t m = i; m < num_tuples; m++){
                        //Insert into the non-empty array.
                        //Try to vectorize here if possible.
                        ht.setBit(tuples[m].key);
                    }
                }
                //logger(INFO, "end setbit: i = %d\n \n",i);
            }
        }
        else{
            //logger(INFO, "num_tuples = %d\n",num_tuples);
            for (uint64_t i = 0; i < num_tuples; i++){
                //Insert into the non-empty array.
                ht.setBit(tuples[i].key);
            }
        }

        ht.computePartPopCount(part, (uint32_t)dst[0][part]);
        
        if(impl_type == CHTFSIMD || impl_type == CHTB2  || impl_type == CHTOPT){
            
            if(num_tuples < 4){
                for (uint64_t i = 0; i < num_tuples; i++){
                    ht.setTuple(tuples[i]);
                }
                return;
            }
            uint64_t i = 0;
            const __m256i BM_VEC = _mm256_set1_epi32(ht.get_bitMapSize() - 1); //Move to CHTJoin.cpp, but bitMapSize is set private.
            const __m256i BP_VEC = _mm256_set1_epi64x(ht.get_bitsPerBucket()-1); //Move to CHTJoin.cpp
            
            for ( i = 0; i < num_tuples-4; i+=4){
                //Insert into the non-empty array.
                //Try to vectorize here if possible.
                ht.setTuples_simd(tuples+i, BM_VEC, BP_VEC);
            }

            //Scalar loop.
            if(i + 4 == num_tuples){
                ht.setTuples_simd(tuples+i, BM_VEC, BP_VEC);
            }
            else{
                for (uint64_t m = i; m < num_tuples; m++){
                //Insert into the non-empty array.
                //Try to vectorize here if possible.
                ht.setTuple(tuples[m]);
                }
            }   
        }
        else if(impl_type == CHTB2S){

            if(num_tuples < 8){
                for (uint64_t i = 0; i < num_tuples; i++){
                    ht.setTuple(tuples[i]);
                }
                return;
            }

            uint64_t i = 0;
            const __m256i BM_VEC = _mm256_set1_epi32(ht.get_bitMapSize() - 1); //Move to CHTJoin.cpp, but bitMapSize is set private.
            const __m256i BP_VEC = _mm256_set1_epi32(ht.get_bitsPerBucket()-1); //Move to CHTJoin.cpp
            //logger(INFO, "part no. : %d, num_tuples = %d\n",part, num_tuples);

            for ( i = 0; i < num_tuples-8; i+=8){
                //Insert into the non-empty array.
                //Try to vectorize here if possible.
                ht.setTuples_simd2(tuples+i, BM_VEC, BP_VEC);
                }
            //Scalar loop.
            if(i + 8 == num_tuples){
                ht.setTuples_simd2(tuples+i, BM_VEC, BP_VEC);
            }
            else{
                for (uint64_t m = i; m < num_tuples; m++){
                //Insert into the non-empty array.
                //Try to vectorize here if possible.
                ht.setTuple(tuples[m]);
                }
            }
        }
        else{
            for (uint64_t i = 0; i < num_tuples; ++i){
                //Insert into the non-empty array.
                //Try to vectorize here if possible.
                ht.setTuple(tuples[i]);
            }
        }
        
    }

    void probe(int threadID, relation_t chunkS)
    {
        uint64_t match = 0;
        uint64_t checksum = 0;
 		tuple_t *tupleS = chunkS.tuples;
		const size_t numS = chunkS.num_tuples;
		uint64_t batchStartUpperBound = numS - PROBE_BATCH_SIZE; //Type casted from const size_t.
        //logger(INFO,"Before probe: match: %d, checksum: %d, threadID: %d, batchStartUpperBound = %d \n", match, checksum, threadID, batchStartUpperBound);
        
        if( impl_type == CHTFSIMD || impl_type == CHTPR1  || impl_type == CHTOPT || impl_type == CHTOPT2 || impl_type == CHTOPT3){            
            //Try to vectorize here if possible.
            __m256i counter_vec = _mm256_setzero_si256();
            for (size_t i = 0; i <= batchStartUpperBound; i += PROBE_BATCH_SIZE)
            {
                counter_vec = ht.batch_probe_simd1(tupleS + i, match, checksum, counter_vec);
            }

            //Adding them all together.
            match += -( _mm256_extract_epi32(counter_vec, 0) + _mm256_extract_epi32(counter_vec, 2) 
            + _mm256_extract_epi32(counter_vec, 4) + _mm256_extract_epi32(counter_vec, 6)  );
        }
        else if(impl_type == CHTPR2){            
            __m256i counter_vec = _mm256_setzero_si256();
            const __m256i BM_VEC = _mm256_set1_epi32(ht.get_bitMapSize() - 1); //Move to CHTJoin.cpp, but bitMapSize is set private.
            const __m256i BP_VEC = _mm256_set1_epi64x(ht.get_bitsPerBucket()-1); //Move to CHTJoin.cpp

            for (size_t i = 0; i <= batchStartUpperBound; i += PROBE_BATCH_SIZE)
            {
                counter_vec = ht.batch_probe_simd2(tupleS + i, match, checksum, counter_vec, BM_VEC, BP_VEC);
            }

            //Adding them all together.
            match += -( _mm256_extract_epi32(counter_vec, 0) + _mm256_extract_epi32(counter_vec, 2) 
            + _mm256_extract_epi32(counter_vec, 4) + _mm256_extract_epi32(counter_vec, 6)  );

        }
        else{ //Original
            for (size_t i = 0; i <= batchStartUpperBound; i += PROBE_BATCH_SIZE)
            {
                ht.batch_probe(tupleS + i, match, checksum);
            }
        }

        //The rest should not be vectorized.
		const size_t leftOver = numS % PROBE_BATCH_SIZE;
		tupleS += numS - leftOver;
        //logger(INFO,"After batch probes: match: %d, checksum: %d, threadID: %d, leftover %d\n", match, checksum, threadID,leftOver);
        //logger(INFO,"match: %d, first_hit_probes = %d, second_hit_probes = %d \n \n",match, ht.first_hit_probes, ht.second_hit_probes);

        //logger(INFO,"leftOver= %d\n",leftOver);
        for (size_t i = 0; i < leftOver; ++i)
        {   
            tuple_t *foundTuple = ht.probe(tupleS[i].key);
            if (foundTuple)
            {
                match++;
                //checksum += foundTuple->payload + tupleS[i].payload;

            }
        }
        matches[threadID] = match;

        checksums[threadID] = checksum;
        //logger(INFO,"After probes: match: %d, checksum: %d, threadID: %d\n", match, checksum, threadID);
    }

public:

    CHTJoin_simd(int _nthreads, int _npart, relation_t *_relR, relation_t *_relS, tuple_t * _partBuffer, int _impl_type) : nthreads(_nthreads),
                                                     npart(_npart),impl_type( _impl_type), relR(_relR) , relS(_relS),partBuffer(_partBuffer),
                                                     partitions(new CHTPartitionQueue(npart)),
                                                     ht(relR->num_tuples, nthreads, npart),
                                                     MASK(npart-1),SHIFT(__builtin_ctz(Utils::nextPowerOfTwo(relR->num_tuples)) - __builtin_ctz(npart)),
                                                     matches(new uint64_t[nthreads]), checksums(new uint64_t[nthreads] ) 
    {
        barrier = new PThreadLockCVBarrier(nthreads);
        hist = (uint64_t **) malloc(nthreads * sizeof(uint64_t *));
        dst = (uint64_t **) malloc(nthreads * sizeof(uint64_t *));

        //logger(INFO,"In create CHTJoin_simd, in CHTJoin.hpp: impl_type=%d.\n",impl_type);
    }


    void join(int threadID) {
        //Init the CHT.
        ht.init(threadID);

        hist[threadID] = (uint64_t*) calloc(npart, sizeof(uint64_t));
        dst[threadID] = (uint64_t*) calloc(npart, sizeof(uint64_t));
        //logger(INFO,"In create join() beginning, in CHTJoin.hpp: impl_type=%d.\n",impl_type);
        barrier->Arrive();

        if (threadID == 0) {
            ocall_get_system_micros(&timers.start);
            ocall_startTimer(&timers.total);
            ocall_startTimer(&timers.timer1);
#ifdef PCM_COUNT
            ocall_set_system_counter_state("Start Partition");
#endif
        }

        // partition phase
        uint64_t numRPerThread = relR->num_tuples / nthreads;
        relation_t myChunkR;
        myChunkR.tuples = relR->tuples + numRPerThread * threadID;
        myChunkR.num_tuples = threadID == nthreads - 1 ? (uint32_t)(relR->num_tuples - numRPerThread * threadID) : (uint32_t)numRPerThread;

        partition(threadID, myChunkR);
        //logger(INFO,"After partition(), in CHTJoin.hpp: impl_type=%d.\n",impl_type);
        barrier->Arrive();

        if (threadID == 0) {
            ocall_get_system_micros(&timers.end);
            ocall_stopTimer(&timers.timer1);
            ocall_startTimer(&timers.timer2);
            partitionTime=timers.end-timers.start;
		}

        // build phase
        int partID;
        while ((partID = partitions->getPartition()) < npart)
        {
            build(threadID, partID);
        }
        //logger(INFO,"After build() , in CHTJoin.hpp: impl_type=%d.\n",impl_type);
        barrier->Arrive();

        if (threadID == 0)
        {
#ifdef PCM_COUNT
            ocall_get_system_counter_state("Partition", 0);
            ocall_set_system_counter_state("Join");
#endif
            ocall_stopTimer(&timers.timer2);
            ocall_startTimer(&timers.timer3);
        }
      //  std::cerr << "build finished" << std::endl;

        // probe phase
        uint64_t numSperThread = relS->num_tuples / nthreads;
        relation_t myChunkS;
        myChunkS.tuples = relS->tuples + numSperThread * threadID;
        myChunkS.num_tuples = threadID == nthreads - 1 ? (uint32_t)(relS->num_tuples - numSperThread * threadID) : (uint32_t)numSperThread;
        probe(threadID, myChunkS);
        //logger(INFO,"In join() ending, in CHTJoin.hpp: impl_type=%d.\n",impl_type);

        barrier->Arrive();

        if (threadID == 0) {
#ifdef PCM_COUNT
            ocall_get_system_counter_state("Join", 0);
#endif
	//		std::cerr << "probe finished" << std::endl;
	        ocall_stopTimer(&timers.timer3);
            ocall_get_system_micros(&timers.end);
            ocall_stopTimer(&timers.total);
            time_usec = timers.end - timers.start;

        }
    }

    join_result_t get_join_result()
    {
        uint64_t m = 0;
        uint64_t c = 0;
        for (int i = 0; i < nthreads; ++i) {
            m += matches[i];
            c += checksums[i];
        }
        join_result_t result;
        result.time_usec = time_usec;
        result.part_usec=partitionTime;
        result.join_usec=time_usec-partitionTime;
        result.matches = m;
        result.checksum = c;
        return result;
    }

    struct timers_t get_timers()
    {
        return timers;
    }

    ~CHTJoin_simd()
    {
		for (int i = 0; i < nthreads; ++i) {
			free(dst[i]);
			free(hist[i]);
		}
		free(dst);
		free(hist);
        delete[] checksums;
        delete[] matches;
        delete 	 barrier;
        delete   partitions;
    }
};
