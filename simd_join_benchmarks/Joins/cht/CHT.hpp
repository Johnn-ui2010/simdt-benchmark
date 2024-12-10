#include <stddef.h>
#include "Barrier.h"
#include "data-types.h"
#include <stdlib.h>
#include <string.h>
#include "Utils.h"
#include <assert.h>
#include <stdexcept>
#include <sstream>
#include <iostream>

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

#define PROBE_BATCH_SIZE 16

#ifndef CORES
#define CORES 8
#endif


/*
    "Implementation of CHT (Concise Hash Table)
    From " R. Barber, G. Lohman, I. Pandis, V. Raman, R. Sidle, G. Attaluri, 
    N. Chainani, S. Lightstone, and D. Sharpe. 2014. Memory-efficient hash joins. 
    Proc. VLDB Endow. 8, 4 (December 2014), 353–364. 
    https://doi.org/10.14778/2735496.2735499 "
*/
class CHT {
private:
    /*  See bitmap description.
        * bits  : actual bitmap
        * count : number of number of 1-bits UP to that word
    */
    struct popCount_t {
        uint32_t bits;
        uint32_t count;
    };

    const int bitsPerBucket=sizeof(popCount_t)*8/2;
    size_t tableSize; //number of table rows (usually of the table R).

    int nThreads;
    int nPartitions;
    size_t bitMapSize; //nextPowerOfTwo(tableSize)
    const size_t partitionSize;
    const int log2PartitionSize;

    /* 
    The bitmap here is an array of 64-bit words aka struct popCount. 
    Only 32 bits depict the actual bitmap (bits).
    The other 32 bits int rest population count represent the number of bits set to ONE in that word.
    */
    popCount_t *bitMap;

    //Array of the non-empty keys and payloads (100% fill factor)
    tuple_t *tupleArray;

    PThreadLockCVBarrier *barrier;

    bool setIfFree(intkey_t pos);
    tuple_t *findTuplePlace(intkey_t key);

    intkey_t nextInPartition(intkey_t pos);

    //SIMD functions
    bool setIfFree_simd(intkey_t pos, intkey_t bits, intkey_t hashBit, intkey_t div);
    intkey_t nextInPartition_simd(intkey_t pos);
    tuple_t *findTuplePlace_simd(intkey_t hash_new, intkey_t hash_new3); //hash_new2 is a helper for hash_new3 (See batch_probe_simd2)

public:

    CHT(size_t tableSize, int nThreads,int nPartitions);
    void computePopCount(int threadId);
    void computePartPopCount(int partID, uint32_t startCount);
    void init(uint64_t threadId);
    void setBit(intkey_t key);
    void setTuple(tuple_t tuple);
    tuple_t* probe(intkey_t key);
	void batch_probe(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum);

    intkey_t get_bitsPerBucket();
    intkey_t get_bitMapSize();
    intkey_t get_partitionSize();

    // TO DO move the desciption to ReadMe.
    void setBits_simd(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC,  __m256i ZERO_VEC);
    void setBits_simd2(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC, __m256i ZERO_VEC);
    void setTuples_simd(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC);
    void setTuples_simd2(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC);
    // SIMD probing1 (Partly SIMD batch_probe, no SIMD for findTuplePlace)
	__m256i batch_probe_simd1(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum, __m256i counter_vec);
    // SIMD probing2 (SIMD batch_probe, SIMD for findTuplePlace)
    __m256i batch_probe_simd2(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum, 
                __m256i counter_vec, __m256i BM_VEC, __m256i BP_VEC);

    //Special count for batch_probe (For debugging only):
    //uint32_t first_hit_probes = 0;
    //uint32_t second_hit_probes = 0;
};



#define _IDHASH_

#if defined(_IDHASH_)
    /** Identity Hashing */
    inline intkey_t hashKey(const intkey_t k) {
        return k;
    }
#elif defined(_FIBHASH_)
	/** Fibonacci Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return (k * 11400714819323198485ull) ;
	}
#elif defined(_CRCHASH_)
	/** CRC Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return _mm_crc32_u64(0, k) ;
	}
#else

/** MurmurHash64A */
inline intkey_t hashKey(intkey_t k) const {
    const intkey_t m = 0xc6a4a7935bd1e995;
    const int r = 47;
    intkey_t h = 0x8445d61a4e774912 ^(8 * m);
    k *= m;
    k ^= k >> r;
    k *= m;
    h ^= k;
    h *= m;
    h ^= h >> r;
    h *= m;
    h ^= h >> r;
    return h ;
}
#endif

//Assign values of class attributes.
CHT::CHT(size_t _tableSize,int _nThreads,int _nPartitions) : tableSize(_tableSize) ,nThreads(_nThreads),nPartitions(_nPartitions),bitMapSize(Utils::nextPowerOfTwo(tableSize)),
                                                          partitionSize((bitMapSize)/(unsigned long) _nPartitions),log2PartitionSize(__builtin_ctz((unsigned int)partitionSize)){

    assert(Utils::nextPowerOfTwo(_tableSize)%(unsigned long) bitsPerBucket==0);

    barrier=new PThreadLockCVBarrier(_nThreads);
}

void CHT::computePopCount(int threadId) {
    //First Version everything done by thread 0:
    if (threadId==0) {
        uint32_t count=0;
        for (size_t i = 0; i < bitMapSize/bitsPerBucket; ++i) {
            bitMap[i].count=count;
            count+=__builtin_popcount(bitMap[i].bits);
        }
    }

    barrier->Arrive();
}

//Not to vectorize, because of the dependency with the element before.
void CHT::computePartPopCount(int partID, uint32_t startCount) {
    uint32_t count = startCount;
    size_t startBuckets = (partitionSize / bitsPerBucket) * partID;
    size_t endBuckets = startBuckets + (partitionSize / bitsPerBucket);
    //logger(INFO, "startBuckets: %d, endBuckets: %d\n",startBuckets, endBuckets);

    //This loop has only 8 iterations.
    for (size_t i = startBuckets; i < endBuckets; ++i) {
        bitMap[i].count = count;
        count += __builtin_popcount(bitMap[i].bits);
    }
}

intkey_t CHT::get_bitsPerBucket(){
    return bitsPerBucket;
}

intkey_t CHT::get_bitMapSize(){
    return bitMapSize;
}  

intkey_t CHT::get_partitionSize(){
    return partitionSize;
}  

//Init the CHT. Allocate the spaces.
void CHT::init(uint64_t threadId) {

    uint64_t memChunkSize = 1024*2;
    uint64_t numChunksTable = (tableSize * sizeof(tuple_t)+memChunkSize-1) / memChunkSize;
    uint64_t numChunksBitmap = ((bitMapSize/bitsPerBucket)* sizeof(popCount_t)+memChunkSize-1) / memChunkSize;

		if (threadId == 0) {
            tupleArray = static_cast<tuple_t*>(malloc(tableSize * sizeof(tuple_t)+memChunkSize));
            bitMap = static_cast<popCount_t *>(malloc((bitMapSize/bitsPerBucket)* sizeof(popCount_t)+memChunkSize));
		}
		barrier->Arrive();

    //Table should be larger than bitmap therefore use same initorder as for table
		intkey_t * initOrderTable = Utils::generateShuffledNumbers((intkey_t)numChunksTable, (intkey_t)numChunksTable, 19650218ULL);
        intkey_t * initOrderBitmap = Utils::generateShuffledNumbers((intkey_t)numChunksBitmap, (intkey_t)numChunksBitmap, 19650218ULL);

		for (uint64_t i = threadId; i < numChunksTable; i += (unsigned long)nThreads) {
			memset(((char*)tupleArray) + (initOrderTable[i] * memChunkSize), 0, memChunkSize);
		}

    for (uint64_t i = threadId; i < numChunksBitmap; i += (unsigned long)nThreads) {
        memset(((char*)bitMap) + (initOrderBitmap[i] * memChunkSize), 0, memChunkSize);
    }
		delete [] initOrderTable;
    delete [] initOrderBitmap;

}

inline bool CHT::setIfFree(intkey_t pos) {
    //
    uint32_t bits=bitMap[pos>> /*(int)log2(bitsPerBucket)*/ 5].bits;
    uint32_t hashBit=(1<<(pos&(bitsPerBucket-1)));
    if ((bits & hashBit) ==0) {
        bitMap[pos/bitsPerBucket].bits|=hashBit;
        return true;
    }
    else
        return false;
}

//Set bit in the bitmap.
inline void CHT::setBit(intkey_t key) {
    intkey_t hash = (intkey_t) (hashKey(key) & (bitMapSize-1));
    //logger(INFO, "key = %d, " ,key);
    //logger(INFO, "correct hash: %d\n" ,hash);
    if(!setIfFree(hash)&&!setIfFree(nextInPartition(hash))) {
//		std::cerr << "setBit " << key << " failed" << std::endl;
//		logger(DBG, "setBit %d failed", key);
//        throw std::runtime_error("TODO insert into Overflow HT should not happen with dense keys!");
//        ocall_throw("TODO insert into Overflow HT should not happen with dense keys!");
	}

}

inline intkey_t CHT::nextInPartition(intkey_t pos) {
    return (intkey_t)((pos&(~(partitionSize-1))) | ((pos+1)&(partitionSize-1)));
}

//Return the tuple place in the non empty CHT array.
inline tuple_t* CHT::findTuplePlace(intkey_t key) {
    intkey_t hash = (intkey_t)(hashKey(key) & (bitMapSize-1));

    return tupleArray+bitMap[hash>>/* (int)log2(bitsPerBucket)*/ 5].count+
            __builtin_popcount(bitMap[hash>>/* (int)log2(bitsPerBucket)*/ 5].bits&
                    ~((~0)<<((hash&(bitsPerBucket-1))))) ;
}

// ---
//Insert into the non-empty array.popcount
inline void CHT::setTuple(tuple_t tuple) {
    tuple_t *toInsert=findTuplePlace(tuple.key);
    //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
    if (toInsert->key==0){ 
        *toInsert=tuple;
        //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
        //logger(INFO, "I'm here inside if statement. key=0\n");
    }
    else {
        //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
        //logger(INFO, "I'm here inside else statement. \n");
        toInsert=tupleArray+ nextInPartition((intkey_t)(toInsert-tupleArray));
        if (toInsert->key==0)
            *toInsert=tuple;

        else {

//            std::stringstream stringstream("TODO insert Tuple into Overflow HT should not happen with dense keys!");
//            stringstream<<" Key: "<<tuple.key<<" toInsert: "<<toInsert->key<<", "<<toInsert->payload;
//            throw std::runtime_error(stringstream.str());
//            logger(ERROR, "Key: %u to Insert: %u, %u", tuple.key, toInsert->key, toInsert->payload);
//            ocall_throw("TODO insert Tuple into Overflow HT should not happen with dense keys!");
        }
    }

}

inline tuple_t *CHT::probe(intkey_t key) {
    tuple_t *toReturn =findTuplePlace(key);
    if (toReturn->key==key){
        //logger(INFO, "Probe: First try correct, %d \n", toReturn->key);
        return toReturn;
    }
    else if ((++toReturn)->key==key){
        //logger(INFO, "Scalar probe: Second try correct, %d \n", toReturn->key);
        return toReturn;
    }
    else {
//		std::cout << "probe key " << key << " failed" << std::endl;
//		throw std::runtime_error("TODO lookup Tuple into Overflow HT should not happen with dense keys!");
		logger(DBG, "Probe key = %d", key);
		ocall_throw("Probe key failed");
	}
}

void CHT::batch_probe(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum) {
	tuple_t *tupleBatch[PROBE_BATCH_SIZE];
	for (int i = 0; i < PROBE_BATCH_SIZE; ++i) {
		tupleBatch[i] = findTuplePlace(probeTuples[i].key);
	}
	for (int i = 0; i < PROBE_BATCH_SIZE; ++i) {
		if (tupleBatch[i]->key == probeTuples[i].key) {
			matches++;
            //first_hit_probes++; //TO comment out
			//checksum += tupleBatch[i]->payload + probeTuples[i].payload;
            //logger(INFO,"In batch_probe probes, First try correct: checksum: %d\n", checksum);
        } else if ((++(tupleBatch[i]))->key == probeTuples[i].key) {
			matches++;
            //second_hit_probes++; //TO comment out
			//checksum += tupleBatch[i]->payload + probeTuples[i].payload;
            //logger(INFO,"In batch_probe probes, Second try correct: checksum: %d\n", checksum);
        } else {

//			std::cout << "batch probe key " << probeTuples[i].key << " " << tupleBatch[i]->key << " failed" << std::endl;
//		   	throw std::runtime_error("TODO lookup Tuple into Overflow HT should not happen with dense keys!");
//		   	logger(DBG, "Batch probe key %u %u failed", probeTuples[i].key, tupleBatch[i]->key);
//		   	ocall_throw("Batch probe key failed");
		}
	}
}

//SIMD helper functions.
//not used. It is actually integrated into setBits_simd().
inline bool CHT::setIfFree_simd(intkey_t pos, intkey_t bits, intkey_t hashBit, intkey_t div) {
    //
        if ((bits & hashBit) ==0) { //hash_numbers_new
            bitMap[div].bits|=hashBit;
        return true;
        }
        return false;
    
}

//Set bit in the bitmap.
//Very optimistic vectorization. We hope that !setIfFree((intkey_t) hash_numbers[j]) only is needed.
//nextPartition will be rarely triggered.
inline void CHT::setBits_simd( tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC, __m256i ZERO_VEC ) {

    //Equivalent to:  hash = (intkey_t) (hashKey(key) & (bitMapSize-1));
    __m256i hash_vec = _mm256_loadu_si256((__m256i const *)(tuples));
    hash_vec = _mm256_and_si256( hash_vec, BM_VEC);
    uint64_t hash_numbers[4];
    _mm256_storeu_pd(( double *) (hash_numbers), _mm256_castsi256_pd (hash_vec) ); 

    //Vectorize !setIfFree((intkey_t) hash_numbers[j])

    //Equivalent to // (pos >> 5)
    //Equivalent to //hash/bitsPerBucket  //bitsPerBucket = 32 = 2^(5), that means shifting right by 5.
    __m256i hash_vec_new = _mm256_srli_epi64( hash_vec, 5);
    uint64_t hash_numbers_new[4];
    _mm256_storeu_pd(( double *) (hash_numbers_new), _mm256_castsi256_pd (hash_vec_new) );

    //Equivalent to //(pos&(bitsPerBucket-1))
    __m256i keyvector = _mm256_and_si256(hash_vec, BP_VEC);
    uint64_t hash_numbers_new2[4];
    _mm256_storeu_pd(( double *) (hash_numbers_new2), _mm256_castsi256_pd (keyvector) ); 
    
    uint64_t bits[4];
    uint64_t hashBit[4];
    for(int j = 0; j < 4; j++){
        bits[j]=bitMap[hash_numbers_new[j]].bits; 
        hashBit[j]=(1<<  hash_numbers_new2[j] );
    }

    //logger(INFO, "\n");
    //Equivalent to //(bits & hashBit) == 0
    __m256i bits_vec = _mm256_loadu_si256((__m256i const *)(bits));
    __m256i hashBit_vec = _mm256_loadu_si256((__m256i const *)(hashBit));
    
    bits_vec = _mm256_and_si256(bits_vec, hashBit_vec);
    bits_vec = _mm256_cmpeq_epi64(bits_vec,ZERO_VEC);

    uint64_t hash_numbers_new3[4];
    _mm256_storeu_pd(( double *) (hash_numbers_new3), _mm256_castsi256_pd (bits_vec) ); 

    for(int j = 0; j < 4; j++){
        bool set1 = false;
        if (hash_numbers_new3[j]) { //hash_numbers_new
            bitMap[hash_numbers_new[j]].bits|=hashBit[j];
            //logger(INFO,"Yes, correct. hash_numbers_new4[j] = %d \n", hash_numbers_new4[j]);
        }
        else {
            set1 = true;
        }
        if( (set1) &&!setIfFree( nextInPartition((intkey_t) hash_numbers[j]))) {

//		std::cerr << "setBit " << key << " failed" << std::endl;
//		logger(DBG, "setBit %d failed", key);
//        throw std::runtime_error("TODO insert into Overflow HT should not happen with dense keys!");
//        ocall_throw("TODO insert into Overflow HT should not happen with dense keys!");
	    }
    }
        //logger(INFO, "\n");
    
}

inline void CHT::setBits_simd2( tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC, __m256i ZERO_VEC ) {

    //Equivalent to:  hash = (intkey_t) (hashKey(key) & (bitMapSize-1));
    __m256i hash_vec = _mm256_setr_epi32(tuples[0].key,tuples[1].key,tuples[2].key,tuples[3].key,
                                        tuples[4].key, tuples[5].key, tuples[6].key, tuples[7].key );
    
    //Alternative, a litte slower.
    //__m256i vindex = _mm256_set_epi32(0,1,2,3,4,5,6,7);
        //int * tuples_int = (int *) tuples;
    //__m256i hash_vec = _mm256_i32gather_epi32(( int *)(tuples), vindex, 8 ); //tuple + i*2, loading every i * 2 tuples. -> It works.

    hash_vec = _mm256_and_si256( hash_vec, BM_VEC);
    uint32_t hash_numbers[8];
    _mm256_storeu_ps(( float *) (hash_numbers), _mm256_castsi256_ps (hash_vec) ); 

    //Vectorize !setIfFree((intkey_t) hash_numbers[j])

    //Equivalent to // (pos >> 5)
    //Equivalent to //hash/bitsPerBucket  //bitsPerBucket = 32 = 2^(5), that means shifting right by 5.
    __m256i hash_vec_new = _mm256_srli_epi32( hash_vec, 5);
    uint32_t hash_numbers_new[8];
    _mm256_storeu_ps(( float *) (hash_numbers_new), _mm256_castsi256_ps (hash_vec_new) );

    //Equivalent to //(pos&(bitsPerBucket-1))
    __m256i keyvector = _mm256_and_si256(hash_vec, BP_VEC);
    uint32_t hash_numbers_new2[8];
    _mm256_storeu_ps(( float *) (hash_numbers_new2), _mm256_castsi256_ps (keyvector) ); 
    
    uint32_t bits[8];
    uint32_t hashBit[8];
    for(int j = 0; j < 8; j++){
        bits[j]=bitMap[hash_numbers_new[j]].bits; 
        hashBit[j]=(1<<  hash_numbers_new2[j] );
        
    }

    //logger(INFO, "\n");
    //Equivalent to //(bits & hashBit) == 0
    __m256i bits_vec = _mm256_loadu_si256((__m256i const *)(bits));
    __m256i hashBit_vec = _mm256_loadu_si256((__m256i const *)(hashBit));
    
    bits_vec = _mm256_and_si256(bits_vec, hashBit_vec);
    bits_vec = _mm256_cmpeq_epi32(bits_vec,ZERO_VEC);

    uint32_t hash_numbers_new3[8];
    _mm256_storeu_ps(( float *) (hash_numbers_new3), _mm256_castsi256_ps (bits_vec) ); 

    for(int j = 0; j < 8; j++){
        bool set1 = false;
        if (hash_numbers_new3[j]) { //hash_numbers_new
            bitMap[hash_numbers_new[j]].bits|=hashBit[j];
            //logger(INFO,"Yes, correct. hash_numbers_new4[j] = %d \n", hash_numbers_new4[j]);
        }
        else {
            set1 = true;
        }
        if( (set1) &&!setIfFree( nextInPartition((intkey_t) hash_numbers[j]))) {

//		std::cerr << "setBit " << key << " failed" << std::endl;
//		logger(DBG, "setBit %d failed", key);
//        throw std::runtime_error("TODO insert into Overflow HT should not happen with dense keys!");
//        ocall_throw("TODO insert into Overflow HT should not happen with dense keys!");
	    }
    }
        //logger(INFO, "\n");
    
}

//Not used.
inline intkey_t CHT::nextInPartition_simd(intkey_t pos) {
    return (intkey_t)((pos&(~(partitionSize-1))) | ((pos+1)&(partitionSize-1)));
}

void CHT::setTuples_simd(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC){
    //Find tuple place part.
    //Similar to batch_probe_simd2().

    //Equivalent to // intkey_t hash = (intkey_t)(hashKey(probeTuples[i].key) & (bitMapSize-1));
    __m256i hash_vec = _mm256_loadu_si256((__m256i const *)(tuples));
    hash_vec = _mm256_and_si256( hash_vec, BM_VEC);
    //uint64_t hash_numbers[4];
    //_mm256_storeu_pd(( double *) (hash_numbers), _mm256_castsi256_pd (hash_vec) ); 

    //Equivalent to // hash>>/* (int)log2(bitsPerBucket)*/ 5
    __m256i hash_vec_new = _mm256_srli_epi64( hash_vec, 5);
    uint64_t hash_numbers_new[4];
    _mm256_storeu_pd(( double *) (hash_numbers_new), _mm256_castsi256_pd (hash_vec_new) );

    //Equivalent to //(hash&(bitsPerBucket-1))
    __m256i keyvector = _mm256_and_si256(hash_vec, BP_VEC);
    uint64_t hash_numbers_new2[4];
    _mm256_storeu_pd(( double *) (hash_numbers_new2), _mm256_castsi256_pd (keyvector) ); 

    //Equivalent to // bitMap[hash_new].bits&  // &~ -> Try to vectorize later.
        //~((~0)<<(hash_new2))
    uint64_t hash_numbers_new3[4];
    for(int i = 0; i < 4; i++){
        //logger(INFO, "tuples[j].key = %d, hash_numbers: %d \n" ,tuples[i].key, hash_numbers_new[i]);

        hash_numbers_new3[i] = (~0) << hash_numbers_new2[i]; 
    }
    
    __m256i A = _mm256_setr_epi64x(bitMap[hash_numbers_new[0]].bits, bitMap[hash_numbers_new[1]].bits, 
    bitMap[hash_numbers_new[2]].bits, bitMap[hash_numbers_new[3]].bits);
    __m256i B = _mm256_loadu_si256((__m256i const *)(hash_numbers_new3 + 0));
    A = _mm256_andnot_si256(B,A); 
    _mm256_storeu_pd(( double *) (hash_numbers_new3), _mm256_castsi256_pd (A) );

    tuple_t *tupleBatch[4]; // Array of ToInserts.

    //logger(INFO, "Problem after tupleBatch in build\n");
    //Problem: Somewhere in tupleBatch[j]=&tuples[j]; or for ( i = 0; i < num_tuples-4; i+=4)
    //Insert part
    for(int j = 0; j < 4; j++){
    //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));

        tupleBatch[j] = findTuplePlace_simd( (uint32_t) hash_numbers_new[j], (uint32_t) hash_numbers_new3[j]);
        if (tupleBatch[j]->key==0){ //Probably, here we will have some problems with with input key = 0.
            *tupleBatch[j]=tuples[j];
            //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(tupleBatch[j]-tupleArray)));
        }
        else {
            //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
            //logger(INFO, "I'm here inside else statement. \n");
            tupleBatch[j]=tupleArray+ nextInPartition((intkey_t)(tupleBatch[j]-tupleArray));
            if (tupleBatch[j]->key==0)
                *tupleBatch[j]=tuples[j];

            else {
    //            std::stringstream stringstream("TODO insert Tuple into Overflow HT should not happen with dense keys!");
    //            stringstream<<" Key: "<<tuple.key<<" toInsert: "<<toInsert->key<<", "<<toInsert->payload;
    //            throw std::runtime_error(stringstream.str());
    //            logger(ERROR, "Key: %u to Insert: %u, %u", tuple.key, toInsert->key, toInsert->payload);
    //            ocall_throw("TODO insert Tuple into Overflow HT should not happen with dense keys!");
            }
        }
        
    }
    //logger(INFO, "Problem after Insert in build\n");
}

void CHT::setTuples_simd2(tuple_t *tuples, __m256i BM_VEC, __m256i BP_VEC){
    //Find tuple place part.
    //Similar to batch_probe_simd2().

    //Equivalent to // intkey_t hash = (intkey_t)(hashKey(probeTuples[i].key) & (bitMapSize-1));
    __m256i hash_vec = _mm256_setr_epi32(tuples[0].key,tuples[1].key,tuples[2].key,tuples[3].key,
                                    tuples[4].key, tuples[5].key, tuples[6].key, tuples[7].key );
    
    hash_vec = _mm256_and_si256( hash_vec, BM_VEC);
    //uint32_t hash_numbers[8];
    //_mm256_storeu_ps(( float *) (hash_numbers), _mm256_castsi256_ps (hash_vec) ); 

    //Equivalent to // hash>>/* (int)log2(bitsPerBucket)*/ 5
    __m256i hash_vec_new = _mm256_srli_epi32( hash_vec, 5);
    uint32_t hash_numbers_new[8];
    _mm256_storeu_ps(( float *) (hash_numbers_new), _mm256_castsi256_ps (hash_vec_new) );

    //Equivalent to //(hash&(bitsPerBucket-1))
    __m256i keyvector = _mm256_and_si256(hash_vec, BP_VEC);
    uint32_t hash_numbers_new2[8];
    _mm256_storeu_ps(( float *) (hash_numbers_new2), _mm256_castsi256_ps (keyvector) ); 

    //Equivalent to // bitMap[hash_new].bits&  // &~ -> Try to vectorize later.
        //~((~0)<<(hash_new2))
    uint32_t hash_numbers_new3[8];
    for(int i = 0; i < 8; i++){
        //logger(INFO, "tuples[j].key = %d, hash_numbers: %d \n" ,tuples[i].key, hash_numbers_new[i]);

        hash_numbers_new3[i] = (~0) << hash_numbers_new2[i]; 
    }
    
    __m256i A = _mm256_setr_epi32(bitMap[hash_numbers_new[0]].bits, bitMap[hash_numbers_new[1]].bits, 
                            bitMap[hash_numbers_new[2]].bits, bitMap[hash_numbers_new[3]].bits,
                            bitMap[hash_numbers_new[4]].bits, bitMap[hash_numbers_new[5]].bits,
                            bitMap[hash_numbers_new[6]].bits, bitMap[hash_numbers_new[7]].bits);

    __m256i B = _mm256_loadu_si256((__m256i const *)(hash_numbers_new3 + 0));
    A = _mm256_andnot_si256(B,A); 
    _mm256_storeu_ps(( float *) (hash_numbers_new3), _mm256_castsi256_ps (A) );

    tuple_t *tupleBatch[8]; // Array of ToInserts.

    //logger(INFO, "Problem after tupleBatch in build\n");
    //Insert part
    for(int j = 0; j < 8; j++){
    //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
        tupleBatch[j] = findTuplePlace_simd( (uint32_t) hash_numbers_new[j], (uint32_t) hash_numbers_new3[j]);
        tuple_t *tuple_tmp = findTuplePlace( (uint32_t) tuples[j].key);
        if(tupleBatch[j]->key != tuple_tmp->key){
            logger(INFO, "In setTuple: tupleBatch[j] = %d,tuple->key = %d\n", tupleBatch[j]->key,tuple_tmp->key );
        }
        
        if (tupleBatch[j]->key==0){ //Probably, here we will have some problems with with input key = 0.
            *tupleBatch[j]=tuples[j];
            //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(tupleBatch[j]-tupleArray)));
        }
        else {
            //logger(INFO, "In setTuple: nextInPart: %d\n",nextInPartition((intkey_t)(toInsert-tupleArray)));
            //logger(INFO, "I'm here inside else statement. \n");
            tupleBatch[j]=tupleArray+ nextInPartition((intkey_t)(tupleBatch[j]-tupleArray));
            if (tupleBatch[j]->key==0)
                *tupleBatch[j]=tuples[j];

            else {
    //            std::stringstream stringstream("TODO insert Tuple into Overflow HT should not happen with dense keys!");
    //            stringstream<<" Key: "<<tuple.key<<" toInsert: "<<toInsert->key<<", "<<toInsert->payload;
    //            throw std::runtime_error(stringstream.str());
    //            logger(ERROR, "Key: %u to Insert: %u, %u", tuple.key, toInsert->key, toInsert->payload);
    //            ocall_throw("TODO insert Tuple into Overflow HT should not happen with dense keys!");
            }
        }
        
    }
    //logger(INFO, "Problem after Insert in build\n");
}


__m256i CHT::batch_probe_simd1(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum, __m256i counter_vec){
	tuple_t *tupleBatch[PROBE_BATCH_SIZE]; //Not stored contiguously. tupleBatch stores pointer of tuple arrays, the tuples.
    //tuple_t tupleBatch[PROBE_BATCH_SIZE];
	for (int i = 0; i < PROBE_BATCH_SIZE; ++i) {
		tupleBatch[i] = findTuplePlace(probeTuples[i].key);
        //logger(INFO, "keys: tupleBatch[%d] = %d, Probetuples[%d] = %d \n", i, tupleBatch[i]->key,i, probeTuples[i].key);
	}

    // Comparing the vector.
    
    //Only the first case covered (Very optimistic): if (tupleBatch[i]->key == probeTuples[i].key) 
    for (int i = 0; i < PROBE_BATCH_SIZE; i += 4) {
        
        __m256i keyvector = _mm256_setr_epi32(tupleBatch[i]->key,-1, tupleBatch[i+1]->key,-1, tupleBatch[i+2]->key,-1, tupleBatch[i+3]->key,-1 );
        //__m256i keyvector = _mm256_loadu_si256(reinterpret_cast<__m256i*>(tupleBatch + i));
        __m256i probe_vector = _mm256_loadu_si256((__m256i const *)(probeTuples + i));
        
        keyvector = _mm256_cmpeq_epi32(keyvector, probe_vector);
        counter_vec = _mm256_add_epi32(keyvector, counter_vec);
    }

    return counter_vec;
    
    //Case not covered: else if ((++(tupleBatch[i]))->key == probeTuples[i].key) 

}

//SIMD return the tuple place in the non empty CHT array.
inline tuple_t* CHT::findTuplePlace_simd(intkey_t hash_new, intkey_t hash_new3) {
    
    //In normal version equivalent to //
    //return tupleArray+bitMap[hash>>/* (int)log2(bitsPerBucket)*/ 5].count+
    //        __builtin_popcount(bitMap[hash>>/* (int)log2(bitsPerBucket)*/ 5].bits&
    //                ~((~0)<<((hash&(bitsPerBucket-1))))) ;
    //hash_new2 =hash&(bitsPerBucket-1);

    return tupleArray+bitMap[hash_new].count+
            __builtin_popcount(hash_new3);
}

__m256i CHT::batch_probe_simd2(tuple_t *probeTuples, uint64_t &matches, uint64_t &checksum, __m256i counter_vec, __m256i BM_VEC, __m256i BP_VEC){
	tuple_t *tupleBatch[PROBE_BATCH_SIZE];
    
    for (int i = 0; i < PROBE_BATCH_SIZE; i+=4) {

        //Equivalent to // intkey_t hash = (intkey_t)(hashKey(probeTuples[i].key) & (bitMapSize-1));
        __m256i hash_vec = _mm256_loadu_si256((__m256i const *)(probeTuples + i));
        hash_vec = _mm256_and_si256( hash_vec, BM_VEC);
        //uint64_t hash_numbers[4];
        //_mm256_storeu_pd(( double *) (hash_numbers), _mm256_castsi256_pd (hash_vec) ); 

        //Equivalent to // hash>>/* (int)log2(bitsPerBucket)*/ 5
        __m256i hash_vec_new = _mm256_srli_epi64( hash_vec, 5);
        uint64_t hash_numbers_new[4];
        _mm256_storeu_pd(( double *) (hash_numbers_new), _mm256_castsi256_pd (hash_vec_new) );

        //Equivalent to //(hash&(bitsPerBucket-1))
        __m256i keyvector = _mm256_and_si256(hash_vec, BP_VEC);
        uint64_t hash_numbers_new2[4];
        _mm256_storeu_pd(( double *) (hash_numbers_new2), _mm256_castsi256_pd (keyvector) ); 

        //Equivalent to ((~0)<<(hash_new2))
        uint64_t hash_numbers_new3[4];
        for(int i = 0; i < 4; i++){
            hash_numbers_new3[i] = (~0) << hash_numbers_new2[i]; 
        }
        //Equivalent to // bitMap[hash_new].bits&  // &~ -> Try to vectorize later.
            //~((~0)<<(hash_new2))
        __m256i A = _mm256_setr_epi64x(bitMap[hash_numbers_new[0]].bits, bitMap[hash_numbers_new[1]].bits, 
                    bitMap[hash_numbers_new[2]].bits, bitMap[hash_numbers_new[3]].bits);
        __m256i B = _mm256_loadu_si256((__m256i const *)(hash_numbers_new3 + 0));
        A = _mm256_andnot_si256(B,A); 
        _mm256_storeu_pd(( double *) (hash_numbers_new3), _mm256_castsi256_pd (A) );

        for(int j = 0; j < 4; j++){
            tupleBatch[i+j] = findTuplePlace_simd( (uint32_t) hash_numbers_new[j], (uint32_t) hash_numbers_new3[j]); 
        }
		
        //logger(INFO, "keys: tupleBatch[%d] = %d, Probetuples[%d] = %d \n", i, tupleBatch[i]->key,i, probeTuples[i].key);
	}

    // Comparing the vector.
    
    //Only the first case covered (Very optimistic): if (tupleBatch[i]->key == probeTuples[i].key) 
    for (int i = 0; i < PROBE_BATCH_SIZE; i += 4) {
        
        __m256i keyvector = _mm256_setr_epi32(tupleBatch[i]->key,-1, tupleBatch[i+1]->key,-1, tupleBatch[i+2]->key,-1, tupleBatch[i+3]->key,-1 );
        __m256i probe_vector = _mm256_loadu_si256((__m256i const *)(probeTuples + i));
        
        keyvector = _mm256_cmpeq_epi32(keyvector, probe_vector);
        counter_vec = _mm256_add_epi32(keyvector, counter_vec);
    }

    return counter_vec;
    
    //Case not covered: else if ((++(tupleBatch[i]))->key == probeTuples[i].key) 

}

