

#include <string.h>
#include "bs.h"
// #include <arm_neon.h>

#include <stdio.h>
#include <inttypes.h> // for PRId64 macro

#define bs2le(x) (x)
#define bs2be(x) (x)
static inline __m512i linear_transform(__m512i x, __m512i map) ;
void bs_transpose(word_t * blocks, word_t width_to_adjacent_block)
{
    word_t transpose[BLOCK_SIZE];
    memset(transpose, 0, sizeof(transpose));
    bs_transpose_dst(transpose,blocks, width_to_adjacent_block);

    int sizeof_transpose = sizeof(transpose);
    memmove(blocks,transpose,sizeof(transpose));

}


__m512i TOWER_TO_AES_MAP;
__m512i AES_TO_TOWER_MAP;
__m512i const_16_512b;

///// TEST LOADS/////
__m512i TEMP1;
__m512i TEMP2;
//// END TEST LOADS////
// Define the constants in the form of vectors
void initialize_maps() {
    // Initialize the constants at runtime
    TOWER_TO_AES_MAP = _mm512_set1_epi64(0x31506aea964e983e);
    AES_TO_TOWER_MAP = _mm512_set1_epi64(0xd1e8863e72a2700c);
    
    const_16_512b = linear_transform(_mm512_set1_epi8(16), TOWER_TO_AES_MAP);

    // remove 2 lines below
    TEMP1 = _mm512_set1_epi64(0x31506aea964e983e);
    TEMP2 = _mm512_set1_epi64(0xd1e8863e72a2700c);


}

// Helper functions

// Updated linear transformation function based on the Rust definition
static inline __m512i linear_transform(__m512i x, __m512i map) {
    // __m512i map_vector = _mm512_set1_epi64(map);
    return _mm512_gf2p8affine_epi64_epi8(x, map, 0); // The constant 0 is the 'b' value, which remains zero.
}

// GF(2^8) multiplication using the GFNI instruction
static inline __m512i gf2p8mul_epi8(__m512i lhs, __m512i rhs) {
    return _mm512_gf2p8mul_epi8(lhs, rhs);
}

// Function that simulates the multiplication in GF(2^8) with optional transformation maps
__m512i gfni_mul(__m512i lhs, __m512i rhs) {
    
    // Apply the TO_AES_MAP transformation if it's not the identity
    // if (TOWER_TO_AES_MAP != 0) {
      lhs = linear_transform(lhs, TOWER_TO_AES_MAP);
       rhs = linear_transform(rhs, TOWER_TO_AES_MAP);
    // }


    // Perform the GF(2^8) multiplication using the GFNI instruction
    __m512i prod_gfni = gf2p8mul_epi8(lhs, rhs);

    // Apply the FROM_AES_MAP transformation if it's not the identity
    // if (AES_TO_TOWER_MAP != 0) {
    prod_gfni = linear_transform(prod_gfni, AES_TO_TOWER_MAP);
    // }

    return prod_gfni;
}

__m512i gfni_mul_const_16(__m512i rhs) {
    // Apply the TO_AES_MAP transformation if it's not the identity
    // if (TOWER_TO_AES_MAP != 0) {
       __m512i rhs_transformed = linear_transform(rhs, TOWER_TO_AES_MAP);
    // }


    // Perform the GF(2^8) multiplication using the GFNI instruction
    __m512i prod_gfni = gf2p8mul_epi8(const_16_512b, rhs_transformed);

    // Apply the FROM_AES_MAP transformation if it's not the identity
    // if (AES_TO_TOWER_MAP != 0) {
      __m512i  prod_tower = linear_transform(prod_gfni, AES_TO_TOWER_MAP);
    // }

    return prod_tower;
}


inline void multiply_512b_using_gfni(__m512i *lhs, __m512i *rhs,__m512i* result ){
    __m512i lhs_512 = _mm512_load_si512(lhs);
    __m512i rhs_512 = _mm512_load_si512(rhs);
    *result = gfni_mul(lhs_512, rhs_512);
}

inline void multiply_constant_512b_using_gfni(__m512i *rhs, __m512i* result ){
    __m512i rhs_512 = _mm512_load_si512(rhs);
    *result = gfni_mul_const_16( rhs_512);
}

#define NUM_INPUTS 64        
#define BYTES_IN_128BIT 16   
#define BYTES_IN_512BIT 64
#define SLICED_OUTPUTS 16


//////////////////////////////////////////////////////
// Byte Slicing: Convert 64 x 128-bit inputs to 16 rows of 512 bits
//////////////////////////////////////////////////////
inline void byte_slice(uint128_t input[NUM_INPUTS], uint64_t output[SLICED_OUTPUTS]) {
    // Cast the input to uint8_t* for easier access to bytes
    uint8_t* input_bytes = (uint8_t*) input;
    uint8_t* output_bytes = (uint8_t*) output;

    // Loop over each byte index (0 to 15 for each 128-bit number)
    for (int byte_index = 0; byte_index < BYTES_IN_128BIT; byte_index++) {
        // Loop over each of the original 128-bit numbers (8 total inputs)
        for (int i = 0; i < NUM_INPUTS; i++) {
            // Place each byte from the input array into the corresponding 64-bit output row
            // the print the index being accesses
            // int index_to = byte_index * NUM_INPUTS + i;
            // int index_from = i * BYTES_IN_128BIT + byte_index;
            output_bytes[byte_index * NUM_INPUTS + i] = input_bytes[i * BYTES_IN_128BIT + byte_index];
        }
    }

}

//////////////////////////////////////////////////////
// Un-byte Slicing:
//////////////////////////////////////////////////////
inline void un_byte_slice(uint64_t input[SLICED_OUTPUTS], uint128_t output[NUM_INPUTS]) {
    // Cast the input to uint8_t* for easier access to bytes
    uint8_t* input_bytes = (uint8_t*) input;
    uint8_t* output_bytes = (uint8_t*) output;

    // Loop over each byte index (0 to 15 for each 128-bit number)
    for (int byte_index = 0; byte_index < BYTES_IN_128BIT; byte_index++) {
        // Loop over each of the 8 128-bit numbers
        for (int i = 0; i < NUM_INPUTS; i++) {
            // Reconstruct the original bytes into the 128-bit numbers
            output_bytes[i * BYTES_IN_128BIT + byte_index] = input_bytes[byte_index * NUM_INPUTS + i];
        }
    }
}
//////////////////////////////////////////////////////
// Byte Slicing: Convert 64 x 128-bit inputs to 16 rows of 512 bits
//////////////////////////////////////////////////////
void byte_slice_avx(__m512i *input, __m512i *output) {
    // Cast the input and output arrays to uint8_t* for byte-wise access
    uint8_t* input_bytes = (uint8_t*)input;
    uint8_t* output_bytes = (uint8_t*)output;

    // Loop over each byte index (0 to 15)
    for (int byte_index = 0; byte_index < BYTES_IN_128BIT; byte_index++) {
        // Collect the bytes from each input vector at the current byte position
        uint8_t temp_bytes[64];

        for (int i = 0; i < NUM_INPUTS; i++) {
            temp_bytes[i] = input_bytes[i * BYTES_IN_128BIT + byte_index];
        }

        // Load the 64 bytes into a 512-bit vector
        __m512i vec = _mm512_loadu_si512((__m512i*)temp_bytes);

        // Store the vector into the output array
        _mm512_store_si512((__m512i*)&output_bytes[byte_index * NUM_INPUTS], vec);
    }

}

inline __m512i wrapper_mm512_xor_si512(__m512i input_1, __m512i input_2){
    // __m512i random;
    // return random;  
    __m512i lhs_512 = _mm512_load_si512(&input_1);
    __m512i rhs_512 = _mm512_load_si512(&input_2);
    return _mm512_xor_si512(lhs_512, rhs_512);
}

void byte_slice_avx_2_mul(__m512i *input_1, __m512i *input_2, __m512i *slice_output_1, __m512i *slice_output_2, __m512i *mul_level_0, __m512i *xor_level_0)  {
    // Cast the input and output arrays to uint8_t* for byte-wise access
    uint8_t* input_bytes_1 = (uint8_t*)input_1;
    uint8_t* input_bytes_2 = (uint8_t*)input_2;
    uint8_t* output_bytes_1 = (uint8_t*)slice_output_1;
    uint8_t* output_bytes_2 = (uint8_t*)slice_output_2;
    uint8_t* mul_bytes = (uint8_t*)mul_level_0;
    uint8_t* xor_bytes = (uint8_t*)xor_level_0;


    // Loop over each byte index (0 to 15)
    for (int byte_index = 0; byte_index < BYTES_IN_128BIT; byte_index++) {
        // Collect the bytes from each input vector at the current byte position
        uint8_t temp_bytes_1[64];
        uint8_t temp_bytes_2[64];
        // __m512i prev_mul_vec = _mm512_set1_epi8(0);

        for (int i = 0; i < NUM_INPUTS; i++) {
            temp_bytes_1[i] = input_bytes_1[i * BYTES_IN_128BIT + byte_index];
            temp_bytes_2[i] = input_bytes_2[i * BYTES_IN_128BIT + byte_index];
        }

        // Load the 64 bytes into a 512-bit vector
        __m512i vec_1 = _mm512_loadu_si512((__m512i*)temp_bytes_1);
        __m512i vec_2 = _mm512_loadu_si512((__m512i*)temp_bytes_2);
        __m512i mul_vec = gfni_mul(vec_1, vec_2);

        // Store the vector into the output array
        _mm512_store_si512((__m512i*)&output_bytes_1[byte_index * NUM_INPUTS], vec_1);
        _mm512_store_si512((__m512i*)&output_bytes_2[byte_index * NUM_INPUTS], vec_2);
        _mm512_store_si512((__m512i*)&mul_bytes[byte_index * NUM_INPUTS], mul_vec);

    }

}


//////////////////////////////////////////////////////
// Un-byte Slicing: Convert 16 rows of 512 bits back to 64 x 128-bit inputs
//////////////////////////////////////////////////////
inline void un_byte_slice_avx(uint64_t input[SLICED_OUTPUTS], uint128_t output[NUM_INPUTS]) {
    uint8_t (*input_matrix)[16] = (uint8_t (*)[16])input;
    uint8_t (*output_matrix)[64] = (uint8_t (*)[64])output;

    for (int i = 0; i < 16; i++) {
        // For each byte position i (from 0 to 15)
        for (int j = 0; j < 64; j += 64 / 4) {
            // Process 16 inputs at a time
            __m512i row = _mm512_set_epi8(
                input_matrix[j + 15][i], input_matrix[j + 14][i], input_matrix[j + 13][i], input_matrix[j + 12][i],
                input_matrix[j + 11][i], input_matrix[j + 10][i], input_matrix[j + 9][i], input_matrix[j + 8][i],
                input_matrix[j + 7][i], input_matrix[j + 6][i], input_matrix[j + 5][i], input_matrix[j + 4][i],
                input_matrix[j + 3][i], input_matrix[j + 2][i], input_matrix[j + 1][i], input_matrix[j + 0][i],
                input_matrix[j + 15][i], input_matrix[j + 14][i], input_matrix[j + 13][i], input_matrix[j + 12][i],
                input_matrix[j + 11][i], input_matrix[j + 10][i], input_matrix[j + 9][i], input_matrix[j + 8][i],
                input_matrix[j + 7][i], input_matrix[j + 6][i], input_matrix[j + 5][i], input_matrix[j + 4][i],
                input_matrix[j + 3][i], input_matrix[j + 2][i], input_matrix[j + 1][i], input_matrix[j + 0][i],
                input_matrix[j + 15][i], input_matrix[j + 14][i], input_matrix[j + 13][i], input_matrix[j + 12][i],
                input_matrix[j + 11][i], input_matrix[j + 10][i], input_matrix[j + 9][i], input_matrix[j + 8][i],
                input_matrix[j + 7][i], input_matrix[j + 6][i], input_matrix[j + 5][i], input_matrix[j + 4][i],
                input_matrix[j + 3][i], input_matrix[j + 2][i], input_matrix[j + 1][i], input_matrix[j + 0][i],
                input_matrix[j + 15][i], input_matrix[j + 14][i], input_matrix[j + 13][i], input_matrix[j + 12][i],
                input_matrix[j + 11][i], input_matrix[j + 10][i], input_matrix[j + 9][i], input_matrix[j + 8][i],
                input_matrix[j + 7][i], input_matrix[j + 6][i], input_matrix[j + 5][i], input_matrix[j + 4][i],
                input_matrix[j + 3][i], input_matrix[j + 2][i], input_matrix[j + 1][i], input_matrix[j + 0][i]
            );

            // Store the row into the output matrix
            _mm512_storeu_si512((__m512i*)&output_matrix[i][j], row);
        }
    }
}

// since all the input is sequential we need to find the next block from the adjacent data block in the sequetial input. 
// for example if every data point is onnly one block deep. then width_to_adjacent_block = 1. if every data point is 2 blocks deep then width_to_adjacent_block = 2.
void bs_transpose_dst(word_t * transpose, word_t * blocks, word_t width_to_adjacent_block)
{
    word_t i,k;
    word_t w;
    for(k=0; k < WORD_SIZE; k++)
    {
        word_t bitpos = ONE << k;
        for (i=0; i < WORDS_PER_BLOCK; i++)
        {
            w = bs2le(blocks[k * WORDS_PER_BLOCK * width_to_adjacent_block + i]);
            word_t offset = i << MUL_SHIFT;

#ifndef UNROLL_TRANSPOSE
            word_t j;
            for(j=0; j < WORD_SIZE; j++)
            {
                // TODO make const time
                transpose[offset + j] |= (w & (ONE << j)) ? bitpos : 0;
            }
#else

            transpose[(offset)+ 0 ] |= (w & (ONE << 0 )) ? (bitpos) : 0;
            transpose[(offset)+ 1 ] |= (w & (ONE << 1 )) ? (bitpos) : 0;
            transpose[(offset)+ 2 ] |= (w & (ONE << 2 )) ? (bitpos) : 0;
            transpose[(offset)+ 3 ] |= (w & (ONE << 3 )) ? (bitpos) : 0;
            transpose[(offset)+ 4 ] |= (w & (ONE << 4 )) ? (bitpos) : 0;
            transpose[(offset)+ 5 ] |= (w & (ONE << 5 )) ? (bitpos) : 0;
            transpose[(offset)+ 6 ] |= (w & (ONE << 6 )) ? (bitpos) : 0;
            transpose[(offset)+ 7 ] |= (w & (ONE << 7 )) ? (bitpos) : 0;
#if WORD_SIZE > 8
            transpose[(offset)+ 8 ] |= (w & (ONE << 8 )) ? (bitpos) : 0;
            transpose[(offset)+ 9 ] |= (w & (ONE << 9 )) ? (bitpos) : 0;
            transpose[(offset)+ 10] |= (w & (ONE << 10)) ? (bitpos) : 0;
            transpose[(offset)+ 11] |= (w & (ONE << 11)) ? (bitpos) : 0;
            transpose[(offset)+ 12] |= (w & (ONE << 12)) ? (bitpos) : 0;
            transpose[(offset)+ 13] |= (w & (ONE << 13)) ? (bitpos) : 0;
            transpose[(offset)+ 14] |= (w & (ONE << 14)) ? (bitpos) : 0;
            transpose[(offset)+ 15] |= (w & (ONE << 15)) ? (bitpos) : 0;
#endif
#if WORD_SIZE > 16
            transpose[(offset)+ 16] |= (w & (ONE << 16)) ? (bitpos) : 0;
            transpose[(offset)+ 17] |= (w & (ONE << 17)) ? (bitpos) : 0;
            transpose[(offset)+ 18] |= (w & (ONE << 18)) ? (bitpos) : 0;
            transpose[(offset)+ 19] |= (w & (ONE << 19)) ? (bitpos) : 0;
            transpose[(offset)+ 20] |= (w & (ONE << 20)) ? (bitpos) : 0;
            transpose[(offset)+ 21] |= (w & (ONE << 21)) ? (bitpos) : 0;
            transpose[(offset)+ 22] |= (w & (ONE << 22)) ? (bitpos) : 0;
            transpose[(offset)+ 23] |= (w & (ONE << 23)) ? (bitpos) : 0;
            transpose[(offset)+ 24] |= (w & (ONE << 24)) ? (bitpos) : 0;
            transpose[(offset)+ 25] |= (w & (ONE << 25)) ? (bitpos) : 0;
            transpose[(offset)+ 26] |= (w & (ONE << 26)) ? (bitpos) : 0;
            transpose[(offset)+ 27] |= (w & (ONE << 27)) ? (bitpos) : 0;
            transpose[(offset)+ 28] |= (w & (ONE << 28)) ? (bitpos) : 0;
            transpose[(offset)+ 29] |= (w & (ONE << 29)) ? (bitpos) : 0;
            transpose[(offset)+ 30] |= (w & (ONE << 30)) ? (bitpos) : 0;
            transpose[(offset)+ 31] |= (w & (ONE << 31)) ? (bitpos) : 0;
#endif
#if WORD_SIZE > 32
            transpose[(offset)+ 32] |= (w & (ONE << 32)) ? (bitpos) : 0;
            transpose[(offset)+ 33] |= (w & (ONE << 33)) ? (bitpos) : 0;
            transpose[(offset)+ 34] |= (w & (ONE << 34)) ? (bitpos) : 0;
            transpose[(offset)+ 35] |= (w & (ONE << 35)) ? (bitpos) : 0;
            transpose[(offset)+ 36] |= (w & (ONE << 36)) ? (bitpos) : 0;
            transpose[(offset)+ 37] |= (w & (ONE << 37)) ? (bitpos) : 0;
            transpose[(offset)+ 38] |= (w & (ONE << 38)) ? (bitpos) : 0;
            transpose[(offset)+ 39] |= (w & (ONE << 39)) ? (bitpos) : 0;
            transpose[(offset)+ 40] |= (w & (ONE << 40)) ? (bitpos) : 0;
            transpose[(offset)+ 41] |= (w & (ONE << 41)) ? (bitpos) : 0;
            transpose[(offset)+ 42] |= (w & (ONE << 42)) ? (bitpos) : 0;
            transpose[(offset)+ 43] |= (w & (ONE << 43)) ? (bitpos) : 0;
            transpose[(offset)+ 44] |= (w & (ONE << 44)) ? (bitpos) : 0;
            transpose[(offset)+ 45] |= (w & (ONE << 45)) ? (bitpos) : 0;
            transpose[(offset)+ 46] |= (w & (ONE << 46)) ? (bitpos) : 0;
            transpose[(offset)+ 47] |= (w & (ONE << 47)) ? (bitpos) : 0;
            transpose[(offset)+ 48] |= (w & (ONE << 48)) ? (bitpos) : 0;
            transpose[(offset)+ 49] |= (w & (ONE << 49)) ? (bitpos) : 0;
            transpose[(offset)+ 50] |= (w & (ONE << 50)) ? (bitpos) : 0;
            transpose[(offset)+ 51] |= (w & (ONE << 51)) ? (bitpos) : 0;
            transpose[(offset)+ 52] |= (w & (ONE << 52)) ? (bitpos) : 0;
            transpose[(offset)+ 53] |= (w & (ONE << 53)) ? (bitpos) : 0;
            transpose[(offset)+ 54] |= (w & (ONE << 54)) ? (bitpos) : 0;
            transpose[(offset)+ 55] |= (w & (ONE << 55)) ? (bitpos) : 0;
            transpose[(offset)+ 56] |= (w & (ONE << 56)) ? (bitpos) : 0;
            transpose[(offset)+ 57] |= (w & (ONE << 57)) ? (bitpos) : 0;
            transpose[(offset)+ 58] |= (w & (ONE << 58)) ? (bitpos) : 0;
            transpose[(offset)+ 59] |= (w & (ONE << 59)) ? (bitpos) : 0;
            transpose[(offset)+ 60] |= (w & (ONE << 60)) ? (bitpos) : 0;
            transpose[(offset)+ 61] |= (w & (ONE << 61)) ? (bitpos) : 0;
            transpose[(offset)+ 62] |= (w & (ONE << 62)) ? (bitpos) : 0;
            transpose[(offset)+ 63] |= (w & (ONE << 63)) ? (bitpos) : 0;
#endif
#endif
                // constant time:
                //transpose[(i<<MUL_SHIFT)+ j] |= (((int64_t)((w & (ONE << j)) << (WORD_SIZE-1-j)))>>(WORD_SIZE-1)) & (ONE<<k);
        }
    }
}

// width_to_adjacent_block should be the same it was transposed with
void bs_transpose_rev(word_t * blocks, word_t width_to_adjacent_block)
{
    word_t i,k;
    word_t w;
    word_t transpose[BLOCK_SIZE];
    memset(transpose, 0, sizeof(transpose));
    for(k=0; k < BLOCK_SIZE; k++)
    {
        w = blocks[k];
        word_t bitpos = bs2be(ONE << (k % WORD_SIZE));
        word_t offset = k / WORD_SIZE;
#ifndef UNROLL_TRANSPOSE
        word_t j;
        for(j=0; j < WORD_SIZE; j++)
        {
            word_t bit = (w & (ONE << j)) ? (ONE << (k % WORD_SIZE)) : 0;
            transpose[j * WORDS_PER_BLOCK * width_to_adjacent_block + (offset)] |= bit;
        }
#else
        transpose[0  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 0 )) ? bitpos : 0;
        transpose[1  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 1 )) ? bitpos : 0;
        transpose[2  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 2 )) ? bitpos : 0;
        transpose[3  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 3 )) ? bitpos : 0;
        transpose[4  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 4 )) ? bitpos : 0;
        transpose[5  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 5 )) ? bitpos : 0;
        transpose[6  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 6 )) ? bitpos : 0;
        transpose[7  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 7 )) ? bitpos : 0;
#if WORD_SIZE > 8
        transpose[8  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 8 )) ? bitpos : 0;
        transpose[9  * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 9 )) ? bitpos : 0;
        transpose[10 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 10)) ? bitpos : 0;
        transpose[11 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 11)) ? bitpos : 0;
        transpose[12 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 12)) ? bitpos : 0;
        transpose[13 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 13)) ? bitpos : 0;
        transpose[14 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 14)) ? bitpos : 0;
        transpose[15 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 15)) ? bitpos : 0;
#endif
#if WORD_SIZE > 16
        transpose[16 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 16)) ? bitpos : 0;
        transpose[17 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 17)) ? bitpos : 0;
        transpose[18 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 18)) ? bitpos : 0;
        transpose[19 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 19)) ? bitpos : 0;
        transpose[20 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 20)) ? bitpos : 0;
        transpose[21 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 21)) ? bitpos : 0;
        transpose[22 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 22)) ? bitpos : 0;
        transpose[23 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 23)) ? bitpos : 0;
        transpose[24 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 24)) ? bitpos : 0;
        transpose[25 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 25)) ? bitpos : 0;
        transpose[26 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 26)) ? bitpos : 0;
        transpose[27 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 27)) ? bitpos : 0;
        transpose[28 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 28)) ? bitpos : 0;
        transpose[29 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 29)) ? bitpos : 0;
        transpose[30 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 30)) ? bitpos : 0;
        transpose[31 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 31)) ? bitpos : 0;
#endif
#if WORD_SIZE > 32
        transpose[32 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 32)) ? bitpos : 0;
        transpose[33 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 33)) ? bitpos : 0;
        transpose[34 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 34)) ? bitpos : 0;
        transpose[35 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 35)) ? bitpos : 0;
        transpose[36 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 36)) ? bitpos : 0;
        transpose[37 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 37)) ? bitpos : 0;
        transpose[38 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 38)) ? bitpos : 0;
        transpose[39 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 39)) ? bitpos : 0;
        transpose[40 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 40)) ? bitpos : 0;
        transpose[41 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 41)) ? bitpos : 0;
        transpose[42 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 42)) ? bitpos : 0;
        transpose[43 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 43)) ? bitpos : 0;
        transpose[44 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 44)) ? bitpos : 0;
        transpose[45 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 45)) ? bitpos : 0;
        transpose[46 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 46)) ? bitpos : 0;
        transpose[47 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 47)) ? bitpos : 0;
        transpose[48 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 48)) ? bitpos : 0;
        transpose[49 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 49)) ? bitpos : 0;
        transpose[50 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 50)) ? bitpos : 0;
        transpose[51 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 51)) ? bitpos : 0;
        transpose[52 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 52)) ? bitpos : 0;
        transpose[53 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 53)) ? bitpos : 0;
        transpose[54 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 54)) ? bitpos : 0;
        transpose[55 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 55)) ? bitpos : 0;
        transpose[56 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 56)) ? bitpos : 0;
        transpose[57 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 57)) ? bitpos : 0;
        transpose[58 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 58)) ? bitpos : 0;
        transpose[59 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 59)) ? bitpos : 0;
        transpose[60 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 60)) ? bitpos : 0;
        transpose[61 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 61)) ? bitpos : 0;
        transpose[62 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 62)) ? bitpos : 0;
        transpose[63 * WORDS_PER_BLOCK + (offset )] |= (w & (ONE << 63)) ? bitpos : 0;
#endif
#endif
    }
    memmove(blocks,transpose,sizeof(transpose));
// /    memcpy(blocks,transpose,sizeof(transpose));
}


void print_word_t_var(word_t var[8]) {
    printf("\n");
    for(int i = 0; i < 8; i++) {
        printf("%lu ", var[i]);
    }
    printf("\n");
}


void print_word_in_hex_and_binary(word_t word) {

    printf("Hex: %" PRIx64 "\n", word);
    for (int i = 63; i >= 0; i--) {
        printf("%llu", (word >> i) & 1);
    }
    printf("\n");
}


