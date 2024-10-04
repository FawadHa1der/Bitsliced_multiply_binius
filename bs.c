

#include <string.h>
#include "bs.h"
// #include <arm_neon.h>

#include <stdio.h>
#include <inttypes.h> // for PRId64 macro

#define bs2le(x) (x)
#define bs2be(x) (x)

void bs_transpose(word_t * blocks, word_t width_to_adjacent_block)
{
    word_t transpose[BLOCK_SIZE];
    memset(transpose, 0, sizeof(transpose));
    bs_transpose_dst(transpose,blocks, width_to_adjacent_block);

    int sizeof_transpose = sizeof(transpose);
    memmove(blocks,transpose,sizeof(transpose));


}

uint8_t multiply_8b_using_log_table(
    uint8_t lhs, uint8_t rhs,
    const uint8_t log_table[256],
    const uint8_t exp_table[256]
);

void multiply_constant_64b_using_log_table(
     uint64_t *rhs, uint64_t* result) {

    // uint8_t *lhs_bytes = (uint8_t *)lhs;
    uint8_t *rhs_bytes = (uint8_t *)rhs;
    uint8_t *result_bytes = (uint8_t *)result;
const uint8_t EXP_TABLE[256] = {
    0x1,  0x13, 0x43, 0x66, 0xab, 0x8c, 0x60, 0xc6, 0x91, 0xca, 0x59, 0xb2, 0x6a, 0x63, 0xf4, 0x53,
    0x17, 0x0f, 0xfa, 0xba, 0xee, 0x87, 0xd6, 0xe0, 0x6e, 0x2f, 0x68, 0x42, 0x75, 0xe8, 0xea, 0xcb,
    0x4a, 0xf1, 0x0c, 0xc8, 0x78, 0x33, 0xd1, 0x9e, 0x30, 0xe3, 0x5c, 0xed, 0xb5, 0x14, 0x3d, 0x38,
    0x67, 0xb8, 0xcf, 0x06, 0x6d, 0x1d, 0xaa, 0x9f, 0x23, 0xa0, 0x3a, 0x46, 0x39, 0x74, 0xfb, 0xa9,
    0xad, 0xe1, 0x7d, 0x6c, 0x0e, 0xe9, 0xf9, 0x88, 0x2c, 0x5a, 0x80, 0xa8, 0xbe, 0xa2, 0x1b, 0xc7,
    0x82, 0x89, 0x3f, 0x19, 0xe6, 0x03, 0x32, 0xc2, 0xdd, 0x56, 0x48, 0xd0, 0x8d, 0x73, 0x85, 0xf7,
    0x61, 0xd5, 0xd2, 0xac, 0xf2, 0x3e, 0x0a, 0xa5, 0x65, 0x99, 0x4e, 0xbd, 0x90, 0xd9, 0x1a, 0xd4,
    0xc1, 0xef, 0x94, 0x95, 0x86, 0xc5, 0xa3, 0x08, 0x84, 0xe4, 0x22, 0xb3, 0x79, 0x20, 0x92, 0xf8,
    0x9b, 0x6f, 0x3c, 0x2b, 0x24, 0xde, 0x64, 0x8a, 0xd,  0xdb, 0x3b, 0x55, 0x7a, 0x12, 0x50, 0x25,
    0xcd, 0x27, 0xec, 0xa6, 0x57, 0x5b, 0x93, 0xeb, 0xd8, 0x09, 0x97, 0xa7, 0x44, 0x18, 0xf5, 0x40,
    0x54, 0x69, 0x51, 0x36, 0x8e, 0x41, 0x47, 0x2a, 0x37, 0x9d, 0x02, 0x21, 0x81, 0xbb, 0xfd, 0xc4,
    0xb0, 0x4b, 0xe2, 0x4f, 0xae, 0xd3, 0xbf, 0xb1, 0x58, 0xa1, 0x29, 0x05, 0x5f, 0xdf, 0x77, 0xc9,
    0x6b, 0x70, 0xb7, 0x35, 0xbc, 0x83, 0x9a, 0x7c, 0x7f, 0x4d, 0x8f, 0x52, 0x04, 0x4c, 0x9c, 0x11,
    0x62, 0xe7, 0x10, 0x71, 0xa4, 0x76, 0xda, 0x28, 0x16, 0x1c, 0xb9, 0xdc, 0x45, 0x0b, 0xb6, 0x26,
    0xff, 0xe5, 0x31, 0xf0, 0x1f, 0x8b, 0x1e, 0x98, 0x5d, 0xfe, 0xf6, 0x72, 0x96, 0xb4, 0x07, 0x7e,
    0x5e, 0xcc, 0x34, 0xaf, 0xc0, 0xfc, 0xd7, 0xf3, 0x2d, 0x49, 0xc3, 0xce, 0x15, 0x2e, 0x7b, 0x00,
};

const uint8_t LOG_TABLE[256] = {
    0x00, 0x00, 0xaa, 0x55, 0xcc, 0xbb, 0x33, 0xee, 0x77, 0x99, 0x66, 0xdd, 0x22, 0x88, 0x44, 0x11,
    0xd2, 0xcf, 0x8d, 0x01, 0x2d, 0xfc, 0xd8, 0x10, 0x9d, 0x53, 0x6e, 0x4e, 0xd9, 0x35, 0xe6, 0xe4,
    0x7d, 0xab, 0x7a, 0x38, 0x84, 0x8f, 0xdf, 0x91, 0xd7, 0xba, 0xa7, 0x83, 0x48, 0xf8, 0xfd, 0x19,
    0x28, 0xe2, 0x56, 0x25, 0xf2, 0xc3, 0xa3, 0xa8, 0x2f, 0x3c, 0x3a, 0x8a, 0x82, 0x2e, 0x65, 0x52,
    0x9f, 0xa5, 0x1b, 0x02, 0x9c, 0xdc, 0x3b, 0xa6, 0x5a, 0xf9, 0x20, 0xb1, 0xcd, 0xc9, 0x6a, 0xb3,
    0x8e, 0xa2, 0xcb, 0x0f, 0xa0, 0x8b, 0x59, 0x94, 0xb8, 0x0a, 0x49, 0x95, 0x2a, 0xe8, 0xf0, 0xbc,
    0x06, 0x60, 0xd0, 0x0d, 0x86, 0x68, 0x03, 0x30, 0x1a, 0xa1, 0x0c, 0xc0, 0x43, 0x34, 0x18, 0x81,
    0xc1, 0xd3, 0xeb, 0x5d, 0x3d, 0x1c, 0xd5, 0xbe, 0x24, 0x7c, 0x8c, 0xfe, 0xc7, 0x42, 0xef, 0xc8,
    0x4a, 0xac, 0x50, 0xc5, 0x78, 0x5e, 0x74, 0x15, 0x47, 0x51, 0x87, 0xe5, 0x05, 0x5c, 0xa4, 0xca,
    0x6c, 0x08, 0x7e, 0x96, 0x72, 0x73, 0xec, 0x9a, 0xe7, 0x69, 0xc6, 0x80, 0xce, 0xa9, 0x27, 0x37,
    0x39, 0xb9, 0x4d, 0x76, 0xd4, 0x67, 0x93, 0x9b, 0x4b, 0x3f, 0x36, 0x04, 0x63, 0x40, 0xb4, 0xf3,
    0xb0, 0xb7, 0x0b, 0x7b, 0xed, 0x2c, 0xde, 0xc2, 0x31, 0xda, 0x13, 0xad, 0xc4, 0x6b, 0x4c, 0xb6,
    0xf4, 0x70, 0x57, 0xfa, 0xaf, 0x75, 0x07, 0x4f, 0x23, 0xbf, 0x09, 0x1f, 0xf1, 0x90, 0xfb, 0x32,
    0x5b, 0x26, 0x62, 0xb5, 0x6f, 0x61, 0x16, 0xf6, 0x98, 0x6d, 0xd6, 0x89, 0xdb, 0x58, 0x85, 0xbd,
    0x17, 0x41, 0xb2, 0x29, 0x79, 0xe1, 0x54, 0xd1, 0x1d, 0x45, 0x1e, 0x97, 0x92, 0x2b, 0x14, 0x71,
    0xe3, 0x21, 0x64, 0xf7, 0x0e, 0x9e, 0xea, 0x5f, 0x7f, 0x46, 0x12, 0x3e, 0xf5, 0xae, 0xe9, 0xe0,
};
    
    // Loop through each byte of lhs and rhs
    for (int i = 0; i < 8; i++) {
        
        result_bytes[i] = multiply_8b_using_log_table(16, rhs_bytes[i], LOG_TABLE, EXP_TABLE);
        
        // // Combine the result back into a uint64_t
        // result_bytes |= ((uint64_t)byte_result << (i * 8));
    }

}


// Wrapper function for uint64_t inputs
void multiply_64b_using_log_table(
    uint64_t *lhs, uint64_t *rhs, uint64_t* result) {

    uint8_t *lhs_bytes = (uint8_t *)lhs;
    uint8_t *rhs_bytes = (uint8_t *)rhs;
    uint8_t *result_bytes = (uint8_t *)result;
const uint8_t EXP_TABLE[256] = {
    0x1,  0x13, 0x43, 0x66, 0xab, 0x8c, 0x60, 0xc6, 0x91, 0xca, 0x59, 0xb2, 0x6a, 0x63, 0xf4, 0x53,
    0x17, 0x0f, 0xfa, 0xba, 0xee, 0x87, 0xd6, 0xe0, 0x6e, 0x2f, 0x68, 0x42, 0x75, 0xe8, 0xea, 0xcb,
    0x4a, 0xf1, 0x0c, 0xc8, 0x78, 0x33, 0xd1, 0x9e, 0x30, 0xe3, 0x5c, 0xed, 0xb5, 0x14, 0x3d, 0x38,
    0x67, 0xb8, 0xcf, 0x06, 0x6d, 0x1d, 0xaa, 0x9f, 0x23, 0xa0, 0x3a, 0x46, 0x39, 0x74, 0xfb, 0xa9,
    0xad, 0xe1, 0x7d, 0x6c, 0x0e, 0xe9, 0xf9, 0x88, 0x2c, 0x5a, 0x80, 0xa8, 0xbe, 0xa2, 0x1b, 0xc7,
    0x82, 0x89, 0x3f, 0x19, 0xe6, 0x03, 0x32, 0xc2, 0xdd, 0x56, 0x48, 0xd0, 0x8d, 0x73, 0x85, 0xf7,
    0x61, 0xd5, 0xd2, 0xac, 0xf2, 0x3e, 0x0a, 0xa5, 0x65, 0x99, 0x4e, 0xbd, 0x90, 0xd9, 0x1a, 0xd4,
    0xc1, 0xef, 0x94, 0x95, 0x86, 0xc5, 0xa3, 0x08, 0x84, 0xe4, 0x22, 0xb3, 0x79, 0x20, 0x92, 0xf8,
    0x9b, 0x6f, 0x3c, 0x2b, 0x24, 0xde, 0x64, 0x8a, 0xd,  0xdb, 0x3b, 0x55, 0x7a, 0x12, 0x50, 0x25,
    0xcd, 0x27, 0xec, 0xa6, 0x57, 0x5b, 0x93, 0xeb, 0xd8, 0x09, 0x97, 0xa7, 0x44, 0x18, 0xf5, 0x40,
    0x54, 0x69, 0x51, 0x36, 0x8e, 0x41, 0x47, 0x2a, 0x37, 0x9d, 0x02, 0x21, 0x81, 0xbb, 0xfd, 0xc4,
    0xb0, 0x4b, 0xe2, 0x4f, 0xae, 0xd3, 0xbf, 0xb1, 0x58, 0xa1, 0x29, 0x05, 0x5f, 0xdf, 0x77, 0xc9,
    0x6b, 0x70, 0xb7, 0x35, 0xbc, 0x83, 0x9a, 0x7c, 0x7f, 0x4d, 0x8f, 0x52, 0x04, 0x4c, 0x9c, 0x11,
    0x62, 0xe7, 0x10, 0x71, 0xa4, 0x76, 0xda, 0x28, 0x16, 0x1c, 0xb9, 0xdc, 0x45, 0x0b, 0xb6, 0x26,
    0xff, 0xe5, 0x31, 0xf0, 0x1f, 0x8b, 0x1e, 0x98, 0x5d, 0xfe, 0xf6, 0x72, 0x96, 0xb4, 0x07, 0x7e,
    0x5e, 0xcc, 0x34, 0xaf, 0xc0, 0xfc, 0xd7, 0xf3, 0x2d, 0x49, 0xc3, 0xce, 0x15, 0x2e, 0x7b, 0x00,
};

const uint8_t LOG_TABLE[256] = {
    0x00, 0x00, 0xaa, 0x55, 0xcc, 0xbb, 0x33, 0xee, 0x77, 0x99, 0x66, 0xdd, 0x22, 0x88, 0x44, 0x11,
    0xd2, 0xcf, 0x8d, 0x01, 0x2d, 0xfc, 0xd8, 0x10, 0x9d, 0x53, 0x6e, 0x4e, 0xd9, 0x35, 0xe6, 0xe4,
    0x7d, 0xab, 0x7a, 0x38, 0x84, 0x8f, 0xdf, 0x91, 0xd7, 0xba, 0xa7, 0x83, 0x48, 0xf8, 0xfd, 0x19,
    0x28, 0xe2, 0x56, 0x25, 0xf2, 0xc3, 0xa3, 0xa8, 0x2f, 0x3c, 0x3a, 0x8a, 0x82, 0x2e, 0x65, 0x52,
    0x9f, 0xa5, 0x1b, 0x02, 0x9c, 0xdc, 0x3b, 0xa6, 0x5a, 0xf9, 0x20, 0xb1, 0xcd, 0xc9, 0x6a, 0xb3,
    0x8e, 0xa2, 0xcb, 0x0f, 0xa0, 0x8b, 0x59, 0x94, 0xb8, 0x0a, 0x49, 0x95, 0x2a, 0xe8, 0xf0, 0xbc,
    0x06, 0x60, 0xd0, 0x0d, 0x86, 0x68, 0x03, 0x30, 0x1a, 0xa1, 0x0c, 0xc0, 0x43, 0x34, 0x18, 0x81,
    0xc1, 0xd3, 0xeb, 0x5d, 0x3d, 0x1c, 0xd5, 0xbe, 0x24, 0x7c, 0x8c, 0xfe, 0xc7, 0x42, 0xef, 0xc8,
    0x4a, 0xac, 0x50, 0xc5, 0x78, 0x5e, 0x74, 0x15, 0x47, 0x51, 0x87, 0xe5, 0x05, 0x5c, 0xa4, 0xca,
    0x6c, 0x08, 0x7e, 0x96, 0x72, 0x73, 0xec, 0x9a, 0xe7, 0x69, 0xc6, 0x80, 0xce, 0xa9, 0x27, 0x37,
    0x39, 0xb9, 0x4d, 0x76, 0xd4, 0x67, 0x93, 0x9b, 0x4b, 0x3f, 0x36, 0x04, 0x63, 0x40, 0xb4, 0xf3,
    0xb0, 0xb7, 0x0b, 0x7b, 0xed, 0x2c, 0xde, 0xc2, 0x31, 0xda, 0x13, 0xad, 0xc4, 0x6b, 0x4c, 0xb6,
    0xf4, 0x70, 0x57, 0xfa, 0xaf, 0x75, 0x07, 0x4f, 0x23, 0xbf, 0x09, 0x1f, 0xf1, 0x90, 0xfb, 0x32,
    0x5b, 0x26, 0x62, 0xb5, 0x6f, 0x61, 0x16, 0xf6, 0x98, 0x6d, 0xd6, 0x89, 0xdb, 0x58, 0x85, 0xbd,
    0x17, 0x41, 0xb2, 0x29, 0x79, 0xe1, 0x54, 0xd1, 0x1d, 0x45, 0x1e, 0x97, 0x92, 0x2b, 0x14, 0x71,
    0xe3, 0x21, 0x64, 0xf7, 0x0e, 0x9e, 0xea, 0x5f, 0x7f, 0x46, 0x12, 0x3e, 0xf5, 0xae, 0xe9, 0xe0,
};
    
    // Loop through each byte of lhs and rhs
    for (int i = 0; i < 8; i++) {
        
        result_bytes[i] = multiply_8b_using_log_table(lhs_bytes[i], rhs_bytes[i], LOG_TABLE, EXP_TABLE);
        
        // // Combine the result back into a uint64_t
        // result_bytes |= ((uint64_t)byte_result << (i * 8));
    }

}

uint8_t multiply_8b_using_log_table(
    uint8_t lhs, uint8_t rhs,
    const uint8_t log_table[256],
    const uint8_t exp_table[256]
) {
    uint8_t result = 0;
    
    if (lhs != 0 && rhs != 0) {
        size_t log_table_index = log_table[lhs] + log_table[rhs];

        if (log_table_index > 254) {
            log_table_index -= 255;
        }
        result = exp_table[log_table_index];
        printf("table look up lhs: %d, rhs: %d, result :%d  \n", lhs, rhs, result);
    }
    
    return result;
}

#define NUM_INPUTS 8        // Number of 128-bit numbers
#define BYTES_IN_128BIT 16   // 16 bytes in a 128-bit number
#define SLICED_OUTPUTS 16


//////////////////////////////////////////////////////
// Byte Slicing: Convert 8 x 128-bit inputs to 16 rows of 64 bits
//////////////////////////////////////////////////////
void byte_slice(uint128_t input[NUM_INPUTS], uint64_t output[SLICED_OUTPUTS]) {
    // Cast the input to uint8_t* for easier access to bytes
    uint8_t* input_bytes = (uint8_t*) input;
    uint8_t* output_bytes = (uint8_t*) output;

    // Loop over each byte index (0 to 15 for each 128-bit number)
    for (int byte_index = 0; byte_index < BYTES_IN_128BIT; byte_index++) {
        // Loop over each of the original 128-bit numbers (8 total inputs)
        for (int i = 0; i < NUM_INPUTS; i++) {
            // Place each byte from the input array into the corresponding 64-bit output row
            // the print the index being accesses
            int index_to = byte_index * NUM_INPUTS + i;
            int index_from = i * BYTES_IN_128BIT + byte_index;
            printf("index to: %d, from: %d\n", index_to, index_from);
            output_bytes[byte_index * NUM_INPUTS + i] = input_bytes[i * BYTES_IN_128BIT + byte_index];
        }
    }
}

//////////////////////////////////////////////////////
// Un-byte Slicing: Convert 16 rows of 64 bits back to 8 x 128-bit inputs
//////////////////////////////////////////////////////
void un_byte_slice(uint64_t input[SLICED_OUTPUTS], uint128_t output[NUM_INPUTS]) {
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


