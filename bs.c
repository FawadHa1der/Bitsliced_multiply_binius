

#include <string.h>
#include "bs.h"
// #include <arm_neon.h>

#include <stdio.h>
#include <inttypes.h> // for PRId64 macro
#include <arm_neon.h>

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


M128 lookup_16x8b(const uint8_t table[256], M128 x) {
    // Load the table as 4 sets of 4 16-byte vectors
    uint8x16x4_t tbl0 = vld1q_u8_x4(&table[0]);   // First 64 bytes (0-63)
    uint8x16x4_t tbl1 = vld1q_u8_x4(&table[64]);  // Second 64 bytes (64-127)
    uint8x16x4_t tbl2 = vld1q_u8_x4(&table[128]); // Third 64 bytes (128-191)
    uint8x16x4_t tbl3 = vld1q_u8_x4(&table[192]); // Fourth 64 bytes (192-255)

    // Perform table lookups
    uint8x16_t y0 = vqtbl4q_u8(tbl0, x);                         // Lookup in the first 64 bytes
    uint8x16_t y1 = vqtbl4q_u8(tbl1, veorq_u8(x, vdupq_n_u8(0x40))); // Lookup in the second 64 bytes
    uint8x16_t y2 = vqtbl4q_u8(tbl2, veorq_u8(x, vdupq_n_u8(0x80))); // Lookup in the third 64 bytes
    uint8x16_t y3 = vqtbl4q_u8(tbl3, veorq_u8(x, vdupq_n_u8(0xC0))); // Lookup in the fourth 64 bytes

    // Combine the results using XOR
    uint8x16_t result = veorq_u8(veorq_u8(y0, y1), veorq_u8(y2, y3));

    return result; // Return the combined result
}

uint8_t multiply_8b_using_log_table(
    uint8_t lhs, uint8_t rhs,
    const uint8_t log_table[256],
    const uint8_t exp_table[256]
);

uint8_t multiply_constant_8b_using_table(
    uint8_t rhs,
    const uint8_t alpha_table[256]
);
    const uint8_t BINARY_TOWER_8B_MUL_ALPHA_MAP  [256] = {
    0x00, 0x10, 0x20, 0x30, 0x40, 0x50, 0x60, 0x70, 0x80, 0x90, 0xa0, 0xb0, 0xc0, 0xd0, 0xe0, 0xf0,
    0x41, 0x51, 0x61, 0x71, 0x01, 0x11, 0x21, 0x31, 0xc1, 0xd1, 0xe1, 0xf1, 0x81, 0x91, 0xa1, 0xb1,
    0x82, 0x92, 0xa2, 0xb2, 0xc2, 0xd2, 0xe2, 0xf2, 0x02, 0x12, 0x22, 0x32, 0x42, 0x52, 0x62, 0x72,
    0xc3, 0xd3, 0xe3, 0xf3, 0x83, 0x93, 0xa3, 0xb3, 0x43, 0x53, 0x63, 0x73, 0x03, 0x13, 0x23, 0x33,
    0x94, 0x84, 0xb4, 0xa4, 0xd4, 0xc4, 0xf4, 0xe4, 0x14, 0x04, 0x34, 0x24, 0x54, 0x44, 0x74, 0x64,
    0xd5, 0xc5, 0xf5, 0xe5, 0x95, 0x85, 0xb5, 0xa5, 0x55, 0x45, 0x75, 0x65, 0x15, 0x05, 0x35, 0x25,
    0x16, 0x06, 0x36, 0x26, 0x56, 0x46, 0x76, 0x66, 0x96, 0x86, 0xb6, 0xa6, 0xd6, 0xc6, 0xf6, 0xe6,
    0x57, 0x47, 0x77, 0x67, 0x17, 0x07, 0x37, 0x27, 0xd7, 0xc7, 0xf7, 0xe7, 0x97, 0x87, 0xb7, 0xa7,
    0xe8, 0xf8, 0xc8, 0xd8, 0xa8, 0xb8, 0x88, 0x98, 0x68, 0x78, 0x48, 0x58, 0x28, 0x38, 0x08, 0x18,
    0xa9, 0xb9, 0x89, 0x99, 0xe9, 0xf9, 0xc9, 0xd9, 0x29, 0x39, 0x09, 0x19, 0x69, 0x79, 0x49, 0x59,
    0x6a, 0x7a, 0x4a, 0x5a, 0x2a, 0x3a, 0x0a, 0x1a, 0xea, 0xfa, 0xca, 0xda, 0xaa, 0xba, 0x8a, 0x9a,
    0x2b, 0x3b, 0x0b, 0x1b, 0x6b, 0x7b, 0x4b, 0x5b, 0xab, 0xbb, 0x8b, 0x9b, 0xeb, 0xfb, 0xcb, 0xdb,
    0x7c, 0x6c, 0x5c, 0x4c, 0x3c, 0x2c, 0x1c, 0x0c, 0xfc, 0xec, 0xdc, 0xcc, 0xbc, 0xac, 0x9c, 0x8c,
    0x3d, 0x2d, 0x1d, 0x0d, 0x7d, 0x6d, 0x5d, 0x4d, 0xbd, 0xad, 0x9d, 0x8d, 0xfd, 0xed, 0xdd, 0xcd,
    0xfe, 0xee, 0xde, 0xce, 0xbe, 0xae, 0x9e, 0x8e, 0x7e, 0x6e, 0x5e, 0x4e, 0x3e, 0x2e, 0x1e, 0x0e,
    0xbf, 0xaf, 0x9f, 0x8f, 0xff, 0xef, 0xdf, 0xcf, 0x3f, 0x2f, 0x1f, 0x0f, 0x7f, 0x6f, 0x5f, 0x4f,
    };

void multiply_constant_128b_using_table(
     M128 *rhs, M128* result) {
    

    // uint8_t *lhs_bytes = (uint8_t *)lhs;
    // uint8_t *rhs_bytes = (uint8_t *)rhs;
    // uint8_t *result_bytes = (uint8_t *)result;
    *result = lookup_16x8b(BINARY_TOWER_8B_MUL_ALPHA_MAP, *rhs);
    
}
const uint8_t EXP_TABLE[256] = {
	0x01, 0x13, 0x43, 0x66, 0xAB, 0x8C, 0x60, 0xC6, 0x91, 0xCA, 0x59, 0xB2, 0x6A, 0x63, 0xF4, 0x53,
	0x17, 0x0F, 0xFA, 0xBA, 0xEE, 0x87, 0xD6, 0xE0, 0x6E, 0x2F, 0x68, 0x42, 0x75, 0xE8, 0xEA, 0xCB,
	0x4A, 0xF1, 0x0C, 0xC8, 0x78, 0x33, 0xD1, 0x9E, 0x30, 0xE3, 0x5C, 0xED, 0xB5, 0x14, 0x3D, 0x38,
	0x67, 0xB8, 0xCF, 0x06, 0x6D, 0x1D, 0xAA, 0x9F, 0x23, 0xA0, 0x3A, 0x46, 0x39, 0x74, 0xFB, 0xA9,
	0xAD, 0xE1, 0x7D, 0x6C, 0x0E, 0xE9, 0xF9, 0x88, 0x2C, 0x5A, 0x80, 0xA8, 0xBE, 0xA2, 0x1B, 0xC7,
	0x82, 0x89, 0x3F, 0x19, 0xE6, 0x03, 0x32, 0xC2, 0xDD, 0x56, 0x48, 0xD0, 0x8D, 0x73, 0x85, 0xF7,
	0x61, 0xD5, 0xD2, 0xAC, 0xF2, 0x3E, 0x0A, 0xA5, 0x65, 0x99, 0x4E, 0xBD, 0x90, 0xD9, 0x1A, 0xD4,
	0xC1, 0xEF, 0x94, 0x95, 0x86, 0xC5, 0xA3, 0x08, 0x84, 0xE4, 0x22, 0xB3, 0x79, 0x20, 0x92, 0xF8,
	0x9B, 0x6F, 0x3C, 0x2B, 0x24, 0xDE, 0x64, 0x8A, 0x0D, 0xDB, 0x3B, 0x55, 0x7A, 0x12, 0x50, 0x25,
	0xCD, 0x27, 0xEC, 0xA6, 0x57, 0x5B, 0x93, 0xEB, 0xD8, 0x09, 0x97, 0xA7, 0x44, 0x18, 0xF5, 0x40,
	0x54, 0x69, 0x51, 0x36, 0x8E, 0x41, 0x47, 0x2A, 0x37, 0x9D, 0x02, 0x21, 0x81, 0xBB, 0xFD, 0xC4,
	0xB0, 0x4B, 0xE2, 0x4F, 0xAE, 0xD3, 0xBF, 0xB1, 0x58, 0xA1, 0x29, 0x05, 0x5F, 0xDF, 0x77, 0xC9,
	0x6B, 0x70, 0xB7, 0x35, 0xBC, 0x83, 0x9A, 0x7C, 0x7F, 0x4D, 0x8F, 0x52, 0x04, 0x4C, 0x9C, 0x11,
	0x62, 0xE7, 0x10, 0x71, 0xA4, 0x76, 0xDA, 0x28, 0x16, 0x1C, 0xB9, 0xDC, 0x45, 0x0B, 0xB6, 0x26,
	0xFF, 0xE5, 0x31, 0xF0, 0x1F, 0x8B, 0x1E, 0x98, 0x5D, 0xFE, 0xF6, 0x72, 0x96, 0xB4, 0x07, 0x7E,
	0x5E, 0xCC, 0x34, 0xAF, 0xC0, 0xFC, 0xD7, 0xF3, 0x2D, 0x49, 0xC3, 0xCE, 0x15, 0x2E, 0x7B, 0x00,
};

const uint8_t LOG_TABLE[256] = {
	0x00, 0x00, 0xAA, 0x55, 0xCC, 0xBB, 0x33, 0xEE, 0x77, 0x99, 0x66, 0xDD, 0x22, 0x88, 0x44, 0x11,
	0xD2, 0xCF, 0x8D, 0x01, 0x2D, 0xFC, 0xD8, 0x10, 0x9D, 0x53, 0x6E, 0x4E, 0xD9, 0x35, 0xE6, 0xE4,
	0x7D, 0xAB, 0x7A, 0x38, 0x84, 0x8F, 0xDF, 0x91, 0xD7, 0xBA, 0xA7, 0x83, 0x48, 0xF8, 0xFD, 0x19,
	0x28, 0xE2, 0x56, 0x25, 0xF2, 0xC3, 0xA3, 0xA8, 0x2F, 0x3C, 0x3A, 0x8A, 0x82, 0x2E, 0x65, 0x52,
	0x9F, 0xA5, 0x1B, 0x02, 0x9C, 0xDC, 0x3B, 0xA6, 0x5A, 0xF9, 0x20, 0xB1, 0xCD, 0xC9, 0x6A, 0xB3,
	0x8E, 0xA2, 0xCB, 0x0F, 0xA0, 0x8B, 0x59, 0x94, 0xB8, 0x0A, 0x49, 0x95, 0x2A, 0xE8, 0xF0, 0xBC,
	0x06, 0x60, 0xD0, 0x0D, 0x86, 0x68, 0x03, 0x30, 0x1A, 0xA1, 0x0C, 0xC0, 0x43, 0x34, 0x18, 0x81,
	0xC1, 0xD3, 0xEB, 0x5D, 0x3D, 0x1C, 0xD5, 0xBE, 0x24, 0x7C, 0x8C, 0xFE, 0xC7, 0x42, 0xEF, 0xC8,
	0x4A, 0xAC, 0x50, 0xC5, 0x78, 0x5E, 0x74, 0x15, 0x47, 0x51, 0x87, 0xE5, 0x05, 0x5C, 0xA4, 0xCA,
	0x6C, 0x08, 0x7E, 0x96, 0x72, 0x73, 0xEC, 0x9A, 0xE7, 0x69, 0xC6, 0x80, 0xCE, 0xA9, 0x27, 0x37,
	0x39, 0xB9, 0x4D, 0x76, 0xD4, 0x67, 0x93, 0x9B, 0x4B, 0x3F, 0x36, 0x04, 0x63, 0x40, 0xB4, 0xF3,
	0xB0, 0xB7, 0x0B, 0x7B, 0xED, 0x2C, 0xDE, 0xC2, 0x31, 0xDA, 0x13, 0xAD, 0xC4, 0x6B, 0x4C, 0xB6,
	0xF4, 0x70, 0x57, 0xFA, 0xAF, 0x75, 0x07, 0x4F, 0x23, 0xBF, 0x09, 0x1F, 0xF1, 0x90, 0xFB, 0x32,
	0x5B, 0x26, 0x62, 0xB5, 0x6F, 0x61, 0x16, 0xF6, 0x98, 0x6D, 0xD6, 0x89, 0xDB, 0x58, 0x85, 0xBD,
	0x17, 0x41, 0xB2, 0x29, 0x79, 0xE1, 0x54, 0xD1, 0x1D, 0x45, 0x1E, 0x97, 0x92, 0x2B, 0x14, 0x71,
	0xE3, 0x21, 0x64, 0xF7, 0x0E, 0x9E, 0xEA, 0x5F, 0x7F, 0x46, 0x12, 0x3E, 0xF5, 0xAE, 0xE9, 0xE0,
};


M128 packed_tower_16x8b_multiply(M128 a, M128 b) {
    // Look up the logarithms of a and b

	// let loga = lookup_16x8b(TOWER_LOG_LOOKUP_TABLE, a).into();
	// let logb = lookup_16x8b(TOWER_LOG_LOOKUP_TABLE, b).into();
	// let logc = unsafe {
	// 	let sum = vaddq_u8(loga, logb);
	// 	let overflow = vcgtq_u8(loga, sum);
	// 	vsubq_u8(sum, overflow)
	// };
	// let c = lookup_16x8b(TOWER_EXP_LOOKUP_TABLE, logc.into()).into();
	// unsafe {
	// 	let a_or_b_is_0 = vorrq_u8(vceqzq_u8(a.into()), vceqzq_u8(b.into()));
	// 	vandq_u8(c, veorq_u8(a_or_b_is_0, M128::fill_with_bit(1).into()))
	// }
	// .into()

    M128 loga = lookup_16x8b(LOG_TABLE, a);
    M128 logb = lookup_16x8b(LOG_TABLE, b);

    // Sum the logarithms and handle overflow using NEON intrinsics
    uint8x16_t sum = vaddq_u8(loga, logb);          // Add logarithms
    uint8x16_t overflow = vcgtq_u8(loga, sum);      // Overflow detection
    uint8x16_t logc = vsubq_u8(sum, overflow);      // Subtract overflow

    // Look up the exponentiation result
    M128 c = lookup_16x8b(EXP_TABLE, logc);

    // Handle case where either a or b is zero
    // uint8x16_t a_is_zero = vceqq_u8(a, vdupq_n_u8(0)); // a == 0 ?
    // uint8x16_t b_is_zero = vceqq_u8(b, vdupq_n_u8(0)); // b == 0 ?
    uint8x16_t a_or_b_is_zero = vorrq_u8(vceqzq_u8(a), vceqzq_u8(b)); // a == 0 || b == 0

    // XOR the result with the condition (if a or b is 0, result should be 0)
    uint8x16_t final_result = vandq_u8(c, veorq_u8(a_or_b_is_zero, vdupq_n_u8(0xFF)));

    // Return the final result
    return final_result;
}
// Wrapper function for uint64_t inputs
void multiply_128b_using_log_table(
    uint8x16_t *lhs, uint8x16_t *rhs, uint8x16_t* result) {
    *result = packed_tower_16x8b_multiply(*lhs, *rhs);
   return;
}




// bsically at bit level we always multiply by 16 in the algorithm and us ethe look up table to get the result.
uint8_t multiply_constant_8b_using_table(
    uint8_t rhs,
    const uint8_t alpha_table[256]
) {
    uint8_t result = 0;
    result = alpha_table[rhs];

    
    // if (lhs != 0 && rhs != 0) {
    //     size_t log_table_index = log_table[lhs] + log_table[rhs];

    //     if (log_table_index > 254) {
    //         log_table_index -= 255;
    //     }
    //     result = exp_table[log_table_index];
    //  //   printf("table look up lhs: %d, rhs: %d, result :%d  \n", lhs, rhs, result);
    // }
    
    return result;
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
     //   printf("table look up lhs: %d, rhs: %d, result :%d  \n", lhs, rhs, result);
    }
    
    return result;
}

#define NUM_INPUTS 16        // Number of 128-bit numbers
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
            // int index_to = byte_index * NUM_INPUTS + i;
            // int index_from = i * BYTES_IN_128BIT + byte_index;
            // printf("index to: %d, from: %d\n", index_to, index_from);
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


