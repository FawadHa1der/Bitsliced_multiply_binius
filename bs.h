
#ifndef _BS_H_
#define _BS_H_

#include <stdint.h>
#include <immintrin.h>

#define BLOCK_SIZE          128
#define BLOCK_SIZE_BYTES    (BLOCK_SIZE / 8)
// #define KEY_SCHEDULE_SIZE   176
#define WORD_SIZE           64
#define WORD_SIZE_BYTES     64

#define BS_BLOCK_SIZE       (BLOCK_SIZE * WORD_SIZE / 8)
#define WORDS_PER_BLOCK     (BLOCK_SIZE / WORD_SIZE)

#if (WORD_SIZE==64)
    typedef uint64_t    word_t;
    #define ONE         1ULL
    #define MUL_SHIFT   6
    #define WFMT        "lx"
    #define WPAD        "016"
    #define __builtin_bswap_wordsize(x) __builtin_bswap64(x)
#elif (WORD_SIZE==32)
    typedef uint32_t    word_t;
    #define ONE         1UL
    #define MUL_SHIFT   5
    #define WFMT        "x"
    #define WPAD        "08"
    #define __builtin_bswap_wordsize(x) __builtin_bswap32(x)
#elif (WORD_SIZE==16)
    typedef uint16_t    word_t;
    #define ONE         1
    #define MUL_SHIFT   4
    #define WFMT        "hx"
    #define WPAD        "04"
    #define __builtin_bswap_wordsize(x) __builtin_bswap16(x)
#elif (WORD_SIZE==8)
    typedef uint8_t     word_t;
    #define ONE         1
    #define MUL_SHIFT   3
    #define WFMT        "hhx"
    #define WPAD        "02"
    #define __builtin_bswap_wordsize(x) (x)
#else
#error "invalid word size"
#endif

#define UNROLL_TRANSPOSE 1
void bs_transpose(word_t * blocks, word_t width_to_adjacent_block);
void bs_transpose_rev(word_t * blocks, word_t width_to_adjacent_block);
void bs_transpose_dst(word_t * transpose, word_t * blocks, word_t width_to_adjacent_block);

// Assuming uint128_t is represented as two uint64_t for low and high parts
typedef struct {
    uint64_t low;
    uint64_t high;
} uint128_t;

void byte_slice(uint128_t *input, uint64_t *output);
void byte_slice_avx(__m512i *input, __m512i *output);
void byte_slice_avx_2_mul(__m512i *input_1, __m512i *input_2, __m512i *slice_output_1, __m512i *slice_output_2, __m512i *mul_level_0, __m512i *xor_level_0 );
__m512i wrapper_mm512_xor_si512(__m512i input_1, __m512i input_2);
void un_byte_slice(uint64_t* input, uint128_t *output);
void un_byte_slice_avx(uint64_t* input, uint128_t *output);
void initialize_maps();
__m512i gfni_mul(__m512i lhs, __m512i rhs);
void multiply_512b_using_gfni(__m512i *lhs, __m512i *rhs,__m512i* result );
void multiply_constant_512b_using_gfni(__m512i *rhs, __m512i* result );
#endif
