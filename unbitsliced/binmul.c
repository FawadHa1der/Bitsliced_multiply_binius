#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
//neon
#include <arm_neon.h>
typedef struct {
    uint64_t low;
    uint64_t high;
} uint128_t;

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
    }
       // printf("table look up lhs: %d, rhs: %d, result :%d  \n", lhs, rhs, result);

    return result;
}


uint64_t binmul64(uint64_t v1, uint64_t v2, uint32_t length, bool is_constant);

int calls = 0;
// Function to multiply two 128-bit integers using binary multiplication based on binius tower construction
uint128_t binmul128(uint128_t v1, uint128_t v2, uint32_t length) {
    // Print the inputs and length for debugging
    calls++;
    uint32_t halflen = length / 2;
    uint32_t quarterlen = length / 4;

    uint64_t halfmask =0;
    if (halflen < 64) {
        halfmask = (1ULL << halflen) - 1;
    } else if (halflen == 64) {
        halfmask = ~0ULL;  // Equivalent to 0xFFFFFFFFFFFFFFFF
    }
    else {
        halfmask = ~0ULL;  // Equivalent to 0xFFFFFFFFFFFFFFFF
        halfmask = (1ULL << (halflen - 64)) - 1;
    }
    uint64_t L1, R1, L2, R2;

    if (length == 128) {
        L1 = v1.low;
        R1 = v1.high;
        L2 = v2.low;
        R2 = v2.high;
    }

    if (L1 == 0 && R1 == 1 ) {
        uint64_t outR_input = 1ULL << quarterlen;
        uint128_t outR;
        outR.high = 0;
        outR.low = binmul64(outR_input,  R2, halflen, false);
        outR.low ^= L2;
        uint128_t ret_value = {(R2 ^ (outR.low << halflen)), 0};
        return ret_value;
    }

    uint128_t L1L2;
    L1L2.high = 0;
    L1L2.low = binmul64(L1, L2, halflen, false);
    uint128_t R1R2;
    R1R2.high = 0;
    R1R2.low = binmul64(R1, R2, halflen, false);
    uint64_t R1R2_high_input = (1ULL << quarterlen);
    uint128_t R1R2_high;
       R1R2_high.high = 0;
    R1R2_high.low = binmul64(R1R2_high_input, R1R2.low, halflen, false);
    

    uint64_t Z3_input_v1 = L1 ^ R1;
    uint64_t Z3_input_v2 = L2 ^ R2;

    uint128_t Z3;
    Z3.high = 0;
    Z3.low = binmul64( Z3_input_v1, Z3_input_v2, halflen, false);

    uint128_t result;
    if (length >= 128) {
        result = (uint128_t) {
            L1L2.low ^ R1R2.low, 
            Z3.low ^ L1L2.low ^ R1R2.low ^ R1R2_high.low
        };
    } else {
        result = (uint128_t) {
            L1L2.low ^ R1R2.low ^ ((Z3.low ^ L1L2.low ^ R1R2.low ^ R1R2_high.low) << halflen), 0
        };
    }

    // printf("binmul128 result:  ");
    // printf("%016llx%016llx\n ", result.high, result.low);

    return result;
}


uint64_t binmul64(uint64_t v1, uint64_t v2, uint32_t length, bool is_constant) {


    if (v1 < 2 || v2 < 2)  {
        uint64_t result = v1 * v2;
        return result;
    }

    if (length == 8){

        uint64_t result = multiply_8b_using_log_table(v1, v2, LOG_TABLE, EXP_TABLE);        
        return result;
    }

    uint32_t halflen = length / 2;
    uint32_t quarterlen = length / 4;

    uint64_t halfmask =0;
    halfmask = (1ULL << halflen) - 1;

    uint64_t L1, R1, L2, R2;

    L1 = v1 & halfmask;
    R1 = v1 >> halflen;

    L2 = v2 & halfmask;
    R2 = v2 >> halflen;

    if (L1 == 0 && R1 == 1) {
        uint64_t outR_input = 1ULL << quarterlen;
        uint64_t outR = binmul64(outR_input,  R2, halflen, true);
        outR ^= L2;
        uint64_t ret_value = (R2 ^ (outR << halflen));
        return ret_value;
    }

    uint64_t L1L2 = binmul64(L1, L2, halflen, false);
    uint64_t R1R2 = binmul64(R1, R2, halflen, false);

    uint64_t R1R2_high_input = (1ULL << quarterlen);
    uint64_t R1R2_high = binmul64(R1R2_high_input, R1R2, halflen, true);
    

    uint64_t Z3_input_v1 = L1 ^ R1;
    uint64_t Z3_input_v2 = L2 ^ R2;

    uint64_t Z3 = binmul64( Z3_input_v1, Z3_input_v2, halflen, false);

    // print all the values for L1L2, R1R2, R1R2_high, Z3 in one line as hex


    // printf("L1L2: %04llx, R1R2: %04llx, R1R2_high: %04llx, Z3: %04llx\n", L1L2, R1R2, R1R2_high, Z3);
    uint64_t upper_result =  (Z3 ^ L1L2 ^ R1R2 ^ R1R2_high) ;

    // printf("upper_result: %04llx\n", upper_result);

    uint64_t result = (uint64_t) L1L2 ^ R1R2 ^ ((upper_result) << halflen);
    return result;
}


// -----------------------------------------------------------------------------
//  Build the 128×128 constant‐times matrix for C.
//    cols[j] = C * (1 << j),  j = 0..127
// -----------------------------------------------------------------------------
void build_matrix128( uint128_t C, uint128_t cols[128] )
{
    for ( int j = 0; j < 128; ++j )
    {
        uint128_t E;
        if ( j < 64 )
        {
            E.low  = (uint64_t)1 << j;
            E.high = 0;
        }
        else
        {
            E.low  = 0;
            E.high = (uint64_t)1 << (j - 64);
        }
        cols[j] = binmul128( C, E, 128 );
    }
}


// Transpose the 128×128 bit–matrix in “cols” into “rows”.
// rows[i][j] = cols[j][i]
static inline void transpose128( const uint128_t cols[128],
                                 uint128_t       rows[128] )
{
    for ( int i = 0; i < 128; ++i )
    {
        uint64_t lo = 0, hi = 0;

        // Build row i, bit by bit:
        for ( int j = 0; j < 128; ++j )
        {
            // extract bit _i_ from column j
            unsigned b;
            if ( i < 64 )
                b = (cols[j].low  >> i) & 1;
            else
                b = (cols[j].high >> (i - 64)) & 1;

            // scatter it into row[i] at position j
            if ( j < 64 )
                lo |= (uint64_t)b << j;
            else
                hi |= (uint64_t)b << (j - 64);
        }

        rows[i].low  = lo;
        rows[i].high = hi;
    }
}
// -----------------------------------------------------------------------------
//  Multiply X by C via the precomputed columns:
//    result = XOR_{j : bit j of X is 1} cols[j].
// -----------------------------------------------------------------------------
uint128_t mul_via_matrix( const uint128_t cols[128], uint128_t X )
{
    // 1) Initialize the accumulator to zero
    //
    uint128_t result = { 0, 0 };

    // 2) As long as X has any 1-bit left…
    //
    while ( X.low  != 0  ||  X.high != 0 )
    {
        unsigned idx;

        // 3) Find the index of the *lowest* set bit in X
        //
        //    - __builtin_ctzll(v) returns “count trailing zeros” in a 64-bit word.
        //      If v =   0b...0101000,  ctzll(v) = 3 (the 0-based position of that single 1).
        //
        if ( X.low != 0 )
        {
            // 3a) If any bit in the low half is 1, pick that first
            idx = __builtin_ctzll( X.low );
        }
        else
        {
            // 3b) Otherwise look in the high half, but add 64 to index bits 64–127
            idx = 64U + __builtin_ctzll( X.high );
        }

        // 4) XOR in the precomputed column for that bit-position
        //
        //    “cols[idx]” is the 128-bit vector = C · E_idx,
        //    so XOR’ing it accumulates C·(sum of all chosen E_idx) = C·X.
        //
        result.low  ^= cols[idx].low;
        result.high ^= cols[idx].high;

        // 5) Clear that lowest set bit so we’ll move on to the next one
        //
        if ( idx < 64 )
            X.low  &= X.low  - 1;  // trick: v & (v−1) clears the least significant 1
        else
            X.high &= X.high - 1;
    }

    // 6) When X==0, we’ve XOR’d in every column whose bit was 1 → final product
    return result;
}

uint128_t mul_via_matrix_rows( const uint128_t rows[128], uint128_t X )
{
  uint128_t result = { 0, 0 };

  for ( int i = 0; i < 128; ++i )
  {
    // 1) mask off only the X-bits this row cares about
    uint64_t lo = rows[i].low  & X.low;
    uint64_t hi = rows[i].high & X.high;

    // 2) compute parity of those 128 bits
    //    __builtin_parityll returns popcount(x)&1
    unsigned bit = __builtin_parityll( lo ) ^ __builtin_parityll( hi );

    // 3) scatter that 1-bit into the correct position of result
    if ( bit )
    {
      if ( i < 64 )
        result.low  |= ((uint64_t)1) << i;
      else
        result.high |= ((uint64_t)1) << (i - 64);
    }
  }

  return result;
}


void mul_via_matrix_rows_bitsliced( const uint128_t rows[128],
                                    const uint128_t X[128],
                                          uint128_t out[128] )
{
  // For each output bit i = 0..127, compute the GF(2) dot-product:
  //    out[i] = ⊕_{j | rows[i][j] == 1} X[j]
  for ( int i = 0; i < 128; ++i )
  {
    // start accumulator at zero
    uint128_t acc = {0,0};

    // for columns 0..63, test rows[i].low
    uint64_t row_lo = rows[i].low;
    while ( row_lo )
    {
      // peel off lowest set bit
      int j = __builtin_ctzll( row_lo );
      row_lo &= row_lo - 1;
      // XOR in that slice
      acc.low  ^= X[j].low;
      acc.high ^= X[j].high;
    }

    // for columns 64..127, test rows[i].high
    uint64_t row_hi = rows[i].high;
    while ( row_hi )
    {
      int k = __builtin_ctzll( row_hi );
      row_hi &= row_hi - 1;
      int j = 64 + k;
      acc.low  ^= X[j].low;
      acc.high ^= X[j].high;
    }

    out[i] = acc;
  }
}

// -----------------------------------------------------------------------------
//  Generate a “random” 128‐bit value by combining rand() calls
// -----------------------------------------------------------------------------
static uint128_t random_u128( void )
{
    uint128_t r;
    r.low  = ((uint64_t)rand() << 32) ^ ((uint64_t)rand() << 16) ^ rand();
    r.high = ((uint64_t)rand() << 32) ^ ((uint64_t)rand() << 16) ^ rand();
    return r;
}

// -----------------------------------------------------------------------------
//  Test that matrix128 + mul_via_matrix matches binmul128 for C
// -----------------------------------------------------------------------------
void test_matrix128( uint128_t C )
{
    uint128_t cols[128], rows[128];
    build_matrix128( C, cols );
    transpose128( cols, rows );

    // 1) Check each basis vector
    for ( int j = 0; j < 128; ++j )
    {
        uint128_t E = {0,0};
        if ( j < 64 )   E.low  = (uint64_t)1 << j;
        else            E.high = (uint64_t)1 << (j - 64);
        uint128_t expect = binmul128( C, E, 128 );
        // uint128_t actual = mul_via_matrix( cols, E );
        uint128_t actual = mul_via_matrix_rows( rows, E );
        assert( expect.low  == actual.low  );
        assert( expect.high == actual.high );
    }

    // 2) Check random vectors
    srand( 42 );
    for ( int i = 0; i < 1000; ++i )
    {
        uint128_t X = random_u128();
        uint128_t expect = binmul128( C, X, 128 );
        uint128_t actual = mul_via_matrix( cols, X );
        assert( expect.low  == actual.low  );
        assert( expect.high == actual.high );
    }

    printf( "✅  matrix128 build & test passed for given constant C\n" );
}

// -----------------------------------------------------------
//  2) Parallel run‐time multiply: bitsliced X[0..127] → Z[0..127].
//     - X[j].bit[i] is bit-j of the i-th 128-bit input.
//     - Z[k].bit[i] = bit-k of (C * (input_i)).
// 
//  Formula:
//    for each output bit k and each parallel index i:
//      Z[k][i] = ⨁_{j: M[k][j]=1} X[j][i]
//
//  Here we loop k then j, XOR’ing the entire 128-bit slice X[j] 
//  into Z[k] whenever the matrix bit M[k][j] is 1.
// -----------------------------------------------------------
// void mul_via_matrix_bitsliced( const uint128_t cols[128],
//                                 const uint128_t X[128],
//                                       uint128_t Z[128] )
// {
//     // 1) Initialize the accumulator to zero
//     //
//     uint128_t result = { 0, 0 };

//     // 2) As long as X has any 1-bit left…
//     //
//     while ( X.low  != 0  ||  X.high != 0 )
//     {
//         unsigned idx;

//         // 3) Find the index of the *lowest* set bit in X
//         //
//         //    - __builtin_ctzll(v) returns “count trailing zeros” in a 64-bit word.
//         //      If v =   0b...0101000,  ctzll(v) = 3 (the 0-based position of that single 1).
//         //
//         if ( X.low != 0 )
//         {
//             // 3a) If any bit in the low half is 1, pick that first
//             idx = __builtin_ctzll( X.low );
//         }
//         else
//         {
//             // 3b) Otherwise look in the high half, but add 64 to index bits 64–127
//             idx = 64U + __builtin_ctzll( X.high );
//         }

//         // 4) XOR in the precomputed column for that bit-position
//         //
//         //    “cols[idx]” is the 128-bit vector = C · E_idx,
//         //    so XOR’ing it accumulates C·(sum of all chosen E_idx) = C·X.
//         //
//         result.low  ^= cols[idx].low;
//         result.high ^= cols[idx].high;

//         // 5) Clear that lowest set bit so we’ll move on to the next one
//         //
//         if ( idx < 64 )
//             X.low  &= X.low  - 1;  // trick: v & (v−1) clears the least significant 1
//         else
//             X.high &= X.high - 1;
//     }

//     // 6) When X==0, we’ve XOR’d in every column whose bit was 1 → final product
//     return result;
// }

void pre_calculate_lookup_table(uint128_t input[8], uint128_t output[256]) {
    // populate all the prev and idx values into arrays and print them after the loop
static uint8_t prev_arr[256] = { 0x00, 0x00, 0x00, 0x02, 0x00, 0x04, 0x04, 0x06, 0x00, 0x08, 0x08, 0x0a, 0x08, 0x0c, 0x0c, 0x0e, 0x00, 0x10, 0x10, 0x12, 0x10, 0x14, 0x14, 0x16, 0x10, 0x18, 0x18, 0x1a, 0x18, 0x1c, 0x1c, 0x1e, 0x00, 0x20, 0x20, 0x22, 0x20, 0x24, 0x24, 0x26, 0x20, 0x28, 0x28, 0x2a, 0x28, 0x2c, 0x2c, 0x2e, 0x20, 0x30, 0x30, 0x32, 0x30, 0x34, 0x34, 0x36, 0x30, 0x38, 0x38, 0x3a, 0x38, 0x3c, 0x3c, 0x3e, 0x00, 0x40, 0x40, 0x42, 0x40, 0x44, 0x44, 0x46, 0x40, 0x48, 0x48, 0x4a, 0x48, 0x4c, 0x4c, 0x4e, 0x40, 0x50, 0x50, 0x52, 0x50, 0x54, 0x54, 0x56, 0x50, 0x58, 0x58, 0x5a, 0x58, 0x5c, 0x5c, 0x5e, 0x40, 0x60, 0x60, 0x62, 0x60, 0x64, 0x64, 0x66, 0x60, 0x68, 0x68, 0x6a, 0x68, 0x6c, 0x6c, 0x6e, 0x60, 0x70, 0x70, 0x72, 0x70, 0x74, 0x74, 0x76, 0x70, 0x78, 0x78, 0x7a, 0x78, 0x7c, 0x7c, 0x7e, 0x00, 0x80, 0x80, 0x82, 0x80, 0x84, 0x84, 0x86, 0x80, 0x88, 0x88, 0x8a, 0x88, 0x8c, 0x8c, 0x8e, 0x80, 0x90, 0x90, 0x92, 0x90, 0x94, 0x94, 0x96, 0x90, 0x98, 0x98, 0x9a, 0x98, 0x9c, 0x9c, 0x9e, 0x80, 0xa0, 0xa0, 0xa2, 0xa0, 0xa4, 0xa4, 0xa6, 0xa0, 0xa8, 0xa8, 0xaa, 0xa8, 0xac, 0xac, 0xae, 0xa0, 0xb0, 0xb0, 0xb2, 0xb0, 0xb4, 0xb4, 0xb6, 0xb0, 0xb8, 0xb8, 0xba, 0xb8, 0xbc, 0xbc, 0xbe, 0x80, 0xc0, 0xc0, 0xc2, 0xc0, 0xc4, 0xc4, 0xc6, 0xc0, 0xc8, 0xc8, 0xca, 0xc8, 0xcc, 0xcc, 0xce, 0xc0, 0xd0, 0xd0, 0xd2, 0xd0, 0xd4, 0xd4, 0xd6, 0xd0, 0xd8, 0xd8, 0xda, 0xd8, 0xdc, 0xdc, 0xde, 0xc0, 0xe0, 0xe0, 0xe2, 0xe0, 0xe4, 0xe4, 0xe6, 0xe0, 0xe8, 0xe8, 0xea, 0xe8, 0xec, 0xec, 0xee, 0xe0, 0xf0, 0xf0, 0xf2, 0xf0, 0xf4, 0xf4, 0xf6, 0xf0, 0xf8, 0xf8, 0xfa, 0xf8, 0xfc, 0xfc, 0xfe, };
static uint8_t idx_arr[256] = { 0x00, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x07, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, };
    for (uint16_t s = 1; s < 256; ++s)
    {
        // uint8_t lsb  = s & -s;                  // isolate lsb
        // uint8_t prev = s ^ lsb;                 // smaller subset
        // uint8_t idx  = __builtin_ctz(lsb);      // position 0…7

        // prev_arr[s] = prev;
        // idx_arr[s]  = idx;
        uint8_t  prev = prev_arr[s] ; // store the previous value
        uint8_t  idx = idx_arr[s];    // store the index value

        
        output[s].low = output[prev].low ^ input[idx].low;                // exactly ONE XOR
        output[s].high = output[prev].high ^ input[idx].high;
    }
    // print our the prev and idx arrays so that I can easily just copy it into as a static array in C code. It should should be syntactically correct C code.
    printf("static uint8_t prev[256] = { ");
    for (int i = 0; i < 256; ++i) {
        printf("0x%02x, ", prev_arr[i]);
    }
    printf("};\n");
    printf("static uint8_t idx[256] = { ");
    for (int i = 0; i < 256; ++i) {
        printf("0x%02x, ", idx_arr[i]);
    }
    printf("};\n");
}

void pre_calculate_lookup_table_128bit(const uint32_t input[8], uint32_t output[256]) {

	static const uint8_t prev_arr[256] = { 0x00, 0x00, 0x00, 0x02, 0x00, 0x04, 0x04, 0x06, 0x00, 0x08, 0x08, 0x0a, 0x08, 0x0c, 0x0c, 0x0e, 0x00, 0x10, 0x10, 0x12, 0x10, 0x14, 0x14, 0x16, 0x10, 0x18, 0x18, 0x1a, 0x18, 0x1c, 0x1c, 0x1e, 0x00, 0x20, 0x20, 0x22, 0x20, 0x24, 0x24, 0x26, 0x20, 0x28, 0x28, 0x2a, 0x28, 0x2c, 0x2c, 0x2e, 0x20, 0x30, 0x30, 0x32, 0x30, 0x34, 0x34, 0x36, 0x30, 0x38, 0x38, 0x3a, 0x38, 0x3c, 0x3c, 0x3e, 0x00, 0x40, 0x40, 0x42, 0x40, 0x44, 0x44, 0x46, 0x40, 0x48, 0x48, 0x4a, 0x48, 0x4c, 0x4c, 0x4e, 0x40, 0x50, 0x50, 0x52, 0x50, 0x54, 0x54, 0x56, 0x50, 0x58, 0x58, 0x5a, 0x58, 0x5c, 0x5c, 0x5e, 0x40, 0x60, 0x60, 0x62, 0x60, 0x64, 0x64, 0x66, 0x60, 0x68, 0x68, 0x6a, 0x68, 0x6c, 0x6c, 0x6e, 0x60, 0x70, 0x70, 0x72, 0x70, 0x74, 0x74, 0x76, 0x70, 0x78, 0x78, 0x7a, 0x78, 0x7c, 0x7c, 0x7e, 0x00, 0x80, 0x80, 0x82, 0x80, 0x84, 0x84, 0x86, 0x80, 0x88, 0x88, 0x8a, 0x88, 0x8c, 0x8c, 0x8e, 0x80, 0x90, 0x90, 0x92, 0x90, 0x94, 0x94, 0x96, 0x90, 0x98, 0x98, 0x9a, 0x98, 0x9c, 0x9c, 0x9e, 0x80, 0xa0, 0xa0, 0xa2, 0xa0, 0xa4, 0xa4, 0xa6, 0xa0, 0xa8, 0xa8, 0xaa, 0xa8, 0xac, 0xac, 0xae, 0xa0, 0xb0, 0xb0, 0xb2, 0xb0, 0xb4, 0xb4, 0xb6, 0xb0, 0xb8, 0xb8, 0xba, 0xb8, 0xbc, 0xbc, 0xbe, 0x80, 0xc0, 0xc0, 0xc2, 0xc0, 0xc4, 0xc4, 0xc6, 0xc0, 0xc8, 0xc8, 0xca, 0xc8, 0xcc, 0xcc, 0xce, 0xc0, 0xd0, 0xd0, 0xd2, 0xd0, 0xd4, 0xd4, 0xd6, 0xd0, 0xd8, 0xd8, 0xda, 0xd8, 0xdc, 0xdc, 0xde, 0xc0, 0xe0, 0xe0, 0xe2, 0xe0, 0xe4, 0xe4, 0xe6, 0xe0, 0xe8, 0xe8, 0xea, 0xe8, 0xec, 0xec, 0xee, 0xe0, 0xf0, 0xf0, 0xf2, 0xf0, 0xf4, 0xf4, 0xf6, 0xf0, 0xf8, 0xf8, 0xfa, 0xf8, 0xfc, 0xfc, 0xfe, };
	static const uint8_t idx_arr[256] = { 0x00, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x07, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x06, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x05, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x04, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, 0x03, 0x00, 0x01, 0x00, 0x02, 0x00, 0x01, 0x00, };

	#pragma unroll
    for (uint16_t s = 1; s < 256; ++s)
    {
        // uint8_t lsb  = s & -s;                  // isolate lsb
        // uint8_t prev = s ^ lsb;                 // smaller subset
        // uint8_t idx  = __builtin_ctz(lsb);      // position 0…7
        uint8_t  prev = prev_arr[s] ; // store the previous value
        uint8_t  idx = idx_arr[s];    // store the index value

        output[s] = output[prev] ^ input[idx];                // exactly ONE XOR
 	}
}


void gpu_mul_via_matrix_bitsliced_four_russians_128bit( const uint128_t rows[128],
                                const uint32_t X[128],
                                      uint32_t Z[128] )
{
    // 1) Initialize the accumulator to zero
    //

    uint32_t lookup[256];
    memset (lookup, 0, sizeof(lookup));

    const int BYTE_SIZE = 8;
    const int OUTER_LOOP = 128 / BYTE_SIZE; // 16 * 8 = 128 bits

    for (int i = 0; i < OUTER_LOOP; ++i)
    {
        pre_calculate_lookup_table_128bit(&X[i * BYTE_SIZE], lookup);
        
        for (int j = 0; j < 128; ++j)
        {
            uint8_t* curren_lookup_rows_bytes = (uint8_t*)&rows[j];
            uint8_t idx = curren_lookup_rows_bytes[i];
            Z[j] ^= lookup[idx];
        }
    }
}



void mul_via_matrix_bitsliced_four_russians_method( const uint128_t rows[128],
                                const uint128_t X[128],
                                      uint128_t Z[128] )
{
    // 1) Initialize the accumulator to zero
    //
    uint128_t lookup[256];
    memset (lookup, 0, sizeof(lookup));

    const int BYTE_SIZE = 8;
    const int OUTER_LOOP = 128 / BYTE_SIZE; // 16 * 8 = 128 bits
    memset (Z, 0, sizeof(uint128_t) * 128);

    for (int i = 0; i < OUTER_LOOP; ++i)
    {
        pre_calculate_lookup_table(&X[i * BYTE_SIZE], lookup);
        
        for (int j = 0; j < 128; ++j)
        {
            uint8_t* curren_lookup_rows_bytes = &rows[j];
            uint8_t idx = curren_lookup_rows_bytes[i];
            Z[j].low ^= lookup[idx].low;
            Z[j].high ^= lookup[idx].high;
        }
    }
}

void mul_via_matrix_bitsliced_simple( const uint128_t rows[128],
                                const uint128_t X[128],
                                      uint128_t Z[128] )
{
    // 1) Initialize the accumulator to zero
    //
    uint128_t lookup[256];
    memset (lookup, 0, sizeof(lookup));

    const int BYTE_SIZE = 8;
    const int OUTER_LOOP = 128 / BYTE_SIZE; // 16 * 8 = 128 bits
    const int INNER_LOOP = 128 / BYTE_SIZE; // 32 * 8 = 256 bits
    // uint128_t Z [128 ];
    // memset (Z, 0, sizeof(Z));
    memset (Z, 0, sizeof(uint128_t) * 128);

    for (int i = 0; i < 16; ++i)
    {
        // pre_calculate_lookup_table(&X[i * BYTE_SIZE], lookup);
        for (int j = 0; j < 128; ++j)
        {
            uint8_t* current_lookup_rows_bytes = &rows[j];
            current_lookup_rows_bytes += i * BYTE_SIZE; 
            for (int k = 0; k < BYTE_SIZE; ++k)
            {
                // for every bit in current_lookup_rows_bytes[0]
                bool bit = (current_lookup_rows_bytes[0] >> k) & 1;
                if (bit)
                {
                    Z[j].low ^= X[ (i * BYTE_SIZE) + k].low;
                    Z[j].high ^= X[ (i * BYTE_SIZE) + k].high;
                }
            }
        }
    }
}



// -----------------------------------------------------------
//  Helpers to pack and unpack bitsliced representations:
// -----------------------------------------------------------
// static uint128_t make_u128( uint64_t hi, uint64_t lo )
// {
//     return (uint128_t){ lo, hi };
// }

// void pack_bitsliced_32bit( const uint32_t IN[128], uint32_t X[128] )
// {
//   for ( int j = 0; j < 128; ++j )
//   {
//     uint64_t lo = 0, hi = 0;

//     // for every input i, extract bit-j of IN[i] and scatter it
//     for ( int i = 0; i < 128; ++i )
//     {
//       unsigned b;
//       if ( j < 64 )
//         b = (IN[i].low  >> j) & 1;            // jth bit lives in .low
//       else
//         b = (IN[i].high >> (j - 64)) & 1;     // jth bit lives in .high

//       if ( i < 64 )
//         lo |= (uint64_t)b << i;              // bit-i of the slice.low
//       else
//         hi |= (uint64_t)b << (i - 64);       // bit-(i-64) of slice.high
//     }

//     X[j].low  = lo;
//     X[j].high = hi;
//   }
// }

#define BITSLICING_BITS_WIDTH 128
#define INTS_PER_UNBITSLICED_VALUE 4

 void transpose32(uint32_t A[32]) {
    int j, k;
    uint32_t m, t;

    m = 0x0000FFFF;
    for (j = 16; j != 0; j = j >> 1, m = m ^ (m << j)) {
        for (k = 0; k < 32; k = (k + j + 1) & ~j) {
            t = ((A[k] >> j) ^ (A[k + j])) & m;
            A[k] = A[k] ^ (t << j);
            A[k + j] = A[k + j] ^ (t);
        }
    }
}

void bitslice_transpose(uint32_t arr_bitsliced[BITSLICING_BITS_WIDTH]) {
		uint32_t tmp[BITSLICING_BITS_WIDTH];  // arr_bitsliced should also be of this size

		memcpy(tmp, arr_bitsliced, BITSLICING_BITS_WIDTH * sizeof(uint32_t));

		for (int i = 0; i < BITSLICING_BITS_WIDTH; ++i) {
			int idx_of_square_transpose = i % INTS_PER_UNBITSLICED_VALUE;
			int idx_within_square_transpose = i / INTS_PER_UNBITSLICED_VALUE;
			int unbitsliced_origin_of_chunk = 32 * idx_of_square_transpose + idx_within_square_transpose;
			arr_bitsliced[unbitsliced_origin_of_chunk] = tmp[i];
		}

		for (int square_chunk = 0; square_chunk < INTS_PER_UNBITSLICED_VALUE; ++square_chunk) {
			transpose32(arr_bitsliced + 32 * square_chunk);
		}
	}

 void bitslice_untranspose(uint32_t arr_bitsliced[BITSLICING_BITS_WIDTH]) {
    uint32_t tmp[BITSLICING_BITS_WIDTH];  // arr_bitsliced should also be of this size

    memcpy(tmp, arr_bitsliced, BITSLICING_BITS_WIDTH * sizeof(uint32_t));

    for (int square_chunk = 0; square_chunk < INTS_PER_UNBITSLICED_VALUE; ++square_chunk) {
        transpose32(tmp + 32 * square_chunk);
    }

    for (int i = 0; i < BITSLICING_BITS_WIDTH; ++i) {
        int chunk_of_number_idx = i / 32;
        int number_idx = i % 32;
        int unbitsliced_destination_of_chunk = INTS_PER_UNBITSLICED_VALUE * number_idx + chunk_of_number_idx;
        arr_bitsliced[unbitsliced_destination_of_chunk] = tmp[i];
    }
}


// ─── PACK: IN[0..127] → bitslices X[0..127] ─────────────────────────────
void pack_bitsliced( const uint128_t IN[128], uint128_t X[128] )
{
  for ( int j = 0; j < 128; ++j )
  {
    uint64_t lo = 0, hi = 0;

    // for every input i, extract bit-j of IN[i] and scatter it
    for ( int i = 0; i < 128; ++i )
    {
      unsigned b;
      if ( j < 64 )
        b = (IN[i].low  >> j) & 1;            // jth bit lives in .low
      else
        b = (IN[i].high >> (j - 64)) & 1;     // jth bit lives in .high

      if ( i < 64 )
        lo |= (uint64_t)b << i;              // bit-i of the slice.low
      else
        hi |= (uint64_t)b << (i - 64);       // bit-(i-64) of slice.high
    }

    X[j].low  = lo;
    X[j].high = hi;
  }
}

// ─── UNPACK: bitslices Z[0..127] → OUT[0..127] ─────────────────────────
void unpack_bitsliced( const uint128_t Z[128], uint128_t OUT[128] )
{
  for ( int i = 0; i < 128; ++i )
  {
    uint64_t lo = 0, hi = 0;

    // for each slice j, pull out bit-i and deposit into OUT[i]
    for ( int j = 0; j < 128; ++j )
    {
      unsigned b;
      if ( i < 64 )
        b = (Z[j].low  >> i) & 1;            // slice.low bit-i
      else
        b = (Z[j].high >> (i - 64)) & 1;     // slice.high bit-(i-64)

      if ( j < 64 )
        lo |= (uint64_t)b << j;              // jth bit of OUT[i].low
      else
        hi |= (uint64_t)b << (j - 64);       // (j-64)th bit of OUT[i].high
    }

    OUT[i].low  = lo;
    OUT[i].high = hi;
  }
}

void mul_via_matrix_bitsliced_four_russians_method_cols(
  const uint128_t cols[128],   // columns of the 128×128 bit-matrix
  const uint128_t X   [128],   // bitsliced input
        uint128_t Z   [128] )  // bitsliced output (zeroed by this call)
{
  const int BYTE_SIZE  = 8;
  const int OUTER_LOOP = 128 / BYTE_SIZE;  // =16

  // zero the outputs
//   for ( int j = 0; j < 128; ++j )
//     Z[j] = {0,0};

  // for each 8-bit “chunk” of X
  for ( int i = 0; i < OUTER_LOOP; ++i )
  {
    // build the 256-entry lookup table for X[i*8..i*8+7]
    uint128_t lookup[256];
    memset( lookup, 0, sizeof lookup );
    pre_calculate_lookup_table( &X[i * BYTE_SIZE], lookup );

    // for each output row j=0..127, gather its 8 bits from cols[ i*8+t ], t=0..7
    for ( int j = 0; j < 128; ++j )
    {
      uint8_t idx = 0;
      // build the byte idx whose bit-t is the j’th bit of column (i*8 + t)
      for ( int t = 0; t < BYTE_SIZE; ++t )
      {
        int c = i * BYTE_SIZE + t;
        bool bit = ( j < 64
                   ? (cols[c].low  >> j) & 1u
                   : (cols[c].high >> (j - 64)) & 1u );
        idx |= bit << t;
      }
      // XOR in that precomputed partial product
      Z[j].low  ^= lookup[idx].low;
      Z[j].high ^= lookup[idx].high;
    }
  }
}

static inline uint64_t rand64(void) {
  return ((uint64_t)rand() << 33) ^ ((uint64_t)rand() << 17) ^ (uint64_t)rand();
}

static inline __uint128_t rand_u128(void) {
  __uint128_t lo = (__uint128_t)rand64();
  __uint128_t hi = (__uint128_t)rand64();
  return (hi << 64) | lo;
}

static inline void u128_to_bytes_le(__uint128_t x, uint8_t b[16]) {
  // Little-endian byte order: b[0] = bits 7..0, ..., b[15] = bits 127..120
  for (int i = 0; i < 16; ++i) b[i] = (uint8_t)(x >> (8*i));
}

static inline __uint128_t make_u128(uint64_t hi, uint64_t lo) {
  return ((__uint128_t)hi << 64) | (__uint128_t)lo;
}

static inline void split_u128(__uint128_t x, uint64_t *hi, uint64_t *lo) {
  *lo = (uint64_t)x;
  *hi = (uint64_t)(x >> 64);
}

static inline uint64x2_t load_u128_as_u64x2(const __uint128_t *src) {
  uint64_t lo = (uint64_t)(*src);
  uint64_t hi = (uint64_t)((*src) >> 64);
  uint64x2_t v = vdupq_n_u64(0);
  v = vsetq_lane_u64(lo, v, 0);
  v = vsetq_lane_u64(hi, v, 1);
  return v;
}

static void fill_random_cols(__uint128_t cols[128]) {
  for (int j = 0; j < 128; ++j) cols[j] = rand_u128();
}

static double ms_since(struct timespec a, struct timespec b) {
  return (b.tv_sec - a.tv_sec) * 1e3 + (b.tv_nsec - a.tv_nsec) / 1e6;
}



void build_byte_tables_from_cols(const __uint128_t cols[128], __uint128_t T[16][256])
{
  for (int pos = 0; pos < 16; ++pos) {
    T[pos][0] = 0;
    for (int v = 1; v < 256; ++v) {
      int lsb = v & -v;
      int bit = __builtin_ctz(lsb);               // 0..7
      T[pos][v] = T[pos][v ^ lsb] ^ cols[pos*8 + bit];
    }
  }
}

static __uint128_t mul_via_matrix_cols_scalar_simple(const __uint128_t cols[128],
                                              __uint128_t X) {
  __uint128_t y = 0;
  for (int j = 0; j < 128; ++j) {
    if ((X >> j) & 1) y ^= cols[j];
  }
  return y;
}


__uint128_t mul_const_neon_bytes(const __uint128_t T[16][256], __uint128_t X)
{
  uint8_t xb[16];
  u128_to_bytes_le(X, xb); // may be we can remove this?

  uint64x2_t acc = vdupq_n_u64(0);
  for (int pos = 0; pos < 16; ++pos) {
    uint64x2_t t = vld1q_u64((const uint64_t*)&T[pos][ xb[pos] ]);
    acc = veorq_u64(acc, t);
  }

  uint64_t out64[2];
  vst1q_u64(out64, acc);
  __uint128_t y = (__uint128_t)out64[0] | ((__uint128_t)out64[1] << 64);
  return y;
}

// -----------------------------------------------------------
//  3) Test end‐to‐end on 128 random inputs:
//     compare each OUT[i] against gf128_mul(C, IN[i]).
// -----------------------------------------------------------
void test_bitsliced( uint128_t C )
{
    uint32_t IN[128];
    for ( int i = 0; i < 128; ++i )
    {
        IN[i] = ((uint64_t)rand() << 32) ^ ((uint64_t)rand() << 16) ^ rand();
    }

    struct timespec t0, t1;
    clock_gettime( CLOCK_MONOTONIC, &t0 );

    //  3a) Build cols[]
    uint128_t cols[128];
    build_matrix128( C, cols );
    //  3b) Transpose cols[] into rows[]
    uint128_t rows[128];
    transpose128( cols, rows );

    // // lets transpose back and check if it is the same
    // uint128_t cols_check[128];
    // transpose128( rows, cols_check );
    // for ( int i = 0; i < 128; ++i )
    // {
    //     assert( cols[i].low  == cols_check[i].low );
    //     assert( cols[i].high == cols_check[i].high );
    // }

    //  3b) Make 128 random test‐vectors IN[i]

    //  3c) Pack into bitsliced X[0..127]
    uint32_t X[128], Z[128], OUT[128], IN_UNBITSLICED[128];
    // memset( X, 0, sizeof(X) );
    memset( Z, 0, sizeof(Z) );
    memset( OUT, 0, sizeof(OUT) );
    // memcpy ( OUT, IN, sizeof(OUT) );
    memcpy ( IN_UNBITSLICED, IN, sizeof(IN_UNBITSLICED) );
    bitslice_transpose( IN ); // Transpose the input to bitsliced format

    // pack_bitsliced( IN, X );

    //lets do a sanity check on the packed bitsliced X and unbitslice to see if it matches IN
    // bitslice_untranspose( IN );
    // for ( int i = 0; i < 128; ++i )
    // {
    //     assert( IN[i]  == OUT[i] );
    // }

    //  3d) Multiply all 128 in parallel
    // mul_via_matrix_bitsliced_four_russians_method( rows, X, Z );
    gpu_mul_via_matrix_bitsliced_four_russians_128bit( rows, IN, OUT ); 
    // mul_via_matrix_bitsliced_four_russians_method_cols( cols, X, Z );
    //  3e) Unpack back to OUT[i]
    bitslice_untranspose( OUT );
    clock_gettime( CLOCK_MONOTONIC, &t1 );
    double elapsed_ms = (t1.tv_sec  - t0.tv_sec ) * 1e3
                      + (t1.tv_nsec - t0.tv_nsec) / 1e6;
    printf( "[benchmark] multiply+unpack took %.3f ms\n", elapsed_ms );

    //  3f) Check each one
    uint128_t *original_unbitsliced = (uint128_t *)IN_UNBITSLICED;
    uint128_t *actual128arr = (uint128_t *)OUT;
    for ( int i = 0; i < 32; ++i )
    {
        uint128_t expect = binmul128( C, original_unbitsliced[i], 128 );
        // printf("IN[%d] = %016llx%016llx, OUT[%d] = %016llx%016llx, expect = %016llx%016llx\n",
        //         i, IN[i].high, IN[i].low, i, OUT[i].high, OUT[i].low,
        //         expect.high, expect.low);
        // // pritn index
        // printf("IN[%d] = %016llx%016llx, OUT[%d] = %016llx%016llx\n",
        //         i, IN[i].high, IN[i].low, i, OUT[i].high, OUT[i].low);

        assert( expect.low  == actual128arr[i].low  );
        assert( expect.high == actual128arr[i].high );
    }
    printf( "✅ bitsliced test passed for constant C\n" );
}

void test_128(){
    // Example 128-bit numbers split into high and low 64-bit parts
    uint128_t v1 = {14143994781733811029ULL, 669260594276348690ULL};  // Example 128-bit number
    v1.low = 15143994781733811029ULL;
    v1.high = 669260594276348690ULL;
    uint128_t v2 = {15875069739565888632ULL, 5354084802999887300ULL};  // Example 128-bit number
    v2.low = 15875069739565888632ULL;
    v2.high = 5354084802999887303ULL;
    // uint128_t result = binmul128(&v1, &v2, 64);
    
        // x[i].low = 14143994781733811029ULL;
        // x[i].high = 669260594276348690ULL;
        // y[i].low = 15875069739565888632ULL;
        // y[i].high = 5354084802999887303ULL;

    // 64 bit numbers
    // uint64_t v1 = 14143994781733811022ULL;  // Example 64-bit number
    // uint64_t v2 = 15875069739565888632ULL;  // Example 64-bit number

    // printf("Input v1 %016llx \n: ", v1);
    // printf("Input v2 %016llx \n: ", v2);
    clock_t start = clock();
    uint128_t result = {0, 0};
    for (int i = 0; i < 64; i++) {
        // uint64_t v1 = i;  // Example 64-bit number
        // uint64_t v2 = i+1;  // Example 64-bit number
        v2.high += 1;
        result = binmul128(v1, v2, 128);
        printf("Combined Result: %016llx%016llx \n", result.high, result.low);
        // printf("", result.high, result.low);
    }
    clock_t end = clock();
    printf("Time taken: %f seconds\n", ((double)(end - start)) / CLOCKS_PER_SEC);

    // printf("C Result: ");
    // printf("%016llx\n ", result);
    // uint128_t combined_result = {result.low, result.high};
    printf("\n");
    // printf("Calls: %d\n", calls);
}

void test_16(){

    // 16 bit numbers
    uint64_t v1 = 5622 ;  // Example 16-bit number
    uint64_t v2 = 7982;  // Example 16-bit number

    // printf("Input v1 %016llx \n: ", v1);
    // printf("Input v2 %016llx \n: ", v2);
    clock_t start = clock();
    uint64_t result = 0;
    // for (int i = 0; i < 128; i++) {
    //     // uint16_t v1 = i;  // Example 16-bit number
    //     // uint16_t v2 = i+1;  // Example 16-bit number
        
    // }
    result = binmul64(v1, v2, 16, false);
    clock_t end = clock();
    printf("Time taken: %f seconds\n", ((double)(end - start)) / CLOCKS_PER_SEC);

    printf("C Result: ");
    printf("%04x\n ", result);
    // uint128_t combined_result = {result.low, result.high};
    // printf("Combined Result: ");
    // printf("%016llx%016llx", combined_result.high, combined_result.low);
    printf("\n");
    // printf("Calls: %d\n", calls);


}




int test_normal_nonbitsliced_matrix_constant_mul(void) {
  srand(42);

  __uint128_t cols[128];
  fill_random_cols(cols);

  // Build the NEON lookup tables
  __uint128_t T[16][256];
  build_byte_tables_from_cols(cols, T);

  // Correctness: compare on many random vectors
  const int N = 5000;
  for (int i = 0; i < N; ++i) {
    __uint128_t X = rand_u128();
    __uint128_t y_ref  = mul_via_matrix_cols_scalar_simple(cols, X);
    __uint128_t y_neon = mul_const_neon_bytes(T, X);
    if (y_ref != y_neon) {
      uint64_t ref_hi, ref_lo, neon_hi, neon_lo, x_hi, x_lo;
      split_u128(y_ref,  &ref_hi,  &ref_lo);
      split_u128(y_neon, &neon_hi, &neon_lo);
      split_u128(X,      &x_hi,    &x_lo);
      printf("Mismatch!\nX=%016llx%016llx\nref =%016llx%016llx\nyNEON=%016llx%016llx\n",
             (unsigned long long)x_hi,   (unsigned long long)x_lo,
             (unsigned long long)ref_hi, (unsigned long long)ref_lo,
             (unsigned long long)neon_hi,(unsigned long long)neon_lo);
      return 1;
    }
  }
  printf("✅ Correctness check passed on %d random vectors.\n", N);

  // Benchmark: scalar vs NEON over M trials
  const int M = 200000;
  struct timespec t0, t1;

  __uint128_t sink = 0; // prevent dead-code elimination

  clock_gettime(CLOCK_MONOTONIC, &t0);
  for (int i = 0; i < M; ++i) {
    __uint128_t X = rand_u128();
    sink ^= mul_via_matrix_cols_scalar_simple(cols, X);
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);
  double ms_scalar = ms_since(t0, t1);

  clock_gettime(CLOCK_MONOTONIC, &t0);
  for (int i = 0; i < M; ++i) {
    __uint128_t X = rand_u128();
    sink ^= mul_const_neon_bytes(T, X);
  }
  clock_gettime(CLOCK_MONOTONIC, &t1);
  double ms_neon = ms_since(t0, t1);

  // Print something involving sink so compiler keeps the loops
  uint64_t s_hi, s_lo; split_u128(sink, &s_hi, &s_lo);
  printf("sink = %016llx%016llx\n", (unsigned long long)s_hi, (unsigned long long)s_lo);

  printf("Scalar baseline:  %8.3f ms total  (%.3f ns/op)\n",
         ms_scalar, (ms_scalar*1e6)/M);
  printf("NEON lookup:      %8.3f ms total  (%.3f ns/op)\n",
         ms_neon,   (ms_neon*1e6)/M);

  return 0;
}
// Main function to test the implementation
int main() {
//    test_16();
//    test_128();


    // uint128_t C = {
    //     .low  = 0xFEDCBA9876543210ULL,
    //     .high = 0x0123456789ABCDEFULL
    // };

    // test_matrix128( C );



    //  Example constant C:
    uint128_t C = { 0xFEDCBA9876543210ULL, 0x0123456789ABCDEFULL };

    //  Run the bitsliced self‐test:
//     srand(42);
//    test_bitsliced( C );
    // Run the nonbitsliced self‐test:
    test_normal_nonbitsliced_matrix_constant_mul();

    return 0;
}
