#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>

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
// Main function to test the implementation
int main() {
//    test_16();
    test_128();
    return 0;
}
