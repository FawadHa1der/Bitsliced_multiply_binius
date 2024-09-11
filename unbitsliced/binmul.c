#include <stdio.h>
#include <stdint.h>
#include <stdbool.h>
#include <time.h>



uint64_t binmul64(uint64_t v1, uint64_t v2, uint32_t length, bool is_constant) {


    if (v1 < 2 || v2 < 2)  {
        uint64_t result = v1 * v2;
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

    uint64_t result = (uint64_t) L1L2 ^ R1R2 ^ ((Z3 ^ L1L2 ^ R1R2 ^ R1R2_high) << halflen);
    return result;
}

// Main function to test the implementation
int main() {
    // Example 128-bit numbers split into high and low 64-bit parts
    // uint128_t v1 = {14143994781733811022ULL, 669260594276348691ULL};  // Example 128-bit number
    // uint128_t v2 = {15875069739565888632ULL, 5354084802999887300ULL};  // Example 128-bit number
    // uint128_t result = binmul128(&v1, &v2, 64);
    

    // 64 bit numbers
    // uint64_t v1 = 14143994781733811022ULL;  // Example 64-bit number
    // uint64_t v2 = 15875069739565888632ULL;  // Example 64-bit number

    // printf("Input v1 %016llx \n: ", v1);
    // printf("Input v2 %016llx \n: ", v2);
    clock_t start = clock();
    uint64_t result = 0;
    for (int i = 0; i < 128; i++) {
        uint64_t v1 = i;  // Example 64-bit number
        uint64_t v2 = i+1;  // Example 64-bit number
        result = binmul64(v1, v2, 64, false);
    }
    clock_t end = clock();
    printf("Time taken: %f seconds\n", ((double)(end - start)) / CLOCKS_PER_SEC);

    // printf("C Result: ");
    // printf("%016llx\n ", result);
    // uint128_t combined_result = {result.low, result.high};
    // printf("Combined Result: ");
    // printf("%016llx%016llx", combined_result.high, combined_result.low);
    printf("\n");
    // printf("Calls: %d\n", calls);
    return 0;
}
