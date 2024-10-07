#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#include "bs_multiply_64.h"
#include "bs_multiply_128.h"
#include "polynomial_byte_slice_mul_2.h" // for 16 bit
#include "polynomial_byte_sliced_mul_16.h" // for 128 bit
#include <time.h>
#include <mach/mach_time.h>




int test_128_byte_slice(){
    uint128_t x[64] = {15143994781733811029ULL, 669260594276348690ULL};
    uint128_t y[64] = {15875069739565888632ULL, 5354084802999887303ULL};
    uint128_t z[64] = {0};
// 14143994781733811029ULL, 669260594276348690ULL
    for (int i = 0; i < 64; i++) {
        x[i].low = 15143994781733811029ULL;
        x[i].high = 669260594276348690ULL;
        y[i].low = 15875069739565888632ULL;
        y[i].high = 5354084802999887303ULL + i + 1;
    }

    uint64_t start = mach_absolute_time();
    byte_slice_transpose_mul_128(x, y, z);
    byte_slice_transpose_mul_128(&x[16], &y[16], &z[16]);
    byte_slice_transpose_mul_128(&x[32], &y[32], &z[32]);
    byte_slice_transpose_mul_128(&x[48], &y[48], &z[48]);
    uint64_t end = mach_absolute_time();

    // Get the timebase info to convert to nanoseconds
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);

    // Calculate the elapsed time in nanoseconds
    uint64_t elapsed = (end - start) * timebase.numer / timebase.denom;

    // Convert nanoseconds to microseconds for display
    double time_taken = (double)elapsed / 1e3;
    printf("Time taken: %f microseconds (%llu nanoseconds)\n", time_taken, elapsed);
    for (int i = 0; i < 16; i++) {
        printf("z[%d] = %016llx%016llx\n", i, z[i].high, z[i].low);
    }

}

int test_16_byte_slice(){
    uint16_t x[8] = {0};
    uint16_t y[8] = {0};
    uint16_t z[8] = {0};
// 14143994781733811029ULL, 669260594276348690ULL
    for (int i = 0; i < 8; i++) {
        x[i] = 5622;
        y[i] = 7982;
    }

    uint64_t start = mach_absolute_time();
    byte_slice_transpose_mul_16(x, y, z);
    uint64_t end = mach_absolute_time();

    // Get the timebase info to convert to nanoseconds
    mach_timebase_info_data_t timebase;
    mach_timebase_info(&timebase);

    // Calculate the elapsed time in nanoseconds
    uint64_t elapsed = (end - start) * timebase.numer / timebase.denom;

    // Convert nanoseconds to microseconds for display
    double time_taken = (double)elapsed / 1e3;
    printf("Time taken: %f microseconds (%llu nanoseconds)\n", time_taken, elapsed);
    for (int i = 0; i < 8; i++) {
        printf("z[%d] = %04x\n", i, z[i]);
    }

}


int main() {
    //test_16_byte_slice();
    test_128_byte_slice();
    return 0;
}