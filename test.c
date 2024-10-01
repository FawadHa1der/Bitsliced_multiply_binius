#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#include "bs_multiply_64.h"
#include "bs_multiply_128.h"
#include "polynomial_byte_sliced_mul_16.h"
#include <time.h>
#include <mach/mach_time.h>




int main() {
    uint128_t x[64] = {0};
    uint128_t y[64] = {0};
    uint128_t z[64] = {0};

    for (int i = 0; i < 64; i++) {
        x[i].low = 14143994781733811022ULL;
        x[i].high = 669260594276348691ULL;
        y[i].low = 15875069739565888632ULL;
        y[i].high = 5354084802999887300ULL;
    }
    uint64_t start = mach_absolute_time();
    byte_slice_transpose_mul_128(x, y, z);
    // transpose_mul(x, y, z);
    // transpose_mul(x, y, z);
    // transpose_mul(x, y, z);
    // bs_transpose(x, 1);
    // bs_transpose(y, 1);
    // bs_multiply_64(x, y, z);
    // bs_multiply_64(&x[64], &y[64], &z[64]);
    // bs_transpose_rev(z, 1);
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
        printf("z[%d] = %016llx%016llx\n", i, z[i].high, z[i].low);
    }

    return 0;
}