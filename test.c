#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#include "bs_multiply_64.h"
#include <time.h>




int main() {
    word_t x[128] = {0};
    word_t y[128] = {0};
    word_t z[128] = {0};
    for (int i = 0; i < 128; i++) {
        x[i] = i;
        y[i] = i+3;
    }
    clock_t start = clock();
    bs_transpose(x, 1);
    bs_transpose(y, 1);
    bs_multiply_64(x, y, z);
    bs_multiply_64(&x[64], &y[64], &z[64]);
    bs_transpose_rev(z, 1);
    clock_t end = clock();

    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;
    printf("Time taken: %f seconds\n", time_taken);


    // for (int i = 0; i < 64; i++) {
    //     printf("%016llx ", z[i]);
    // }

    return 0;
}


