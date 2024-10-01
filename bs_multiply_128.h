#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#ifndef _BS_MULTIPLY_128_H_
#define _BS_MULTIPLY_128_H_

// typedef struct {
//     uint64_t low;
//     uint64_t high;
// } uint128_t;

// // Z is the output and assumed to be ord_length 
// void bs_multiply_128(word_t x[128], word_t y[128], word_t *z);

void transpose_mul(uint64_t *x, uint64_t *y, uint64_t* z);

#endif
