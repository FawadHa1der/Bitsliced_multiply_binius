#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#ifndef _BS_MULTIPLY_H_
#define _BS_MULTIPLY_H_

// Z is the output and assumed to be ord_length 
void bs_multiply_64(word_t x[64], word_t y[64], word_t *z);



#endif
