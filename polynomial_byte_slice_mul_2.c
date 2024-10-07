#include <stdio.h>
#include <assert.h>
#include <stdint.h>
#include "bs.h"
#include "polynomial_byte_slice_mul_2.h"

#define NUM_INPUTS 8        // Number of 16-bit numbers
#define BYTES_IN_16BIT 2   // 2 bytes in a 16-bit number
#define SLICED_OUTPUTS 2

void byte_slice_mul_16(uint64_t x[2], uint64_t y[2], uint64_t* z);

//////////////////////////////////////////////////////
// Byte Slicing: Convert 8 x 16-bit inputs to 2 rows of 64 bits
//////////////////////////////////////////////////////
void byte_slice_16(uint16_t input[NUM_INPUTS], uint64_t output[SLICED_OUTPUTS]) {
    // Cast the input to uint8_t* for easier access to bytes
    uint8_t* input_bytes = (uint8_t*) input;
    uint8_t* output_bytes = (uint8_t*) output;

    // Loop over each byte index (0 to 15 for each 16-bit number)
    for (int byte_index = 0; byte_index < BYTES_IN_16BIT; byte_index++) {
        // Loop over each of the original 16-bit numbers (8 total inputs)
        for (int i = 0; i < NUM_INPUTS; i++) {
            // Place each byte from the input array into the corresponding 64-bit output row
            // the print the index being accesses
            int index_to = byte_index * NUM_INPUTS + i;
            int index_from = i * BYTES_IN_16BIT + byte_index;
            printf("index to: %d, from: %d\n", index_to, index_from);
            output_bytes[byte_index * NUM_INPUTS + i] = input_bytes[i * BYTES_IN_16BIT + byte_index];
        }
    }
}

//////////////////////////////////////////////////////
// Un-byte Slicing: Convert 16 rows of 64 bits back to 8 x 16-bit inputs
//////////////////////////////////////////////////////
void un_byte_slice_16(uint64_t input[SLICED_OUTPUTS], uint16_t output[NUM_INPUTS]) {
    // Cast the input to uint8_t* for easier access to bytes
    uint8_t* input_bytes = (uint8_t*) input;
    uint8_t* output_bytes = (uint8_t*) output;

    // Loop over each byte index (0 to 15 for each 16-bit number)
    for (int byte_index = 0; byte_index < BYTES_IN_16BIT; byte_index++) {
        // Loop over each of the 8 16-bit numbers
        for (int i = 0; i < NUM_INPUTS; i++) {
            // Reconstruct the original bytes into the 16-bit numbers
            output_bytes[i * BYTES_IN_16BIT + byte_index] = input_bytes[byte_index * NUM_INPUTS + i];
        }
    }
}

void byte_slice_transpose_mul_16(uint16_t x[8], uint16_t y[8], uint16_t* z){
    uint64_t x_transposed[2], test_x_transposed[2];
    uint64_t y_transposed[2];
    uint64_t z_transposed[2];
    
    byte_slice_16(x, x_transposed);
    // un_byte_slice_16(x_transposed, test_x_transposed);
    // uint64_t *test_x = x;
    // // test both are equal
    // for (int i = 0; i < 2; i++) {
    //     assert(test_x[i] == test_x_transposed[i]);
    // }
    byte_slice_16(y, y_transposed);
    byte_slice_mul_16(x_transposed, y_transposed, z_transposed);
    un_byte_slice_16(z_transposed, z);

}


void byte_slice_mul_16(uint64_t x[2], uint64_t y[2], uint64_t* z)
{
  uint64_t n5 , n6 , n7 , n8 , n9 , n10 , n11 , n12, n6_high_r1r2 ;
  //n5 = x[0] & y[0] ;
  multiply_128b_using_log_table( &x[0], &y[0],&n5);

  //n6 = x[1] & y[1] ;
  multiply_128b_using_log_table( &x[1],&y[1], &n6);

  n10 = n6 ^ n5 ;
  multiply_constant_128b_using_table(&n6,&n6_high_r1r2 ); 


  n7 = x[1] ^ x[0] ;
  n8 = y[1] ^ y[0] ;
  //n9 = n7 & n8 ;
  multiply_128b_using_log_table( &n7,&n8, &n9);
  
  n11 = n10 ^ n9 ;
  n12 = n11 ^ n6_high_r1r2 ;
  printf("L1L2: %016llx, R1R2: %016llx, R1R2_high: %016llx, Z3: %016llx\n", n5, n6, n6_high_r1r2, n9);
  uint64_t temp_result = n11 ^ n6_high_r1r2 ;
    printf("temp_result: %016llx\n", temp_result);
  z[0] = n10 ;
  z[1] = n12  ;
  // lets print the individual components making up n12
//   printf("L1L2: %016llx, R1R2: %016llx, L1L2: %016llx \n", n5,n6,   );
  printf("z[0]: %016llx, z[1]: %016llx\n", z[0], z[1]);
}



// module top( x0 , x1 , x2 , x3 , y0 , y1 );
//   input x0 , x1 , x2 , x3 ;
//   output y0 , y1 ;
// endmodule
