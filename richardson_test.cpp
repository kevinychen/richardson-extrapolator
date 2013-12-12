/*
 * richardson_test.cpp
 *
 * Unit test for richardson.cpp.
 */

#include "richardson.h"

#include <gmp.h>

using namespace richardson;


mpf_t ONE;

/* Sample SequenceFuncs. */
void reciprocal(index_t index, mpf_t result) {
    // Calculate 1/index
    mpf_init_set_d(result, index);
    mpf_div(result, ONE, result);
}

void approxPi(index_t index, mpf_t result) {
    // Approximate PI according to the first [index] terms of the sequence
    //
    // PI / 4 = 1 - 1/3 + 1/5 - 1/7 + ...
    mpf_set_d(result, 0);
    mpf_t subtotal;
    mpf_init_set_d(subtotal, 0);
    for (index_t i = 0; i < index; i++) {
        mpf_set_d(subtotal, (i % 2 == 0 ? 1 : -1) * (2 * i + 1));
        mpf_div(subtotal, ONE, subtotal);
        mpf_add(result, result, subtotal);
    }
}


int main() {
    // Initialization code
    mpf_init_set_d(ONE, 1.0);

    // Tests
    {
        const index_t num_samples = 10;
        const index_t start_index = 1;
        mpf_t ans;
        mpf_init(ans);
        extrapolate(num_samples, start_index, &reciprocal, ans);
        mpf_out_str(NULL, 10, 10, ans);
        // TODO test equality with actual answer.
    }
}

