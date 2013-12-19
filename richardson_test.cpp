/*
 * richardson_test.cpp
 *
 * Unit test for richardson.cpp.
 */

#include "richardson.h"

#include <gmp.h>

using namespace richardson;


mpf_t ZERO;
mpf_t ONE;
mpf_t PI;

/* Sample SequenceFuncs. */
void reciprocal(index_t index, mpf_t result) {
    // Calculate 1/index
    mpf_init_set_d(result, index);
    mpf_div(result, ONE, result);
}

void onePlusReciprocal(index_t index, mpf_t result) {
    // Calculate 1 + 1/index
    mpf_init_set_d(result, index);
    mpf_div(result, ONE, result);
    mpf_add(result, ONE, result);
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

// Test Richardson extrapolating the function f with the specified number
// of samples and the specified start index, and ensure that it matches
// the given expected value to num_precision decimal places.
#define TEST(num_samples, start_index, f, expected, num_precision) \
{ \
    mpf_t error; \
    mpf_init_set_d(error, 0.1); \
    mpf_pow_ui(error, error, num_precision); \
 \
    mpf_t ans; \
    mpf_init(ans); \
    extrapolate(num_samples, start_index, &f, ans); \
 \
    mpf_t diff; \
    mpf_init(diff); \
    mpf_sub(diff, ans, expected); \
    mpf_abs(diff, diff); \
    if (mpf_cmp(diff, error) > 0) { \
        printf("%s test failed. Expected ", #f); \
        mpf_out_str(NULL, 10, num_precision, expected); \
        printf(", but got "); \
        mpf_out_str(NULL, 10, num_precision, ans); \
        printf("\n"); \
        numFailures++; \
    } \
}

int main() {
    // Initialization code
    mpf_init_set_d(ZERO, 0.0);
    mpf_init_set_d(ONE, 1.0);
    mpf_init_set_d(PI, 3.1415926535897932384626433);

    int numFailures = 0;

    TEST(10, 1, reciprocal, ZERO, 10);
    TEST(10, 1, onePlusReciprocal, ONE, 10);
    TEST(20, 1, approxPi, PI, 18);

    if (numFailures) {
        printf("Found %d failures.\n", numFailures);
        return 1;
    } else {
        printf("Success!\n");
        return 0;
    }
}

