/*
 * richardson_test.cpp
 *
 * Unit test for richardson.cpp.
 */

#include "richardson.h"

#include <gmp.h>

using namespace richardson;


/* Sample SequenceFuncs. */
void reciprocal(index_t index, mpf_t result) {
    // Calculate 1/index
    mpf_set_d(result, index);
    mpf_ui_div(result, 1, result);
}

void onePlusReciprocal(index_t index, mpf_t result) {
    // Calculate 1 + 1/index
    mpf_set_d(result, index);
    mpf_ui_div(result, 1, result);
    mpf_add_ui(result, result, 1);
}

void approxPi(index_t index, mpf_t result) {
    // Approximate PI according to the first [index] terms of the sequence
    //
    // PI / 4 = 1 - 1/3 + 1/5 - 1/7 + ...
    mpf_set_d(result, 0);
    mpf_t subtotal;
    mpf_init2(subtotal, mpf_get_prec(result));
    for (index_t i = 0; i < index; i++) {
        mpf_set_d(subtotal, 2 * i + 1);
        mpf_ui_div(subtotal, 1, subtotal);
        if (i % 2 == 0)
            mpf_add(result, result, subtotal);
        else
            mpf_sub(result, result, subtotal);
    }
    mpf_mul_ui(result, result, 4);
}

// Test Richardson extrapolating the function f with the specified number
// of samples at the specified precision and the specified start index, and
// ensure that it matches the given expected value to within 2^-diff_precision.
#define TEST(num_samples, start_index, f, expected, precision, diff_precision) \
{ \
    mpf_t error; \
    mpf_init_set_d(error, 0.5); \
    mpf_pow_ui(error, error, diff_precision); \
 \
    mpf_t ans; \
    mpf_init2(ans, precision); \
    extrapolate(num_samples, start_index, &f, ans); \
 \
    mpf_t diff; \
    mpf_init2(diff, precision); \
    mpf_sub(diff, ans, expected); \
    mpf_abs(diff, diff); \
    if (mpf_cmp(diff, error) > 0) { \
        printf("%s test failed. Expected ", #f); \
        mpf_out_str(NULL, 10, precision, expected); \
        printf(", but got "); \
        mpf_out_str(NULL, 10, precision, ans); \
        printf("\n"); \
        numFailures++; \
    } \
}

int main() {
    // Initialization code
    mpf_t ZERO;
    mpf_init_set_d(ZERO, 0.0);
    mpf_t ONE;
    mpf_init_set_d(ONE, 1.0);
    mpf_t PI;
    mpf_init2(PI, 128);
    mpf_set_str(PI, "3.1415926535897932384626433832795028841971", 10);

    int numFailures = 0;

    TEST(10, 1, reciprocal, ZERO, 64, 10);
    TEST(10, 1, onePlusReciprocal, ONE, 64, 10);
    TEST(16, 1, approxPi, PI, 128, 80);

    if (numFailures) {
        printf("Found %d failures.\n", numFailures);
        return 1;
    } else {
        printf("Success!\n");
        return 0;
    }
}

