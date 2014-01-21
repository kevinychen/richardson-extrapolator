/*
 * richardson.cpp
 *
 * Library functions for Richardson Extrapolation.
 */

#include "richardson.h"

#include <cstdlib>
#include <gmp.h>

namespace richardson {

void extrapolate(index_t num_samples, index_t start_index,
        SequenceFunc f, mpf_t ans) {
    // Calculate the desired samples, then pass to the other extrapolate
    // function.
    mpf_t *samples = (mpf_t*)malloc(sizeof(mpf_t) * num_samples);
    mpf_t *lim = samples + num_samples;
    mp_bitcnt_t precision = mpf_get_prec(ans);
    for (mpf_t *ptr = samples; ptr < lim; ptr++, start_index <<= 1) {
        mpf_init2(*ptr, precision);
        f(start_index, *ptr);
    }
    extrapolate(num_samples, samples, ans);

    // Clean
    for (mpf_t *ptr = samples; ptr < lim; ptr++)
        mpf_clear(*ptr);
    free(samples);
}

void extrapolate(index_t num_samples, mpf_t *samples, mpf_t ans) {
    // The Richardson extrapolation recursive formula is
    //
    // A_n+1(x) = (2^n A_n(2x) - A_n(x)) / (2^n - 1)

    mpf_t mult;  // mult = 2^n
    mpf_init_set_d(mult, 1);

    mpf_t denom;  // denom = 1 / (mult - 1)
    mp_bitcnt_t precision = mpf_get_prec(ans);
    mpf_init2(denom, precision);

    mpf_t *end = samples + num_samples;
    for (mpf_t *lim = samples; lim < end; lim++) {  // lim - samples = n
        mpf_mul_ui(mult, mult, 2);
        mpf_set(denom, mult);
        mpf_sub_ui(denom, denom, 1);
        mpf_ui_div(denom, 1, denom);
        // evaluate all extrapolations
        for (mpf_t *ptr = end; --ptr > lim; ) {
            // A_n+1(x) = (mult A_n(2x) - A_n(x)) * denom
            mpf_mul(*ptr, *ptr, mult);
            mpf_sub(*ptr, *ptr, *(ptr - 1));
            mpf_mul(*ptr, *ptr, denom);
        }
    }
    mpf_set(ans, *(end - 1));  // move to ans

    // Clean
    mpf_clear(mult);
    mpf_clear(denom);
}

}

