/*
 * richardson.cpp
 *
 * Library functions for Richardson Extrapolation.
 *
 * Usage: ./richardson [filename] -p [precision]
 *   the file given by [filename] should contain k+1 (the number of samples)
 *     lines, with a single floating point number on each line, representing
 *     the samples a_n, a_2n, a_4n, a_8n, ... a_(2^k)n.
 *   the precision is the number of significant bits that the Richardson
 *     extrapolation is calculated with.
 */

#include "richardson.h"

#include <cstdlib>
#include <cstring>
#include <fstream>
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

using namespace richardson;

int main(int argc, char *argv[]) {
    const mp_bitcnt_t DEFAULT_PRECISION = 64;

    // Parse arguments
    char *filename = NULL;
    mp_bitcnt_t precision = -1;
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-p") == 0)
            precision = atoi(argv[++i]);
        else if (argv[i][0] != '-' && filename == NULL)
            filename = argv[i];
    }

    // Validate arguments
    if (filename == NULL) {
        printf(
"Usage: ./richardson [filename] -p [precision]\n"
"  the file given by [filename] should contain k+1 (the number of samples)\n"
"    lines, with a single floating point number on each line, representing\n"
"    the samples a_n, a_2n, a_4n, a_8n, ... a_(2^k)n.\n"
"  the precision is the number of significant bits that the Richardson\n"
"    extrapolation is calculated with.\n"
        );
        return 1;
    }
    if (precision == -1) {
        printf("Using default precision of %lu\n", DEFAULT_PRECISION);
        precision = DEFAULT_PRECISION;
    }
    if (precision <= 0) {
        printf("Error: precision must be a positive integer\n");
        return 1;
    }

    // Count number of samples
    std::ifstream in;
    in.open(filename);
    index_t num_samples = 0;
    std::string s;
    while (in >> s)
        num_samples++;
    in.close();

    // Validate number of samples
    if (num_samples < 1 || num_samples >= 8 * sizeof(index_t)) {
        printf("Error: number of samples must be an integer in [1, %lu)\n",
                8 * sizeof(index_t));
        return 1;
    }
    printf("Found %u samples:\n", num_samples);

    // Initialize answer mpf_t
    mpf_t ans;
    mpf_init2(ans, precision);

    // Read in samples
    mpf_t *samples = (mpf_t*)malloc(sizeof(mpf_t) * num_samples);
    in.open(filename);
    for (index_t index = 0; index < num_samples; index++) {
        in >> s;
        mpf_init2(samples[index], precision);
        mpf_set_str(samples[index], s.c_str(), 10);
        mpf_out_str(NULL, 10, precision, samples[index]);
        printf("\n");
    }
    in.close();

    extrapolate(num_samples, samples, ans);

    // Print answer
    printf("Extrapolation: ");
    mpf_out_str(NULL, 10, precision, ans);
    printf("\n");

    for (index_t index = 0; index < num_samples; index++)
        mpf_clear(samples[index]);
    free(samples);
    mpf_clear(ans);

    return 0;
}
