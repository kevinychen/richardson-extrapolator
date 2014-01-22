/*
 * binary.cpp
 *
 * Main for Richardson Extrapolation on a data file.
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
