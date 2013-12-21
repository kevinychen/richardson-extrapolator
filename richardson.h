/*
 * richardson.h
 *
 * Library functions for Richardson Extrapolation.
 *
 * The library contains two main functions, both named extrapolate(). The
 * first one (with four arguments) can be used with a general function that
 * can calculate the nth term of a sequence. The second one (with three
 * arguments) should be used when computing an arbitrary term is expensive,
 * and it is more efficient to calculate all n terms of a sequence in one
 * run, then passing the necessary terms directly to extrapolate().
 *
 * Example usage:
 *
 * // Define SequenceFunc describing the sequence 1/n
 * void f(index_t index, mpf_t result) {
 *   mpf_set_d(result, index);
 *   mpf_ui_div(result, 1, result);
 * }
 *
 * const index_t num_samples = 10;  // must be < 32
 * const index_t start_index = 1;
 * mpf_t ans;
 * mpf_init(ans);
 *
 * // Richardson extrapolate from the samples a1, a2, a4, ... a512.
 * extrapolate(num_samples, start_index, &f, ans);
 *
 * // Alternatively, we can compute only the required samples, then directly
 * // feed them into the three argument variant of extrapolate().
 * mpf_t samples[num_samples];
 * for (index_t i = 0; i < num_samples; i++)
 *   f(start_index << i, samples + i);
 * extrapolate(num_samples, samples, ans);
 *
 * // do something with ans
 *
 * // Destroy the mpf_t used to store the answer
 * mpf_clear(ans);
 */

#ifndef RICHARDSON_H_
#define RICHARDSON_H_

#include <gmp.h>
#include <stdint.h>

namespace richardson {

/*
 * index_t is the basic nonnegative integer type, representing indices of a
 * sequence starting from 0:
 *
 * a0, a1, a2, ...
 */
typedef uint32_t index_t;

/*
 * SequenceFunc is the basic function type, representing a sequence. It takes
 * an index as input and updates the mpf_t with the value of the corresponding
 * term of the sequence.
 */
typedef void (*SequenceFunc)(index_t, mpf_t);


/*
 * extrapolate() takes as input a sequence, which is described by a SequenceFunc.
 * A SequenceFunc takes an index n and outputs the nth element of the sequence
 * (0-indexed). extrapolate() will only select a few samples from the sequence
 * to perform the Richardson extrapolation. Specifically, given a sequence
 *
 * a_0, a_1, a_2, a_3, ...
 *
 * the function will Richardson extrapolate according to the samples
 *
 * a_n, a_2n, a_4n, a_8n, ... a_(2^k)n,
 *
 * where n and k are provided parameters.
 *
 * Parameters:
 *   num_samples: the number of samples to take from the sequence, or k+1
 *     in the above example.
 *   start_index: the index of the first sample to take, or n in the above
 *     example.
 *   f: a SequenceFunc that takes an index n and outputs the nth element
 *     of the sequence (0-indexed).
 *   ans: an initialized mpf_t where the answer will be stored upon
 *     successful completion of the function.
 */
void extrapolate(index_t num_samples, index_t start_index, SequenceFunc f, mpf_t ans);


/*
 * A variant of the extrapolate() function. This function takes already
 * computed samples from a sequence, specifically the terms
 *
 * a_n, a_2n, a_4n, a_8n, ... a_(2^k)n,
 *
 * and Richardson extrapolates the limit of the sequence.
 *
 * Parameters:
 *   num_samples: the number of samples taken, or k+1 in the above example.
 *   samples: an array consisting of the samples in the order shown above,
 *     i.e., samples[i] = a_(2^i)n.
 *   ans: an initialized mpf_t where the answer will be stored upon
 *     successful completion of the function.
 */
void extrapolate(index_t num_samples, mpf_t *samples, mpf_t ans);

}

#endif  // RICHARDSON_H_
