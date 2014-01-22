Richardson Extrapolator
=======================

A C++ library that performs numerically stable Richardson Extrapolation on arbitrary size input.

The library requires the GMP library, which on Ubuntu can be installed as follows:

        $ sudo apt-get install libgmp3-dev


To run via command line:
------------------------

1. Build the binary. This requires g++.

        $ make

2. This will create a binary called richardson. Entering

        $ ./richardson

   will show the usage:

        Usage: ./richardson [filename] -p [precision]
          the file given by [filename] should contain k+1 (the number of samples)
            lines, with a single floating point number on each line, representing
            the samples a_n, a_2n, a_4n, a_8n, ... a_(2^k)n.
          the precision is the number of significant bits that the Richardson
            extrapolation is calculated with.

3. As an example, we can perform Richardson Extrapolation on the sample input file, using a precision of 128:

        $ ./richardson sample.txt -p 128


To use as a library:
--------------------

The library contains two main functions, both named extrapolate(). The
first one (with four arguments) can be used with a general function that
can calculate the nth term of a sequence. The second one (with three
arguments) should be used when computing an arbitrary term is expensive,
and it is more efficient to calculate all n terms of a sequence in one
run, then passing the necessary terms directly to extrapolate().

Example usage:

        // Define SequenceFunc describing the sequence 1/n
        void f(index_t index, mpf_t result) {
          mpf_set_d(result, index);
          mpf_ui_div(result, 1, result);
        }
        
        const index_t num_samples = 10;  // must be < 32
        const index_t start_index = 1;
        mpf_t ans;
        mpf_init(ans);
        
        // Richardson extrapolate from the samples a1, a2, a4, ... a512.
        extrapolate(num_samples, start_index, &f, ans);
        
        // Alternatively, we can compute only the required samples, then directly
        // feed them into the three argument variant of extrapolate().
        mpf_t samples[num_samples];
        for (index_t index = 0; index < num_samples; index++) {
          mpf_init(samples[index]);
          f(start_index << index, samples[index]);
        }
        extrapolate(num_samples, samples, ans);
        
        // do something with ans
        
        // Destroy the mpf_t used to store the samples and answer
        for (index_t index = 0; index < num_samples; index++)
          mpf_clear(samples[index]);
        mpf_clear(ans);


To run tests:
-------------

Run the `make test` command. This requires g++.

