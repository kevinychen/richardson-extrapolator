# Makefile
#
# run:
#   make  (to create the binary called richardson to perform Richardson
#          Extrapolation via command line)
#   make test  (to run the tests)
#   make clean

all:
	g++ binary.cpp richardson.cpp -lgmp -o richardson

test:
	g++ richardson_test.cpp richardson.cpp -lgmp -o richardson_test
	./richardson_test
	rm richardson_test

clean:
	rm richardson -f
