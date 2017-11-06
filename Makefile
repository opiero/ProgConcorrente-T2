all:
	mpicc gaussjordan.c -fopenmp -O3 -o gaussjordan.out

clean:
	rm *.out
	
run:
	mpirun ./gaussjordan.out

test:
	mpicc gaussjordan.c -Wall -Wextra -g -fopenmp -o gaussjordan.out
	mpirun ./gaussjordan.out