all:
	mpicc gaussjordan.c -Wall -Wextra -fopenmp -O3 -o gaussjordan.out

clean:
	rm *.out
	
run:
	mpirun --hostfile hosts -np 12 ./gaussjordan.out 8

test:
	mpicc gaussjordan.c -Wall -Wextra -fopenmp -o gaussjordan.out -g
	mpirun --hostfile hosts -np 12 ./gaussjordan.out 8