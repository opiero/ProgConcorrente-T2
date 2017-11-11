all:
	mpicc gaussjordan.c -Wall -Wextra -fopenmp -O3 -o gaussjordan.out

clean:
	rm *.out
	
run:
	mpirun ./gaussjordan.out -np 8 8 time.txt

test:
	mpicc gaussjordan.c -Wall -Wextra -fopenmp -o gaussjordan.out -g
	mpirun ./gaussjordan.out -np 8 8 time.txt

val:
	mpicc gaussjordan.c -Wall -Wextra -fopenmp -o gaussjordan.out -g
	valgrind --leak-check=full --track-origins=yes mpirun ./gaussjordan.out -np 8 8 time.txt