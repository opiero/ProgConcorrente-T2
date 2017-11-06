#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

#define MATRIX_FILE "matriz.txt"
#define VECTOR_FILE "vetor.txt"
#define OUTPUT_FILE "resultado.txt"

typedef double element;



element ** read_mat(char * filename, int * row, int * col) {
	(*row) = (*col) = 0;
	element ** mat = NULL;

	return mat;
}


int main (int argc, char * argv[]) {

	int row, col;
	element ** mat = read_mat(MATRIX_FILE, &row, &col);

	int my_rank = -1;
	int num_proc;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);



	MPI_Finalize();

	return EXIT_SUCCESS;
}
