/*
TRABALHO 2 DE PROGRAMACAO CONCORRENTE

Alunos: 	Bernardo Barcellos de Castro Cunha,
			Eduardo Santos Carlos de Souza,
			Gustavo Cabral de Barros,
			Piero Lima Capelo

NUSP:		9293380,
			9293481,
			9293028,
			9293115

Professor: 	Paulo Sergio Lopes de Souza
*/

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

void delete_mat (element ** mat, int row) {
	int i;

	for (i = 0; i < row; i++)
		free(element[i]);

	free(element);
}

int main (int argc, char * argv[]) {

	int row, col;
	element ** mat = read_mat(MATRIX_FILE, &row, &col);

	int order;
	int my_rank = -1;
	int num_proc;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &num_proc);

	order = atoi(argv[1]);

	int job = row/order;

	for (i = 0; i < row; i++) {

		int pivo = mat[i][i];

		for (j = my_rank*job; j < (my_rank+1) * job; j++) {

			//zerar

		}

	}


	MPI_Finalize();

	delete_mat(mat, row);

	return EXIT_SUCCESS;
}
