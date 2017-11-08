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
#include <string.h>
#include <omp.h>
#include <mpi.h>

#define MATRIX_FILE "matriz.txt"
#define VECTOR_FILE "vetor.txt"
#define OUTPUT_FILE "resultado.txt"

#define ALLOC_INIT_SIZE 8 //Nao pode ser < 0
#define DELIMITERS " \t\n\0,;"

typedef double element;



element * read_vect(char * filename, int * row){
	(*row) = 0;
	int alloced_row = ALLOC_INIT_SIZE;
	element * vect = (element*) malloc(alloced_row*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * vect_file = fopen(filename, "r");
	while(getline(&line, &linesize, vect_file) != -1){
		while ((*row) >= alloced_row){
			alloced_row*=2;
			vect = (element*) realloc(vect, alloced_row*sizeof(element));
		}

		sscanf(line, "%lf", vect + ((*row)++));
		free(line);
		line = NULL;
	}
	free(line);
	fclose(vect_file);

	return vect;
}

void print_vect(element * vect, int row){
	int i;
	for (i=0; i<row; i++) printf("%.3lf\n", vect[i]);
}



element ** read_mat(char * filename, int * row, int * col) {
	(*row) = (*col) = 0;
	int alloced_col = ALLOC_INIT_SIZE, alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	(*mat) = (element*) malloc(alloced_col*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * mat_file = fopen(filename, "r");
	getline(&line, &linesize, mat_file);
	char * aux_str = strtok(line, DELIMITERS);
	while(aux_str != NULL){
		while ((*col) >= alloced_col){
			alloced_col*=2;
			(*mat) = (element*) realloc((*mat), alloced_col*sizeof(element));
		}

		sscanf(aux_str, "%lf", (*mat) + ((*col)++));
		aux_str = strtok(NULL, DELIMITERS);
	}
	if (alloced_col < (*col)+1) (*mat) = (element*) realloc((*mat), ((*col)+1)*sizeof(element));
	free(line);

	line = NULL;
	(*row) = 1;
	while (getline(&line, &linesize, mat_file) != -1){
		while ((*row) >= alloced_row){
			alloced_row*=2;
			mat = (element**) realloc(mat, alloced_row*sizeof(element*));
		}
		mat[(*row)] = (element*) malloc(((*col)+1)*sizeof(element));

		int cur_col = 0;
		aux_str = strtok(line, DELIMITERS);
		while(aux_str != NULL) {
			sscanf(aux_str, "%lf", mat[(*row)] + (cur_col++));
			aux_str = strtok(NULL, DELIMITERS);
		}

		(*row)++;
		free(line);
		line = NULL;
	}
	free(line);
	fclose(mat_file);

	return mat;
}

void print_mat(element ** mat, int row, int col) {
	int i, j;
	for (i=0; i<row; i++){
		for (j=0; j<col; j++) printf("%.3lf ", mat[i][j]);
		printf("\n");
	}
}

void free_mat(element ** mat, int row) {
	int i;
	for (i=0; i<row; i++) free(mat[i]);
	free(mat);
}



void append_col(element ** mat, int mrow, int * mcol, element * vect) {
	int i;
	for (i=0; i<mrow; i++) mat[i][(*mcol)] = vect[i];
	(*mcol)++;
}



int main (int argc, char * argv[]) {
	int vrow, mrow, mcol;
	element ** mat = read_mat(MATRIX_FILE, &mrow, &mcol);
	element * vect = read_vect(VECTOR_FILE, &vrow);
	append_col(mat, mrow, &mcol, vect);
	free(vect);

	print_vect(vect, vrow);
	print_mat(mat, mrow, mcol);
/*
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

	MPI_Finalize();*/

	free_mat(mat, mrow);
	return EXIT_SUCCESS;
}
