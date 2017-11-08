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



element * read_vect(char * filename, int * nrows){
	(*nrows) = 0;
	int alloced_row = ALLOC_INIT_SIZE;
	element * vect = (element*) malloc(alloced_row*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * vect_file = fopen(filename, "r");
	while(getline(&line, &linesize, vect_file) != -1){
		while ((*nrows) >= alloced_row){
			alloced_row*=2;
			vect = (element*) realloc(vect, alloced_row*sizeof(element));
		}

		sscanf(line, "%lf", vect + ((*nrows)++));
		free(line);
		line = NULL;
	}
	free(line);
	fclose(vect_file);

	return vect;
}

void print_vect(element * vect, int nrows){
	int i;
	for (i=0; i<nrows; i++) printf("%.3lf\n", vect[i]);
}



element ** read_mat(char * filename, int * nrows, int * ncols) {
	(*nrows) = (*ncols) = 0;
	int alloced_col = ALLOC_INIT_SIZE, alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	(*mat) = (element*) malloc(alloced_col*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * mat_file = fopen(filename, "r");
	getline(&line, &linesize, mat_file);
	char * aux_str = strtok(line, DELIMITERS);
	while(aux_str != NULL){
		while ((*ncols) >= alloced_col){
			alloced_col*=2;
			(*mat) = (element*) realloc((*mat), alloced_col*sizeof(element));
		}

		sscanf(aux_str, "%lf", (*mat) + ((*ncols)++));
		aux_str = strtok(NULL, DELIMITERS);
	}
	if (alloced_col < (*ncols)+1) (*mat) = (element*) realloc((*mat), ((*ncols)+1)*sizeof(element));
	free(line);

	line = NULL;
	(*nrows) = 1;
	while (getline(&line, &linesize, mat_file) != -1){
		while ((*nrows) >= alloced_row){
			alloced_row*=2;
			mat = (element**) realloc(mat, alloced_row*sizeof(element*));
		}
		mat[(*nrows)] = (element*) malloc(((*ncols)+1)*sizeof(element));

		int cur_col = 0;
		aux_str = strtok(line, DELIMITERS);
		while(aux_str != NULL) {
			sscanf(aux_str, "%lf", mat[(*nrows)] + (cur_col++));
			aux_str = strtok(NULL, DELIMITERS);
		}

		(*nrows)++;
		free(line);
		line = NULL;
	}
	free(line);
	fclose(mat_file);

	return mat;
}

void print_mat(element ** mat, int nrows, int ncols) {
	int i, j;
	for (i=0; i<nrows; i++){
		for (j=0; j<ncols; j++) printf("%.3lf ", mat[i][j]);
		printf("\n");
	}
}

void free_mat(element ** mat, int nrows) {
	int i;
	for (i=0; i<nrows; i++) free(mat[i]);
	free(mat);
}



void append_col(element ** mat, int nrows, int * ncols, element * vect) {
	int i;
	for (i=0; i<nrows; i++) mat[i][(*ncols)] = vect[i];
	(*ncols)++;
}



int sequential_max_line(element ** mat, int nrows, int col) {
	int i, max_idx = -1;
	element max = -1, mod;
	for (i=0; i<nrows; i++){
		mod = (mat[i][col] >= 0) ? mat[i][col] : -mat[i][col];
		if (mod > max){
			max = mod;
			max_idx = i;
		}
	}

	return max_idx;
}

element * seuquential_gaussjordan(element ** mat, int nrows, int ncols) {
	ncols--;
	int i, j, k, min = (ncols <= nrows) ? ncols : nrows;
	for (j=0; j<min; j++){
		int max_idx = sequential_max_line(mat, nrows, j);
		element * aux_row = mat[j];
		mat[j] = mat[max_idx];
		mat[max_idx] = aux_row;

		for (i=0; i<nrows; i++){
			if (i!=j){
				element rat = mat[i][j]/mat[j][j];
				for (k=0; k<ncols+1; k++) mat[i][k] -= rat*mat[j][k];
			}
		}
	}

	element * res = (element*) malloc(nrows*sizeof(element));
	for (i=0; i<nrows; i++) res[i] = mat[i][ncols]/mat[i][i];

	return res;
}



int main (int argc, char * argv[]) {
	int vnrows, mnrows, mncols;
	element ** mat = read_mat(MATRIX_FILE, &mnrows, &mncols);
	element * vect = read_vect(VECTOR_FILE, &vnrows);
	append_col(mat, mnrows, &mncols, vect);
	free(vect);

	element * res = seuquential_gaussjordan(mat, mnrows, mncols);
	print_vect(res, mnrows);
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

	free_mat(mat, mnrows);
	return EXIT_SUCCESS;
}
