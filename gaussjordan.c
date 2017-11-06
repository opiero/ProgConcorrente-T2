#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

#define MATRIX_FILE "matriz.txt"
#define VECTOR_FILE "vetor.txt"
#define OUTPUT_FILE "resultado.txt"

#define ALLOC_INIT_SIZE 8 //Nao pode ser 0
#define LINE_ENDER 13

typedef double element;



long int get_filesize(FILE * file) {
	long int cur_pos = ftell(file);
	fseek(file, 0, SEEK_END);
	long int filesize = ftell(file);
	fseek(file, cur_pos, SEEK_SET);

	return filesize;
}



element * read_vect(char * filename, int * col){
	(*col) = 0;
	int alloced_col = ALLOC_INIT_SIZE;
	element * vect = (element*) malloc(alloced_col*sizeof(element));

	FILE * vect_file = fopen(filename, "r");
	long int filesize = get_filesize(vect_file);
	while(ftell(vect_file) < filesize){
		if ((*col) == alloced_col){
			alloced_col*=2;
			vect = (element*) realloc(vect, alloced_col*sizeof(element));
		}

		if (fscanf(vect_file, "%lf", vect + (*col)) != 1) printf("Reading error\n");
		fgetc(vect_file); fgetc(vect_file);
		(*col)++;
	}

	fclose(vect_file);
	return vect;
}

void print_vect(element * vect, int col){
	int i;
	for (i=0; i<col; i++) printf("%.3lf\n", vect[i]);
}



element ** read_mat(char * filename, int * row, int * col) {
	(*col) = 0;
	int alloced_col = ALLOC_INIT_SIZE;
	element * aux_vect = (element*) malloc(alloced_col*sizeof(element));

	char cur_char = 0;
	FILE * mat_file = fopen(filename, "r");
	while (cur_char != LINE_ENDER){
		if ((*col) == alloced_col){
			alloced_col*=2;
			aux_vect = (element*) realloc(aux_vect, alloced_col*sizeof(element));
		}

		if (fscanf(mat_file, "%lf", aux_vect + (*col)) != 1) printf("Reading error\n");
		cur_char = fgetc(mat_file);
		(*col)++;
	}
	fgetc(mat_file);
	if ((*col) == alloced_col){
		alloced_col++;
		aux_vect = (element*) realloc(aux_vect, alloced_col*sizeof(element));
	}

	int alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	mat[0] = aux_vect;
	(*row) = 1;
	long int filesize = get_filesize(mat_file);
	while(ftell(mat_file) < filesize){
		if ((*row) == alloced_row){
			alloced_row*=2;
			mat = (element**) realloc(mat, alloced_row*sizeof(element*));
		}

		int i;
		mat[(*row)] = (element*) malloc(((*col)+1)*sizeof(element));
		for (i=0; i<(*col); i++) if (fscanf(mat_file, "%lf", mat[(*row)] + i) != 1) printf("Reading error\n");
		(*row)++;
		fgetc(mat_file); fgetc(mat_file);
	}

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
	int vcol, mrow, mcol;
	element ** mat = read_mat(MATRIX_FILE, &mrow, &mcol);
	element * vect = read_vect(VECTOR_FILE, &vcol);
	append_col(mat, mrow, &mcol, vect);	
	
	free_mat(mat, mrow);
	free(vect);
	return 0;

	/*int row, col;
	element ** mat = read_mat(MATRIX_FILE, &row, &col);

	int my_rank = -1;
	int num_proc;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);



	MPI_Finalize();

	return EXIT_SUCCESS;*/
}
