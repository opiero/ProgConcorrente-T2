#include <stdio.h>
#include <stdlib.h>

#define MATRIX_FILE "matriz.txt"
#define VECTOR_FILE "vetor.txt"
#define OUTPUT_FILE "resultado.txt"

typedef double element;



element ** read_mat(char * filename, int * row, int * col)
{
	(*row) = (*col) = 0;
	element ** mat = NULL;

	return mat;
}

int main()
{
	int row, col;
	element ** mat = read_mat(MATRIX_FILE, &row, &col);

	return 0;
}