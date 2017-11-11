/*
TRABALHO 2 DE PROGRAMACAO CONCORRENTE

Alunos:
Bernardo Barcellos de Castro Cunha	9293380
Eduardo Santos Carlos de Souza	9293481
Gustavo Cabral de Barros	9293028
Piero Lima Capelo	9293115

Professor:
Paulo Sergio Lopes de Souza
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
#define DELIMITERS " \t,;"

typedef double element;
#define ELEMENT_MPI MPI_DOUBLE
#define ELEMENT_READ_MASK "%lf"
#define ELEMENT_PRINT_MASK "%.3lf\n"

#define MOD(n) (((n) >= 0) ? (n) : (-(n)))



element * read_vect(char * filename, int * nrows)
{
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

		sscanf(line, ELEMENT_READ_MASK, vect + ((*nrows)++));
		free(line);
		line = NULL;
	}
	free(line);
	fclose(vect_file);

	return vect;
}

void print_vect(char * filename, element * vect, int nrows)
{
	FILE * out_file;
	if (strcmp("stdout", filename) == 0) out_file = stdout;
	else if (strcmp("stderr", filename) == 0) out_file = stderr;
	else out_file = fopen(filename, "w");

	int i;
	for (i=0; i<nrows; i++) fprintf(out_file, ELEMENT_PRINT_MASK, vect[i]);
	
	if ((out_file != stdout) && (out_file != stderr)) fclose(out_file);
}



element ** read_mat(char * filename, int * nrows, int * ncols)
{
	(*nrows) = (*ncols) = 0;
	int alloced_col = ALLOC_INIT_SIZE, alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	(*mat) = (element*) malloc(alloced_col*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * mat_file = fopen(filename, "r");
	if (getline(&line, &linesize, mat_file) == -1) printf("No lines on file\n");
	char * aux_str = strtok(line, DELIMITERS);
	while(aux_str != NULL){
		while ((*ncols) >= alloced_col){
			alloced_col*=2;
			(*mat) = (element*) realloc((*mat), alloced_col*sizeof(element));
		}

		sscanf(aux_str, ELEMENT_READ_MASK, (*mat) + ((*ncols)++));
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
			sscanf(aux_str, ELEMENT_READ_MASK, mat[(*nrows)] + (cur_col++));
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

void free_mat(void ** mat, int nrows)
{
	int i;
	for (i=0; i<nrows; i++) free(mat[i]);
	free(mat);
}

void append_col(element ** mat, int nrows, int * ncols, element * col)
{
	int i;
	for (i=0; i<nrows; i++) mat[i][(*ncols)] = col[i];
	(*ncols)++;
}



element * sequential_gaussjordan(element ** mat, int nrows, int ncols)
{
	ncols--;
	int i, j, k, max_idx, min = (ncols <= nrows) ? ncols : nrows;
	element max;
	for (j=0; j<min; j++){
		max_idx = -1;
		max = 0;
		for (i=j; i<nrows; i++){
			if (MOD(mat[i][j]) > MOD(max)){
				max = mat[i][j];
				max_idx = i;
			}
		}
		element * aux_row = mat[j];
		mat[j] = mat[max_idx];
		mat[max_idx] = aux_row;

		ncols++;
		for (i=0; i<nrows; i++){
			if (i!=j){
				element rat = mat[i][j]/mat[j][j];
				for (k=0; k<ncols; k++) mat[i][k] -= rat*mat[j][k];
			}
		}
		ncols--;
	}

	element * res = (element*) malloc(nrows*sizeof(element));
	for (i=0; i<nrows; i++) res[i] = mat[i][ncols]/mat[i][i];

	return res;
}



void send_initial_lines(element ** mat, int nrows, int ncols, int nproc, int * init_tag)
{
	MPI_Bcast(&nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int * job_lines_idxs = (int*) malloc((nrows/nproc + nrows%nproc)*sizeof(int));
	int i, j, job_size = nrows/nproc, last_job = job_size;
	for (i=1; i<nproc; i++){
		if (i==nproc-1) job_size += nrows%nproc;
		for (j=last_job; j<last_job+job_size; j++) job_lines_idxs[j-last_job] = j;
		last_job = j;

		MPI_Send(&job_size, 1, MPI_INT, i, (*init_tag), MPI_COMM_WORLD);
		MPI_Send(job_lines_idxs, job_size, MPI_INT, i, (*init_tag)+1, MPI_COMM_WORLD);
		for (j=0; j<job_size; j++) MPI_Send(mat[job_lines_idxs[j]], ncols, ELEMENT_MPI, i, (*init_tag)+2, MPI_COMM_WORLD);
	}
	free(job_lines_idxs);
	(*init_tag)+=2;
}

element ** recv_initial_lines(int * nrows, int * ncols, int * job_size, int ** job_lines_idxs, int * init_tag)
{
	MPI_Status status;
	MPI_Bcast(nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);

	int i;
	MPI_Recv(job_size, 1, MPI_INT, 0, (*init_tag), MPI_COMM_WORLD, &status);
	(*job_lines_idxs) = (int*) malloc((*job_size)*sizeof(int));
	element ** mat = (element**) malloc((*job_size)*sizeof(element*));
	for (i=0; i<(*job_size); i++) mat[i] = (element*) malloc((*ncols)*sizeof(element));
	
	MPI_Recv((*job_lines_idxs), (*job_size), MPI_INT, 0, (*init_tag)+1, MPI_COMM_WORLD, &status);
	for (i=0; i<(*job_size); i++) MPI_Recv(mat[i], (*ncols), ELEMENT_MPI, 0, (*init_tag)+2, MPI_COMM_WORLD, &status);
	
	(*init_tag)+=2;
	return mat;
}

void send_final_answer(element ** mat, int ncols, int job_size, int * job_lines_idxs, int * init_tag)
{
	int i;
	element * res = (element*) malloc(job_size*sizeof(element));
	for (i=0; i<job_size; i++) res[i] = mat[i][ncols];

	MPI_Send(job_lines_idxs, job_size, MPI_INT, 0, (*init_tag), MPI_COMM_WORLD);
	MPI_Send(res, job_size, ELEMENT_MPI, 0, (*init_tag)+1, MPI_COMM_WORLD);
	(*init_tag)+=1;
	free(res);
}

element * recv_final_answer(int nrows, int nproc,  int * init_tag)
{
	int i, j, aux_size = (nrows / nproc);
	int * sizes = (int*) malloc((nproc-1)*sizeof(int));
	for (i=0; i<nproc-1; i++){
		if (i==nproc-2) aux_size += (nrows % nproc);
		sizes[i] = aux_size;
	}
	
	MPI_Status status;
	int ** idxs = (int**) malloc((nproc-1)*sizeof(int*));
	for (i=0; i<nproc-1; i++){
		idxs[i] = (int*) malloc(sizes[i]*sizeof(int));
		MPI_Recv(idxs[i], sizes[i], MPI_INT, i+1, (*init_tag), MPI_COMM_WORLD, &status);
	}

	element * aux_vect = (element*) malloc(aux_size*sizeof(element));
	element * res = (element*) malloc(nrows*sizeof(element));
	for (i=0; i<nproc-1; i++){
		MPI_Recv(aux_vect, sizes[i], ELEMENT_MPI, i+1, (*init_tag)+1, MPI_COMM_WORLD, &status);
		for (j=0; j<sizes[i]; j++) res[idxs[i][j]] = aux_vect[j];
	}
	(*init_tag)+=1;

	free_mat((void**)idxs, nproc-1);
	free(aux_vect);
	free(sizes);
	return res;
}

element * parallel_gaussjordan(element ** mat, int nrows, int ncols, int job_size, int * job_lines_idxs, int nproc, int rank, int * init_tag, int n_threads)
{
	ncols--;
	int i, j, k, min = (ncols <= nrows) ? ncols : nrows;
	int local_max_idx, global_max_idx, global_max_proc;
	int * recv_idx, * aux_max_idx_vect = (int*) malloc(n_threads*sizeof(int));
	element local_max, global_max;
	element * global_max_line = (element*) malloc((ncols+1)*sizeof(element)), * recv_max, * aux_max_vect = (element*) malloc(n_threads*sizeof(element));
	if (rank==0){
		recv_idx = (int*) malloc(nproc*sizeof(int));
		recv_max = (element*) malloc(nproc*sizeof(element));
	}

	omp_set_dynamic(0);
	omp_set_num_threads(n_threads);
	for (j=0; j<min; j++){
		local_max = 0;
		local_max_idx = -1;
		memset(aux_max_vect, 0, n_threads*sizeof(element));
		memset(aux_max_idx_vect, -1, n_threads*sizeof(int));
		#pragma omp parallel for shared(aux_max_vect, aux_max_idx_vect) private(i)
		for (i=0; i<job_size; i++)
		{
			int id = omp_get_thread_num();
			if ((job_lines_idxs[i] >= j) && (MOD(mat[i][j]) > MOD(aux_max_vect[id]))){
				aux_max_vect[id] = mat[i][j];
				aux_max_idx_vect[id] = i;
			}
		}
		for (i=0; i<n_threads; i++){
			if (MOD(aux_max_vect[i]) >= MOD(local_max)){
				local_max = aux_max_vect[i];
				local_max_idx = aux_max_idx_vect[i];
			}
		}

		MPI_Gather(job_lines_idxs + local_max_idx, 1, MPI_INT, recv_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(&local_max, 1, ELEMENT_MPI, recv_max, 1, ELEMENT_MPI, 0, MPI_COMM_WORLD);

		if (rank==0){
			global_max = 0;
			global_max_proc = global_max_idx = -1;
			for (k=0; k<nproc; k++){
				if (MOD(recv_max[k]) > MOD(global_max)) {
					global_max = recv_max[k];
					global_max_idx = recv_idx[k];
					global_max_proc = k;
				}
			}
		}
		
		MPI_Bcast(&global_max, 1, ELEMENT_MPI, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_max_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_max_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

		#pragma omp parallel for private(i)
		for (i=0; i<job_size; i++)
		{
			if (job_lines_idxs[i]==j) job_lines_idxs[i] = global_max_idx;
		}
		if (rank == global_max_proc){
			job_lines_idxs[local_max_idx] = j;
			memcpy(global_max_line, mat[local_max_idx], (ncols+1)*sizeof(element));
		}

		//printf("%d\n", global_max_proc);
		MPI_Bcast(global_max_line, ncols+1, ELEMENT_MPI, global_max_proc, MPI_COMM_WORLD);

		ncols++;
		#pragma omp parallel for private(i)
		for (i=0; i<job_size; i++)
		{
			if (job_lines_idxs[i] != j){
				element rat = mat[i][j]/global_max;
				for (k=0; k<ncols; k++) mat[i][k] -= rat*global_max_line[k];
			}
		}
		ncols--;

		MPI_Barrier(MPI_COMM_WORLD);
	}
	#pragma omp parallel for private(i)
	for (i=0; i<job_size; i++)
	{
		mat[i][ncols] /= mat[i][job_lines_idxs[i]];
	}

	free(aux_max_idx_vect);
	free(aux_max_vect);
	free(global_max_line);
	if (rank==0){
		free(recv_idx);
		free(recv_max);
	}

	if (rank == 0){
		element * res = recv_final_answer(nrows, nproc, init_tag);
		for (i=0; i<job_size; i++) res[job_lines_idxs[i]] = mat[i][ncols];
		return res;
	}
	else{
		send_final_answer(mat, ncols, job_size, job_lines_idxs, init_tag);
		return NULL;
	}
}



int main (int argc, char * argv[])
{
	int nrows, ncols, job_size;
	int * job_lines_idxs;
	element ** mat;

	int rank, msgtag = 0, nproc;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	if (rank == 0) {
		int vnrows;
		element * vect = read_vect(VECTOR_FILE, &vnrows);
		mat = read_mat(MATRIX_FILE, &nrows, &ncols);
		append_col(mat, nrows, &ncols, vect);
		free(vect);

		int i;
		job_size = nrows/nproc;
		job_lines_idxs = (int*) malloc(job_size*sizeof(int));
		for (i=0; i<job_size; i++) job_lines_idxs[i] = i;

		send_initial_lines(mat, nrows, ncols, nproc, &msgtag);
	}
	else mat = recv_initial_lines(&nrows, &ncols, &job_size, &job_lines_idxs, &msgtag);

	element * res = parallel_gaussjordan(mat, nrows, ncols, job_size, job_lines_idxs, nproc, rank, &msgtag, 1);
	if (rank == 0){
		print_vect("stdout", res, nrows);
		free(res);
	}
	
	free(job_lines_idxs);
	free_mat((void**)mat, ((rank==0) ? nrows : job_size));
	MPI_Finalize();
	return 0;
}
