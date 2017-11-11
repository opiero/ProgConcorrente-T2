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



/*
Funcao abre o arquivo verificando se não é stdin, stdout ou stderr
*/
FILE * open_file(char * filename, const char * mode)
{
	FILE * file;
	if (strcmp("stdout", filename) == 0) file = stdout;
	else if (strcmp("stderr", filename) == 0) file = stderr;
	else if (strcmp("stdin", filename) == 0) file = stdin;
	else file = fopen(filename, mode);

	return file;
}

/*
Funcao fecha o arquivo verificando se não é stdin, stdout ou stderr
*/
void close_file(FILE * file)
{
	if ((file != stdout) && (file != stderr) && (file != stdin))
		fclose(file);
}



element * read_vect(char * filename, int * nrows)
{
	(*nrows) = 0;
	int alloced_row = ALLOC_INIT_SIZE;
	element * vect = (element*) malloc(alloced_row*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * vect_file = open_file(filename, "r");
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
	close_file(vect_file);

	return vect;
}

/*
Funcao printa o vetor coluna de elementos no arquivo com nome filename
*/
void print_vect(char * filename, element * vect, int nrows)
{
	int i;
	FILE * out_file = open_file(filename, "w");
	//Imprime cada elemento da matriz como uma linha
	for (i=0; i<nrows; i++) fprintf(out_file, ELEMENT_PRINT_MASK, vect[i]);
	close_file(out_file);
}



element ** read_mat(char * filename, int * nrows, int * ncols)
{
	(*nrows) = (*ncols) = 0;
	int alloced_col = ALLOC_INIT_SIZE, alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	(*mat) = (element*) malloc(alloced_col*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * mat_file = open_file(filename, "r");
	if (getline(&line, &linesize, mat_file) == -1) printf("No lines on file\n");
	char * aux_str = strtok(line, DELIMITERS);
	while(aux_str != NULL){
		while ((*ncols) >= alloced_col){
			alloced_col*=2;
			(*mat) = (element*) realloc((*mat), alloced_col*sizeof(element));
		}

		if (sscanf(aux_str, ELEMENT_READ_MASK, (*mat) + (*ncols)) == 1) (*ncols)++;
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
	close_file(mat_file);

	return mat;
}

/*
Libera da heap uma matriz generica
*/
void free_mat(void ** mat, int nrows)
{
	int i;
	for (i=0; i<nrows; i++) free(mat[i]);
	free(mat);
}

/*
Adiciona uma coluna col na matriz mat, e incremeneta
o numero de colunas da matriz. A mariz ja deve possuir a
memeoria necessaria alocada
*/
void append_col(element ** mat, int nrows, int * ncols, element * col)
{
	int i;
	for (i=0; i<nrows; i++) mat[i][(*ncols)] = col[i];
	(*ncols)++;
}



/*
Realiza o algoritimo de gaussjordan utilizando apenas
uma maquina e uma thread, sequencialemnet. Retorna o vetor
resposta com o valor das variaveis. Aceita matrizes nao-quadradas.
*/
element * sequential_gaussjordan(element ** mat, int nrows, int ncols)
{
	//ncols-1 para desconsiderar a ultima coluna, que e a coluna com os vator resultado
	int i, j, k, max_idx, min = ((ncols-1) <= nrows) ? (ncols-1) : nrows;
	element max, rat;
	element * aux_row;
	for (j=0; j<min; j++){
		/*
		Pega o elemento com maior valor em modulo na coluna para reduzir erro numerico.
		Depois troca a linha desse elemento com a linha da coluna da iteracao atual.
		Desconsidera as linhas das iteracoes anteriores porque elas nao podem ser utilizadas,
		ja que elas mudariam os valores das posicoes ja zeradas.
		*/
		max_idx = -1;
		max = 0;
		for (i=j; i<nrows; i++){
			if (MOD(mat[i][j]) > MOD(max)){
				max = mat[i][j];
				max_idx = i;
			}
		}
		aux_row = mat[j];
		mat[j] = mat[max_idx];
		mat[max_idx] = aux_row;

		/*
		Para todas as linhas menos a da iteracao atual e feita
		a multiplicacao e soma das linhascom o objetivo de zerar
		todos os elementos da coluna menos a diagonal principal
		*/
		for (i=0; i<nrows; i++){
			if (i != j){
				rat = mat[i][j] / mat[j][j];
				for (k=0; k<ncols; k++) mat[i][k] -= rat*mat[j][k];
			}
		}
	}

	/*
	Calcular o vetor resposta dividindo a ultima coluna (coluna b em Ax = b)
	pelos coeficientes na diagonal principal
	*/
	element * res = (element*) malloc(nrows*sizeof(element));
	for (i=0; i<nrows; i++) res[i] = mat[i][ncols-1] / mat[i][i];

	return res;
}



/*
Cada processo atribui quais linhas da matriz ele computara,
baseado em seu rank
*/
void assign_initial_lines(int * nrows, int * ncols, int * job_size, int ** job_lines_idxs, int rank, int nproc)
{
	//Broadcas das dimensoes da matriz
	MPI_Bcast(nrows, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(ncols, 1, MPI_INT, 0, MPI_COMM_WORLD);

	//Atribuicao inicial de quem tera cada linha. Obs: O vetor job_lines_idxs não fica constante
	(*job_size) = (*nrows) / nproc;
	int i, init_pos = rank * (*job_size);
	if (rank == nproc-1) (*job_size) += (*nrows) % nproc;
	(*job_lines_idxs) = (int*) malloc((*job_size)*sizeof(int));
	for (i=init_pos; i<init_pos+(*job_size); i++) (*job_lines_idxs)[i-init_pos] = i;
}

/*
Processo 0 envia as linhas corretas para os outrso processos
*/
void send_initial_lines(element ** mat, int nrows, int ncols, int nproc, int * init_tag)
{	
	//Cada iteracao atribui as linhas para o processo i
	int i, j, job_size = nrows/nproc, last_iter;
	for (i=1; i<nproc; i++){
		last_iter = (i+1)*job_size;
		//Aumentar o numero de linhas para o ultimo processo
		if (i == nproc-1) last_iter += nrows % nproc;
		//Enviar para o processo i todas sua linhas
		for (j=i*job_size; j<last_iter; j++) MPI_Send(mat[j], ncols, ELEMENT_MPI, i, (*init_tag)+1, MPI_COMM_WORLD);
	}

	(*init_tag)+=1;
}

/*
Os processos recebem do processo 0 as linhas que irao computar
*/
element ** recv_initial_lines(int ncols, int job_size, int * init_tag)
{
	int i;
	MPI_Status status;
	element ** mat = (element**) malloc(job_size*sizeof(element*));
	//Aloca cada linha e recebe o vetor do processo 0 
	for (i=0; i<job_size; i++){
		mat[i] = (element*) malloc(ncols*sizeof(element));
		MPI_Recv(mat[i], ncols, ELEMENT_MPI, 0, (*init_tag)+1, MPI_COMM_WORLD, &status);
	}
	
	(*init_tag)+=1;
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
	else {
		recv_idx = NULL;
		recv_max = NULL;
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



/*
Imprime o tempo decorrido do processo no arquivo filename
*/
void print_time(char * filename, double elap_time)
{
	FILE * time_file = open_file(filename, "a");
	fprintf(time_file, "%lf\n", elap_time);
	close_file(time_file);
}



int main (int argc, char * argv[])
{
	//Variaveis necesarias para todos os processos
	int nrows = 0, ncols = 0, job_size = 0, n_threads = 1;
	int * job_lines_idxs = NULL;
	element * res = NULL;
	element ** mat = NULL;
	double elap_time = 0;

	//Variaveis para o OpenMPI necesarias para todos os processos
	int rank = 0, msgtag = 0, nproc = 1;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);
	if (rank == 0) {
		//Leitura do vetor e da matriz apenas no processo de rank 0
		int vnrows;
		element * vect = read_vect(VECTOR_FILE, &vnrows);
		mat = read_mat(MATRIX_FILE, &nrows, &ncols);
		//Adiciona o vetor lido como a ultima coluna da matriz para facilitar as operacoes
		append_col(mat, nrows, &ncols, vect);
		free(vect);
	}
	
	if (nproc == 1){
		if (rank==0){
			elap_time = omp_get_wtime();
			res = sequential_gaussjordan(mat, nrows, ncols);
			elap_time = omp_get_wtime() - elap_time;
			print_vect(OUTPUT_FILE, res, nrows);
			print_time(argv[2], elap_time);
			free(res);
			free_mat((void**)mat, nrows);
		}
	}
	else{
		//Definicao do numero de threads e armazenamento do tempo inicial
		n_threads = atoi(argv[1]);
		if (rank==0) elap_time = omp_get_wtime();

		//Atribuicao das linhas para cada processo
		assign_initial_lines(&nrows, &ncols, &job_size, &job_lines_idxs, rank, nproc);
		//Passagem das linhas das matrizes para os procesos apropriados
		if (rank == 0) send_initial_lines(mat, nrows, ncols, nproc, &msgtag);
		else mat = recv_initial_lines(ncols, job_size, &msgtag);

		/*
		Cada processo realiza suas operacoes do gauss-jordan, e para 
		o processo 0 e retornado o vetor com as respostas
		*/
		res = parallel_gaussjordan(mat, nrows, ncols, job_size, job_lines_idxs, nproc, rank, &msgtag, n_threads);
		if (rank == 0){
			elap_time = omp_get_wtime() - elap_time;
			print_vect(OUTPUT_FILE, res, nrows);
			print_time(argv[2], elap_time);
			free(res);
		}
		
		free(job_lines_idxs);
		free_mat((void**)mat, ((rank==0) ? nrows : job_size));
	}

	MPI_Finalize();
	return 0;
}
