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
#define TIME_FILE "time.txt"

//Nao pode ser <= 0
#define ALLOC_INIT_SIZE 8
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



/*
Le o vetor b de Ax = b do arquivo filename
*/
element * read_vect(char * filename, int * nrows)
{
	/*
	Alloca o vetor com um tamanho inicial. O tamanho do vetor cresce 
	dinamicamente e exponencialmente para reduzir o numero de reallocs
	*/
	(*nrows) = 0;
	int alloced_row = ALLOC_INIT_SIZE;
	element * vect = (element*) malloc(alloced_row*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * vect_file = open_file(filename, "r");
	//Enquanto ainda existir uma nova linha
	while(getline(&line, &linesize, vect_file) != -1){
		//Ignora linhas vazias
		if (strcmp(line, "\n") != 0){
			//Incrementar o tamanho alocado do vetor caso necessario
			while ((*nrows) >= alloced_row){
				// *=2 ao inves de ++ para reduzir o numero de reallocs
				alloced_row*=2;
				vect = (element*) realloc(vect, alloced_row*sizeof(element));
			}

			//Incrementa o numero de elementos do vetor caso haja um elemento na linha
			if (sscanf(line, ELEMENT_READ_MASK, vect + (*nrows)) == 1) (*nrows)++;
			free(line);
			//line=NULL e necessario para o getline()
			line = NULL;
		}
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



/*
Le a matriz A de Ax = b do arquivo filename
*/
element ** read_mat(char * filename, int * nrows, int * ncols)
{
	/*
	Alloca a matriz com um tamanho inicial. O tamanho da matriz cresce 
	dinamicamente e exponencialmente para reduzir o numero de reallocs
	*/
	(*nrows) = (*ncols) = 0;
	int alloced_col = ALLOC_INIT_SIZE, alloced_row = ALLOC_INIT_SIZE;
	element ** mat = (element**) malloc(alloced_row*sizeof(element*));
	(*mat) = (element*) malloc(alloced_col*sizeof(element));

	size_t linesize = 0;
	char * line = NULL;
	FILE * mat_file = open_file(filename, "r");
	if (getline(&line, &linesize, mat_file) == -1) printf("No lines on file\n");
	char * aux_str = strtok(line, DELIMITERS);
	//Leitura da primeira linha para definir o numero de colunas da matriz
	while(aux_str != NULL){
		//Aumenta o numero de colunas exponencialmente se necessario.
		while ((*ncols) >= alloced_col){
			alloced_col*=2;
			(*mat) = (element*) realloc((*mat), alloced_col*sizeof(element));
		}

		//Le no numero
		if (sscanf(aux_str, ELEMENT_READ_MASK, (*mat) + (*ncols)) == 1) (*ncols)++;
		//Pega o proximo numero na linha
		aux_str = strtok(NULL, DELIMITERS);
	}
	//Alloca 1 coluna a mais para inserir o vetor b posteriormente
	if (alloced_col < (*ncols)+1) (*mat) = (element*) realloc((*mat), ((*ncols)+1)*sizeof(element));
	free(line);

	line = NULL;
	(*nrows) = 1;
	//Leitura do resto da matriz
	while (getline(&line, &linesize, mat_file) != -1){
		//Ignora linhas vazias
		if (strcmp(line, "\n") != 0){
			//Aumenta o numero de linhas exponencialmente se necessario.
			while ((*nrows) >= alloced_row){
				alloced_row*=2;
				mat = (element**) realloc(mat, alloced_row*sizeof(element*));
			}
			//Alloca 1 coluna a mais para inserir o vetor b posteriormente
			mat[(*nrows)] = (element*) malloc(((*ncols)+1)*sizeof(element));

			int cur_col = 0;
			aux_str = strtok(line, DELIMITERS);
			while(aux_str != NULL) {
				//Le o numero atual e pega o proximo
				if (sscanf(aux_str, ELEMENT_READ_MASK, mat[(*nrows)] + cur_col) == 1) cur_col++;
				aux_str = strtok(NULL, DELIMITERS);
			}

			(*nrows)++;
			free(line);
			//line=NULL e necessario para o getline()
			line = NULL;
		}
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



/*
Depois que terminou a computacao, os processos
mandam ao root o valor da coluna de resposta deles, junto
com o indice que aquela resposta pertence na matriz original 
*/
void send_final_answer(element ** mat, int ncols, int job_size, int * job_lines_idxs, int * init_tag)
{
	//Criacao de um vetor auxiliar para mandar com apenas 1 send
	int i;
	element * res = (element*) malloc(job_size*sizeof(element));
	for (i=0; i<job_size; i++) res[i] = mat[i][ncols-1];

	//Envio dos indices na matriz original e da resposta 
	MPI_Send(job_lines_idxs, job_size, MPI_INT, 0, (*init_tag), MPI_COMM_WORLD);
	MPI_Send(res, job_size, ELEMENT_MPI, 0, (*init_tag)+1, MPI_COMM_WORLD);
	(*init_tag)+=1;
	free(res);
}

/*
Depois que terminou a computacao, o processo root
recebe o valor da coluna de resposta dos outros processos,
junto com o indice que aquela resposta pertence na matriz original 
*/
element * recv_final_answer(int nrows, int nproc,  int * init_tag)
{
	//Recebe de cada processo quais indices da matriz original ele tem
	MPI_Status status;
	int i, j, cur_size = (nrows / nproc);	
	int ** idxs = (int**) malloc((nproc-1)*sizeof(int*));
	for (i=1; i<nproc; i++){
		//Caso especial para o ultimo processo
		if (i==nproc-1) cur_size += nrows % nproc;
		idxs[i-1] = (int*) malloc(cur_size*sizeof(int));
		MPI_Recv(idxs[i-1], cur_size, MPI_INT, i, (*init_tag), MPI_COMM_WORLD, &status);
	}

	/*
	Recebe de cada processo o pedaco da resposta que ele tem, e
	armazena na resposta final ja na posicao correta
	*/
	cur_size = (nrows / nproc);	
	element * aux_vect = (element*) malloc((cur_size + (nrows % nproc))*sizeof(element));
	element * res = (element*) malloc(nrows*sizeof(element));
	for (i=1; i<nproc; i++){
		//Caso especial para o ultimo processo
		if (i==nproc-1) cur_size += nrows % nproc;
		MPI_Recv(aux_vect, cur_size, ELEMENT_MPI, i, (*init_tag)+1, MPI_COMM_WORLD, &status);
		//Armazena na resposta final na posicao correta, dada pelo vetor de indices rebebido
		for (j=0; j<cur_size; j++) res[idxs[i-1][j]] = aux_vect[j];
	}
	(*init_tag)+=1;

	free_mat((void**)idxs, nproc-1);
	free(aux_vect);
	return res;
}



/*
Realiza o algoritimo em paralelo com nproc processos
e nthreads threads por processo
*/
element * parallel_gaussjordan(element ** mat, int nrows, int ncols, int job_size, int * job_lines_idxs, int nproc, int rank, int * init_tag, int nthreads)
{
	/*
	Declaracao das variaveis. min recebe a menor dimensao para realizar o problema para matrizes nao quadradas
	ncols-1 pois precisa desconsiderar a ultima coluna com o vetor b de Ax = b
	*/
	int i, j, k, id, min = ((ncols-1) <= nrows) ? (ncols-1) : nrows;
	int local_max_idx, global_max_idx, global_max_proc;
	int * recv_idx = NULL, * aux_max_idx_vect = NULL;
	element local_max, global_max, rat;
	element * global_max_line = NULL, * recv_max = NULL, * aux_max_vect = NULL;
	
	//Allocacao dos vetores
	global_max_line = (element*) malloc(ncols*sizeof(element));
	aux_max_vect = (element*) malloc(nthreads*sizeof(element));
	aux_max_idx_vect = (int*) malloc(nthreads*sizeof(int));
	if (rank==0){
		recv_idx = (int*) malloc(nproc*sizeof(int));
		recv_max = (element*) malloc(nproc*sizeof(element));
	}

	//Especificao o numero de threads
	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
	//Iterar nas colunas da matriz
	for (j=0; j<min; j++){
		/*
		Encontra o maximo valor e indice daquela coluna para cada processo
		utilizando diversas threads. Reduce de omp nao funciona com double, entao
		foi utilizado um vetor auxiliar para armazenar o maximo de cada threads,
		e depois foi feito o maximo desse vetor
		*/
		local_max = 0;
		local_max_idx = -1;
		memset(aux_max_vect, 0, nthreads*sizeof(element));
		memset(aux_max_idx_vect, -1, nthreads*sizeof(int));
		//Maximo por threads
		#pragma omp parallel for shared(aux_max_vect, aux_max_idx_vect) private(i, id)
		for (i=0; i<job_size; i++)
		{
			id = omp_get_thread_num();
			if ((job_lines_idxs[i] >= j) && (MOD(mat[i][j]) > MOD(aux_max_vect[id]))){
				aux_max_vect[id] = mat[i][j];
				aux_max_idx_vect[id] = i;
			}
		}
		//Maximo do processo
		for (i=0; i<nthreads; i++){
			if (MOD(aux_max_vect[i]) >= MOD(local_max)){
				local_max = aux_max_vect[i];
				local_max_idx = aux_max_idx_vect[i];
			}
		}

		/*
		O processo 0 recebe o maximo e o indice do maximo de cada um dos outros processos,
		e calcula o maximo global da matriz. Nao e necessario usar multithread ou reduction
		pois o numero de processos é pequeno, entao nao haveria ganho significativo em desempenho
		*/
		MPI_Gather(job_lines_idxs + local_max_idx, 1, MPI_INT, recv_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Gather(&local_max, 1, ELEMENT_MPI, recv_max, 1, ELEMENT_MPI, 0, MPI_COMM_WORLD);

		//Maximo da matriz
		if (rank==0){
			global_max = 0;
			global_max_proc = global_max_idx = -1;
			for (i=0; i<nproc; i++){
				if (MOD(recv_max[i]) > MOD(global_max)) {
					global_max = recv_max[i];
					global_max_idx = recv_idx[i];
					global_max_proc = i;
				}
			}
		}
		
		//Depois que o processo 0 calculou o maximo, ele faz um broadcast de quem e esse maximo
		MPI_Bcast(&global_max, 1, ELEMENT_MPI, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_max_idx, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&global_max_proc, 1, MPI_INT, 0, MPI_COMM_WORLD);

		/*
		Troca de linhas, onde o processo que tem a linha de mesmo indice que a coluna da
		iteracao atual troca de posicao com a linha do maximo atual, e o processo que tem a
		linha do maximo faz o contrario. Essa troca e feita simplesmente alterando o
		vetor que mapeia a linha do processo pra linha da matriz original
		*/
		#pragma omp parallel for private(i)
		for (i=0; i<job_size; i++)
		{
			if (job_lines_idxs[i]==j) job_lines_idxs[i] = global_max_idx;
		}
		if (rank == global_max_proc){
			job_lines_idxs[local_max_idx] = j;
			//Armazenamento da linha maxima para o broadcast dela
			memcpy(global_max_line, mat[local_max_idx], ncols*sizeof(element));
		}

		//Broadcast da linha pivo da iteracao atual
		MPI_Bcast(global_max_line, ncols, ELEMENT_MPI, global_max_proc, MPI_COMM_WORLD);

		//Multiplicacao e soma das linhas para zerar a coluna, menos a diagonal principal
		#pragma omp parallel for private(i, rat)
		for (i=0; i<job_size; i++)
		{
			if (job_lines_idxs[i] != j){
				rat = mat[i][j]/global_max;
				for (k=0; k<ncols; k++) mat[i][k] -= rat*global_max_line[k];
			}
		}

		//Sincronizacao entre os processos
		MPI_Barrier(MPI_COMM_WORLD);
	}
	//Divisao da coluna resposta pela diagonal principal para obter a resposta final
	#pragma omp parallel for private(i)
	for (i=0; i<job_size; i++)
	{
		mat[i][ncols-1] /= mat[i][job_lines_idxs[i]];
	}

	free(aux_max_idx_vect);
	free(aux_max_vect);
	free(global_max_line);
	if (rank==0){
		free(recv_idx);
		free(recv_max);
	}

	/*
	Depois de ter sido computado os resultados, os processos enviam ao processo
	0 o pedaco da resposta que eles tem
	*/
	if (rank == 0){
		element * res = recv_final_answer(nrows, nproc, init_tag);
		for (i=0; i<job_size; i++) res[job_lines_idxs[i]] = mat[i][ncols-1];
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
	int nrows = 0, ncols = 0, job_size = 0, nthreads = 1;
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
	
	//Definicao do numero de threads e armazenamento do tempo inicial
	if (argc >= 2) nthreads = atoi(argv[1]);
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
	res = parallel_gaussjordan(mat, nrows, ncols, job_size, job_lines_idxs, nproc, rank, &msgtag, nthreads);
	if (rank == 0){
		elap_time = omp_get_wtime() - elap_time;
		print_vect(OUTPUT_FILE, res, nrows);
		print_time(TIME_FILE, elap_time);
		free(res);
	}
	
	free(job_lines_idxs);
	free_mat((void**)mat, ((rank==0) ? nrows : job_size));

	MPI_Finalize();
	return 0;
}
