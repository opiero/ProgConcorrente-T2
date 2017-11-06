#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <mpi.h>

int main (int argc, char * argv[]) {

  int my_rank = -1;
  int num_proc;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);

  


  MPI_Finalize();

  return EXIT_SUCCESS;
}
