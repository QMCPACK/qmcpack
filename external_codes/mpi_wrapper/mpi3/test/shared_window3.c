#if COMPILATION_INSTRUCTIONS
mpicc -Wfatal-errors $0 -o $0x.x &&  mpirun `#srun -p pdebug` -n 2 $0x.x $@ && rm -f $0x.x; exit
#endif

#include<mpi.h>
#include<assert.h>
#include<stdio.h>

int main(int argc, char* argv[]){
  MPI_Init(&argc, &argv);
  MPI_Comm node;
  int ss = MPI_Comm_split_type(MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &node);
  assert(ss == MPI_SUCCESS);
  
  int rank = -1;
  MPI_Comm_rank(node, &rank);


  #define N 65536
  
  MPI_Win win[N];				     
  void* base_ptr[N];
  int s[N];
  MPI_Aint size[N];
  int a[N];
  char* ptr[N];

  for(int i = 0; i != N; ++i){
    base_ptr[i] = 0;
    printf("before alloc shared # %i\n", i); fflush(stdout);
	int sam = MPI_Alloc_mem((rank==0?99:0)*sizeof(char), MPI_INFO_NULL, &base_ptr[i]);
    assert(sam == MPI_SUCCESS);
//    s[i] = MPI_Win_create(base_ptr[i], (rank==0?99:0)*sizeof(char), 1, MPI_INFO_NULL, node, &win[i]);
    s[i] = MPI_Win_allocate_shared((rank==0?99:0)*sizeof(char), 1, MPI_INFO_NULL, node, &base_ptr[i], &win[i]);
    assert(s[i] == MPI_SUCCESS);
    if(rank == 0) assert(base_ptr[i] != 0);
    size[i] = -1;		
    a[i] = -1;
    ptr[i] = 0;
//  MPI_Win_get_attr(&win, MPI_WIN_BASE, &ptr[i], 0);
    MPI_Win_shared_query(win[i], 0, &size[i], &a[i], &ptr[i]);
    assert(&ptr[i] != 0);	
    printf("alloc # %i ok\n", i); fflush(stdout);		
  }

  MPI_Barrier(node);
  printf("barrier at %i\n", rank); fflush(stdout);

  for(int i = 0; i != N; ++i){
    if(rank == 0) ptr[i][3] = 'z'; 
    MPI_Barrier(node);
    assert(ptr[i][3] == 'z');
  }

  for(int i = 0; i != N; ++i){
    if(win[i] != MPI_WIN_NULL) MPI_Win_free(&win[i]);
  }

  MPI_Finalize();
  return 0;
}

