
#include "include/JabberWocky.h"
#include "include/process_file.h"
#include "include/bp_open.h"

#define MAX_FILE 100

void usage(){
  printf("Usage: JabberWocky <dir>\n");
  exit(0);
}

int main(int argc, char ** argv){
  DIR           *pDir;
  struct dirent *pDirent;
  int total_file;
  
  //process parameter
  if(argc < 2) 
    usage();

  //initialize MPI
  int rank, size;
  MPI_Init (&argc, &argv);  /* starts MPI */
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);  /* get current process id */
  MPI_Comm_size (MPI_COMM_WORLD, &size);  /* get number of processes */

  if(rank ==0){
    pDir = opendir (argv[1]);
    if (pDir == NULL) {
      printf ("Cannot open directory '%s'\n", argv[1]);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }

    //gather all files
    char names[MAX_FILE][200];
    int i = 0;
    while ((pDirent = readdir(pDir)) != NULL) {
      sprintf(names[i], "%s", pDirent->d_name);
      i++;
    }
    total_file = i;
    closedir (pDir);


    //process each file
    struct BP_FILE * fh_v[MAX_FILE];
    for(i = 0; i<total_file; i++){
      if(strcmp(names[i], ".")==0 || strcmp(names[i], "..")==0) continue;
      struct BP_FILE * fh;
      fh_v[i] = fh;
      process_file(argv[1], names[i], fh_v[i]);
    }
  }

  //finalize MPI
  MPI_Finalize();

  return 0;
}
