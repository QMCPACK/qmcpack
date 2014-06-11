
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
      struct BP_FILE *fh = (struct BP_FILE *)malloc(sizeof(struct BP_FILE));
      fh_v[i] = fh;
      fh_v[i]->b = (struct adios_bp_buffer_struct_v1 *)malloc(sizeof(struct adios_bp_buffer_struct_v1));
      fh_v[i]->pgs_root = (struct bp_index_pg_struct_v1 *)malloc(sizeof(struct bp_index_pg_struct_v1));
      fh_v[i]->vars_root = (struct adios_index_var_struct_v1 *)malloc(sizeof(struct adios_index_var_struct_v1));
      fh_v[i]->attrs_root = (struct adios_index_attribute_struct_v1 *)malloc(sizeof(struct adios_index_attribute_struct_v1));
      fh_v[i]->gvar_h = (struct BP_GROUP_VAR *)malloc(sizeof(struct BP_GROUP_VAR));
      fh_v[i]->gattr_h = (struct BP_GROUP_ATTR *)malloc(sizeof(struct BP_GROUP_ATTR));
      process_file(argv[1], names[i], fh_v[i]);
    }
  }

  //take down the fh_v structure

  //finalize MPI
  MPI_Finalize();

  return 0;
}
