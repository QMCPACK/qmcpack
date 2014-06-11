
#include "include/process_file.h"

void process_file(char *dir, char *filename, struct BP_FILE * fh){
  char fullname[200];
  sprintf(fullname, "%s%s", dir, filename);
  printf("%s\n", fullname);
  bp_open (fullname, MPI_COMM_WORLD, fh);
}

