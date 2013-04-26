#include "DataSporkAnalysis.h"

int main(int ac, char* av[])
{
  DataSporkAnalysis ds;
  int dummy = ds.getOptions(ac,av);
  ds.execute();
  return dummy;
}
