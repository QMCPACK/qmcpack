#ifdef HAVE_ADIOS
#include <Configuration.h>
#include <adios_read.h>
#include "Particle/Walker.h"
#include "Particle/ParticleSet.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

class AdiosWalkerInput
{
public :

  std::string FileRoot;
  typedef ParticleSet::Walker_t::ParticlePos_t R_t;
  int particle_num;
  MCWalkerConfiguration targetW;
  Communicate* myComm;

  AdiosWalkerInput(MCWalkerConfiguration& w, Communicate* c);
  ~AdiosWalkerInput();

  bool put(std::vector<xmlNodePtr>& wset);

  void read(int nprocs, ADIOS_FILE* file_handle, int& walker_win, vector<int>& nw);
  void append_walkers(R_t& walker_buff,int, int);
  void setMCWalker(vector<int> nw);

  string getFileRoot();
};
}

#endif


















