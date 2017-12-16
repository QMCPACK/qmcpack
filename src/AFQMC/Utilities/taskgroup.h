
#ifndef AFQMC_TASK_GROUP_H
#define AFQMC_TASK_GROUP_H

#include<vector>
#include<string>
#include<map>
#include<ctime>
#include<sys/time.h>
#include<cstdlib>
#include<ctype.h>
#include<algorithm>
#include<iostream>
#include<ostream>
#include <mpi.h>

#include"AFQMC/config.h"
#include<Message/MPIObjectBase.h>
#include<Message/CommOperators.h>

namespace qmcplusplus
{

namespace afqmc
{

// sets up communicators and task groups
// Various divisions are setup:
//   1. head_of_nodes: used for all shared memory setups
//   2. breaks global comm into groups of ncores_per_TG x nnodes_per_TG
//      and sets up appropriate communicators     
//   Right now does not allow communications outside the TG.
//   This option must be enabled in order to implement algorithms that
//   need to calculate properties that involve all walkers, e.g. pure estimators    
class TaskGroup: public MPIObjectBase, public AFQMCInfo {

  public:

  TaskGroup(Communicate *c, std::string name):MPIObjectBase(c),tgname(name),initialized(false),
     verbose(true)
  {}
  ~TaskGroup() {};

  void setBuffer(SPComplexSMVector* buf) { commBuff = buf; }

  // right now using std::vector and std::string to make the initial implementatino
  //   easier, but this is not efficient and can lead to memory fragmentation for large 
  //   processor counts (e.g. > 10k)
  bool setup(int ncore, int nnode, bool print=false) { 
 
    verbose = print;
    ncores_per_TG = ncore;
    nnodes_per_TG = nnode;
   
    app_log()<<std::endl
             <<"**************************************************************\n"
             <<" Setting up Task Group: " <<name <<std::endl; 
    // do setup based on hostname and rank on myComm  
    std::vector<char> names;
    if(myComm->rank()==0) names.resize(myComm->size()*HOST_NAME_MAX);
    std::vector<char> myname(HOST_NAME_MAX);
    // get hostname
    gethostname(myname.data(),HOST_NAME_MAX);

    myComm->gather(myname,names,0);    
    std::vector<int> data(4);
 
    if(myComm->rank() == 0) {
      // check for consistency of split 
      std::vector< hostinfo > node_map;
      node_map.push_back( hostinfo(myname.data(),0) );
      for(int i=1; i<myComm->size(); i++) {
        int k = look_for_match(node_map,names.data()+i*HOST_NAME_MAX); 
        if( k >= 0 ) { 
          node_map[k].cnt++; 
        } else {
          node_map.push_back( hostinfo( names.data()+i*HOST_NAME_MAX,0) );
        } 
      }

      if( node_map.size()%nnodes_per_TG != 0  ) {
        std::cerr<<"Found " <<node_map.size() <<" nodes. " <<std::endl;
        APP_ABORT(" Error in TaskGroup setup(): Number of nodes found is not divisible by requested number of nodes per Task Group. \n");
      }
      int cpn, ntg;   
      cpn = node_map[0].cnt;
      for(int i=1; i<node_map.size(); i++) {
        if(node_map[i].cnt != cpn) {
          app_error()<<" Error: Inconsistent number of cores in node: " <<i <<std::endl;
          app_error()<<" Expected: " <<cpn+1  <<"  Found: " <<node_map[i].cnt+1 <<std::endl;   
          APP_ABORT(" Error in TaskGroup::setup(): Found inconsistent number of cores in nodes. All nodes must have the same number of cores. \n\n\n");
        }
      } 
      cpn++; 
      if( cpn%ncores_per_TG != 0  ) {
        std::cerr<<"Found " <<cpn <<" cores per node. " <<std::endl;
        APP_ABORT(" Error in TaskGroup setup(): Number of cores per node found is not divisible by requested number of cores in Task Group. \n");
      }
      ntg = (node_map.size()/nnodes_per_TG) * (cpn/ncores_per_TG);
      app_log()<<" Found: " <<node_map.size() <<" nodes, each with: " <<cpn <<" cores. " <<std::endl; 
      app_log()<<" Task Group named: " <<tgname <<" will be split in " <<ntg <<" groups. \n"
               <<" Each group contains " <<nnodes_per_TG <<" nodes * " <<ncores_per_TG <<" cores/node " <<std::endl;  

      // reset cnter
      for(int i=0; i<node_map.size(); i++) node_map[i].cnt = 0; 
      node_map[0].cnt++;
     
      // assign keys
      data[0] = 0;                 // node number: number of the node the current task belongs to 
      data[1] = 0;                 // core rank 
      data[2] = node_map.size();   // number of nodes
      data[3] = cpn;               // number of cores per node
      for(int i=1; i<myComm->size(); i++) {
        int k = look_for_match(node_map,names.data()+i*HOST_NAME_MAX);
        if( k >= 0 ) {
          data[0] = k;
          data[1] = node_map[k].cnt++;
          myComm->send(data.data(),4,i,1010,myComm->getMPI());
        } else {
          APP_ABORT(" Error in TaskGroup::setup(): This should not happen. \n\n\n ");
        }
      }      

      data[0] = data[1] = 0;

    } else {
      // receive key
      MPI_Status st;  
      myComm->recv(data.data(),4,0,1010,myComm->getMPI(),&st);
    }

    tot_nodes = data[2];
    tot_cores = data[3];
    node_number = data[0];
    core_number = data[1];

    // split communicator 
    nrows = tot_cores/ncores_per_TG;
    ncols = tot_nodes/nnodes_per_TG; 
    mycol = node_number/nnodes_per_TG; 
    node_in_TG = node_number%nnodes_per_TG; 
    myrow = core_number/ncores_per_TG; 
    TG_number = mycol + ncols*myrow; 
    number_of_TGs = nrows*ncols;
    myComm->split_comm(TG_number,MPI_COMM_TG);
    MPI_Comm_rank(MPI_COMM_TG,&TG_rank);
    MPI_Comm_size(MPI_COMM_TG,&TG_nproc);
    // assign a unique number to each local group 
    int TG_number_local;
//    TG_number_local = TG_number*nnodes_per_TG + node_in_TG; 
//    myComm->split_comm(TG_number_local,MPI_COMM_TG_LOCAL);
    TG_root = false;
    if(TG_rank==0) TG_root = true;
    core_rank = core_number%ncores_per_TG;
    core_root = (core_rank==0);

    // setup list of roots in a TG, which are the only ones who communicate
    ranks_of_core_roots.reserve(nnodes_per_TG);
    next_core_root = prev_core_root = -1; 
    std::vector<int> tmp(2);
    position_in_ranks_of_core_roots=-1;
    for(int i=0; i<TG_nproc; i++) {
      if( TG_rank == i ) tmp[0] = core_rank;  
      myComm->bcast(tmp.data(),1,i,MPI_COMM_TG); 
      if(tmp[0]==0) ranks_of_core_roots.push_back(i);
      if( core_root && TG_rank == i ) position_in_ranks_of_core_roots = ranks_of_core_roots.size()-1; 
    }
    if(core_root) {
     // std::cout<<myComm->rank() <<" " <<position_in_ranks_of_core_roots <<" " <<ranks_of_core_roots[position_in_ranks_of_core_roots] <<std::endl;
      if( position_in_ranks_of_core_roots < 0 || position_in_ranks_of_core_roots >= nnodes_per_TG ) {
       std::cerr<<" TaskGroup::setup() position_in_ranks_of_core_roots: " <<position_in_ranks_of_core_roots <<std::endl;
        APP_ABORT(" Logic error in TaskGroup::setup(). \n\n\n ");
      }
      if( ranks_of_core_roots[position_in_ranks_of_core_roots] != TG_rank ) {
        std::cerr<<" TaskGroup::setup() ranks_of_core_roots[position_in_ranks_of_core_roots]: " <<ranks_of_core_roots[position_in_ranks_of_core_roots] <<std::endl; 
        APP_ABORT(" Logic error in TaskGroup::setup(). \n\n\n ");
      }
      if( position_in_ranks_of_core_roots == 0 ) {
        next_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots+1];
        prev_core_root = ranks_of_core_roots[nnodes_per_TG-1];
      } else if( position_in_ranks_of_core_roots == nnodes_per_TG-1 ) {
        next_core_root = ranks_of_core_roots[0];
        prev_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots-1];
      } else {
        next_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots+1];
        prev_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots-1];
      }
    }
    app_log()<<"**************************************************************\n";
    initialized=true;
    return true;
  }

  // sets up new TG with global information from previously defined TG 
  bool quick_setup(int ncore, int nnode, int node_number_, int core_number_, int tot_nodes_, int tot_cores_ , bool print=true ) {

  tot_nodes = tot_nodes_;
  tot_cores = tot_cores_;
  node_number = node_number_;
  core_number = core_number_;

//    if(!initialized) {
//      app_error()<<" Error: Call to TaskGroup::quick_setup in uninitialized state. \n";
//      return false; 
//    }

    verbose = print;
    ncores_per_TG = ncore;
    nnodes_per_TG = nnode;

    app_log()<<std::endl
             <<"**************************************************************\n"
             <<" Setting up Task Group: " <<name <<std::endl;

    if( tot_nodes%nnodes_per_TG != 0  ) {
      APP_ABORT(" Error in TaskGroup::quick_setup(): Number of nodes is not divisible by requested number of nodes per Task Group. \n");
      return false;
    }
    if( tot_cores%ncores_per_TG != 0  ) {
      APP_ABORT(" Error in TaskGroup::quick_setup(): Number of cores per node is not divisible by requested number of cores in Task Group. \n");
      return false;
    }

    // split communicator 
    nrows = tot_cores/ncores_per_TG;
    ncols = tot_nodes/nnodes_per_TG; 
    mycol = node_number/nnodes_per_TG; 
    node_in_TG = node_number%nnodes_per_TG; 
    myrow = core_number/ncores_per_TG; 
    TG_number = mycol + ncols*myrow; 
    number_of_TGs = nrows*ncols;
    myComm->split_comm(TG_number,MPI_COMM_TG);
    MPI_Comm_rank(MPI_COMM_TG,&TG_rank);
    MPI_Comm_size(MPI_COMM_TG,&TG_nproc);
    // assign a unique number to each local group 
    int TG_number_local;
//    TG_number_local = TG_number*nnodes_per_TG + node_in_TG; 
//    myComm->split_comm(TG_number_local,MPI_COMM_TG_LOCAL);
    TG_root = false;
    if(TG_rank==0) TG_root = true;
    core_rank = core_number%ncores_per_TG;
    core_root = (core_rank==0);

    if(verbose) {
      app_log()<<" System contains " <<tot_nodes <<" nodes, each with: " <<tot_cores <<" cores. " <<std::endl;
      app_log()<<" Task Group named: " <<tgname <<" will be split in " <<nrows*ncols <<" groups. \n"
               <<" Each group contains " <<nnodes_per_TG <<" nodes * " <<ncores_per_TG <<" cores/node " <<std::endl;
    }

    // setup list of roots in a TG, which are the only ones who communicate
    ranks_of_core_roots.reserve(nnodes_per_TG);
    next_core_root = prev_core_root = -1; 
    std::vector<int> tmp(2);
    position_in_ranks_of_core_roots=-1;
    for(int i=0; i<TG_nproc; i++) {
      if( TG_rank == i ) tmp[0] = core_rank;  
      myComm->bcast(tmp.data(),1,i,MPI_COMM_TG); 
      if(tmp[0]==0) ranks_of_core_roots.push_back(i);
      if( core_root && TG_rank == i ) position_in_ranks_of_core_roots = ranks_of_core_roots.size()-1; 
    }
    if(core_root) {
     // std::cout<<myComm->rank() <<" " <<position_in_ranks_of_core_roots <<" " <<ranks_of_core_roots[position_in_ranks_of_core_roots] <<std::endl;
      if( position_in_ranks_of_core_roots < 0 || position_in_ranks_of_core_roots >= nnodes_per_TG ) {
       std::cerr<<" TaskGroup::setup() position_in_ranks_of_core_roots: " <<position_in_ranks_of_core_roots <<std::endl;
        APP_ABORT(" Logic error in TaskGroup::setup(). \n\n\n ");
      }
      if( ranks_of_core_roots[position_in_ranks_of_core_roots] != TG_rank ) {
        std::cerr<<" TaskGroup::setup() ranks_of_core_roots[position_in_ranks_of_core_roots]: " <<ranks_of_core_roots[position_in_ranks_of_core_roots] <<std::endl; 
        APP_ABORT(" Logic error in TaskGroup::setup(). \n\n\n ");
      }
      if( position_in_ranks_of_core_roots == 0 ) {
        next_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots+1];
        prev_core_root = ranks_of_core_roots[nnodes_per_TG-1];
      } else if( position_in_ranks_of_core_roots == nnodes_per_TG-1 ) {
        next_core_root = ranks_of_core_roots[0];
        prev_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots-1];
      } else {
        next_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots+1];
        prev_core_root = ranks_of_core_roots[position_in_ranks_of_core_roots-1];
      }
    }
    app_log()<<"**************************************************************" <<std::endl;
    initialized=true;
    return true;
  }

  void set_min_max(int min, int max) {
    min_index=min;
    max_index=max;
  }

  void get_min_max(int& min, int& max) const {
    min=min_index;
    max=max_index;
  }

  // over full TG using mpi communicator 
  void barrier() {
    MPI_Barrier(MPI_COMM_TG);
  }

  // over local node using boost sync 
  void local_barrier() {
    //commBuff->barrier();
    MPI_Barrier(MPI_COMM_TG_LOCAL);
  }

  MPI_Comm getTGCOMM() const { return MPI_COMM_TG; }

  MPI_Comm getTGCommLocal() const { return MPI_COMM_TG_LOCAL; }

  void setTGCommLocal(MPI_Comm cm) { MPI_COMM_TG_LOCAL = cm; }

  MPI_Comm getNodeCommLocal() const { return MPI_COMM_NODE_LOCAL; }

  void setNodeCommLocal(MPI_Comm cm) { MPI_COMM_NODE_LOCAL = cm; }

  MPI_Comm getHeadOfNodesComm() const { return MPI_COMM_HEAD_OF_NODES;}

  void setHeadOfNodesComm(MPI_Comm cm) {MPI_COMM_HEAD_OF_NODES = cm;}

  void allgather_TG(std::vector<int>& l, std::vector<int>& g) {
    myComm->allgather(l,g,l.size(),MPI_COMM_TG);
  }  

  // size is in units of ComplexType and represents (walker_size)*(number_of_walkers) 
  void resize_buffer(int& size) 
  {
    std::vector<int> sz(1);
    sz[0]=size;    
    myComm->gmax(sz,MPI_COMM_TG);
    size = sz[0];
    // reset SM is necessary 
    commBuff->resize(size);
  }

  // on entry, nblock has the number of blocks that should be sent 
  // on return, nblock has the number of blocks received
  void rotate_buffer(int& nblock, int block_size)
  {
    int n0 = nblock;
    if(commBuff->size() < nblock*block_size) {
      APP_ABORT(" Error in TaskGroup::rotate_buffer(). Buffer size is too small. \n\n\n ");
    }
    commBuff->barrier();
    if(core_root) {
      local_buffer.resize(commBuff->size());  // this guarantees that I'll be able to receive any message
      if(nnodes_per_TG%2 != 0) {
        APP_ABORT("Error: TaskGroup::rotate_buffer curently implemented for an even number on nodes per task group. Aborting!!! \n\n\n");
      }
      // simple algorithm for now, make this efficient later
// this can be made much faster and efficient 
      if(position_in_ranks_of_core_roots%2==0) {
        MPI_Status status;
        myComm->send(commBuff->values(),n0*block_size,next_core_root,1001,MPI_COMM_TG);
        myComm->recv(local_buffer.data(),local_buffer.size(),prev_core_root,1002,MPI_COMM_TG,&status);
        // assuming doubles for now, FIX FIX FIX
        MPI_Get_count(&status,MPI_DOUBLE,&nblock);
        nblock = nblock/2/block_size; // since I'm communicating std::complex
        std::copy(local_buffer.begin(),local_buffer.begin()+nblock*block_size,commBuff->begin());
      } else {
        MPI_Status status;
        myComm->recv(local_buffer.data(),local_buffer.size(),prev_core_root,1001,MPI_COMM_TG,&status);
        // assuming doubles for now, FIX FIX FIX
        MPI_Get_count(&status,MPI_DOUBLE,&nblock);
        nblock = nblock/(block_size*(sizeof(SPComplexType)/sizeof(double))); 
        myComm->send(commBuff->values(),n0*block_size,next_core_root,1002,MPI_COMM_TG);
        std::copy(local_buffer.begin(),local_buffer.begin()+nblock*block_size,commBuff->begin());
      }
    } 
    commBuff->share(&nblock,1,core_root);
  } 

  int getTotalNodes() const { return tot_nodes; }

  int getTotalCores() const { return tot_cores; }

  int getNodeID() const { return node_number; }

  int getCoreID() const { return core_number; }

  int getCoreRank() const { return core_rank; }

  int getLocalNodeNumber() const { return node_in_TG; }

  int getTGNumber() const { return TG_number; }

  int getNumberOfTGs() const { return number_of_TGs; }

  int getTGRank() const { return TG_rank; }

  int getTGSize() const { return TG_nproc; }

  int getNCoresPerTG() const { return ncores_per_TG; }

  int getNNodesPerTG() const { return nnodes_per_TG; }
 
  void getRanksOfRoots(std::vector<int>& ranks, int& pos ) const { 
    ranks=ranks_of_core_roots; 
    pos=position_in_ranks_of_core_roots; 
  }

  void getSetupInfo(std::vector<int>& data) const
  {  
    data.resize(5);
    data[0]=node_number;
    data[1]=core_number;
    data[2]=tot_nodes;
    data[3]=tot_cores;
    data[4]=ncores_per_TG; 
  }
 
  // must be setup externally to be able to reuse between different TG 
  SPComplexSMVector* commBuff;  

  std::string tgname;

  bool verbose;
  bool initialized;
  
  int node_number, core_number, tot_nodes, tot_cores;
  int TG_number; 
  int number_of_TGs; 
  // TGs are defined in a 2-D framwork. Rows correspond to different groups in a node 
  // Cols correspond to division of nodes into groups. Every MPI task belongs to a specific
  // TG given by the myrow and mycol integer.      
  int myrow, mycol, nrows, ncols, node_in_TG;
  bool TG_root;            // over full TG
  int TG_rank, TG_nproc;   // over full TG, notice that nproc = ncores_per_TG * nnodes_per_TG 
  bool core_root;          // over local node
  int core_rank;           // over local node 
  std::vector<int> ranks_of_core_roots;
  int position_in_ranks_of_core_roots; 
  int next_core_root, prev_core_root; // only meaningful at core_root processes  
  MPI_Comm MPI_COMM_TG;   // Communicator over all cores in a given TG 
  MPI_Comm MPI_COMM_TG_LOCAL;   // Communicator over all cores in a given TG that reside in the given node 
  MPI_Comm MPI_COMM_NODE_LOCAL; // Communicator over all cores of a node. Must be created externally. Same above
  MPI_Comm MPI_COMM_HEAD_OF_NODES;  // deceiving name for historical reasons, this is a split of COMM_WORLD over core_number. 
  std::vector<SPComplexType> local_buffer;

  int ncores_per_TG;  // total number of cores in all nodes must be a multiple 
  int nnodes_per_TG;  // total number of nodes in communicator must be a multiple  
  int min_index, max_index;  

  struct hostinfo { 
    char id[HOST_NAME_MAX];
    int cnt; 
    hostinfo( char* name, int n) {
      std::memcpy(id,name,HOST_NAME_MAX*sizeof(char));
      cnt=n;
    } 
  };
  int look_for_match(std::vector<hostinfo>& v, char* n) {
    for(int i=0; i<v.size(); i++) {
      if( std::strcmp(v[i].id,n)==0 ) 
        return i;
    }
    return -1;
  }

};

}
 
}


#endif
