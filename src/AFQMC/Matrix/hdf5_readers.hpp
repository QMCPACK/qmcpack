//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
////
//// File developed by: 
////
//// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////////////


#include<cassert>
#include<complex>
#include<cstdlib>
#include<algorithm>
#include<utility>
#include<vector>
#include<numeric>
#if defined(USE_MPI)
#include<mpi.h>
#endif

#include "Configuration.h"
#include "Utilities/UtilityFunctions.h"
#include "io/hdf_archive.h"
#include "Message/CommOperators.h"

#include<AFQMC/config.0.h>

#define MORE  1555
#define LAST  2777

namespace qmcplusplus
{

namespace afqmc
{


// look for a dataset call ..._descriptor, which will contain integer types 
// describing the properties of the datasets, e.g.
// - transposed?
// - sorted?
// - etc.
//
// If this dataset is not found, then assume generic
//
// SpMatrix_Partition both partitions the matrix into local segments belonging to a 
// given node, as well as applies cutoffs to values.
//
// For good efficiency, it is possible to partition the local block 
// (obtained through split.local(i,j,val))
// among equivalent nnodes (same node id within TG and using split.local_split(i,j,val)) 
// and read/sort only this local segment.
// After the call to compress (which only does it over the local segment),
// assemble the full local partition with a call to Allgather with the corresponding
// nodes in COMM_WORLD.
//   -- minimizes sorting
//
// Algorith 1:
//  1. head_of_nodes read a block from the hdf5 in parallel and sequentially
//  2. loop over head_of_nodes, each one 
//      i. bcasts his block to the group
//      ii. from each bcast, adds to the local list if split.local(i,j)==true
//  *** overlap the calls to bcast (with Ibcast) and the block processing
//  3. after read, compress
//  4. after compress, gather segments from other equivalent nodes 
//  5. If memory is not a problem, further speedup can be obtained by having
//     other cores also read/bcast/process blocks.   
//       --  possible downside:  1) memory, 2) I/O bottleneck (e.g. BGQ) 
template<class SpMatrix,
         class SpMatrix_Partition,
         class task_group
         >
inline bool multiple_reader_hdf5_SpMat(SpMatrix& SpM, SpMatrix_Partition& split, hdf_archive& dump, std::string name, int nblk, task_group& TG, int n_working_cores)
{

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID() , coreid = TG.getCoreID();
  int node_rank;

  // process hdf5 file
  if( coreid < n_working_cores ) {

    bool need_lock = n_working_cores>1;
    MPI_Comm_rank(TG.getHeadOfNodesComm(),&node_rank); // this should be equal to nodeid   
    int fct = 1;
#if defined(QMC_COMPLEX)
    fct=2;
#endif

    std::vector<int> block_size(nblk);
    if(!dump.read(block_size,name+std::string("_block_sizes"))) {
      app_error()<<" Error in multiple_reader_hdf5_SpMat: Problems reading ***_block_sizes dataset. \n";
      return false;
    }

    int ntmax = *std::max_element(block_size.begin(), block_size.end());
    std::vector<IndexType> ivec, ivec2;
    std::vector<ValueType> vvec, vvec2;
    ivec.reserve(2*ntmax);
    vvec.reserve(ntmax);    
    ivec2.reserve(2*ntmax);
    vvec2.reserve(ntmax);    

    std::vector<int> pool_dist;
    // divide blocks over n_working_cores groups
    FairDivide(nblk,n_working_cores,pool_dist);

    int first_local_block = pool_dist[coreid];
    int last_local_block = pool_dist[coreid+1];
    int nbtot = last_local_block-first_local_block;
    int niter = nbtot/nnodes + std::min(nbtot%nnodes,1);

    for(int iter=0; iter<niter; iter++) { 
      int myblock_number = first_local_block + nnodes*iter + nodeid; 
      int first_block = first_local_block + nnodes*iter; 
      int last_block = std::min(first_block+nnodes,last_local_block);
      if(myblock_number < last_local_block) {
        ivec.resize(2*block_size[myblock_number]);
        vvec.resize(block_size[myblock_number]);
        if(!dump.read(ivec,name+std::string("_index_")+std::to_string(myblock_number))) {
          app_error()<<" Error in multiple_reader_hdf5_SpMat: Problems reading ***_index_" <<myblock_number <<" dataset. \n";
          return false;
        }
        if(!dump.read(vvec,name+std::string("_vals_")+std::to_string(myblock_number))) {
          app_error()<<" Error in multiple_reader_hdf5_SpMat: Problems reading ***_vals_" <<myblock_number <<" dataset. \n";
          return false;
        }
      }    
      for(int k=first_block,ipr=0; k<last_block; k++,ipr++) {
        ivec2.resize(2*block_size[k]);
        vvec2.resize(block_size[k]);
        if(ipr==node_rank) {
          assert(myblock_number==k);
          std::copy(ivec.begin(),ivec.end(),ivec2.begin());
          MPI_Bcast(ivec2.data(), ivec2.size(), MPI_INT, ipr, TG.getHeadOfNodesComm() );
          std::copy(vvec.begin(),vvec.end(),vvec2.begin());
          MPI_Bcast(vvec2.data(), fct*vvec2.size(), MPI_DOUBLE, ipr, TG.getHeadOfNodesComm() );
        } else {
          MPI_Bcast(ivec2.data(), ivec2.size(), MPI_INT, ipr, TG.getHeadOfNodesComm() );
          MPI_Bcast(vvec2.data(),fct*vvec2.size(), MPI_DOUBLE, ipr, TG.getHeadOfNodesComm() );
        }

        for(int ik=0, ikend=block_size[k]; ik<ikend; ik++) 
          if(split.local(ivec2[2*ik],ivec2[2*ik+1], vvec2[ik])) 
            SpM.add(ivec2[2*ik],ivec2[2*ik+1], static_cast<SPValueType>(vvec2[ik]),need_lock); 

      }
    }

  } 

  // compress matrix
  SpM.compress(TG.getNodeCommLocal());

  // gather segment 
  //std::vector<int> counts, displ;
  //if(coreid == 0) {
  //  counts.resize(nnodes);
  //  displ.resize(nnodes+1);
  //}
 
  return true;
}

template<class SpMatrix,
         class SpMatrix_Partition,
         class task_group
         >
inline bool single_reader_hdf5_SpMat(SpMatrix& SpM, SpMatrix_Partition& split, hdf_archive& dump, std::string name, int nblk, task_group& TG)
{
  return true;
}


template<class SpMatrix_Partition,
         class IVec,
         class task_group
         >
inline bool multiple_reader_count_entries(hdf_archive& dump, SpMatrix_Partition& split, std::string name, int nblk, bool byRow, IVec& counts, task_group& TG, int n_working_cores)
{
  double cut = split.getCutoff();
  int dim;
  if(byRow) 
    std::tie(dim,std::ignore) = split.getDimensions();
  else 
    std::tie(std::ignore,dim) = split.getDimensions();
  counts.resize(dim);
  std::fill(counts.begin(),counts.end(),0); 
  assert(sizeof(counts[0]) == sizeof(int));

  int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID() , coreid = TG.getCoreID();
  int node_rank;

  // process hdf5 file
  if( coreid < n_working_cores ) {

    std::vector<int> block_size(nblk);
    if(!dump.read(block_size,name+std::string("_block_sizes"))) {
      app_error()<<" Error in multiple_reader_count_entries: Problems reading ***_block_sizes dataset. \n";
      return false;
    }

    int ntmax = *std::max_element(block_size.begin(), block_size.end());
    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(2*ntmax);
    vvec.reserve(ntmax);
  
    std::vector<int> pool_dist;
    // divide blocks over n_working_cores groups
    FairDivide(nblk,n_working_cores,pool_dist); 

    int first_local_block = pool_dist[coreid];
    int last_local_block = pool_dist[coreid+1];

    for(int ib=first_local_block; ib<last_local_block; ib++) {
      if( ib%nnodes != nodeid ) continue;
      ivec.resize(2*block_size[ib]);
      vvec.resize(block_size[ib]);
      if(!dump.read(ivec,name+std::string("_index_")+std::to_string(ib))) {
        app_error()<<" Error in multiple_reader_count_entries: Problems reading ***_index_" <<ib <<" dataset. \n";
        return false;
      }
      if(!dump.read(vvec,name+std::string("_vals_")+std::to_string(ib))) {
        app_error()<<" Error in multiple_reader_count_entries: Problems reading ***_vals_" <<ib <<" dataset. \n";
        return false;
      } 
      if(byRow) { 
        for(int ik=0, ikend=block_size[ib]; ik<ikend; ik++) {
          if(std::abs(vvec[ik]) > cut) {
            assert(ivec[2*ik] < dim);
            counts[ivec[2*ik]]++;
          }
        }
      } else {
        for(int ik=0, ikend=block_size[ib]; ik<ikend; ik++) {
          if(std::abs(vvec[ik]) > cut) {
            assert(ivec[2*ik+1] < dim);
            counts[ivec[2*ik+1]]++;
          }
        }
      }
    }
  }

  IVec a(counts);
  //  FIX FIX FIX: use image communicator, not COMM_WORLD
  MPI_Allreduce(a.data(),counts.data(),counts.size(),MPI_INT,MPI_SUM,MPI_COMM_WORLD);
  return true;
}

template<class SpMatrix_Partition,
         class IVec,
         class task_group
         >
inline bool single_reader_count_entries(hdf_archive& dump, SpMatrix_Partition& split, std::string name, int nblk, bool byRow, IVec& counts, task_group& TG)
{
  if(TG.getCoreID() != 0)
    return true;

  int dim;
  if(byRow) 
    std::tie(dim,std::ignore) = split.getDimensions();
  else 
    std::tie(std::ignore,dim) = split.getDimensions();
  double cut = split.getCutoff();
  counts.resize(dim);
  std::fill(counts.begin(),counts.end(),0); 

  std::vector<int> block_size(nblk);
  if(!dump.read(block_size,name+std::string("_block_sizes"))) {
    app_error()<<" Error in single_reader_count_entries: Problems reading ***_block_sizes dataset. \n";
    return false;
  }

  int ntmax = *std::max_element(block_size.begin(), block_size.end());
  std::vector<IndexType> ivec;
  std::vector<ValueType> vvec;
  ivec.reserve(2*ntmax);
  vvec.reserve(ntmax);

  for(int ib=0; ib<nblk; ib++) { 

    ivec.resize(2*block_size[ib]);
    vvec.resize(block_size[ib]);
    if(!dump.read(ivec,name+std::string("_index_")+std::to_string(ib))) {
      app_error()<<" Error in single_reader_count_entries: Problems reading ***_index_" <<ib <<" dataset. \n";
      return false;
    }
    if(!dump.read(vvec,name+std::string("_vals_")+std::to_string(ib))) {
      app_error()<<" Error in single_reader_count_entries: Problems reading ***_vals_" <<ib <<" dataset. \n";
      return false;
    }
    if(byRow) {
      for(int ik=0, ikend=block_size[ib]; ik<ikend; ik++) {
        if(std::abs(vvec[ik]) > cut) {
          assert(ivec[2*ik] < dim);
          counts[ivec[2*ik]]++;
        }
      }
    } else {
      for(int ik=0, ikend=block_size[ib]; ik<ikend; ik++) {
        if(std::abs(vvec[ik]) > cut) {
          assert(ivec[2*ik+1] < dim);
          counts[ivec[2*ik+1]]++;
        }
      }
    }
  }
 
  return true;
}

template<class SpMatrix,
         class SpMatrix_Partition,
         class task_group
         >
inline bool read_hdf5_SpMat(SpMatrix& SpM, SpMatrix_Partition& split, hdf_archive& dump, std::string name, int nblk, task_group& TG, bool parallel_read=true, int n_working_cores=1)
{
  assert(n_working_cores>=1);
  if(parallel_read)
    return multiple_reader_hdf5_SpMat(SpM,split,dump,name,nblk,TG,n_working_cores);
  else
    return single_reader_hdf5_SpMat(SpM,split,dump,name,nblk,TG);
}

template<class SpMatrix_Partition,
         class IVec,
         class task_group
         >
inline bool count_entries_hdf5_SpMat(hdf_archive& dump, SpMatrix_Partition& split, std::string name, int nblk, bool byRow, IVec& counts, task_group& TG, bool parallel_read=true, int n_working_cores=1)
{
  assert(n_working_cores>=1);
  if(parallel_read)
    return multiple_reader_count_entries(dump,split,name,nblk,byRow,counts,TG,n_working_cores);
  else
    return single_reader_count_entries(dump,split,name,nblk,byRow,counts,TG);
}

template<class SpMatrix,
         class task_group 
         >
inline std::tuple<int,int> write_hdf5_SpMat(SpMatrix& SpM, hdf_archive& dump, std::string name, int nterms_per_blk, task_group& TG)
{
  std::size_t nblocks=0,ntot=0;
  if(TG.getTGNumber() > 0)
    return std::make_tuple(nblocks,ntot);
  int nnodes_per_TG = TG.getNNodesPerTG();
  int node_number = TG.getLocalNodeNumber();
  int core_rank = TG.getCoreRank();
  if(core_rank==0 && node_number==0) {

    // get ranks of node heads on TG
    std::vector<int> ranks;
    int pos=0;
    if(nnodes_per_TG>1) TG.getRanksOfRoots(ranks,pos);
    else {
      ranks.push_back(0);
    }

    std::size_t sz = SpM.size();
    int nb_loc = int( sz/std::size_t(nterms_per_blk) + std::min(sz%std::size_t(nterms_per_blk),std::size_t(1)));
    if(nnodes_per_TG>1) {
      int nt;
      MPI_Reduce(&nb_loc,&nt,1,MPI_INT,MPI_SUM,0,TG.getTGCOMM());
      nblocks = std::size_t(nt);
    } else {
      nblocks = nb_loc; 
    }

    std::vector<int> bsize;
    bsize.reserve(nblocks);

    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(2*nterms_per_blk);
    vvec.reserve(nterms_per_blk);

    typename SpMatrix::int_iterator col = SpM.cols_begin();
    typename SpMatrix::int_iterator row = SpM.rows_begin();
    typename SpMatrix::iterator val = SpM.vals_begin();

    // do mine first
    int cnt=0, nterms; 
    IndexType nblk=0; 
    for(int i=0; i<nb_loc; i++) {
      if( (sz-std::size_t(cnt)) < std::size_t(nterms_per_blk) )
        nterms = int(sz-std::size_t(cnt));
      else
        nterms = nterms_per_blk;
      ivec.clear();
      vvec.clear();
      for(int k=0; k<nterms; k++, cnt++, row++, col++, val++) {
        ivec.push_back(*row);
        ivec.push_back(*col);
        vvec.push_back(static_cast<ValueType>(*val));  
      }
      bsize.push_back(nterms);
      dump.write(ivec,std::string("Spvn_index_")+std::to_string(nblk));
      dump.write(vvec,std::string("Spvn_vals_")+std::to_string(nblk));
      nblk++;
    } 

    ntot = sz;
    MPI_Status st;
    for(int rk=1; rk<nnodes_per_TG; rk++) {
     
      bool done = false;
      do { 
        // resize to receive
        ivec.resize(2*nterms_per_blk);
        vvec.resize(nterms_per_blk);
#if defined(QMC_COMPLEX)
        MPI_Recv(vvec.data(),2*nterms_per_blk,MPI_DOUBLE,ranks[rk],MPI_ANY_TAG,TG.getTGCOMM(),&st);
#else
        MPI_Recv(vvec.data(),nterms_per_blk,MPI_DOUBLE,ranks[rk],MPI_ANY_TAG,TG.getTGCOMM(),&st);
#endif
        MPI_Recv(ivec.data(),2*nterms_per_blk,MPI_INT,ranks[rk],MPI_ANY_TAG,TG.getTGCOMM(),&st);

        int nt_;
        MPI_Get_count(&st,MPI_INT,&nt_);
        nterms = std::size_t(nt_/2);
        if(st.MPI_TAG == MORE) {
          assert(nterms == nterms_per_blk);
        } else if(st.MPI_TAG == LAST) {
          assert(nterms <= nterms_per_blk);
          done=true;
        } else {
          APP_ABORT(" Error: This should not happen. \n\n\n");
        }

        ntot += std::size_t(nterms); 
        bsize.push_back(nterms);
        // resize to write
        ivec.resize(2*nterms);
        vvec.resize(nterms);
        dump.write(ivec,std::string("Spvn_index_")+std::to_string(nblk));
        dump.write(vvec,std::string("Spvn_vals_")+std::to_string(nblk));
        nblk++;
      } while(!done);  

    }

    dump.write(bsize,"Spvn_block_sizes");
  
  } else if(nnodes_per_TG>1 && core_rank == 0) {

    std::size_t sz = SpM.size();
    int nt,nb_loc = int(sz/nterms_per_blk + std::min(sz%nterms_per_blk,std::size_t(1))); 
    MPI_Reduce(&nb_loc,&nt,1,MPI_INT,MPI_SUM,0,TG.getTGCOMM());

    std::vector<IndexType> ivec;
    std::vector<ValueType> vvec;
    ivec.reserve(2*nterms_per_blk);
    vvec.reserve(nterms_per_blk);

    typename SpMatrix::int_iterator col = SpM.cols_begin();
    typename SpMatrix::int_iterator row = SpM.rows_begin();
    typename SpMatrix::iterator val = SpM.vals_begin();

    int cnt=0, nterms;
    IndexType nblk=0;
    for(int i=0; i<nb_loc; i++) {
      if( (sz-std::size_t(cnt)) < std::size_t(nterms_per_blk) )
        nterms = int(sz-std::size_t(cnt));
      else
        nterms = nterms_per_blk;
      ivec.clear();
      vvec.clear();
      for(int k=0; k<nterms; k++, cnt++, row++, col++, val++) {
        ivec.push_back(*row);
        ivec.push_back(*col);
        vvec.push_back(static_cast<ValueType>(*val));
      }
#if defined(QMC_COMPLEX)
      MPI_Send(vvec.data(),2*nterms,MPI_DOUBLE,0, (i==nb_loc-1)?(LAST):(MORE) ,TG.getTGCOMM());
#else
      MPI_Send(vvec.data(),nterms,MPI_DOUBLE,0, (i==nb_loc-1)?(LAST):(MORE) ,TG.getTGCOMM());
#endif
      MPI_Send(ivec.data(),2*nterms,MPI_INT,0, (i==nb_loc-1)?(LAST):(MORE) ,TG.getTGCOMM());
    }

  } else if(nnodes_per_TG>1) {

    int nb_loc=0,nt;
    MPI_Reduce(&nb_loc,&nt,1,MPI_INT,MPI_SUM,0,TG.getTGCOMM());

  } 
  return std::make_tuple(nblocks,ntot);

} 

}

}
