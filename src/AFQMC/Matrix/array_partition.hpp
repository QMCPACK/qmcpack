//////////////////////////////////////////////////////////////////////////////////////
////// This file is distributed under the University of Illinois/NCSA Open Source License.
////// See LICENSE file in top directory for details.
//////
////// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//////
////// File developed by: 
//////
////// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////////

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
#include "Utilities/UtilityFunctions.h"

#include "AFQMC/Utilities/Utils.h"
#include<AFQMC/config.0.h>

namespace qmcplusplus
{

namespace afqmc
{

/*
 * object that encapsulates the partitioning of a 2-D array (matrix) 
 */
template<class task_group,
         typename IType,
         typename DType
         >
struct simple_matrix_partition
{
  simple_matrix_partition(IType nr, IType nc, DType rc=0):cut(std::abs(rc)),r0(0),r1(nr),c0(0),c1(nc),nrows(nr),ncols(nc) {}

  ~simple_matrix_partition() {}

  template<typename DType2>
  inline bool local(const IType i, const IType j, const DType2 v) const { 
    return (i>=r0) && (i<r1) &&
           (j>=c0) && (j<c1) &&
           (std::abs(v) > cut);
  }  

  // this breaks a matrix dimension over TG.nnodes, assuming homogeneous blocks 
  inline void partition(const task_group& TG, bool byRow, std::vector<IType>& sets)
  {
    int nnodes = TG.getNNodesPerTG();
    IType cnt=0;
    int nblk; 
    if(byRow)
      nblk = nrows;
    else
      nblk = ncols;
    assert(nblk >= nnodes);
    sets.resize(nnodes+1);
    sets[0] = 0;
    if(nnodes > 1) {
      FairDivide(nblk,nnodes,sets);
      int node_number = TG.getLocalNodeNumber(); 
      if(byRow) {
        r0=sets[node_number];
        r1=sets[node_number+1];
      } else {
        c0=sets[node_number];
        c1=sets[node_number+1];
      }
    } else {
      if(byRow) {
        r0=0;
        r1=nrows;
        sets[1]=nrows;
      } else {
        c0=0;
        c1=ncols;
        sets[1]=ncols;
      }
    }    
  }

  // this breaks a matrix dimension over TG.nnodes 
  inline void partition(const task_group& TG, bool byRow, const std::vector<IType>& counts, std::vector<IType>& sets)
  {
    int nnodes = TG.getNNodesPerTG();
    int nblk = counts.size();
    IType cnt=0;
    assert(nblk >= nnodes);
    if(byRow) 
      assert(nblk == nrows);
    else
      assert(nblk == ncols);
    sets.resize(nnodes+1);
    sets[0] = 0;
    if(nnodes > 1) {
      std::vector<IType> nv(counts.size()+1);
      nv[0]=0;
      typename std::vector<IType>::iterator itn=nv.begin()+1;
      for(typename std::vector<IType>::const_iterator itc=counts.begin(), ite=counts.end(); itc!=ite; itc++,itn++) {
        cnt+=*(itc);
        (*itn)=cnt;
      }
      balance_partition_ordered_set(counts.size(),nv.data(),sets);
      int node_number = TG.getLocalNodeNumber(); 
      if(byRow) {
        r0=sets[node_number];
        r1=sets[node_number+1];
      } else {
        c0=sets[node_number];
        c1=sets[node_number+1];
      }
    } else {
      if(byRow) {
        r0=0;
        r1=nrows;
        sets[1]=nrows;
      } else {
        c0=0;
        c1=ncols;
        sets[1]=ncols;
      }
    }    
  }

  // this breaks a matrix dimension over TG.
  inline void partition_over_TGs(const task_group& TG, bool byRow, const std::vector<IType>& counts, std::vector<IType>& sets)
  {
    int ngrps = TG.getNumberOfTGs();
    int nblk = counts.size();
    IType cnt=0;
    assert(nblk >= ngrps);
    if(byRow)
      assert(nblk == nrows);
    else
      assert(nblk == ncols);
    sets.resize(ngrps+1);
    sets[0] = 0;
    if(ngrps > 1) {
      std::vector<IType> nv(counts.size()+1);
      nv[0]=0;
      typename std::vector<IType>::iterator itn=nv.begin()+1;
      for(typename std::vector<IType>::const_iterator itc=counts.begin(), ite=counts.end(); itc!=ite; itc++,itn++) {
        cnt+=*(itc);
        (*itn)=cnt;
      }
      balance_partition_ordered_set(counts.size(),nv.data(),sets);
      int node_number = TG.getTGNumber();
      if(byRow) {
        r0=sets[node_number];
        r1=sets[node_number+1];
      } else {
        c0=sets[node_number];
        c1=sets[node_number+1];
      }
    } else {
      if(byRow) {
        r0=0;
        r1=nrows;
        sets[1]=nrows;
      } else {
        c0=0;
        c1=ncols;
        sets[1]=ncols;
      }
    }
  }

  // this breaks a local segment over TG.ncores_per_TG, assumes homogeneous blocks 
  inline void sub_partition(const task_group& TG, bool byRow, std::vector<IType>& sets) 
  {
    int ncores = TG.getNCoresPerTG();
    int nblk; 
    if(byRow)
      nblk = r1-r0;
    else
      nblk = c1-c0;
    assert(nblk >= ncores);
    sets.resize(ncores+1);
    sets[0] = 0;
    if(ncores > 1) {
      FairDivide(nblk,ncores,sets);
      int core_rank = TG.getCoreRank();  
      if(byRow) {
        sr0=r0+sets[core_rank];
        sr1=r0+sets[core_rank+1];
      } else {
        sc0=c0+sets[core_rank];
        sc1=c0+sets[core_rank+1];
      }
    } else {
      if(byRow) {
        sr0=r0;
        sr1=r1;
        sets[1]=r1;
      } else {
        sc0=c0;
        sc1=c1;
        sets[1]=c1;
      }
    }

  }

  // this breaks a local segment over TG.ncores_per_TG 
  inline void sub_partition(const task_group& TG, bool byRow, const std::vector<IType>& counts, std::vector<IType>& sets) 
  {
    int ncores = TG.getNCoresPerTG();
    int nblk = counts.size();
    IType cnt=0;
    assert(nblk >= ncores);
    if(byRow)
      assert(nblk == r1-r0);
    else
      assert(nblk == c1-c0);
    sets.resize(ncores+1);
    sets[0] = 0;
    if(ncores > 1) {
      std::vector<IType> nv(counts.size()+1);
      nv[0]=0;
      typename std::vector<IType>::iterator itn=nv.begin()+1;
      for(typename std::vector<IType>::const_iterator itc=counts.begin(), ite=counts.end(); itc!=ite; itc++,itn++) {
        cnt+=*(itc);
        (*itn)=cnt;
      }
      balance_partition_ordered_set(counts.size(),nv.data(),sets);
      int core_rank = TG.getCoreRank();  
      if(byRow) {
        sr0=r0+sets[core_rank];
        sr1=r0+sets[core_rank+1];
      } else {
        sc0=c0+sets[core_rank];
        sc1=c0+sets[core_rank+1];
      }
    } else {
      if(byRow) {
        sr0=r0;
        sr1=r1;
        sets[1]=r1;
      } else {
        sc0=c0;
        sc1=c1;
        sets[1]=c1;
      }
    }

  }

  std::tuple<int,int,int,int> getLocalPartition() const {
    return std::make_tuple(r0,r1,c0,c1);
  }

  std::tuple<int,int,int,int> getSubPartition() const {
    return std::make_tuple(sr0,sr1,sc0,sc1);
  }

  DType getCutoff() const {
    return cut; 
  }

  std::tuple<int,int> getDimensions() const {
    return std::make_tuple(nrows,ncols);
  }
  
  private:

  DType cut;
  IType nrows, ncols;
  IType r0, r1;  // lower and upper bond of row segment
  IType c0, c1;  // lower and upper bound of column segment 
  IType sr0, sr1;  // lower and upper bond of row sub-partition (relative to 0) 
  IType sc0, sc1;  // lower and upper bound of column sub-partition (relative to 0)

};
  

}

}

