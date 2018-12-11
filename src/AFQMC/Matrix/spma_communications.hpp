////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory 
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_SPMA_COMMUNICATIONS_HPP
#define QMCPLUSPLUS_AFQMC_SPMA_COMMUNICATIONS_HPP

#include "mpi.h"
#include "AFQMC/config.h"

namespace qmcplusplus
{
 
namespace afqmc
{

template<class task_group>
void redistribute_sparse_matrix_size(std::vector<std::size_t>& nnz_per_row, task_group& TG)
{

  int node_in_TG = TG.getLocalNodeNumber();
  int nodeid = TG.getNodeID();
  MPI_Comm comm;
  MPI_Comm_split(TG.Cores().impl_,node_in_TG,nodeid,&comm);

  std::vector<std::size_t> n = nnz_per_row;
  MPI_Allreduce(n.data(),nnz_per_row.data(),n.size(),MPI_UNSIGNED_LONG,MPI_SUM,comm);

  MPI_Comm_free(&comm);
}
 
template<class SpMatrix,
          class taskgroup>
inline void redistribute_sparse_matrix(taskgroup& TG, SpMatrix& SpMat)
{

  using T = typename SpMatrix::value_type;

  TG.global_barrier(); 

  if(TG.getNNodesPerTG() > 1) {

    APP_ABORT(" Error: Distributed version of redistribute_sparse_matrix not yet implemented. \n");

    int nnodes = TG.getTotalNodes(), nodeid = TG.getNodeID();
    int ncores = TG.getTotalCores(), coreid = TG.getCoreID();
    int node_in_TG = TG.getLocalNodeNumber();
    MPI_Comm comm;
    MPI_Comm_split(TG.Cores().impl_,node_in_TG,nodeid,&comm);
    int grpid,ngrps;
    MPI_Comm_size(comm,&ngrps);
    MPI_Comm_rank(comm,&grpid);


    int nmax = 1000000;
    std::vector<int> cols, index;
    std::vector<T> vals;
    cols.reserve(nmax);
    vals.reserve(nmax);

    int nr = SpMat.rows();
    int nr0, nr1;
    std::tie(nr0, nr1) = FairDivideBoundary(coreid, nr, ncores);
    index.resize(nr1-nr0+1);

    typename SpMatrix::int_iterator nnz = SpMat.nnz_per_row_begin() + nr0;

    // record initial number of elements in each row
    std::vector<int> initial_sizes(nr1-nr0);
    for(int r=0; r<nr1-nr0; r++)
      initial_sizes[r] = *(nnz+r);

    SpMat.barrier();
    // loop over blocks, copy data to buffer, bcast, add to matrix
    for(int grp = 0; grp < ngrps; grp++) {

      int nterms, nblk, cnt=0;
      if( grp == grpid ) nterms = std::accumulate(initial_sizes.begin(),initial_sizes.end(),int(0));
      MPI_Bcast(&nterms,1,MPI_INT,grp,comm);
      nblk = nterms/nmax + ((nterms%nmax!=0)?(1):(0));

      for(int blk=0, nitot=0; blk<nblk; blk++) {
        int ni = std::min(nitot+nmax,nterms);

        if(grp == grpid) {
          cols.clear();
          vals.clear();
          // must recalculate this every time in case rowIndex changes when adding elements below
          typename SpMatrix::indx_iterator SpMat_index = SpMat.rowIndex_begin() + nr0;
          index[0] = 0;
          for(int r=0, nsum=0; r<nr1-nr0; r++) {

            if(nsum + initial_sizes[r] <= nitot) {
              // skip row, already sent
              nsum += initial_sizes[r];
              index[r+1] = index[r]; // should be zero
            } else {

              int n0 = std::max(0,nitot-nsum);
              int nn = std::min(initial_sizes[r],n0+ni-index[r]);
              std::size_t rI0 = *(SpMat.row_index(r)) + static_cast<std::size_t>(n0);
              std::size_t rIN = *(SpMat.row_index(r)) + static_cast<std::size_t>(nn);

              index[r+1] = index[r]+nn;

              // dump data
              for(std::size_t i=rI0; i<rIN; i++) {
                vals.push_back(*(SpMat.values(i)));
                cols.push_back(*(SpMat.column_data(i)));
                nsum++;
              }

            }
            if(nsum == nitot+ni) break;
            if(nsum > nitot+ni)  {
              app_error()<<" Bug in redistribute_sparse_matrix. " <<std::endl;
              APP_ABORT(" Bug in redistribute_sparse_matrix. ");
            }

          }

          if(vals.size() != ni || cols.size() != ni) {
            app_error()<<" Bug in redistribute_sparse_matrix (1). " <<std::endl;
            APP_ABORT(" Bug in redistribute_sparse_matrix (1). ");
          }

        } else {
          cols.resize(ni);
          vals.resize(ni);
        }

        nitot += ni;

        MPI_Bcast(vals.data(),vals.size()*sizeof(T),MPI_CHAR,grp,comm);
        MPI_Bcast(cols.data(),cols.size(),MPI_INT,grp,comm);
        MPI_Bcast(index.data(),index.size(),MPI_INT,grp,comm);

        if(grp != grpid) {
          for(int r=0; r<nr1-nr0; r++)
            for(int n=index[r]; n<index[r+1]; n++)
              SpMat.add(r+nr0,cols[n],vals[n],false);  // no need for locks, row segments are non-overlapping
        }
      }
    }
    MPI_Comm_free(&comm);

  } else if(TG.getCoreID() == 0) {

    int rk = TG.getNodeID(), npr = TG.getTotalNodes(); 
    long ptr,n0;
    n0 = SpMat.size(); // my number of terms, always from zero to n0
    ptr = n0; // position to copy elements to 
    std::vector<long> size(npr);
    TG.Cores().all_gather_value(n0,size.begin()); 
    long ntot = 0;
    for(int i=0; i<npr; i++) ntot+=size[i];

std::cout<<TG.Cores().rank() <<" " <<n0 <<" " <<ntot <<std::endl;
TG.Cores().barrier();

    if(ntot > SpMat.capacity()) {
      app_error()<<" Problems gathering hamiltonian. Capacity of std::vector is not sufficient: " <<ntot <<" " <<SpMat.capacity() <<" \n";
      APP_ABORT("Problems gathering hamiltonian. Capacity of std::vector is not sufficient. \n");
    }

    TG.Cores().barrier();  
    SpMat.resize_serial(ntot);
    TG.Cores().barrier();  

    for(int i=0; i<npr; i++) {
      if(i==rk) { // I send
        MPI_Bcast(SpMat.row_data(),n0,MPI_INT,i,TG.Cores().impl_);
        MPI_Bcast(SpMat.column_data(),n0,MPI_INT,i,TG.Cores().impl_);
        MPI_Bcast(SpMat.values(),n0*sizeof(T),MPI_CHAR,i,TG.Cores().impl_);
//        TG.Cores().broadcast_n(SpMat.rows_begin(),n0,i);
//        TG.Cores().broadcast_n(SpMat.cols_begin(),n0,i);
//        TG.Cores().broadcast_n(SpMat.vals_begin(),n0,i);
      } else { // I reveive
        MPI_Bcast(SpMat.row_data()+ptr,size[i],MPI_INT,i,TG.Cores().impl_);
        MPI_Bcast(SpMat.column_data()+ptr,size[i],MPI_INT,i,TG.Cores().impl_);
        MPI_Bcast(SpMat.values()+ptr,size[i]*sizeof(T),MPI_CHAR,i,TG.Cores().impl_);
//        TG.Cores().broadcast_n(SpMat.rows_begin()+ptr,size[i],i);
//        TG.Cores().broadcast_n(SpMat.cols_begin()+ptr,size[i],i);
//        TG.Cores().broadcast_n(SpMat.vals_begin()+ptr,size[i],i);
        ptr+=size[i];
      }
    }
  }
  TG.global_barrier(); 
}

}

}

#endif
