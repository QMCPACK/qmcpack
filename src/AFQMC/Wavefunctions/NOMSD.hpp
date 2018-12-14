//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_AFQMC_NOMSD_H
#define QMCPLUSPLUS_AFQMC_NOMSD_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/config.h"
#include "mpi3/shm/mutex.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/mpi3_SHMBuffer.hpp"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"


namespace qmcplusplus
{

namespace afqmc
{

/*
 * Class that implements a multi-Slater determinant trial wave-function.
 * Single determinant wfns are also allowed. 
 * No relation between different determinants in the expansion is assumed.
 * Designed for non-orthogonal MSD expansions. 
 * For particle-hole orthogonal MSD wfns, use FastMSD.
 */  
class NOMSD: public AFQMCInfo
{

  using CVector = boost::multi_array<ComplexType,1>;  
  using CMatrix = boost::multi_array<ComplexType,2>;  
  using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;  
  using shared_mutex = boost::mpi3::shm::mutex;  

  public:

    NOMSD(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, HamiltonianOperations&& hop_, 
          std::vector<ComplexType>&& ci_, std::vector<PsiT_Matrix>&& orbs_, 
          WALKER_TYPES wlk, ValueType nce, int targetNW=1):
                AFQMCInfo(info),TG(tg_),
                SDetOp(((wlk!=2)?(NMO):(2*NMO)),((wlk!=2)?(NAEA):(NAEA+NAEB))),
                HamOp(std::move(hop_)),ci(std::move(ci_)),OrbMats(std::move(orbs_)),
                walker_type(wlk),NuclearCoulombEnergy(nce),
                shmbuff_for_E(nullptr),
                mutex(std::make_unique<shared_mutex>(TG.TG_local())),
                last_number_extra_tasks(-1),last_task_index(-1),
                //local_group_comm(nullptr),
                //local_group_comm(std::make_unique<shared_communicator>(TG.TG_local().split(0))),
                local_group_comm(),
                shmbuff_for_G(nullptr),
                req_Gsend(MPI_REQUEST_NULL),
                req_Grecv(MPI_REQUEST_NULL),
                req_SMsend(MPI_REQUEST_NULL),
                req_SMrecv(MPI_REQUEST_NULL)
    {
      compact_G_for_vbias = (ci.size()==1); // this should be input, since it is determined by HOps 
      transposed_G_for_vbias_ = HamOp.transposed_G_for_vbias();  
      transposed_G_for_E_ = HamOp.transposed_G_for_E();  
      transposed_vHS_ = HamOp.transposed_vHS();  

      excitedState = false;  
      std::string excited_file("");  
      int i_=-1,a_=-1;  
      ParameterSet m_param;
      m_param.add(excited_file,"excited","std::string");
      // generalize this to multi-particle excitations, how do I read a list of integers???
      m_param.add(i_,"i","int");
      m_param.add(a_,"a","int");
      m_param.put(cur);

      if(excited_file != "" && 
         i_ >= 0 &&
         a_ >= 0) {
        if(i_ < NMO && a_ < NMO) { 
          if(i_ >= NAEA || a_ < NAEA)
            APP_ABORT(" Errors: Inconsistent excited orbitals for alpha electrons. \n");
          excitedState=true;
          maxOccupExtendedMat = {a_,NAEB};
          numExcitations = {1,0};
          excitations.push_back({i_,a_});
        } else if(i_ >= NMO && a_ >= NMO) {
          if(i_ >= NMO+NAEB || a_ < NMO+NAEB)
            APP_ABORT(" Errors: Inconsistent excited orbitals for beta electrons. \n");
          excitedState=true;
          maxOccupExtendedMat = {NAEA,a_-NMO};
          numExcitations = {0,1};
          excitations.push_back({i_-NMO,a_-NMO});
        } else {
          APP_ABORT(" Errors: Inconsistent excited orbitals. \n");
        }    
        readWfn(excited_file,excitedOrbMat,NMO,maxOccupExtendedMat.first,maxOccupExtendedMat.second);
      }    
    }

    ~NOMSD() {
        if(req_SMrecv!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_SMrecv);
        if(req_SMsend!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_SMsend);
        if(req_Grecv!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_Grecv);
        if(req_Gsend!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_Gsend);
    }

    NOMSD(NOMSD const& other) = delete;
    NOMSD& operator=(NOMSD const& other) = delete;
    NOMSD(NOMSD&& other) = default;
    NOMSD& operator=(NOMSD&& other) = default;

    int local_number_of_cholesky_vectors() const 
    { return HamOp.local_number_of_cholesky_vectors(); }
    int global_number_of_cholesky_vectors() const 
    { return HamOp.global_number_of_cholesky_vectors(); }
    bool distribution_over_cholesky_vectors() const 
    { return HamOp.distribution_over_cholesky_vectors(); }

    int size_of_G_for_vbias() const 
    {  return dm_size(!compact_G_for_vbias);  }

    bool transposed_G_for_vbias() const { return transposed_G_for_vbias_; }
    bool transposed_G_for_E() const { return transposed_G_for_E_; }
    bool transposed_vHS() const { return transposed_vHS_; }

/*
    const std::vector<PsiT_Matrix>& getOrbMat() { return OrbMats; }
    int getOrbSize () { return 2*NMO; }
    const std::vector<ComplexType>& getCiCoeff() { return ci; }
*/

    template<class Vec>
    void vMF(Vec&& v);

    CMatrix getOneBodyPropagatorMatrix(TaskGroup_& TG, CVector const& vMF)
    { return HamOp.getOneBodyPropagatorMatrix(TG,vMF); }

    /*
     * local contribution to vbias for the Green functions in G 
     * G: [size_of_G_for_vbias()][nW]
     * v: [local # Chol. Vectors][nW]
     */
    template<class MatG, class MatA>
    void vbias(const MatG& G, MatA&& v, double a=1.0) {
      if(transposed_G_for_vbias_) {
        assert( G.shape()[0] == v.shape()[1] );
        assert( G.shape()[1] == size_of_G_for_vbias() );
      } else {  
        assert( G.shape()[0] == size_of_G_for_vbias() );
        assert( G.shape()[1] == v.shape()[1] );
      }  
      assert( v.shape()[0] == HamOp.local_number_of_cholesky_vectors());
      if(ci.size()==1) {
        // HamOp expects a compact Gc with alpha/beta components 
        HamOp.vbias(G,std::forward<MatA>(v),a);
      } else {
        if(walker_type == CLOSED )
          HamOp.vbias(G,std::forward<MatA>(v),a); // factor of 2 now in HamOps
        else if(walker_type == NONCOLLINEAR)
          HamOp.vbias(G,std::forward<MatA>(v),a);
        else {
          // HamOp expects either alpha or beta, so must be called twice 
          HamOp.vbias(G[indices[range_t(0,NMO*NMO)][range_t()]]
                ,std::forward<MatA>(v),a,0.0);
          HamOp.vbias(G[indices[range_t(NMO*NMO,2*NMO*NMO)][range_t()]]
                ,std::forward<MatA>(v),a,1.0);
        }
      }  
      TG.local_barrier();    
    }

    /*
     * local contribution to vHS for the Green functions in G 
     * X: [# chol vecs][nW]
     * v: [NMO^2][nW]
     */
    template<class MatX, class MatA>
    void vHS(const MatX& X, MatA&& v, double a=1.0) {
      assert( X.shape()[0] == HamOp.local_number_of_cholesky_vectors() );
      if(transposed_G_for_vbias_)
        assert( X.shape()[1] == v.shape()[0] );
      else    
        assert( X.shape()[1] == v.shape()[1] );
      HamOp.vHS(X,std::forward<MatA>(v),a);
      TG.local_barrier();    
    }

    /*
     * Calculates the local energy and overlaps of all the walkers in the set and stores
     * them in the wset data
     */
    template<class WlkSet>
    void Energy(WlkSet& wset) {
      int nw = wset.size();
      if(ovlp.num_elements() != nw)
        ovlp.resize(extents[nw]);
      if(eloc.shape()[0] != nw || eloc.shape()[1] != 3)
        eloc.resize(extents[nw][3]);
      Energy(wset,eloc,ovlp);
      TG.local_barrier();
      if(TG.getLocalTGRank()==0) {
	int p=0;
	for(typename WlkSet::iterator it=wset.begin(); it!=wset.end(); ++it, ++p) {
	  it->overlap() = ovlp[p];
	  it->E1() = eloc[p][0];		
	  it->EXX() = eloc[p][1];		
	  it->EJ() = eloc[p][2];		
	}
      }  
      TG.local_barrier();
    }

    /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
    template<class WlkSet, class Mat, class TVec> 
    void Energy(const WlkSet& wset, Mat&& E, TVec&& Ov) {
      if(TG.getNNodesPerTG() > 1)
        Energy_distributed(wset,std::forward<Mat>(E),std::forward<TVec>(Ov));
      else
        Energy_shared(wset,std::forward<Mat>(E),std::forward<TVec>(Ov));
    }

    /*
     * Calculates the mixed density matrix for all walkers in the walker set. 
     * Options:
     *  - compact:   If true (default), returns compact form with Dim: [NEL*NMO], 
     *                 otherwise returns full form with Dim: [NMO*NMO]. 
     *  - transpose: If false (default), returns standard form with Dim: [XXX][nW]
     *                 otherwise returns the transpose with Dim: [nW][XXX}
     */
    template<class WlkSet, class MatG>
    void MixedDensityMatrix(const WlkSet& wset, MatG&& G, bool compact=true, bool transpose=false) {
      int nw = wset.size();
      if(ovlp.num_elements() != nw)
        ovlp.resize(extents[nw]);
      MixedDensityMatrix(wset,std::forward<MatG>(G),ovlp,compact,transpose);
    }

    template<class WlkSet, class MatG, class TVec>
    void MixedDensityMatrix(const WlkSet& wset, MatG&& G, TVec&& Ov, bool compact=true, bool transpose=false);

    /*
     * Calculates the mixed density matrix for all walkers in the walker set
     *   with a format consistent with (and expected by) the vbias routine.
     * This is implementation dependent, so this density matrix should ONLY be used
     * in conjunction with vbias. 
     */
    template<class WlkSet, class MatG>
    void MixedDensityMatrix_for_vbias(const WlkSet& wset, MatG&& G) {
      int nw = wset.size();
      if(ovlp.num_elements() != nw)
        ovlp.resize(extents[nw]);	
      MixedDensityMatrix(wset,std::forward<MatG>(G),ovlp,compact_G_for_vbias,transposed_G_for_vbias_);
    }

    /*
     * Calculates the overlaps of all walkers in the set. Returns values in arrays. 
     */
    template<class WlkSet, class TVec>
    void Overlap(const WlkSet& wset, TVec&& Ov);

    /*
     * Calculates the overlaps of all walkers in the set. Updates values in wset. 
     */
    template<class WlkSet>
    void Overlap(WlkSet& wset)
    {
      int nw = wset.size();
      if(ovlp.num_elements() != nw)
        ovlp.resize(extents[nw]);
      Overlap(wset,ovlp);
      TG.local_barrier();
      if(TG.getLocalTGRank()==0) {
        int p=0;
        for(typename WlkSet::iterator it=wset.begin(); it!=wset.end(); ++it, ++p) 
          it->overlap() = ovlp[p];
      }	 
      TG.local_barrier();
    }


    /*
     * Orthogonalizes the Slater matrices of all walkers in the set.  
     * Options:
     *  - bool importanceSamplingt(default=true): use algorithm appropriate for importance sampling. 
     *         This means that the determinant of the R matrix in the QR decomposition is ignored.
     *         If false, add the determinant of R to the weight of the walker. 
     */
    template<class WlkSet>
    void Orthogonalize(WlkSet& wset, bool impSamp); 

    /*
     * Orthogonalizes the Slater matrix of a walker in an excited state calculation.
     */
    template<class Mat>
    void OrthogonalizeExcited(Mat&& A, SpinTypes spin);

  protected: 

    TaskGroup_& TG;
 
    SlaterDetOperations<ComplexType> SDetOp;
  
    HamiltonianOperations HamOp;

    std::vector<ComplexType> ci;

    // eventually switched from CMatrix to SMHSparseMatrix(node)
    std::vector<PsiT_Matrix> OrbMats;

    std::unique_ptr<SHM_Buffer> shmbuff_for_E;

    std::unique_ptr<shared_mutex> mutex;

    // in both cases below: closed_shell=0, UHF/ROHF=1, GHF=2
    WALKER_TYPES walker_type;

    bool compact_G_for_vbias;

    // in the 3 cases, true means [nwalk][...], false means [...][nwalk]
    bool transposed_G_for_vbias_;
    bool transposed_G_for_E_;
    bool transposed_vHS_;

    ValueType NuclearCoulombEnergy; 

    // not elegant, but reasonable for now
    int last_number_extra_tasks;
    int last_task_index;

    // shared_communicator for parallel work within TG_local()
    //std::unique_ptr<shared_communicator> local_group_comm; 
    shared_communicator local_group_comm; 
    std::unique_ptr<SHM_Buffer> shmbuff_for_G;

    // excited states
    bool excitedState;
    std::vector<std::pair<int,int>> excitations;
    boost::multi_array<ComplexType,3> excitedOrbMat; 
    CMatrix extendedMatAlpha;
    CMatrix extendedMatBeta;
    std::pair<int,int> maxOccupExtendedMat;
    std::pair<int,int> numExcitations; 

    // buffers and work arrays, careful here!!!
    CVector ovlp, localGbuff, ovlp2;
    CMatrix eloc, eloc2, eloc3;

    MPI_Request req_Gsend, req_Grecv;
    MPI_Request req_SMsend, req_SMrecv;

    /* 
     * Computes the mixed density matrix of a single given determinant in the trial wave function.
     * Intended to be used in combination with the energy evaluation routine.
     * G and Ov are expected to be in shared memory.
     */
    template<class WlkSet, class MatG, class TVec>
    void MixedDensityMatrix_for_E(const WlkSet& wset, MatG&& G, TVec&& Ov, int nd);

    template<class MatSM, class MatG, class TVec>
    void MixedDensityMatrix_for_E_from_SM(const MatSM& SM, MatG&& G, TVec&& Ov, int nd); 

    /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
    template<class WlkSet, class Mat, class TVec>
    void Energy_shared(const WlkSet& wset, Mat&& E, TVec&& Ov);

    /*
     * Calculates the local energy and overlaps of all the walkers in the set and 
     * returns them in the appropriate data structures
     */
    template<class WlkSet, class Mat, class TVec>
    void Energy_distributed(const WlkSet& wset, Mat&& E, TVec&& Ov) {
      if(ci.size()==1)
        Energy_distributed_singleDet(wset,std::forward<Mat>(E),std::forward<TVec>(Ov));
      else
        Energy_distributed_multiDet(wset,std::forward<Mat>(E),std::forward<TVec>(Ov));
    }

    template<class WlkSet, class Mat, class TVec>
    void Energy_distributed_singleDet(const WlkSet& wset, Mat&& E, TVec&& Ov);

    template<class WlkSet, class Mat, class TVec>
    void Energy_distributed_multiDet(const WlkSet& wset, Mat&& E, TVec&& Ov);

    int dm_size(bool full) const {
      switch(walker_type) {
        case CLOSED: // closed-shell RHF
          return (full)?(NMO*NMO):(NAEA*NMO);
          break;
        case COLLINEAR:
          return (full)?(2*NMO*NMO):((NAEA+NAEB)*NMO);
          break;
        case NONCOLLINEAR:
          return (full)?(4*NMO*NMO):((NAEA+NAEB)*2*NMO);
          break;
        default:
          APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
          return -1;
      }
    }
    // dimensions for each component of the DM. 
    std::pair<int,int> dm_dims(bool full, SpinTypes sp=Alpha) const {
      using arr = std::pair<int,int>;
      switch(walker_type) {
        case CLOSED: // closed-shell RHF
          return (full)?(arr{NMO,NMO}):(arr{NAEA,NMO});
          break;
        case COLLINEAR:
          return (full)?(arr{NMO,NMO}):((sp==Alpha)?(arr{NAEA,NMO}):(arr{NAEB,NMO}));
          break;
        case NONCOLLINEAR:
          return (full)?(arr{2*NMO,2*NMO}):(arr{NAEA+NAEB,2*NMO});
          break;
        default:
          APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
          return arr{-1,-1};
      }
    }

};

}

}

#include "AFQMC/Wavefunctions/NOMSD.icc"

#endif

