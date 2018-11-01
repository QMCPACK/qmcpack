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

#ifndef QMCPLUSPLUS_AFQMC_PHMSD_HPP
#define QMCPLUSPLUS_AFQMC_PHMSD_HPP

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <tuple>

#include "AFQMC/Utilities/readWfn.h"
#include "AFQMC/config.h"
#include "mpi3/shm/mutex.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "AFQMC/Utilities/taskgroup.h"
#include "AFQMC/Matrix/mpi3_SHMBuffer.hpp"
#include "AFQMC/Matrix/array_of_sequences.hpp"

#include "AFQMC/HamiltonianOperations/HamiltonianOperations.hpp"
#include "AFQMC/SlaterDeterminantOperations/SlaterDetOperations.hpp"

#include "AFQMC/Wavefunctions/phmsd_helpers.hpp"
#include "AFQMC/Wavefunctions/Excitations.hpp"


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
 * NOTE: Optimization note: CLOSED and NONCOLLINEAR calculations with a unique reference
 *                          only need a single set of unique overlaps/determinants/energies.
 *                          Fix this to improve performance!!!! 
 * THERE IS A PROBLEM WITH CLOSED SHELL CALCULATIONS!!!!
 * INCONSISTENCY WHEN SPIN DEPENDENT QUANTITIES ARE REUESTED, e.g. G. 
 * SOLUTION: For CLOSED and NONCOLLINEAR, force a single reference PHMSD.
 *           Then, only calculate alpha component of things and assume beta==alpha when needed,
 *           even if this is not true.
 */  
class PHMSD: public AFQMCInfo
{

  using CVector = boost::multi_array<ComplexType,1>;  
  using CMatrix = boost::multi_array<ComplexType,2>;  
  using SHM_Buffer = mpi3_SHMBuffer<ComplexType>;  
  using shared_mutex = boost::mpi3::shm::mutex;  
  using shared_CMatrix = boost::multi::array<ComplexType,2,shared_allocator<ComplexType>>;
  using shared_C3Tensor = boost::multi::array<ComplexType,3,shared_allocator<ComplexType>>;
  using shared_C4Tensor = boost::multi::array<ComplexType,4,shared_allocator<ComplexType>>;
  using index_aos = ma::sparse::array_of_sequences<int,int,
                                                   boost::mpi3::intranode::allocator<int>,
                                                   boost::mpi3::intranode::is_root>;

  public:

    PHMSD(AFQMCInfo& info, xmlNodePtr cur, afqmc::TaskGroup_& tg_, HamiltonianOperations&& hop_, 
          std::map<int,int>&& acta2mo_, std::map<int,int>&& actb2mo_,
          ph_excitations<int,ComplexType>&& abij_, 
          index_aos&& beta_coupled_to_unique_alpha__,
          index_aos&& alpha_coupled_to_unique_beta__,
          std::vector<PsiT_Matrix>&& orbs_, 
          WALKER_TYPES wlk, ValueType nce, int targetNW=1):
                AFQMCInfo(info),TG(tg_),
                SDetOp(((wlk!=2)?(NMO):(2*NMO)),((wlk!=2)?(NAEA):(NAEA+NAEB))),
                HamOp(std::move(hop_)),
                acta2mo(std::move(acta2mo_)),
                actb2mo(std::move(actb2mo_)),
                abij(std::move(abij_)),
                OrbMats(std::move(orbs_)),
                walker_type(wlk),NuclearCoulombEnergy(nce),
                shmbuff_for_E(nullptr),
                mutex(std::make_unique<shared_mutex>(TG.TG_local())),
                last_number_extra_tasks(-1),last_task_index(-1),
                local_group_comm(),
                shmbuff_for_G(nullptr),
                req_Gsend(MPI_REQUEST_NULL),
                req_Grecv(MPI_REQUEST_NULL),
                req_SMsend(MPI_REQUEST_NULL),
                req_SMrecv(MPI_REQUEST_NULL),
                maxnactive(std::max(OrbMats[0].shape()[0],OrbMats.back().shape()[0])),
                max_exct_n(std::max(abij.maximum_excitation_number()[0],
                                    abij.maximum_excitation_number()[1])),
                maxn_unique_confg(    
                    std::max(abij.number_of_unique_excitations()[0],
                             abij.number_of_unique_excitations()[1])),
                unique_overlaps({2,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                unique_Etot({2,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                QQ0inv0({1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                QQ0inv1({1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                GA2D0_shm({1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                GB2D0_shm({1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                local_ov(extents[2][maxn_unique_confg]),
                local_etot(extents[2][maxn_unique_confg]),
                local_QQ0inv0({OrbMats[0].shape()[0],NAEA}),
                local_QQ0inv1({OrbMats.back().shape()[0],NAEB}),
                Qwork(extents[max_exct_n][max_exct_n]),
                Gwork(extents[NAEA][maxnactive]),
                Ovmsd({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                Emsd({1,1,1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                QQ0A({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                QQ0B({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                GrefA({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}),
                GrefB({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                KEright({1,1,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                KEleft({1,1},shared_allocator<ComplexType>{TG.TG_local()}), 
                det_couplings{std::move(beta_coupled_to_unique_alpha__),
                              std::move(alpha_coupled_to_unique_beta__)}
    {
      /* To me, PHMSD is not compatible with walker_type=CLOSED unless
       * the MSD expansion is symmetric with respect to spin. For this, 
       * it is better to write a specialized class that assumes either spin symmetry
       * or e.g. Perfect Pairing.
       */
      if(walker_type == CLOSED) 
        APP_ABORT("Error: PHMSD requires walker_type != CLOSED.\n");

      compact_G_for_vbias = true; 
      transposed_G_for_vbias_ = HamOp.transposed_G_for_vbias();  
      transposed_G_for_E_ = HamOp.transposed_G_for_E();  
      transposed_vHS_ = HamOp.transposed_vHS();  
      fast_ph_energy = HamOp.fast_ph_energy();

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

    ~PHMSD() {
        if(req_SMrecv!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_SMrecv);
        if(req_SMsend!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_SMsend);
        if(req_Grecv!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_Grecv);
        if(req_Gsend!=MPI_REQUEST_NULL)
            MPI_Request_free(&req_Gsend);
    }

    PHMSD(PHMSD const& other) = delete;
    PHMSD& operator=(PHMSD const& other) = delete;
    PHMSD(PHMSD&& other) = default;
    PHMSD& operator=(PHMSD&& other) = default;

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
      assert( v.shape()[0] == HamOp.local_number_of_cholesky_vectors());
      double scl = (walker_type==COLLINEAR)?0.5:1.0;
      if(transposed_G_for_vbias_) {
        assert( G.shape()[0] == v.shape()[1] );
        assert( G.shape()[1] == size_of_G_for_vbias() );
        HamOp.vbias(G[indices[range_t()][range_t(0,OrbMats[0].shape()[0]*NMO)]],
                    std::forward<MatA>(v),scl*a,0.0);
        if(walker_type==COLLINEAR) 
          HamOp.vbias(G[indices[range_t()][range_t(OrbMats[0].shape()[0]*NMO,G.shape()[1])]],
                      std::forward<MatA>(v),scl*a,1.0);
      } else {  
        assert( G.shape()[0] == size_of_G_for_vbias() );
        assert( G.shape()[1] == v.shape()[1] );
        HamOp.vbias(G[indices[range_t(0,OrbMats[0].shape()[0]*NMO)][range_t()]],
                    std::forward<MatA>(v),scl*a,0.0);
        if(walker_type==COLLINEAR) 
          HamOp.vbias(G[indices[range_t(OrbMats[0].shape()[0]*NMO,G.shape()[0])][range_t()]],
                      std::forward<MatA>(v),scl*a,1.0);
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

    std::map<int,int> acta2mo;
    std::map<int,int> actb2mo;

    ph_excitations<int,ComplexType> abij;

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

    // shared memory arrays for temporary calculations
    bool fast_ph_energy;
    size_t maxn_unique_confg; // maximum number of unque configurations 
    size_t maxnactive;   // maximum number of states in active space
    size_t max_exct_n;   // maximum excitation number (number of electrons excited simultaneously)
    // used by OVerlap and MixedDensityMatrix
    shared_CMatrix unique_overlaps;
    shared_CMatrix unique_Etot;
    shared_CMatrix QQ0inv0;  // Q * inv(Q0) 
    shared_CMatrix QQ0inv1;  // Q * inv(Q0) 
    shared_CMatrix GA2D0_shm;  
    shared_CMatrix GB2D0_shm; 
    boost::multi_array<ComplexType,2> local_ov;
    boost::multi_array<ComplexType,2> local_etot;
    boost::multi::array<ComplexType,2> local_QQ0inv0;
    boost::multi::array<ComplexType,2> local_QQ0inv1;
    boost::multi_array<ComplexType,2> Qwork;     
    boost::multi_array<ComplexType,2> Gwork; 
    // used by Energy_shared 
    boost::multi_array<ComplexType,1> wgt; 
    boost::multi_array<ComplexType,1> opSpinEJ; 
    shared_C3Tensor Ovmsd;   // [nspins][maxn_unique_confg][nwalk]
    shared_C4Tensor Emsd;    // [nspins][maxn_unique_confg][nwalk][3]
    shared_C3Tensor QQ0A;    // [nwalk][NAOA][NAEA]
    shared_C3Tensor QQ0B;    // [nwalk][NAOB][NAEB]
    shared_C3Tensor GrefA;     // [nwalk][NAOA][NMO]
    shared_C3Tensor GrefB;     // [nwalk][NAOB][NMO]
    shared_C3Tensor KEright;   
    shared_CMatrix KEleft;     
     

    // array of sequence structure storing the list of connected alpha/beta configurations
    std::array<index_aos,2> det_couplings; 

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
    void Energy_distributed(const WlkSet& wset, Mat&& E, TVec&& Ov); 

    int dm_size(bool full) const {
      switch(walker_type) {
        case CLOSED: // closed-shell RHF
          return (full)?(NMO*NMO):(OrbMats[0].shape()[0]*NMO);
          break;
        case COLLINEAR:
          return (full)?(2*NMO*NMO):((OrbMats[0].shape()[0]+OrbMats.back().shape()[0])*NMO);
          break;
        case NONCOLLINEAR:
          return (full)?(4*NMO*NMO):((OrbMats[0].shape()[0])*2*NMO);
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
          return (full)?(arr{NMO,NMO}):(arr{OrbMats[0].shape()[0],NMO});
          break;
        case COLLINEAR:
          return (full)?(arr{NMO,NMO}):((sp==Alpha)?(arr{OrbMats[0].shape()[0],NMO}):(arr{OrbMats.back().shape()[0],NMO}));
          break;
        case NONCOLLINEAR:
          return (full)?(arr{2*NMO,2*NMO}):(arr{OrbMats[0].shape()[0],2*NMO});
          break;
        default:
          APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
          return arr{-1,-1};
      }
    }
    std::pair<int,int> dm_dims_ref(bool full, SpinTypes sp=Alpha) const {
      using arr = std::pair<int,int>;
      switch(walker_type) {
        case CLOSED: // closed-shell RHF
          return (full)?(arr{NMO,NMO}):(arr{NAEA,NMO});
          break;
        case COLLINEAR:
          return (full)?(arr{NMO,NMO}):((sp==Alpha)?(arr{NAEA,NMO}):(arr{NAEB,NMO}));
          break;
        case NONCOLLINEAR:
          return (full)?(arr{2*NMO,2*NMO}):(arr{NAEA,2*NMO});
          break;
        default:
          APP_ABORT(" Error: Unknown walker_type in dm_size. \n");
          return arr{-1,-1};
      }
    }

};

}

}

#include "AFQMC/Wavefunctions/PHMSD.icc"

#endif

