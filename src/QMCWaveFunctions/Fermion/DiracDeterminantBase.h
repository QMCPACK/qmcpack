//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file DiracDeterminantBaseBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#define QMCPLUSPLUS_DIRACDETERMINANTWITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Utilities/NewTimer.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"

namespace qmcplusplus
  {

  class DiracDeterminantBase: public OrbitalBase
    {
    protected:
      ParticleSet *targetPtcl;
    public:
      bool Optimizable;
      void registerTimers();
      NewTimer UpdateTimer, RatioTimer, InverseTimer, BufferTimer, SPOVTimer, SPOVGLTimer;
      // Optimizable parameters
      opt_variables_type myVars;


      typedef SPOSetBase::IndexVector_t IndexVector_t;
      typedef SPOSetBase::ValueVector_t ValueVector_t;
      typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
      typedef SPOSetBase::GradVector_t  GradVector_t;
      typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
      typedef SPOSetBase::HessMatrix_t  HessMatrix_t;
      typedef SPOSetBase::HessType      HessType;

      /** constructor
       *@param spos the single-particle orbital set
       *@param first index of the first particle
       */
      DiracDeterminantBase(SPOSetBasePtr const &spos, int first=0);

      ///default destructor
      virtual ~DiracDeterminantBase();

      /**copy constructor
       * @param s existing DiracDeterminantBase
       *
       * This constructor makes a shallow copy of Phi.
       * Other data members are allocated properly.
       */
      DiracDeterminantBase(const DiracDeterminantBase& s);

      DiracDeterminantBase& operator=(const DiracDeterminantBase& s);

      ///** return a clone of Phi
      // */
      //SPOSetBasePtr clonePhi() const;

      SPOSetBasePtr getPhi() { return Phi; };

      inline IndexType rows() const { return NumPtcls; }

      inline IndexType cols() const { return NumOrbitals; }

      /** set the index of the first particle in the determinant and reset the size of the determinant
       *@param first index of first particle
       *@param nel number of particles in the determinant
       */
      virtual void set(int first, int nel);
      virtual RealType getAlternatePhaseDiff(){return 0.0;}
      virtual RealType getAlternatePhaseDiff(int iat){return 0.0;}

      ///set BF pointers
      virtual
      void setBF(BackflowTransformation* BFTrans) {}

      ///optimizations  are disabled
      virtual inline void checkInVariables(opt_variables_type& active)
      {
        Phi->checkInVariables(active);
        Phi->checkInVariables(myVars);
      }

      virtual inline void checkOutVariables(const opt_variables_type& active)
      {
	Phi->checkOutVariables(active);
	myVars.clear();
	myVars.insertFrom(Phi->myVars);
	myVars.getIndex(active);
      }
      
      virtual void resetParameters(const opt_variables_type& active)
      {
        Phi->resetParameters(active);
	for(int i=0; i<myVars.size(); ++i) {
	  int ii=myVars.Index[i];
	  if(ii>=0) myVars[i]= active[ii];
	}
      }

      virtual void evaluateDerivatives(ParticleSet& P,
				       const opt_variables_type& active,
				       vector<RealType>& dlogpsi,
				       vector<RealType>& dhpsioverpsi);

      // used by DiracDeterminantWithBackflow 
      virtual void evaluateDerivatives(ParticleSet& P,
                                            const opt_variables_type& active,
                                            int offset,
                                            Matrix<RealType>& dlogpsi,
                                            Array<GradType,3>& dG,
                                            Matrix<RealType>& dL) {}

      inline void reportStatus(ostream& os)
      {
      }
      void resetTargetParticleSet(ParticleSet& P)
      {
        Phi->resetTargetParticleSet(P); targetPtcl = &P;
      }

      ///reset the size: with the number of particles and number of orbtials
      virtual void resize(int nel, int morb);

      virtual RealType registerData(ParticleSet& P, PooledData<RealType>& buf);
      
      virtual void registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType=0);

      virtual void memoryUsage_DataForDerivatives(ParticleSet& P,long& orbs_only, long& orbs, long& invs, long& dets)
      {
      // mmorales: not sure you need to store myL,myG for nonlocal psp optimization 
          orbs_only += NumPtcls*NumOrbitals; 
          orbs += NumPtcls*NumOrbitals + NP*4 + 2; 
      }

      virtual void copyToDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);

      virtual void copyFromDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);
      
      virtual RealType evaluateLogForDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);
      
      virtual RealType evaluateLogFromDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);
      
      virtual RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

      virtual void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

      /** dump the inverse to the buffer
       */
      void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf);

      /** copy the inverse from the buffer
       */
      void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

      /** return the ratio only for the  iat-th partcle move
       * @param P current configuration
       * @param iat the particle thas is being moved
       */
      virtual ValueType ratio(ParticleSet& P, int iat);


      virtual ValueType alternateRatio(ParticleSet& P)
      {
        return 1.0;
      }
      /** return the ratio
       * @param P current configuration
       * @param iat particle whose position is moved
       * @param dG differential Gradients
       * @param dL differential Laplacians
       *
       * Data member *_temp contain the data assuming that the move is accepted
       * and are used to evaluate differential Gradients and Laplacians.
       */
      virtual ValueType ratio(ParticleSet& P, int iat,
                              ParticleSet::ParticleGradient_t& dG,
                              ParticleSet::ParticleLaplacian_t& dL);

      virtual ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat);
      virtual GradType evalGrad(ParticleSet& P, int iat);
      virtual GradType evalGradSource(ParticleSet &P, ParticleSet &source,
                                      int iat);

      virtual GradType evalGradSource
      (ParticleSet& P, ParticleSet& source, int iat,
       TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
       TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

      virtual GradType evalGradSourcep
      (ParticleSet& P, ParticleSet& source, int iat,
       TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
       TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad);

      virtual ValueType logRatio(ParticleSet& P, int iat,
                                 ParticleSet::ParticleGradient_t& dG,
                                 ParticleSet::ParticleLaplacian_t& dL);

      /** move was accepted, update the real container
       */
      virtual void acceptMove(ParticleSet& P, int iat);

      /** move was rejected. copy the real container to the temporary to move on
       */
      virtual void restore(int iat);

      virtual void update(ParticleSet& P,
                          ParticleSet::ParticleGradient_t& dG,
                          ParticleSet::ParticleLaplacian_t& dL,
                          int iat);

      virtual RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);


      ///evaluate log of determinant for a particle set: should not be called
      virtual RealType
      evaluateLog(ParticleSet& P,
                  ParticleSet::ParticleGradient_t& G,
                  ParticleSet::ParticleLaplacian_t& L) ;


      virtual ValueType
      evaluate(ParticleSet& P,
               ParticleSet::ParticleGradient_t& G,
               ParticleSet::ParticleLaplacian_t& L);

      virtual OrbitalBasePtr makeClone(ParticleSet& tqp) const;

      /** cloning function 
       * @param tqp target particleset
       * @param spo spo set
       *
       * This interface is exposed only to SlaterDet and its derived classes
       * can overwrite to clone itself correctly.
       */
      virtual DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;
//       virtual DiracDeterminantBase* makeCopy(ParticleSet& tqp, SPOSetBase* spo) const {return makeCopy(spo); };

      virtual void get_ratios(ParticleSet& P, vector<ValueType>& ratios);
      ///total number of particles
      int NP;
      ///number of single-particle orbitals which belong to this Dirac determinant
      int NumOrbitals;
      ///number of particles which belong to this Dirac determinant
      int NumPtcls;
      ///index of the first particle with respect to the particle set
      int FirstIndex;
      ///index of the last particle with respect to the particle set
      int LastIndex;
      ///index of the particle (or row)
      int WorkingIndex;
      ///a set of single-particle orbitals used to fill in the  values of the matrix
      SPOSetBasePtr Phi;

      /////Current determinant value
      //ValueType CurrentDet;
      /// psiM(j,i) \f$= \psi_j({\bf r}_i)\f$
      ValueMatrix_t psiM, psiM_temp;

      /// temporary container for testing
      ValueMatrix_t psiMinv;

      /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
      GradMatrix_t  dpsiM, dpsiM_temp;

      /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
      ValueMatrix_t d2psiM, d2psiM_temp;

      /// Used for force computations
      GradMatrix_t grad_source_psiM, grad_lapl_source_psiM;
      HessMatrix_t  grad_grad_source_psiM;
      GradMatrix_t phi_alpha_Minv, grad_phi_Minv;
      ValueMatrix_t lapl_phi_Minv;
      HessMatrix_t grad_phi_alpha_Minv;

      /// value of single-particle orbital for particle-by-particle update
      ValueVector_t psiV;
      GradVector_t dpsiV;
      ValueVector_t d2psiV;
      ValueVector_t workV1, workV2;
      GradVector_t workG;

      Vector<ValueType> WorkSpace;
      Vector<IndexType> Pivot;

      ValueType curRatio,cumRatio;
      ValueType *FirstAddressOfG;
      ValueType *LastAddressOfG;
      ValueType *FirstAddressOfdV;
      ValueType *LastAddressOfdV;
      //    double ComputeExtraTerms(int ptcl_gradient, int elDim,int ionDim);
      ParticleSet::ParticleGradient_t myG, myG_temp;
      ParticleSet::ParticleLaplacian_t myL, myL_temp;
//
      virtual inline void setLogEpsilon(ValueType x) { }

#ifdef QMC_CUDA
    /////////////////////////////////////////////////////
    // Functions for vectorized evaluation and updates //
    /////////////////////////////////////////////////////
    virtual void recompute(MCWalkerConfiguration &W, bool firstTime)
      { cerr << "Need specialization of DiracDetermiantBase::recompute.\n";  abort(); }

    virtual void 
    reserve (PointerPool<gpu::device_vector<CudaRealType> > &pool)
      { cerr << "Need specialization of DiracDetermiantBase::reserve.\n";  abort(); }

    virtual void 
    addLog (MCWalkerConfiguration &W, vector<RealType> &logPsi)
      { cerr << "Need specialization of DiracDetermiantBase::addLog.\n";  abort(); }

    virtual void 
    ratio (MCWalkerConfiguration &W, int iat,
	   vector<ValueType> &psi_ratios,	vector<GradType>  &grad,
	   vector<ValueType> &lapl)
      { cerr << "Need specialization of DiracDetermiantBase::ratio.\n";  abort(); }

    virtual void 
    ratio (vector<Walker_t*> &walkers,    vector<int> &iatList,
	   vector<PosType> &rNew, vector<ValueType> &psi_ratios, 
	   vector<GradType>  &grad, vector<ValueType> &lapl)
      { cerr << "Need specialization of DiracDetermiantBase::ratio.\n";  abort(); }

    virtual void 
    addGradient(MCWalkerConfiguration &W, int iat,
    		vector<GradType> &grad)
      { cerr << "Need specialization of DiracDetermiantBase::addGradient.\n";  abort(); }

    virtual void update (vector<Walker_t*> &walkers, int iat)
      { cerr << "Need specialization of DiracDetermiantBase::update.\n";  abort(); }

    virtual void 
    update (const vector<Walker_t*> &walkers, 
	    const vector<int> &iatList) 
      { cerr << "Need specialization of DiracDetermiantBase::update.\n";  abort(); }
    
    void 
    gradLapl (MCWalkerConfiguration &W, GradMatrix_t &grads,
	      ValueMatrix_t &lapl)
      { cerr << "Need specialization of DiracDetermiantBase::gradLapl.\n";  abort(); }

    virtual void 
    NLratios (MCWalkerConfiguration &W,  vector<NLjob> &jobList,
	      vector<PosType> &quadPoints, vector<ValueType> &psi_ratios)
      { cerr << "Need specialization of DiracDetermiantBase::NLratios.\n";  abort(); }
#endif
    };



}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
