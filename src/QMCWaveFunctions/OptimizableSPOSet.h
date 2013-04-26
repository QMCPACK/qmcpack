//////////////////////////////////////////////////////////////////
// (c) Copyright 2010-  by Ken Esler and Jeongnim Kim
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: esler@uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////

#ifndef OPTIMIZABLE_SPO_SET_H
#define OPTIMIZABLE_SPO_SET_H

#include "QMCWaveFunctions/SPOSetBase.h"
#include "Numerics/OptimizableFunctorBase.h"
#include "Optimize/VariableSet.h"

namespace qmcplusplus
{

class OptimizableSPOSet : public SPOSetBase
{
protected:
  typedef optimize::VariableSet opt_variables_type;
  ///typedef for name-value lists
  typedef optimize::VariableSet::variable_map_type variable_map_type;

  // Number of occupied states, number of basis states
  int N, M;

  RealType derivScale;
  RealType thr;

  // If BasisOrbitals==NULL, only GSOrbitals is used and it's evaluate
  // functions should return N+M orbitals.
  SPOSetBase *GSOrbitals, *BasisOrbitals;

  // The first index is the orbital to be optimized, the second is the
  // basis element
  //    ValueMatrix_t OptCoefs;

  // This maps the parameter indices to elements in the C
  // matrix.  Using pointers allows the use of complex C
  // while maintaining real optimizable parameters.
  vector<RealType*> ParamPointers;
  // Maps the parameter index in myVars to an index in the C array
  vector<TinyVector<int,2> > ParamIndex;

  void addParameter (string id, int iorb, int basis);

  ValueVector_t GSVal, BasisVal, GSLapl, BasisLapl;
  GradVector_t  GSGrad, BasisGrad;

  ValueMatrix_t GSValMatrix, BasisValMatrix, GSLaplMatrix, BasisLaplMatrix, GradTmpSrc, GradTmpDest;
  GradMatrix_t  GSGradMatrix, BasisGradMatrix;

  // Cache the positions to avoid recomputing GSVal, BasisVal, etc.
  // if unnecessary.
  vector<PosType> CachedPos;
public:
  ///set of variables to be optimized;  These are mapped to the
  ///C matrix.  Update:  Moved to SPOSetBase
  // opt_variables_type myVars;

  // For each occupied orbital, this lists which of the M
  // basis functions are used to construct the optimal orbital.
  // This must be initialized by the SPOSet builder.  If it is not
  // initialized, we assume that C is fully dense.
  // The first element of the TinyVector is the basis element
  // number.  The second is the corresponding parameter number.
  vector<vector<TinyVector<int,2> > > ActiveBasis;


  OptimizableSPOSet() : N(0), M(0), derivScale(10.0), thr(0.0), GSOrbitals(0), BasisOrbitals(0)
  {
    Optimizable = true;
  }

  OptimizableSPOSet(int num_orbs, SPOSetBase *gsOrbs,
                    SPOSetBase* basisOrbs=0) :
    GSOrbitals(gsOrbs), BasisOrbitals(basisOrbs), derivScale(10.0)
  {
    N = num_orbs;
    setOrbitalSetSize(N);
    if (BasisOrbitals)
      M = BasisOrbitals->getOrbitalSetSize();
    else
      M = GSOrbitals->getOrbitalSetSize() - N;
    resize(N,M);
    Optimizable = true;
  }



  //    bool put(xmlNodePtr cur, SPOPool_t &spo_pool);
  void resetTargetParticleSet(ParticleSet& P);
  void setOrbitalSetSize(int norbs);

  void set_active_basis (vector<vector<int> > &active)
  {
    //ActiveBasis = active;
  }

  // Read coefficients, or create the XML element.
  bool put (xmlNodePtr node, SPOPool_t &spo_pool);

  // This stores new orbital coefficients int C
  void resetParameters(const opt_variables_type& optvars);

  // Evaluate the derivative of the optimized orbitals with
  // respect to the parameters
  void evaluateDerivatives(ParticleSet& P, int iat,
                           const opt_variables_type& active,
                           ValueMatrix_t& d_phi,
                           ValueMatrix_t& d_lapl_phi);

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);


  // Evaluate routines.  These call GSOrbitals->evaluate and possibly
  // BasisOrbitals->evaluate, then does the matrix product with
  // C.
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);
  void evaluate(const ParticleSet& P, const PosType& r, vector<RealType> &psi);
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi,
                GradVector_t& dpsi, ValueVector_t& d2psi);
  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi,
                GradVector_t& dpsi, HessVector_t& d2psi);
  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet,
                            ValueMatrix_t& d2logdet);
  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of OptimizableOrbitalSet::evaluate_notranspose() for grad_grad_logdet. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Need specialization of OptimazableOrbitalSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
  }

  void evaluateBasis (const ParticleSet &P, int first, int last,
                      ValueMatrix_t &basis_val, GradMatrix_t &basis_grad,
                      ValueMatrix_t &basis_lapl);
  void copyParamsFromMatrix (const opt_variables_type& active,
                             const Matrix<RealType> &mat,
                             vector<RealType> &destVec);
  void copyParamsFromMatrix (const opt_variables_type& active,
                             const Matrix<ComplexType> &mat,
                             vector<RealType> &destVec);


  void resize(int n, int m)
  {
    N=n;
    M=m;
    if (BasisOrbitals)
    {
      GSVal.resize(N);
      GSGrad.resize(N);
      GSLapl.resize(N);
      BasisVal.resize(M);
      BasisGrad.resize(M);
      BasisLapl.resize(M);
      GSValMatrix.resize (N,N);
      GSGradMatrix.resize(N,N);
      GSLaplMatrix.resize(N,N);
      BasisValMatrix.resize (M,N);
      BasisGradMatrix.resize(M,N);
      BasisLaplMatrix.resize(M,N);
    }
    else
    {
      GSVal.resize(N+M);
      GSGrad.resize(N+M);
      GSLapl.resize(N+M);
      GSValMatrix.resize(N,N+M);
      GSGradMatrix.resize(N,N+M);
      GSLaplMatrix.resize(N,N+M);
    }
    GradTmpSrc.resize(M,N);
    GradTmpDest.resize(N,N);
    C.resize(N,M);
    ActiveBasis.resize(N);
    BasisSetSize = M;
    ///from inherited class
    t_logpsi.resize(N,N);
  }

  // Make a copy of myself
  SPOSetBase* makeClone() const;
};
}

#endif
