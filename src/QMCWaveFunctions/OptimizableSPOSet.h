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

    ValueVector_t GSVal, BasisVal, GSLapl, BasisLapl;
    GradVector_t  GSGrad, BasisGrad;

    // Cache the positions to avoid recomputing GSVal, BasisVal, etc.
    // if unnecessary.
    vector<PosType> CachedPos;
  public:
    ///set of variables to be optimized;  These are mapped to the
    ///C matrix
    opt_variables_type myVars;

    // For each occupied orbital, this lists which of the M
    // basis functions are used to construct the optimal orbital.
    // This must be initialized by the SPOSet builder.  If it is not
    // initialized, we assume that C is fully dense.
    // The first element of the TinyVector is the basis element
    // number.  The second is the corresponding parameter number.
    vector<vector<TinyVector<int,2> > > ActiveBasis;

    OptimizableSPOSet(int num_orbs, SPOSetBase *gsOrbs, 
		      SPOSetBase* basisOrbs=NULL) :
      GSOrbitals(gsOrbs), BasisOrbitals(basisOrbs)
    {
      N = num_orbs;
      if (BasisOrbitals) {
	M = BasisOrbitals->getOrbitalSetSize();
	GSVal.resize(N);    GSGrad.resize(N);    GSLapl.resize(N);
	BasisVal.resize(M); BasisGrad.resize(M); BasisGrad.resize(N);
      }
      else {
	M = GSOrbitals->getOrbitalSetSize() - N;
	GSVal.resize(N+M);  GSGrad.resize(N+M);  GSLapl.resize(N+M);
      }
    
      C.resize(N,M);
      ActiveBasis.resize(N);
    }

    void set_active_basis (vector<vector<int> > &active)
    { 
      //ActiveBasis = active; 
    }
  
    // Read coefficients, or create the XML element.
    bool put (xmlNodePtr node);

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
    void evaluate(const ParticleSet& P, PosType r, vector<RealType> &psi);
    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, 
		  GradVector_t& dpsi, ValueVector_t& d2psi);
    void evaluate(const ParticleSet& P, int first, int last,
		  ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
		  ValueMatrix_t& d2logdet);
    void evaluate_notranspose(const ParticleSet& P, int first, int last,
			      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
			      ValueMatrix_t& d2logdet);
    void evaluateBasis (const ParticleSet &P, int first, int last,
			ValueMatrix_t &basis_val, GradMatrix_t &basis_grad,
			ValueMatrix_t &basis_lapl);
      

    // Make a copy of myself
    SPOSetBase* makeClone() const;
  };
}

#endif
