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

#ifndef DIRAC_DETERMINANT_OPT_H
#define DIRAC_DETERMINANT_OPT_H

#include "QMCWaveFunctions/Fermion/DiracDeterminantBase.h"

namespace qmcplusplus
{

class DiracDeterminantOpt : public DiracDeterminantBase
{
protected:
  typedef optimize::VariableSet opt_variables_type;
  opt_variables_type myVars;

  // Basis for optimization
//     SPOSetBasePtr Basis;

  // First index is basis element
  // Second index is the orbital number being optimized
//     Array<ValueType,2> ExcitedCoefs;
  int NumOrbitals, NumBasis;
//     SPOSetBasePtr ExcitedStates;

  // Inverse of Aopt -- not transposed as in DiracDeterminantBase
  ValueMatrix_t AoptInv;
  // Basis functions evaluated at all of my electron positions
  // First index is electron, second index is basis index
  ValueMatrix_t BasisVals, BasisLapl;
  GradMatrix_t  BasisGrad, dgrad_dC, G_gamma;
  // Matrix product of Ainv and BasisVals
  ValueMatrix_t dlogdet_dC;
  // Stores the C_{ij} derivative of (\nabla^2 det)/det
  ValueMatrix_t dlapl_dC;
  // Intermediate used to compute the above.
  ValueMatrix_t Gamma, L_gamma;

  // This vector maps variable number into ExcitedCoefs matrix
  vector<TinyVector<int,2> > VarIndex;
  //
  vector<PosType> MyG;

public:
  DiracDeterminantBase* makeCopy(SPOSetBase* spo) const;

  DiracDeterminantOpt(ParticleSet &ptcl, SPOSetBasePtr const &gs_spos, int first=0);
  // This stores new orbital coefficients and updates the
  // inverse matrices.
  void resetParameters(const opt_variables_type& optvars);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           vector<RealType>& dlogpsi,
                           vector<RealType>& dhpsioverpsi);

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
};

}

#endif
