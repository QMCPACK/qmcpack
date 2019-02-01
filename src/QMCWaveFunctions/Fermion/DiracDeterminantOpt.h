//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef DIRAC_DETERMINANT_OPT_H
#define DIRAC_DETERMINANT_OPT_H

#include "QMCWaveFunctions/Fermion/DiracDeterminant.h"

namespace qmcplusplus
{

class DiracDeterminantOpt : public DiracDeterminant
{
protected:
  typedef optimize::VariableSet opt_variables_type;
  opt_variables_type myVars;

  // Basis for optimization
//     SPOSetPtr Basis;

  // First index is basis element
  // Second index is the orbital number being optimized
//     Array<ValueType,2> ExcitedCoefs;
  int NumOrbitals, NumBasis;
//     SPOSetPtr ExcitedStates;

  // Inverse of Aopt -- not transposed as in DiracDeterminant
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
  std::vector<TinyVector<int,2> > VarIndex;
  //
  std::vector<PosType> MyG;

public:
  DiracDeterminant* makeCopy(SPOSet* spo) const;

  DiracDeterminantOpt(ParticleSet &ptcl, SPOSetPtr const &gs_spos, int first=0);
  // This stores new orbital coefficients and updates the
  // inverse matrices.
  void resetParameters(const opt_variables_type& optvars);

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi);

  void checkInVariables(opt_variables_type& active);
  void checkOutVariables(const opt_variables_type& active);
};

}

#endif
