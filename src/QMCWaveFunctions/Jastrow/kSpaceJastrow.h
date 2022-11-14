//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file kSpaceJastrow.h
 * @brief Declaration of Long-range TwoBody Jastrow
 *
 * The initial coefficients are evaluated according to RPA.
 * A set of k-dependent coefficients are generated.
 * Optimization should use rpa_k0, ...rpa_kn as the names of the
 * optimizable parameters.
 */
#ifndef QMCPLUSPLUS_LR_KSPACEJASTROW_H
#define QMCPLUSPLUS_LR_KSPACEJASTROW_H

#include "QMCWaveFunctions/WaveFunctionComponent.h"
#include "OhmmsData/libxmldefs.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "LongRange/LRHandlerBase.h"

namespace qmcplusplus
{
/** Functor which return \f$frac{Rs}{k^2 (k^2+(1/Rs)^2)}\f$
 */
template<typename T>
struct RPA0
{
  T Rs;
  T OneOverRsSq;
  RPA0(T rs = 1) : Rs(rs) { OneOverRsSq = 1.0 / (Rs * Rs); }
  inline T operator()(T kk)
  {
    T k2 = std::sqrt(kk);
    return Rs / (k2 * (k2 + OneOverRsSq));
    //return (-0.5+0.5*std::pow(1.0+12.0*OneOverRs3/kk/kk,0.5));
  }
};


template<typename T>
class kSpaceCoef
{
public:
  // The coefficient
  T cG;
  // The first and last indices of the G-vectors to which this
  // coefficient belongs
  int firstIndex, lastIndex;
  inline void set(std::vector<T>& unpackedCoefs)
  {
    for (int i = firstIndex; i <= lastIndex; i++)
      unpackedCoefs[i] = cG;
  }
};

class kSpaceJastrow : public WaveFunctionComponent, public OptimizableObject
{
public:
  typedef enum
  {
    CRYSTAL,
    ISOTROPIC,
    NOSYMM
  } SymmetryType;

private:
  using ComplexType = std::complex<RealType>;
  ////////////////
  // Basic data //
  ////////////////
  RealType CellVolume, NormConstant;
  int num_elecs;
  int NumSpins;
  int NumIons, NumIonSpecies;

  // Gvectors included in the summation for the
  // one-body and two-body Jastrows, respectively
  std::vector<PosType> OneBodyGvecs;
  std::vector<PosType> TwoBodyGvecs;

  // Vectors of unique coefficients, one for each group of
  // symmetry-equivalent G-vectors
  std::vector<kSpaceCoef<ComplexType>> OneBodySymmCoefs;
  std::vector<kSpaceCoef<RealType>> TwoBodySymmCoefs;

  // Vector of coefficients, one for each included G-vector
  // in OneBodyGvecs/TwoBodyGvecs
  std::vector<ComplexType> OneBodyCoefs;
  std::vector<RealType> TwoBodyCoefs;

  // Stores how we choose to group the G-vectors by symmetry
  SymmetryType OneBodySymmType, TwoBodySymmType;

  // Rho_k data
  // For the ions
  std::vector<ComplexType> Ion_rhoG;
  // For the electrons
  std::vector<ComplexType> OneBody_rhoG, TwoBody_rhoG;
  // Holds the phase for the electrons with size commensurate with
  // OneBodyGvecs, and TwoBodyGvecs, respectively
  std::vector<RealType> OneBodyPhase, TwoBodyPhase;
  //
  std::vector<ComplexType> OneBody_e2iGr, TwoBody_e2iGr_new, TwoBody_e2iGr_old;
  Matrix<ComplexType> Delta_e2iGr;

  // Map of the optimizable variables:
  //std::map<std::string,RealType*> VarMap;

  //////////////////////
  // Member functions //
  //////////////////////

  // Enumerate G-vectors with nonzero structure factors
  void setupGvecs(RealType kcut, std::vector<PosType>& gvecs, bool useStructFact);
  void setupCoefs();

  // Sort the G-vectors into appropriate symmtry groups
  template<typename T>
  void sortGvecs(std::vector<PosType>& gvecs, std::vector<kSpaceCoef<T>>& coefs, SymmetryType symm);

  // Returns true if G1 and G2 are equivalent by crystal symmetry
  bool Equivalent(PosType G1, PosType G2);
  void StructureFactor(PosType G, std::vector<ComplexType>& rho_G);

  const ParticleSet& Ions;
  std::string OneBodyID;
  std::string TwoBodyID;
  double Prefactor;

public:
  kSpaceJastrow(const ParticleSet& ions,
                ParticleSet& elecs,
                SymmetryType oneBodySymm,
                RealType oneBodyCutoff,
                std::string onebodyid,
                bool oneBodySpin,
                SymmetryType twoBodySymm,
                RealType twoBodyCutoff,
                std::string twobodyid,
                bool twoBodySpin);

  kSpaceJastrow(const ParticleSet& ions);

  void setCoefficients(std::vector<RealType>& oneBodyCoefs, std::vector<RealType>& twoBodyCoefs);

  std::string getClassName() const override { return "kSpaceJastrow"; }
  //implement virtual functions for optimizations
  bool isOptimizable() const override { return true; }
  void checkOutVariables(const opt_variables_type& active) override;

  void extractOptimizableObjectRefs(UniqueOptObjRefs& opt_obj_refs) override { opt_obj_refs.push_back(*this); }

  void checkInVariablesExclusive(opt_variables_type& active) final;
  void resetParametersExclusive(const opt_variables_type& active) final;

  LogValueType evaluateLog(const ParticleSet& P,
                           ParticleSet::ParticleGradient& G,
                           ParticleSet::ParticleLaplacian& L) override;

  PsiValueType ratio(ParticleSet& P, int iat) override;

  GradType evalGrad(ParticleSet& P, int iat) override;
  PsiValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat) override;

  void restore(int iat) override;
  void acceptMove(ParticleSet& P, int iat, bool safe_to_delay = false) override;

  // Allocate per-walker data in the PooledData buffer
  void registerData(ParticleSet& P, WFBufferType& buf) override;
  // Walker move has been accepted -- update the buffer
  LogValueType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch = false) override;
  // Pull data from the walker buffer at the beginning of a block of
  // single-particle moves
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf) override;

  ///process input file
  bool put(xmlNodePtr cur);
  ///output jastrow coefficients
  void printOneBody(std::ostream& os);
  void printTwoBody(std::ostream& os);

  // Implements strict weak ordering with respect to the
  // structure factors.  Used to sort the G-vectors according to
  // crystal symmetry
  bool operator()(PosType G1, PosType G2);
  std::unique_ptr<WaveFunctionComponent> makeClone(ParticleSet& tqp) const override;

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& active,
                           Vector<ValueType>& dlogpsi,
                           Vector<ValueType>& dhpsioverpsi) override;

  /** evaluate the ratio
  */
  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios) override;

private:
  void copyFrom(const kSpaceJastrow& old);
  std::vector<int> TwoBodyVarMap;
  std::vector<int> OneBodyVarMap;
};
} // namespace qmcplusplus
#endif
