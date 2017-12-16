//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#ifndef QMCPLUSPLUS_CONVERT2_RADIALGRID_H
#define QMCPLUSPLUS_CONVERT2_RADIALGRID_H

#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCTools/MolecularOrbitalBasis.h"
#include "QMCTools/RGFBuilderBase.h"

namespace qmcplusplus
{

/** derived class from OrbitalBuilderBase
 *
 * Create a basis set of molecular orbital types as defined in MolecularOrbitalBasis
 * with radial wave functions on the radial grids.
 */
class GridMolecularOrbitals: public OrbitalBuilderBase
{

public:

  //@typedef radial grid type
  typedef OneDimGridBase<RealType>                       GridType;
  //@typedef \f$R_{nl}\f$ radial functor type defined on a grid
  typedef OneDimGridFunctor<RealType>                    RadialOrbitalType;
  //@typedef centered orbital type which represents a set of \f$R_{nl}Y_{lm}\f$  for a center
  typedef SphericalOrbitalSet<RadialOrbitalType,GridType> CenteredOrbitalType;
  //@typedef molecuar orbital basis composed of multiple CenteredOrbitalType s
  typedef MolecularOrbitalBasis<CenteredOrbitalType>      BasisSetType;

  /** constructor
   * \param els reference to the electrons
   * \param psi reference to the wavefunction
   * \param ions reference to the ions
   */
  GridMolecularOrbitals(ParticleSet& els, TrialWaveFunction& psi, ParticleSet& ions);

  /** initialize the Antisymmetric wave function for electrons
   *@param cur the current xml node
   *
   */
  bool put(xmlNodePtr cur);

  /** process basis element to build the basis set
   *@param cur the current xml node
   *@return a pointer to the BasisSet
   *
   *This member function is necessary in order to use SDSetBuilderWithBasisSet
   *to initialize fermion wavefunctions with slater  determinants.
   */
  BasisSetType* addBasisSet(xmlNodePtr cur);

private:

  enum {DONOT_EXPAND=0, GAUSSIAN_EXPAND=1, NATURAL_EXPAND};

  ///reference to the ionic system with which the basis set are associated
  ParticleSet& IonSys;

  ///pointer to the BasisSet built by GridMolecularOrbitals
  BasisSetType*      BasisSet;

  ///distance table pointer
  DistanceTableData* d_table;

  ///Current RadiaaGridFunctorBuilder
  RGFBuilderBase* rbuilder;
  ///map for the radial orbitals
  std::map<std::string,int>    RnlID;
  ///map for the centers
  std::map<std::string,int>    CenterID;

  ///map for (n,l,m,s) to its quantum number index
  std::map<std::string,int> nlms_id;

  ///append Ylm channels
  int expandYlm(const std::string& rnl, const QuantumNumberType& nlms, int num,
                CenteredOrbitalType* aos, xmlNodePtr cur1,
                int expandlm=DONOT_EXPAND);

};
}
#endif
