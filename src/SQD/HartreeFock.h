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
    
    


#ifndef OHMMS_ATOMICHF_HARTREEFOCK_H
#define OHMMS_ATOMICHF_HARTREEFOCK_H
#include "SQD/HFConfiguration.h"
#include "SQD/RadialPotentialSet.h"
#include "OhmmsData/libxmldefs.h"
#include <fstream>

namespace ohmmshf
{

/**Hartree-Fock solver
 */
class HartreeFock
{

public:

  typedef SphericalOrbitalTraits::BasisSetType BasisSetType;
  typedef SphericalOrbitalTraits::value_type value_type;
  typedef SphericalOrbitalTraits::RadialOrbital_t RadialOrbital_t;

  HartreeFock(RadialPotentialSet& pot, BasisSetType& psi);

  bool put(xmlNodePtr cur);
  void setRoot(const std::string& aroot);
  bool solve();
  int report();
  inline value_type getE(int i) const
  {
    return eigVal[i];
  }

private:
  ///maximum number of scf iterations
  int maxiter;
  ///number of closed shells
  int num_closed_shells;
  ///tolerance for the eigen solver
  value_type eig_tol;
  ///tolerance for the self-consistent field solver
  value_type scf_tol;
  ///mixing ratio of self-consistent field
  value_type ratio;
  ///root of all the output files
  std::string RootFileName;
  ///name of the log file
  std::string LogFileName;
  std::string AtomName, PotType, GridType;


  xmlNodePtr grid_ptr;
  xmlNodePtr pot_ptr;
  xmlNodePtr orb_ptr;
  ///pointer to the radial grid
  SphericalOrbitalTraits::RadialGrid_t* myGrid;
  ///container for the eigenvalues
  std::vector<value_type> eigVal;
  ///the radial potential
  RadialPotentialSet& Pot;
  ///the radial orbitals
  BasisSetType& Psi;

  bool initGrid();
  bool initHamiltonian();
  bool initOrbitalSet();

  /** exceute HF-SCF calculations
   * @param norb number of orbitals
   *
   * Transform_t handles boundary (cusp) conditions
   */
  template<typename Transform_t> void run(int norb);

};

}

#endif
