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
#include "AtomicHF/HFConfiguration.h"
#include "AtomicHF/RadialPotentialSet.h"
#include "OhmmsData/libxmldefs.h"
#include <libxml++/libxml++.h>
#include <fstream>

namespace ohmmshf
{

/**Hartree-Fock solver
 */
struct HartreeFock
{

  typedef SphericalOrbitalTraits::BasisSetType BasisSetType;
  typedef SphericalOrbitalTraits::value_type value_type;
  typedef SphericalOrbitalTraits::RadialOrbital_t RadialOrbital_t;

  int maxiter;
  value_type eig_tol, scf_tol, ratio;
  std::string LogFileName;
  std::vector<value_type> eigVal;

  RadialPotentialSet& Pot;
  BasisSetType& Psi;

  HartreeFock(RadialPotentialSet& pot,
              BasisSetType& psi,
              const xmlpp::Node* root);

  bool put(const xmlpp::Node* q);

  void setRoot(const std::string& aroot);

  template<class Transform_t> void run(Transform_t* fake, int norb);

  void solve( std::string pottype, std::string gridtype, int norb);
};

bool parseXMLFile(RadialPotentialSet& Pot,
                  SphericalOrbitalTraits::BasisSetType& Psi,
                  std::string& name,std::string& pottype,
                  std::string& gridtype,
                  const xmlpp::Node* root);
}

#endif
