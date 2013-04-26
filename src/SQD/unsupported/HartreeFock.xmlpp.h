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
  string LogFileName;
  std::vector<value_type> eigVal;

  RadialPotentialSet& Pot;
  BasisSetType& Psi;

  HartreeFock(RadialPotentialSet& pot,
              BasisSetType& psi,
              const xmlpp::Node* root);

  bool put(const xmlpp::Node* q);

  void setRoot(const string& aroot);

  template<class Transform_t> void run(Transform_t* fake, int norb);

  void solve(string pottype, string gridtype, int norb);
};

bool parseXMLFile(RadialPotentialSet& Pot,
                  SphericalOrbitalTraits::BasisSetType& Psi,
                  string& name,string& pottype,
                  string& gridtype,
                  const xmlpp::Node* root);
}

#endif
