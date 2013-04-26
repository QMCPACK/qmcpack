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
  void setRoot(const string& aroot);
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
  string RootFileName;
  ///name of the log file
  string LogFileName;
  string AtomName, PotType, GridType;


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
