#ifndef OHMMS_PSEUDOGEN_HARTREEFOCK_H
#define OHMMS_PSEUDOGEN_HARTREEFOCK_H
#include "SQD/HFConfiguration.h"
#include "SQD/RadialPotentialSet.h"
#include "OhmmsData/libxmldefs.h"
#include "Optimize/Minimize.h"
#include "Optimize/VarList.h"
#include <fstream>

namespace ohmmshf
{

/**Hartree-Fock solver
 */
class PseudoGen : public MinimizeFunction
{

public:

  typedef SphericalOrbitalTraits::BasisSetType BasisSetType;
  typedef SphericalOrbitalTraits::value_type value_type;
  typedef SphericalOrbitalTraits::RadialGrid_t RadialGrid_t;
  typedef SphericalOrbitalTraits::RadialOrbital_t RadialOrbital_t;

  PseudoGen(RadialPotentialSet& pot, BasisSetType& psi);

  bool put(xmlNodePtr cur);
  void setRoot(const string& aroot);
  bool solve();
  int report();
  inline value_type getE(int i) const
  {
    return PPeigVal[i];
  }
  ///assign optimization parameter i
  scalar& Params(int i)
  {
    return OptParams[i];
  }
  ///return optimization parameter i
  scalar Params(int i) const
  {
    return OptParams[i];
  }
  ///return the number of optimizable parameters
  int NumParams()
  {
    return OptParams.size();
  }
  void WriteStuff()
  {
    cout << "calling WriteStuff" << endl;
  }

  double Cost();
  bool optimize();
  bool run();

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
  ///conjugate gradient tolerance
  value_type cg_tolerance;
  ///conjugate gradient stepsize
  value_type cg_stepsize;
  ///conjugate gradient epsilon
  value_type cg_epsilon;
  ///root of all the output files
  string RootFileName;
  ///name of the log file
  string LogFileName;
  string AtomName, PotType, GridType;

  xmlNodePtr grid_ptr;
  xmlNodePtr pot_ptr;
  xmlNodePtr orb_ptr;
  xmlNodePtr opt_ptr;
  ///pointer to the radial grid
  SphericalOrbitalTraits::RadialGrid_t* myGrid;

  ///the radial potential
  RadialPotentialSet& Pot;
  ///the radial orbitals
  BasisSetType& Psi;

  ///All-Electron eigenvalues
  std::vector<value_type> AEeigVal;
  ///Pseudo eigenvalues
  std::vector<value_type> PPeigVal;
  ///parameters to be optimized
  vector<scalar> OptParams;
  ///ID tag for each optimizable parameter
  vector<string> IDtag;
  ///All-Electron orbitals
  vector<RadialOrbital_t> AEorbitals;
  ///All-Electron partial norms
  vector<RadialOrbital_t> AEorbitals_norm;
  ///registry for which variables are added for optimization
  VarRegistry<value_type> vreg;

  bool initGrid();
  bool initHamiltonian();
  bool initOrbitalSet();
  bool initOptimizer();
  bool initAEOrbital(const std::string& grpname,const int& index);

  ///matching radius for evaluating partial norms
  value_type rmatch;
  ///weights for eigenvalues and partial norms in the cost function
  value_type weight_eig, weight_norm;
  ///number of times cost function evaluated
  int NumCostCalls;

  ///pointer to HDF5 file for All-Electron orbitals
  hid_t afile;
  ///pointer to the main group of the HDF5 file
  hid_t group_id;

  bool putOptParams();
  bool getOptParams();

  void plot_ascii();
  void plot_siesta_grid();

  template<class Transform_t> double runHF(Transform_t* fake, int norb);
};

}

#endif
