//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "SQD/HFConfiguration.h"
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "SQD/SphericalPotential/StepPotential.h"
#include "SQD/SphericalPotential/RegularLinearTransform.h"
#include "SQD/SphericalPotential/RegularLogTransform.h"
#include "SQD/SphericalPotential/NuclearLinearTransform.h"
#include "SQD/SphericalPotential/NuclearLogTransform.h"
#include "SQD/SphericalPotential/NuclearRelLogTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "SQD/HartreeFock.h"

namespace ohmmshf
{

/** Solve a HF eigen problem for a spherically-symmetric external potential.
 * @param norb the number of eigen vectors to be obtained
 */
template<typename Transform_t>
inline void
HartreeFock::run(int norb)
{
  typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
  value_type Vtotal,KEnew, KEold,E;
  value_type lowerbound, upperbound;
  std::vector<value_type> energy(Pot.size());
  eigVal.resize(norb);
  int iter = 0;
  Vtotal = Pot.evaluate(Psi,energy,norb);
  Pot.mix(0.0);
  KEnew = Pot.calcKE(Psi,0,norb);
  std::string label("spdf");
  std::ofstream log_stream(LogFileName.c_str());
  log_stream.precision(8);
  do
  {
    KEold = KEnew;
    value_type eigsum = 0.0;
    //loop over the orbitals
    for(int ob=0; ob < norb; ob++)
    {
      //set the number of nodes of the eigen vector
      Pot.V[ob].setNumOfNodes(Pot.getNumOfNodes(Psi.N[ob],Psi.L[ob]));
      //calculate the lower and upper bounds for the eigenvalue
      Pot.EnergyBound(lowerbound,upperbound);
      //set up the transformer
      Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob],Psi.CuspParam,
                     Pot.getMass());
      //initialize the numerov solver
      Numerov_t numerov(es,Psi(ob));
      //calculate the eigenvalue and the corresponding orbital
      eigsum += (eigVal[ob] =
                   numerov.solve(lowerbound, upperbound, eig_tol));
      log_stream << Psi.N[ob]<< label[Psi.L[ob]] << '\t'
                 << eigVal[ob] << std::endl;
    }
    log_stream << std::endl;
    //normalize the orbitals
    Psi.normalize(norb);
    //restrict the orbitals
    Psi.applyRestriction(norb);
    //calculate the new kinetic energy
    KEnew = Pot.calcKE(Psi,eigsum,norb);
    //the total energy
    //  E = KEnew + Vtotal;
    //for the new orbitals Psi, calculate the new SCF potentials
    Vtotal = Pot.evaluate(Psi,energy,norb);
    //calculate the total energy
    E = KEnew + Vtotal;
    //restrict the potential
    Pot.applyRestriction(Psi);
    //mix the new SCF potential with the old
    Pot.mix(ratio);
    log_stream.precision(10);
    log_stream << "Iteration #" << iter+1 << std::endl;
    log_stream << "KE    = " << std::setw(15) << KEnew
               << "  PE     = " << std::setw(15) << Vtotal << std::endl;
    log_stream << "PE/KE = " << std::setw(15) << Vtotal/KEnew
               << "  Energy = " << std::setw(15) << E << std::endl;
    log_stream << std::endl;
    iter++;
    //continue the loop until the kinetic energy converges
  }
  while(std::abs(KEnew-KEold)>scf_tol && iter<maxiter);
  log_stream << "V_External = " << energy[0] << std::endl;
  log_stream << "V_Hartree = "  << energy[1] << std::endl;
  log_stream << "V_Exchange = " << energy[2] << std::endl;
  log_stream << "E_tot = " << E << std::endl;
}

/** Instantiate a Transformation function based on the potential and grid type and call run.
 */
bool HartreeFock::solve()
{
  int norb = Psi.size();
  if(PotType == "nuclear")
  {
    if(GridType == "linear")
    {
      run<NuclearLinearTransform<RadialOrbital_t> >(norb);
    }
    else if(GridType == "log")
    {
      run<NuclearLogTransform<RadialOrbital_t> >(norb);
    }
  }
  else if(PotType == "nuclear_scalar_rel")
  {
    if(GridType == "log")
    {
      run<NuclearRelLogTransform<RadialOrbital_t> >(norb);
    }
  }
  else
  {
    if(GridType == "linear")
    {
      run<RegularLinearTransform<RadialOrbital_t> >(norb);
    }
    else if(GridType == "log")
    {
      run<RegularLogTransform<RadialOrbital_t> >(norb);
    }
  }
  return true;
}


}

