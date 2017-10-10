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
    
    
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include "QMCHamiltonians/LocalPPotential.h"
#include "Numerics/OneDimGridBase.h"
#include "Numerics/OneDimGridFunctor.h"
#include "Numerics/OneDimCubicSpline.h"
//#include "Numerics/HDFNumericAttrib.h"
#include "Utilities/SimpleParser.h"
#include "Utilities/OhmmsInfo.h"

namespace qmcplusplus
{

LocalPPotential::RadialPotentialSet::~RadialPotentialSet()
{
  for(int ig=0; ig<grid_m.size(); ig++)
    delete grid_m[ig];
  for(int ip=0; ip<lpp_m.size(); ip++)
    delete lpp_m[ip];
}

/*!
 *\param ions the positions of the ions
 *\param els the positions of the electrons
 *\brief the constructor
 *
 * For each ion-type, an ASCII file "*.psf" must
 be provided.  The "*.psf" must contain two columns,
 the first column being the grid, the second being
 the potential on the grid.
*/

LocalPPotential::LocalPPotential(ParticleSet& ions, ParticleSet& els):
  Centers(ions.GroupID), d_table(NULL)
{
  d_table = DistanceTable::getTable(DistanceTable::add(ions,els));
  std::vector<std::string> vlist;
  const SpeciesSet& Species(ions.getSpeciesSet());
  for(int i=0; i< Species.getTotalNum(); i++)
  {
    std::vector<double> grid_temp, pp_temp;
    std::string species(Species.speciesName[i]);
    std::string fname = species+".psf";
    std::ifstream fin(fname.c_str(),std::ios_base::in);
    if(!fin)
    {
      ERRORMSG("Could not open file " << fname)
      exit(-1);
    }
    XMLReport("Reading a file for the PseudoPotential for species " << species)
    double r, f1;
    while(fin >> r >> f1)
    {
      grid_temp.push_back(r);
      pp_temp.push_back(f1);
    }
    int npts = grid_temp.size();
    //add a new local potential to the list
    PP.push_back(new RadialPotentialSet);
    GridType *agrid = new NumericalGrid<double>(grid_temp);
    LocalPotentialType* app = new OneDimCubicSpline<double>(agrid,pp_temp);
    fin.close();
    //calculate the cusp condition
    int imin = 1;
    RealType yprime_i = ((*app)(imin+1)-(*app)(imin))/app->dr(imin);
    app->spline(imin,yprime_i,app->size()-1,0.0);
    PP[i]->add(agrid,app);
  }
} // Centers

///destructor
LocalPPotential::~LocalPPotential()
{
  for(int pp=0; pp<PP.size(); pp++)
    delete PP[pp];
}
}

