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
#include "Numerics/Clebsch_Gordan.h"
#include "SQD/fillshells.h"
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "SQD/SphericalPotential/StepPotential.h"
#include "SQD/SphericalPotential/SJPseudoPotential.h"
#include "SQD/SphericalPotential/SHEGPotential.h"
#include "SQD/HartreeFock.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"

namespace ohmmshf
{

/** The contructor.
 *@param pot the potential
 * @param psi the wavefunction
 */
HartreeFock::HartreeFock(RadialPotentialSet& pot,
                         SphericalOrbitalTraits::BasisSetType& psi):
  Pot(pot), Psi(psi),
  num_closed_shells(0),
  maxiter(1000), eig_tol(1e-12),
  scf_tol(1e-8), ratio(0.35),
  GridType("none"),
  grid_ptr(NULL), orb_ptr(NULL),
  pot_ptr(NULL), myGrid(0) { }

/** Sets the root for all output files.
 *@param aroot the root for all output files
 */

void HartreeFock::setRoot(const std::string& aroot)
{
  RootFileName = aroot;
  LogFileName = RootFileName + ".log";
}

/** Set the parameters for the eigen solver
 *@param q the current xml node which contains the parameter definitions for the eigen solver
 *@return true if successful
 *
 *Available parameters
 <ul>
 <li> max_iter, the maximum self-consistent iterations, default=1000
 <li> eig_tol, the eigen-value tolerance, default = \f$1 e^{-12}\f$
 <li> en_tol, the tolerance of the self-consistent loops,
 default = \f$1 e^{-8}\f$
 <li> mix_ratio, the mixing ratio of the charge density, default = 0.35
 </ul>
 */
bool HartreeFock::put(xmlNodePtr d_root)
{
  xmlNodePtr cur = NULL;
  //using XPath instead of recursive search
  xmlXPathContextPtr m_context = xmlXPathNewContext(d_root->doc);
  xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//atom",m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    ERRORMSG("Missing Atom information. Exit")
    return false;
  }
  else
  {
    cur = result->nodesetval->nodeTab[0];
    xmlAttrPtr att = cur->properties;
    while(att != NULL)
    {
      std::string aname((const char*)(att->name));
      const char *avalue = (const char*)(att->children->content);
      if(aname == "name")
      {
        AtomName = avalue;
      }
      else if (aname == "num_closed_shells")
      {
        num_closed_shells = atoi(avalue);
      }
      att = att->next;
    }
    XMLReport("Atom name = " << AtomName);
    XMLReport("Number of closed shells = " << num_closed_shells);
    //initialize xmlNode pointers for the grid, wavefunction and potential
    xmlNodePtr cur1 = cur->xmlChildrenNode;
    while(cur1 != NULL)
    {
      std::string cname1((const char*)(cur1->name));
      if(cname1 == "grid")
      {
        grid_ptr = cur1;
      }
      else if(cname1 == "orbitalset")
      {
        orb_ptr = cur1;
      }
      else if(cname1 == "hamiltonian")
      {
        pot_ptr = cur1;
        OhmmsAttributeSet aAttrib;
        aAttrib.add(PotType,"type");
        aAttrib.put(cur1);
      }
      cur1 = cur1->next;
    }
  }
  if(orb_ptr == NULL || pot_ptr == NULL || grid_ptr == NULL)
  {
    ERRORMSG("Missing one of nodes: <grid/>, <orbitalset/> or <hamiltonian/>. Exit");
    return false;
  }
  xmlXPathFreeObject(result);
  //initialize the grid
  bool success = initGrid();
  if(!success)
  {
    ERRORMSG("Failed to create a grid");
    return false;
  }
  Psi.m_grid = myGrid;
  //initialize the wavefunction
  success = initOrbitalSet();
  if(!success)
  {
    ERRORMSG("Failed to create/initialiaze orbitals");
    return false;
  }
  //initialize the internal storage for a potential
  Pot.initialize(Psi);
  //initialize the Hamiltonian
  success = initHamiltonian();
  if(!success)
  {
    ERRORMSG("Failed to create/initialiaze hamiltonian");
    return false;
  }
  //initialize the eigen solver
  result = xmlXPathEvalExpression((const xmlChar*)"//eigensolve",m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    WARNMSG("Using default values for eigen solver");
  }
  else
  {
    ParameterSet params;
    params.add(maxiter,"max_iter","int");
    params.add(eig_tol,"eig_tol","double");
    params.add(scf_tol,"etot_tol","double");
    params.add(ratio,"mix_ratio","double");
    params.put(result->nodesetval->nodeTab[0]);
  }
  xmlXPathFreeObject(result);
  XMLReport("maximum iterations = " << maxiter);
  XMLReport("eigentolerance = " << eig_tol);
  XMLReport("scftolerance = " << scf_tol);
  XMLReport("ratio = " << ratio);
  xmlXPathFreeContext(m_context);
  return true;
}

/** Initialize the radial grid.
 *
 *Available parameters
 <ul>
 <li> type: the grid type (log or linear)
 <li> scale: the scaling factor, default 1.0
 <li> min: the minimum value of the grid, default 0.001
 <li> max: the maximum value of the grid, default 1000.0
 <li> npts: the number of grid points, default 2001
 </ul>
 */

bool HartreeFock::initGrid()
{
  xmlNodePtr cur = grid_ptr;
  double scale = 1.0;
  double rmin = 0.001;
  double rmax = 1000.0;
  int npts = 2001;
  //use common method
  OhmmsAttributeSet aAttrib;
  aAttrib.add(scale,"scale");
  aAttrib.add(GridType,"type");
  aAttrib.add(rmin,"ri");
  aAttrib.add(rmax,"rf");
  aAttrib.add(npts,"npts");
  aAttrib.put(cur);
  //use old method via <parameter/>
  XMLReport("Grid type = " << GridType);
  ParameterSet params;
  params.add(rmin,"min","double");
  params.add(rmax,"max","double");
  params.add(npts,"npts","int");
  params.put(cur->xmlChildrenNode);
  if(GridType == "log")
  {
    myGrid = new LogGrid<double>;
    rmin/=scale;
    rmax/=sqrt(scale);
    myGrid->set(rmin,rmax,npts);
  }
  else if(GridType == "linear")
  {
    myGrid = new LinearGrid<double>;
    rmin/=scale;
    rmax/=sqrt(scale);
    myGrid->set(rmin,rmax,npts);
  }
  else
  {
    ERRORMSG("Grid Type Options: Log or Linear.");
    return false;
  }
  XMLReport("Radial Grid: type " << GridType << " rmin = "
            << rmin << ", rmax = " << rmax << ", npts = " << npts);
  //     << ", dh = " << myGrid->dh() << ", npts = " << npts)
  return myGrid != 0;
}

/** Initialize the wavefunction.
 *
 *Available parameters
 <ul>
 <li> rest_type: the restriction type (spin, spin_space, none)
 </ul>
 */

bool HartreeFock::initOrbitalSet()
{
  //find the restriction type
  xmlNodePtr cur = orb_ptr;
  if(xmlHasProp(cur,(const xmlChar*)"condition"))
  {
    Psi.Restriction = (const char*)(xmlGetProp(cur, (const xmlChar *) "condition"));
    XMLReport("Orbital restriction type = " << Psi.Restriction);
  }
  //fill the closed shells (shell filling depends of potential)
  if(num_closed_shells)
  {
    if(PotType.find("nuclear")<PotType.size())
    {
      FillShellsNucPot(Psi,num_closed_shells);
    }
    if(PotType == "harmonic" || PotType == "step" || PotType == "heg")
    {
      FillShellsHarmPot(Psi,num_closed_shells);
    }
    //if(PotType == "heg") {
    //  FillShellsHarmPot(Psi,num_closed_shells,1);
    //}
  }
  //pass the xmlNode to Psi
  Psi.put(cur);
  app_log() << "Total number of orbitals = " << Psi.size() << std::endl;
  app_log() << "Total number of unique radial orbitals = " <<  Psi.NumUniqueOrb << std::endl;
  app_log() << "(Orbital index, Number of Orbitals, n,l,m,s)"<< std::endl;;
  for(int j=0; j < Psi.size(); j++)
  {
    int id = Psi.ID[j];
    app_log()<< "(" << id << ", " << Psi.IDcount[id]
             << ", " << Psi.N[j]
             << ", " << Psi.L[j]
             << ", " << Psi.M[j]
             << ", " << Psi.S[j]
             << ")" << std::endl;
  }
  //return false if there is no wave functions
  return Psi.size() != 0;
}

/** Initialize the Hamiltonian.
 *
 *Available Hamiltonians
 <ul>
 <li> Nuclear
 <li> Harmonic
 <li> Step
 <li> Pseudo
 <li> Nuclear Scalar Relativistic
 </ul>
 */
bool HartreeFock::initHamiltonian()
{
  //set charge to zero
  Psi.Charge=0;

  if(PotType.find("nuclear")<PotType.size())
  {
    XMLReport("Creating a Nuclear Potential.");
    double Z = Psi.size();
    ParameterSet params;
    params.add(Z,"z","double");
    params.put(pot_ptr);
    XMLReport("Potential Paramter: Z = " << Z);
    Pot.add(new ZOverRPotential(Z), true);
    Psi.CuspParam = Z;
    Psi.Charge=Z;
  } // Z/r potential
  else if(PotType=="heg")
  {
    XMLReport("Creating a Spherical Jellium potential.");
    SHEGPotential* ap=new SHEGPotential(Psi.size());
    ap->put(pot_ptr);
    Pot.add(ap);
    Psi.CuspParam = 0.0;//set cusp to zero
  }//spherical Jellium
  else if(PotType == "harmonic")
  {
    XMLReport("Creating a Harmonic Potential.");
    double Omega=0.5;
    ParameterSet params;
    params.add(Omega,"z","double");
    params.add(Omega,"omega","double");
    params.put(pot_ptr);
    XMLReport("Potential Parameter: Omega = " << Omega);
    Pot.add(new HarmonicPotential(Omega),true);
    Psi.CuspParam = 0.0;
  } //Harmonic
  else if(PotType == "step")
  {
    XMLReport("Creating a Step-Function Potential.");
    StepPotential* apot = new StepPotential;
    apot->put(pot_ptr);
    Pot.add(apot,true);
    Psi.CuspParam = 0.0;
  } //step potential
  else if(PotType == "pseudo")
  {
    XMLReport("Creating a Starkloff-Joannopoulos Pseudo-Potential.");
    double rc, lambda, Zeff;
    ParameterSet params;
    params.add(rc,"rc","double");
    params.add(lambda,"lambda","double");
    params.add(Zeff,"zeff","double");
    params.put(pot_ptr);
    XMLReport("Potential Parameter: Effective Charge = " << Zeff);
    XMLReport("Potential Parameter: Core Radius = " << rc);
    XMLReport("Potential Parameter: Lambda = " << lambda);
    Pot.add(new SJPseudoPotential(Zeff,rc,lambda));
    Psi.CuspParam = 0.0;
    Psi.Charge=Zeff;
  } // if Pseudo
  else
  {
    ERRORMSG("Unknown potential type" << PotType);
    return false;
  }
  //determine maximum angular momentum
  int lmax = 0;
  int norb = Psi.size();
  for(int i=0; i<norb; i++)
    if(Psi.L[i] > lmax)
      lmax = Psi.L[i];
  XMLReport("Maximum Angular Momentum = " << lmax);
  //add Hartree and Exchange common to all
  Clebsch_Gordan* CG_coeff = new Clebsch_Gordan(lmax);
  Pot.add(new HartreePotential(CG_coeff,norb));
  Pot.add(new ExchangePotential(CG_coeff,norb));
  double m=1.0;
  ParameterSet params1;
  params1.add(m,"mass","double");
  params1.put(pot_ptr);
  Pot.setMass(m);
  XMLReport("Effective mass = " << Pot.getMass());
  return true;
}

/** Print the orbitals and their corresponding eigenvalues to a file "AtomName.orb.dat".
 *@return the maximum radius of the orbital set
 *
 *The orbitals can easily be plotted by using gnuplot or xmgrace.
 */
int HartreeFock::report()
{
  std::string fileforplot(RootFileName);
  fileforplot.append(".orb.dat");
  std::ofstream fout(fileforplot.c_str());
  fout << "#Results for " << AtomName << " with " << PotType
       << " potential on " << GridType << " grid." << std::endl;
  fout << "#Eigen values " << std::endl;
  fout.precision(10);
  fout.setf(std::ios::scientific,std::ios::floatfield);
  for(int orb=0; orb<eigVal.size(); orb++)
  {
    fout << "# n=" << Psi.N[orb] << " " <<  " l=" << Psi.L[orb]
         << std::setw(25) << eigVal[orb] << " " << std::endl;
  }
  fout << "#The number of unique radial orbitals " << Psi.NumUniqueOrb << std::endl;
  //find the maximum radius for the orbital set
  int max_rad_all = 0;
  int orbindex = 0;
  for(int orb=0; orb<Psi.NumUniqueOrb; orb++)
  {
    int max_rad = myGrid->size()-1;
    while(std::abs(Psi(orbindex,max_rad)) < 1e-12)
      max_rad--;
    max_rad += 2;
    max_rad_all = std::max(max_rad_all,max_rad);
    orbindex += Psi.IDcount[orb];
  }
  fout.precision(12);
  for(int ig=0; ig<max_rad_all; ig++)
  {
    fout << std::setw(22) << myGrid->r(ig);
    orbindex=0;
    for(int orb=0; orb<Psi.NumUniqueOrb; orb++)
    {
      fout << std::setw(22) << Psi(orbindex,ig);
      orbindex += Psi.IDcount[orb];
    }
    fout << std::endl;
  }
  Psi.print_HDF5(RootFileName,GridType,eigVal);
  Psi.print_basis(AtomName,RootFileName,GridType);
  Pot.getStorage(Psi,RootFileName);
  return max_rad_all;
}
}




