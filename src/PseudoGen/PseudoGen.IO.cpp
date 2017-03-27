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
    
    


#include "SQD/HFConfiguration.h"
#include "Numerics/Clebsch_Gordan.h"
#include "Numerics/HDFNumericAttrib.h"
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "SQD/SphericalPotential/SJPseudoPotential.h"
#include "PseudoGen/PseudoGen.h"

namespace ohmmshf
{

/**
   @param q the current xml node which contains the parameter definitions
   *for the eigen solver
   @return true if successful
   @brief Set the parameters for the eigen solver
   *
   *Available parameters
   <ul>
   <li> max_iter, the maximum self-consistent iterations, default=1000
   <li> eig_tol, the eigen-value tolerance, default = \f$1 e^{-12}\f$
   <li> en_tol, the tolerance of the self-consistent loops,
   *default = \f$1 e^{-8}\f$
   <li> mix_ratio, the mixing ratio of the charge density, default = 0.35
   </ul>
*/
bool PseudoGen::put(xmlNodePtr d_root)
{
  xmlNodePtr cur = NULL;
  //using XPath instead of recursive search
  xmlXPathContextPtr m_context = xmlXPathNewContext(d_root->doc);
  xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//atom",
                             m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
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
      else
        if (aname == "num_closed_shells")
        {
          num_closed_shells = atoi(avalue);
        }
      att = att->next;
    }
    //initialize xmlNode pointers for the grid, wavefunction and potential
    xmlNodePtr cur1 = cur->xmlChildrenNode;
    while(cur1 != NULL)
    {
      std::string cname1((const char*)(cur1->name));
      if(cname1 == "grid")
      {
        grid_ptr = cur1;
      }
      else
        if(cname1 == "orbitalset")
        {
          orb_ptr = cur1;
        }
        else
          if(cname1 == "potential")
          {
            pot_ptr = cur1;
            if(xmlHasProp(cur1,(const xmlChar*)"type"))
            {
              PotType = (const char*)(xmlGetProp(cur1, (const xmlChar *) "type"));
              XMLReport("The type of potential " << PotType)
            }
            else
            {
              ERRORMSG("Potential type is undefined. Exit")
              return false;
            }
          }
      cur1 = cur1->next;
    }
  }
  if(orb_ptr == NULL || pot_ptr == NULL || grid_ptr == NULL)
  {
    ERRORMSG("Missing one of nodes: <grid/>, <orbitalset/> or <potential/>. Exit")
    return false;
  }
  xmlXPathFreeObject(result);
  result = xmlXPathEvalExpression((const xmlChar*)"//pseudosolve",m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    WARNMSG("Using default values for pseudo solver.")
  }
  else
  {
    opt_ptr = result->nodesetval->nodeTab[0];
  }
  //initialize the grid
  bool success = initGrid();
  if(!success)
  {
    ERRORMSG("Failed to create a grid")
    return false;
  }
  Psi.m_grid = myGrid;
  //initialize the wavefunction
  success = initOrbitalSet();
  if(!success)
  {
    ERRORMSG("Failed to create/initialiaze orbitals")
    return false;
  }
  //initialize the wavefunction
  success = initOptimizer();
  if(!success)
  {
    ERRORMSG("Failed to initialize optimizer.")
    return false;
  }
  //initialize the internal storage for a potential
  Pot.initialize(Psi);
  //initialize the Hamiltonian
  success = initHamiltonian();
  if(!success)
  {
    ERRORMSG("Failed to create/initialiaze potential")
    return false;
  }
  //initialize the eigen solver
  result = xmlXPathEvalExpression((const xmlChar*)"//eigensolve",m_context);
  if(xmlXPathNodeSetIsEmpty(result->nodesetval))
  {
    WARNMSG("Using default values for eigen solver")
  }
  else
  {
    cur = result->nodesetval->nodeTab[0]->children;
    while(cur != NULL)
    {
      std::string cname((const char*)(cur->name));
      if(cname == "parameter")
      {
        xmlAttrPtr att = cur->properties;
        while(att != NULL)
        {
          std::string aname((const char*)(att->name));
          std::string vname((const char*)(att->children->content));
          if(aname == "name")
          {
            if(vname == "max_iter")
            {
              putContent(maxiter,cur);
            }
            else
              if(vname == "eig_tol")
              {
                putContent(eig_tol,cur);
              }
              else
                if(vname == "etot_tol")
                {
                  putContent(scf_tol,cur);
                }
                else
                  if(vname == "mix_ratio")
                  {
                    putContent(ratio,cur);
                  }
          }
          att = att->next;
        }
      }
      cur = cur->next;
    }
  }
  xmlXPathFreeObject(result);
  XMLReport("maximum iterations = " << maxiter)
  XMLReport("eigentolerance = " << eig_tol)
  XMLReport("scftolerance = " << scf_tol)
  XMLReport("ratio = " << ratio)
  xmlXPathFreeContext(m_context);
  //set the root file name
  setRoot(AtomName);
  return true;
}

/**
   @brief Initialize the radial grid.
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

bool PseudoGen::initGrid()
{
  xmlNodePtr cur = grid_ptr;
  double scale = 1.0;
  double min = 0.001;
  double max = 1000.0;
  int npts = 2001;
  xmlAttrPtr att = cur->properties;
  while(att != NULL)
  {
    std::string aname((const char*)(att->name));
    const char *avalue = (const char*)(att->children->content);
    if(aname == "type")
    {
      GridType = avalue;
    }
    else
      if(aname == "scale")
      {
        scale = atof(avalue);
      }
    att= att->next;
  }
  XMLReport("Grid type = " << GridType)
  xmlNodePtr cur1 = cur->xmlChildrenNode;
  while(cur1 != NULL)
  {
    if(!xmlStrcmp(cur1->name, (const xmlChar*)"parameter"))
    {
      att = cur1->properties;
      while(att != NULL)
      {
        std::string aname((const char*)(att->name));
        if(aname == "name")
        {
          std::string vname((const char*)(att->children->content));
          if(vname == "min")
          {
            putContent(min,cur1);
          }
          else
            if(vname == "max")
            {
              putContent(max,cur1);
            }
            else
              if(vname == "npts")
              {
                putContent(npts,cur1);
              }
        }
        att=att->next;
      }
    }
    cur1 = cur1->next;
  }
  if(GridType == "log")
  {
    myGrid = new LogGrid<double>;
    min/=scale;
    max/=sqrt(scale);
    myGrid->set(min,max,npts);
  }
  else
    if(GridType == "linear")
    {
      myGrid = new LinearGrid<double>;
      min/=scale;
      max/=sqrt(scale);
      myGrid->set(min,max,npts);
    }
    else
    {
      ERRORMSG("Grid Type Options: Log or Linear.")
      return false;
    }
  XMLReport("Radial Grid: type " << GridType << " rmin = "
            << min << ", rmax = " << max << ", npts = " << npts)
  return myGrid != NULL;
}

/**
   @brief Initialize the wavefunction.
   *
   *Available parameters
   <ul>
   <li> rest_type: the restriction type (spin, spin_space, none)
   </ul>
*/

bool PseudoGen::initOrbitalSet()
{
  std::string afilename;
  //find the restriction type and file name for the All-Electron Orbitals
  xmlNodePtr cur = orb_ptr;
  xmlAttrPtr att = cur->properties;
  while(att != NULL)
  {
    std::string aname((const char*)(att->name));
    const char *avalue = (const char*)(att->children->content);
    if(aname == "rest_type")
    {
      Psi.Restriction = avalue;
    }
    else
      if(aname == "file")
      {
        afilename = avalue;
      }
    att= att->next;
  }
  XMLReport("Orbital restriction type = " << Psi.Restriction)
  XMLReport("Opening file  " << afilename)
  afile = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  group_id = H5Gopen(afile,"radial_basis_states");
  int n;
  int l;
  int m;
  int s;
  std::string grpname;
  value_type occ;
  value_type rm;
  int index = 0;
  int norb = 0;
  xmlNodePtr tcur = cur->xmlChildrenNode;
  while(tcur != NULL)
  {
    if (!(xmlStrcmp(tcur->name, (const xmlChar *) "orbital")))
      norb++;
    tcur = tcur->next;
  }
  XMLReport("The number of orbitals = " << norb)
  AEeigVal.resize(norb);
  PPeigVal.resize(norb);
  for(int ob=0; ob<norb; ob++)
  {
    AEorbitals.push_back(RadialOrbital_t(Psi.m_grid));
    AEorbitals_norm.push_back(RadialOrbital_t(Psi.m_grid));
  }
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    if (!(xmlStrcmp(cur->name, (const xmlChar *) "orbital")))
    {
      LOGMSG("Found Orbital");
      n = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "n"));
      if(n < 0)
      {
        ERRORMSG("Invalid value for n");
        return false;
      }
      l = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "l"));
      if(l < 0)
      {
        ERRORMSG( "Invalid value for l");
        return false;
      }
      m = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "m"));
      if(m < -l && m > 1)
      {
        ERRORMSG( "Invalid value for m");
        return false;
      }
      s = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "s"));
      if(s != 1 && s != -1)
      {
        ERRORMSG( "Invalid value for spin");
        return false;
      }
      occ = atof((const char*)xmlGetProp(cur, (const xmlChar *) "c"));
      if(static_cast<int>(occ) !=0 && static_cast<int>(occ) != 1)
      {
        ERRORMSG( "Invalid value for occupation");
        return false;
      }
      //check to see if orbital with same quantum numbers has
      //already been added
      //  NLMSIndex nlms(n,l,m,s);
      //       NLMS_Map_t::iterator it = OccNo.find(nlms);
      //       if(it == OccNo.end()) {
      // 	add(n,l,m,s,occ);
      // 	LOGMSG("Adding Orbital: n=" << n << ", l=" << l <<
      // 	  ", m=" << m << ", s=" << s << ", occ=" << occ);
      // 	//OccNo[nlmn] = Num-1;
      // 	OccNo[nlms] = 1;
      //       } else {
      // 	ERRORMSG( "Error, orbital " << n << l << m << s
      // 	     << " already occupied");
      // 	return false;
      //       }
      grpname = (const char*)xmlGetProp(cur, (const xmlChar *) "id");
      XMLReport("Adding Pseudo Orbital (" << n << "," << l << ","
                << m << "," << s << ")")
      Psi.add(n,l,m,s,occ);
      if(s == 1)
        Psi.Nup++;
      else
        Psi.Ndown++;
      initAEOrbital(grpname,index);
      index++;
    }
    cur = cur->next;
  }
  H5Gclose(group_id);
  H5Fclose(afile);
  XMLReport("Closing external file")
  LOGMSG("Total number of orbitals = " << Psi.size());
  LOGMSG("(Orbital index, Number of Orbitals)")
  for(int j=0; j < Psi.size(); j++)
  {
    int id = Psi.ID[j];
    LOGMSG("(" << id << ", " << Psi.IDcount[id] << ")");
  }
  //return false if there is no wave functions
  return Psi.size() != 0;
}

bool PseudoGen::initAEOrbital(const std::string& grpname,
                              const int& index)
{
  //open the group containing the proper orbital
  //char grpname[128];
  //sprintf(grpname,"orbital%04d",grpid);
  hid_t group_id_orb = H5Gopen(group_id,grpname.c_str());
  XMLReport("Opening group " << grpname)
  Vector<double> EigVal;
  HDFAttribIO<Vector<double> > EigValHDFIn(EigVal);
  EigValHDFIn.read(group_id_orb,"eigenvalue");
  AEeigVal[index] = EigVal[0];
  XMLReport("Eigenvalue = " << AEeigVal[index])
  //check that the quantum numbers are consistent
  Vector<int> quantum_numbers;
  HDFAttribIO<Vector<int> > QuantumNoHDFIn(quantum_numbers);
  QuantumNoHDFIn.read(group_id_orb,"quantum_numbers");
  XMLReport("Quantum Numbers " << quantum_numbers[0] << " " <<
            quantum_numbers[1] << " " <<  quantum_numbers[2])
  //read the radial orbital data
  Vector<double> rad_orb;
  HDFAttribIO<Vector<double> > RadOrbIn(rad_orb);
  RadOrbIn.read(group_id_orb,"radial_orbital");
  H5Gclose(group_id_orb);
  XMLReport("Reading in radial orbital")
  RadialOrbital_t temp(Psi(0));
  //    std::cout << "index of ae orbital = " << index << std::endl;
  for(int j=0; j<rad_orb.size(); j++)
  {
    AEorbitals[index](j) = rad_orb[j];
    temp(j) = rad_orb[j]*rad_orb[j];
  }
  std::cout << "Checking the normalization: norm = " << setprecision(12)
       <<integrate_RK2_forward(temp,AEorbitals_norm[index]) << std::endl;
  return true;
}

/**
   @brief Initialize the Hamiltonian.
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

bool PseudoGen::initHamiltonian()
{
  //class to generate Clebsh Gordan coeff.
  Clebsch_Gordan* CG_coeff = NULL;
  XMLReport("Creating a Starkloff-Joannopoulos Pseudo-Potential.")
  double rc = 1;
  double lambda = 10.0;
  double Zeff = 4.0;
  xmlNodePtr cur = pot_ptr->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "parameter")
    {
      xmlAttrPtr att = cur->properties;
      while(att != NULL)
      {
        std::string vname((const char*)(att->children->content));
        if(vname == "rc")
        {
          putContent(rc,cur);
        }
        else
          if(vname == "lambda")
          {
            putContent(lambda,cur);
          }
          else
            if(vname == "zeff")
            {
              putContent(Zeff,cur);
            }
        att= att->next;
      }
    }
    cur = cur->next;
  }
  // determine maximum angular momentum
  int lmax = 0;
  int norb = Psi.size();
  for(int i=0; i<norb; i++)
    if(Psi.L[i] > lmax)
      lmax = Psi.L[i];
  lmax++; // increment by 1
  XMLReport("Maximum Angular Momentum = " << lmax)
  CG_coeff = new Clebsch_Gordan(lmax);
  Pot.add(new SJPseudoPotential(vreg,Zeff,rc,lambda));
  Pot.add(new HartreePotential(CG_coeff,norb));
  Pot.add(new ExchangePotential(CG_coeff,norb));
  Psi.CuspParam = 0.0;
  xmlNodePtr cur_sub = pot_ptr->xmlChildrenNode;
  while(cur_sub != NULL)
  {
    std::string pname((const char*)(cur_sub->name));
    if(pname == "parameter")
    {
      xmlAttrPtr att=cur_sub->properties;
      while(att != NULL)
      {
        std::string vname((const char*)(att->children->content));
        if(vname == "mass")
        {
          double m=1.0;
          putContent(m,cur_sub);
          Pot.setMass(m);
          XMLReport("The effective mass is set to " << Pot.getMass())
        }
        att = att->next;
      }
    }
    cur_sub = cur_sub->next;
  }
  return true;
}

bool PseudoGen::initOptimizer()
{
  //initialize the eigen solver
  xmlNodePtr cur = opt_ptr->children;
  while(cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == "parameter")
    {
      xmlAttrPtr att = cur->properties;
      while(att != NULL)
      {
        std::string aname((const char*)(att->name));
        std::string vname((const char*)(att->children->content));
        if(aname == "name")
        {
          if(vname == "cg_tolerance")
          {
            putContent(cg_tolerance,cur);
          }
          else
            if(vname == "cg_epsilon")
            {
              putContent(cg_epsilon,cur);
            }
            else
              if(vname == "weight_norm")
              {
                putContent(weight_norm,cur);
              }
              else
                if(vname == "weight_eig")
                {
                  putContent(weight_eig,cur);
                }
                else
                  if(vname == "cg_stepsize")
                  {
                    putContent(cg_stepsize,cur);
                  }
                  else
                    if(vname == "r_match")
                    {
                      putContent(rmatch,cur);
                    }
        }
        att = att->next;
      }
    }
    cur = cur->next;
  }
  XMLReport("Conjugate Gradient Tolerance = " << cg_tolerance)
  XMLReport("Conjugate Gradient Epsilon = " << cg_epsilon)
  XMLReport("Conjugate Gradient Step Size = " << cg_stepsize)
  XMLReport("Weight for eigenvalues = " << weight_eig)
  XMLReport("Weight for normalization = " << weight_norm)
  XMLReport("Matching radius = " << rmatch)
  return true;
}


/**
   @return the maximum radius of the orbital set
   @brief Print the orbitals and their corresponding eigenvalues
   to a file "AtomName.orb.dat".  The orbitals can easily be
   plotted by using gnuplot or xmgrace.
*/

int PseudoGen::report()
{
  std::string fileforplot(RootFileName);
  fileforplot.append(".orb.dat");
  std::ofstream fout(fileforplot.c_str());
  fout << "#Results for " << AtomName << " with " << PotType
       << " potential on " << GridType << " grid." << std::endl;
  fout << "#Eigen values " << std::endl;
  fout.precision(10);
  fout.setf(std::ios::scientific,std::ios::floatfield);
  for(int orb=0; orb<PPeigVal.size(); orb++)
  {
    fout << "# n=" << Psi.N[orb] << " " <<  " l=" << Psi.L[orb] << std::setw(15)
         << PPeigVal[orb] << " " << std::endl;
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
  Psi.print_HDF5(RootFileName,GridType,PPeigVal);
  Psi.print_basis(AtomName,RootFileName,GridType);
  return max_rad_all;
}
}




