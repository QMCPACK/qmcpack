//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "SQD/HFConfiguration.h"
#include "Numerics/Clebsch_Gordan.h"
#include "SQD/fillshells.h"
#include "SQD/SphericalPotential/ZOverRPotential.h"
#include "SQD/SphericalPotential/HarmonicPotential.h"
#include "SQD/SphericalPotential/StepPotential.h"
#include "SQD/SphericalPotential/SJPseudoPotential.h"
#include "SQD/HartreeFock.h"

namespace ohmmshf {

  /**
     @param pot the potential
     @param psi the wavefunction
     @brief The contructor.
  */

  HartreeFock::HartreeFock(RadialPotentialSet& pot, 
			   SphericalOrbitalTraits::BasisSetType& psi):
    Pot(pot), Psi(psi), 
    num_closed_shells(0),
    maxiter(1000), eig_tol(1e-12),
    scf_tol(1e-8), ratio(0.35), 
    GridType("none"),
    grid_ptr(NULL), orb_ptr(NULL), 
    pot_ptr(NULL), myGrid(NULL) { }

  /**
     @param aroot the root for all output files
     @brief Sets the root for all output files.
  */

  void HartreeFock::setRoot(const string& aroot) {
    RootFileName = aroot;
    LogFileName = RootFileName + ".log";
  }

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
     default = \f$1 e^{-8}\f$
     <li> mix_ratio, the mixing ratio of the charge density, default = 0.35
     </ul>
  */
  bool HartreeFock::put(xmlNodePtr d_root){

    xmlNodePtr cur = NULL;

    //using XPath instead of recursive search
    xmlXPathContextPtr m_context = xmlXPathNewContext(d_root->doc);


    xmlXPathObjectPtr result = xmlXPathEvalExpression((const xmlChar*)"//atom",m_context);
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      ERRORMSG("Missing Atom information. Exit")
	return false;
    } else {
      cur = result->nodesetval->nodeTab[0];
      xmlAttrPtr att = cur->properties;
      while(att != NULL) {
	string aname((const char*)(att->name));
	const char *avalue = (const char*)(att->children->content);
	if(aname == "name") {
	  AtomName = avalue;
	} else if (aname == "num_closed_shells") {
	  num_closed_shells = atoi(avalue);
	}
	att = att->next;
      }

      XMLReport("Atom name = " << AtomName)
	XMLReport("Number of closed shells = " << num_closed_shells)
	//initialize xmlNode pointers for the grid, wavefunction and potential
	xmlNodePtr cur1 = cur->xmlChildrenNode;
      while(cur1 != NULL) {
	string cname1((const char*)(cur1->name));
	if(cname1 == "grid") {
	  grid_ptr = cur1;
	} else if(cname1 == "orbitalset") {
	  orb_ptr = cur1;
	} else if(cname1 == "potential") {
	  pot_ptr = cur1;
	  if(xmlHasProp(cur1,(const xmlChar*)"type")) {
	    PotType = (const char*)(xmlGetProp(cur1, (const xmlChar *) "type"));
	    XMLReport("The type of potential " << PotType)
	      } else {
		ERRORMSG("Potential type is undefined. Exit")
		  return false;
	      }
	}
	cur1 = cur1->next;
      }
    }
    if(orb_ptr == NULL || pot_ptr == NULL || grid_ptr == NULL) {
      ERRORMSG("Missing one of nodes: <grid/>, <orbitalset/> or <potential/>. Exit")
	return false;
    }
    xmlXPathFreeObject(result);

    //initialize the grid
    bool success = initGrid();
    if(!success) {
      ERRORMSG("Failed to create a grid")
	return false;
    }

    Psi.m_grid = myGrid;
    //initialize the wavefunction
    success = initOrbitalSet();
    if(!success) {
      ERRORMSG("Failed to create/initialiaze orbitals")
	return false;
    }

    //initialize the internal storage for a potential
    Pot.initialize(Psi);

    //initialize the Hamiltonian
    success = initHamiltonian();
    if(!success) {
      ERRORMSG("Failed to create/initialiaze potential")
	return false;
    }
   
    //initialize the eigen solver
    result = xmlXPathEvalExpression((const xmlChar*)"//eigensolve",m_context);
    if(xmlXPathNodeSetIsEmpty(result->nodesetval)) {
      WARNMSG("Using default values for eigen solver")
	} else {
	  cur = result->nodesetval->nodeTab[0]->children;
	  while(cur != NULL) {
	    string cname((const char*)(cur->name));
	    if(cname == "parameter") {
	      xmlAttrPtr att = cur->properties;
	      while(att != NULL) {
		string aname((const char*)(att->name));
		string vname((const char*)(att->children->content));
		if(aname == "name") {
		  if(vname == "max_iter"){
		    putContent(maxiter,cur);
		  } else if(vname == "eig_tol"){
		    putContent(eig_tol,cur);
		  } else if(vname == "etot_tol"){
		    putContent(scf_tol,cur);
		  } else if(vname == "mix_ratio"){
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

  bool HartreeFock::initGrid() {

    xmlNodePtr cur = grid_ptr;
    double scale = 1.0;
    double min = 0.001;
    double max = 1000.0;
    int npts = 2001;
    xmlAttrPtr att = cur->properties;
    while(att != NULL) {
      string aname((const char*)(att->name));
      const char *avalue = (const char*)(att->children->content);
      if(aname == "type") {
	GridType = avalue;
      } else if(aname == "scale") {
	scale = atof(avalue);
      }
      att= att->next;
    }

    XMLReport("Grid type = " << GridType)
      xmlNodePtr cur1 = cur->xmlChildrenNode;	  
    while(cur1 != NULL) {
      if(!xmlStrcmp(cur1->name, (const xmlChar*)"parameter")) {
	att = cur1->properties;
	while(att != NULL) {
	  string aname((const char*)(att->name));
	  if(aname == "name") {
	    string vname((const char*)(att->children->content));
	    if(vname == "min") {
	      putContent(min,cur1);
	    } else if(vname == "max") { 
	      putContent(max,cur1);
	    } else if(vname == "npts") {
	      putContent(npts,cur1);
	    }
	  }
	  att=att->next;
	}
      }
      cur1 = cur1->next;
    }
    if(GridType == "log"){
      myGrid = new LogGrid<double>;
      min/=scale;
      max/=sqrt(scale);
      myGrid->set(min,max,npts);
    } else if(GridType == "linear"){
      myGrid = new LinearGrid<double>;
      min/=scale;
      max/=sqrt(scale);
      myGrid->set(min,max,npts);
    } else {
      ERRORMSG("Grid Type Options: Log or Linear.")
	return false;
    }
    XMLReport("Radial Grid: type " << GridType << " rmin = " 
	      << min << ", rmax = " << max << ", npts = " << npts) 
      //     << ", dh = " << myGrid->dh() << ", npts = " << npts) 
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

  bool HartreeFock::initOrbitalSet() {

    //find the restriction type
    xmlNodePtr cur = orb_ptr;
    if(xmlHasProp(cur,(const xmlChar*)"rest_type")) {
      Psi.Restriction = (const char*)(xmlGetProp(cur, (const xmlChar *) "rest_type"));
      XMLReport("Orbital restriction type = " << Psi.Restriction)
	}

    //fill the closed shells (shell filling depends of potential)
    if(num_closed_shells) {
      if(PotType == "harmonic" || PotType == "step"){
	FillShellsHarmPot(Psi,num_closed_shells);
      } else { 
	FillShellsNucPot(Psi,num_closed_shells);
      }
    }
    //pass the xmlNode to Psi
    Psi.put(cur);

    LOGMSG("Total number of orbitals = " << Psi.size()); 

    LOGMSG("(Orbital index, Number of Orbitals)")
      for(int j=0; j < Psi.size(); j++){
	int id = Psi.ID[j];
	LOGMSG("(" << id << ", " << Psi.IDcount[id] << ")");
      }

    //return false if there is no wave functions
    return Psi.size() != 0;
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

  bool HartreeFock::initHamiltonian() {

    //class to generate Clebsh Gordan coeff.
    Clebsch_Gordan* CG_coeff = NULL;

    if(PotType == "nuclear"){
      XMLReport("Creating a Nuclear Potential.")
	double Z = Psi.size();

      xmlNodePtr cur = pot_ptr->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "parameter") {
	  xmlAttrPtr att = cur->properties;
	  while(att != NULL) {
	    string vname((const char*)(att->children->content));
	    if(vname == "z") {
	      putContent(Z,cur);
	    }
	    att= att->next;
	  }
	}
	cur = cur->next;
      }
      XMLReport("Potential Paramter: Z = " << Z)
	// determine maximum angular momentum
	int lmax = 0;
      for(int i=0; i<Psi.L.size(); i++)
	if(Psi.L[i] > lmax) lmax = Psi.L[i]; 
      XMLReport("Maximum Angular Momentum = " << lmax)
	CG_coeff = new Clebsch_Gordan(lmax);

      Pot.add(new ZOverRPotential(Z), true);
      Pot.add(new HartreePotential(CG_coeff));
      Pot.add(new ExchangePotential(CG_coeff));
      Psi.CuspParam = Z;
    } // if Nuclear
    else if(PotType == "nuclear_scalar_rel"){
      XMLReport("Creating a Nuclear Scalar Relativistic Potential.")
	double Z = Psi.size();

      xmlNodePtr cur = pot_ptr->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "parameter") {
	  xmlAttrPtr att = cur->properties;
	  while(att != NULL) {
	    string vname((const char*)(att->children->content));
	    if(vname == "z") {
	      putContent(Z,cur);
	    }
	    att= att->next;
	  }
	}
	cur = cur->next;
      }
      XMLReport("Potential Paramter: Z = " << Z)
	int lmax = 0;
      for(int i=0; i<Psi.L.size(); i++)
	if(Psi.L[i] > lmax) lmax = Psi.L[i]; 
      XMLReport("Maximum Angular Momentum = " << lmax)
	CG_coeff = new Clebsch_Gordan(lmax);
      
      Pot.add(new ZOverRPotential(Z),true);
      Pot.add(new HartreePotential(CG_coeff));
      Pot.add(new ExchangePotential(CG_coeff));
      Psi.CuspParam = Z;
    } // if Nuclear Scalar Relativistic
    else if(PotType == "harmonic"){
      XMLReport("Creating a Harmonic Potential.")
	double Omega=0.5;

      xmlNodePtr cur = pot_ptr->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "parameter") {
	  xmlAttrPtr att = cur->properties;
	  while(att != NULL) {
	    string vname((const char*)(att->children->content));
	    if(vname == "omega" || vname == "z") {
	      putContent(Omega,cur);
	    }
	    att= att->next;
	  }
	}
	cur = cur->next;
      }
      XMLReport("Potential Parameter: Omega = " << Omega)
	//determine maximum angular momentum
	int lmax = 0;
      for(int i=0; i<Psi.L.size(); i++)
	if(Psi.L[i] > lmax) lmax = Psi.L[i]; 
      lmax++; // increment by 1
      XMLReport("Maximum Angular Momentum = " << lmax)
	CG_coeff = new Clebsch_Gordan(lmax);

      Pot.add(new HarmonicPotential(Omega),true);
      Pot.add(new HartreePotential(CG_coeff));
      Pot.add(new ExchangePotential(CG_coeff));
      Psi.CuspParam = 0.0;
    } // if Harmonic
    else if(PotType == "step"){
      XMLReport("Creating a Step-Function Potential.")
	//determine maximum angular momentum
	int lmax = 0;
      for(int i=0; i<Psi.L.size(); i++)
	if(Psi.L[i] > lmax) lmax = Psi.L[i];
      lmax++; // increment by 1
      XMLReport("Maximum Angular Momentum = " << lmax)
	CG_coeff = new Clebsch_Gordan(lmax);

      StepPotential* apot = new StepPotential;
      apot->put(pot_ptr);
      Pot.add(apot,true);
      Pot.add(new HartreePotential(CG_coeff));
      Pot.add(new ExchangePotential(CG_coeff));
      Psi.CuspParam = 0.0;
    } // if Step
    else if(PotType == "pseudo"){
      XMLReport("Creating a Starkloff-Joannopoulos Pseudo-Potential.")
	double rc, lambda, Zeff;

      xmlNodePtr cur = pot_ptr->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "parameter") {
	  xmlAttrPtr att = cur->properties;
	  while(att != NULL) {
	    string vname((const char*)(att->children->content));
	    if(vname == "rc") {
	      putContent(rc,cur);
	    } else if(vname == "lambda") {
	      putContent(lambda,cur);
	    } else if(vname == "zeff") {
	      putContent(Zeff,cur);
	    }
	    att= att->next;
	  }
	}
	cur = cur->next;
      }
      XMLReport("Potential Parameter: Effective Charge = " << Zeff)
	XMLReport("Potential Parameter: Core Radius = " << rc)
	XMLReport("Potential Parameter: Lambda = " << lambda)
	// determine maximum angular momentum
	int lmax = 0;
      for(int i=0; i<Psi.L.size(); i++)
	if(Psi.L[i] > lmax) lmax = Psi.L[i]; 
      lmax++; // increment by 1
      XMLReport("Maximum Angular Momentum = " << lmax)
	CG_coeff = new Clebsch_Gordan(lmax);

      Pot.add(new SJPseudoPotential(Zeff,rc,lambda));
      Pot.add(new HartreePotential(CG_coeff));
      Pot.add(new ExchangePotential(CG_coeff));
      Psi.CuspParam = 0.0;
    } // if Pseudo
    else {
      ERRORMSG("Unknown potential type" << PotType)
	return false;
    }


    xmlNodePtr cur_sub = pot_ptr->xmlChildrenNode;
    while(cur_sub != NULL) {
      string pname((const char*)(cur_sub->name));
      if(pname == "parameter") {
	xmlAttrPtr att=cur_sub->properties;
	while(att != NULL) {
	  string vname((const char*)(att->children->content));
	  if(vname == "mass") {
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

  /**
     @return the maximum radius of the orbital set
     @brief Print the orbitals and their corresponding eigenvalues
     to a file "AtomName.orb.dat".  The orbitals can easily be
     plotted by using gnuplot or xmgrace.
  */

  int HartreeFock::report() {

    string fileforplot(AtomName);
    fileforplot.append(".orb.dat");
    ofstream fout(fileforplot.c_str());

    fout << "#Results for " << AtomName << " with " << PotType << " potential on " 
	 << GridType << " grid." <<endl;
    fout << "#Eigen values " << endl;
    for(int orb=0; orb<eigVal.size(); orb++) {
      fout << "# n=" << Psi.N[orb] << " " <<  " l=" << Psi.L[orb] << setw(15) 
	   << eigVal[orb] << " " << endl;
    }

    fout << "#The number of unique radial orbitals " << Psi.NumUniqueOrb << endl;
    //find the maximum radius for the orbital set
    int max_rad_all = 0;
    int orbindex = 0;
    for(int orb=0; orb<Psi.NumUniqueOrb; orb++){
      int max_rad = myGrid->size()-1;
      while(fabs(Psi(orbindex,max_rad)) < 1e-12) max_rad--;   
      max_rad += 2;
      max_rad_all = std::max(max_rad_all,max_rad);
      orbindex += Psi.IDcount[orb];  
    }

    for(int ig=0; ig<max_rad_all; ig++) {
      fout << setw(15) << myGrid->r(ig);
      orbindex=0;
      for(int orb=0; orb<Psi.NumUniqueOrb; orb++){
	fout << setw(15) << Psi(orbindex,ig);
	orbindex += Psi.IDcount[orb];  
      }
      fout << endl;
    }
    Psi.print_HDF5(RootFileName,GridType,eigVal);
    Psi.print_basis(AtomName,RootFileName,GridType);
    return max_rad_all;
  }
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

  
  
