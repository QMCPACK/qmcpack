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
#include "AtomicHF/HFConfiguration.h"
#include "AtomicHF/Clebsch_Gordan.h"
#include "PseudoGen/PseudoGen.h"
#include "AtomicHF/fillshells.h"
#include "Numerics/HDFNumericAttrib.h"
#include "OhmmsData/OhmmsElementBase.h"
#include "Utilities/OhmmsSpecies.h"
#include <fstream>
#include <libxml/xpath.h>

namespace ohmmshf {

  bool parseXMLFile(HFAtomicOrbitals& Psi, string& gridtype,
		    const xmlpp::Node* root) {
    using namespace xmlpp;
    //  string gridtype;
    OneDimGridBase<double>* grid;
    NodeSet gset = root->find("./Grid");
    if(gset.empty()){
      ERRORMSG("No grid information")
	return false;
    } else {
      Element* p = dynamic_cast<Element*>(gset[0]);
      Element::AttributeList atts = p->get_attributes();
      Element::AttributeList::iterator it = atts.begin();
	
      double scale = 32.0;
      double min = 0.00001;
      double max = 1000.0;
      int npts = 2001;
      while(it != atts.end()) {  
	const string& aname = (*it)->get_name();
	if(aname == "type") {
	  gridtype = (*it)->get_value();
	} else if(aname == "scale") {
	  scale = atof((*it)->get_value().c_str());
	}
	it++;
      }
	
	XMLReport("Grid type = " << gridtype)
	  
	  NodeSet gpset = gset[0]->find("./Parameter");
	if(gpset.empty()){
	  WARNMSG("Using default grid values")
	    } else {
	      for(int i=0; i<gpset.size(); i++){
		Element* cur = dynamic_cast<Element*>(gpset[i]);
		Attribute* att = cur->get_attribute("name");
		if(att) {
		  const string& aname = att->get_value();
		  xmlNode* curc = cur->cobj();
		  if(aname == "min") putContent(min,curc);
		  else if(aname == "max") putContent(max,curc);
		  else if(aname == "npts") putContent(npts,curc);
		}
	      }
	    }

	
	if(gridtype == "log"){
	  grid = new LogGrid<double>;
	  min/=scale;
	  max/=sqrt(scale);
	  grid->set(min,max,npts);
	  XMLReport("rmin = " << min << ", rmax = " << max 
		    << ", dh = " << grid->dh() << ", npts = " << npts) 
	    } else if(gridtype == "linear"){
	      grid = new LinearGrid<double>;
	      min/=scale;
	      max/=sqrt(scale);
	      grid->set(min,max,npts);
	      XMLReport("rmin = " << min << ", rmax = " << max 
			<< ", dh = " << grid->dh() << ", npts = " << npts) 
		} else {
		  ERRORMSG("Grid Type Options: Log or Linear.")
		    return false;
		}
	
    }
      // pass grid to wavefunction
    Psi.m_grid = grid;
    NodeSet oset = root->find("./OrbitalSet");
    if(oset.empty()){
      ERRORMSG("Must specify OrbitalSet")
	return false;
    } else {
      Element* cur = dynamic_cast<Element*>(oset[0]);
      Attribute* att = cur->get_attribute("rest_type");
      if(att) Psi.Restriction = att->get_value();
      XMLReport("Orbital restriction type = " << Psi.Restriction)
	}
  }
}

int main(int argc, char **argv) {
  using namespace ohmmshf;
  using namespace xmlpp;
   DomParser parser;
  parser.parse_file(argv[1]);
  Node* root
   = parser.get_document()->get_root_node(); //deleted by DomParser.
   string gridtype;
//   double scale = 32.0;
//   double min = 0.00001;
//   double max = 1000.0;
//   int npts = 2001;
//  OneDimGridBase<double>* grid;
//  grid = new LogGrid<double>;
//  min/=scale;
//  max/=sqrt(scale);
//  grid->set(min,max,npts);
  RadialPotentialSet Pot;
  HFAtomicOrbitals Psi;
  parseXMLFile(Psi,gridtype,root);
  //  Psi.m_grid = grid;
  //  Psi.Restriction = "spin+space";
  int lmax = 1;
  XMLReport("Maximum Angular Momentum = " << lmax)
    Clebsch_Gordan* CG_coeff;
  CG_coeff = new Clebsch_Gordan(lmax);
  
  Pot.add(new HartreePotential(CG_coeff));
  Pot.add(new ExchangePotential(CG_coeff));

  XMLReport("Adding Pseudo 1s down orbital.")
  Psi.add(1,0,0,-1,1.0);
  XMLReport("Adding Pseudo 1s up orbital.")
  Psi.add(1,0,0,1,1.0);
  XMLReport("Adding Pseudo 2px up orbital.")
  Psi.add(2,1,1,1,1.0);
  XMLReport("Adding Pseudo 2pz up orbital.")
  Psi.add(2,1,0,1,1.0);
  Pot.initialize(Psi);
  //  Psi.CuspParam = 0.0;
  // Psi.EigenBoundParam = -8.0;
  PseudoHF HFSolver(Pot,Psi,root);
  HFSolver.run();
  string element = "Ge";
  string fname = element+".pp.h5";
  hid_t afile = H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,
			  H5P_DEFAULT,H5P_DEFAULT);

  int orbindex = 0;
  
  hid_t group_id = H5Gcreate(afile,element.c_str(),0);
 
  hid_t group_id_pot = H5Gcreate(group_id,"Potential",0);
  Vector<double> cusp_temp(1);
  cusp_temp(0) = Psi.CuspParam;
  HDFAttribIO<Vector<double> > CuspHDFOut(cusp_temp);
  CuspHDFOut.write(group_id_pot,"CuspParam");
  Vector<double> pot_temp(1);
  pot_temp(0) = Psi.EigenBoundParam;
  HDFAttribIO<Vector<double> > PotHDFOut(pot_temp);
  PotHDFOut.write(group_id_pot,"PotParam");
  H5Gclose(group_id_pot);     

  gridtype+="grid";
  hid_t group_id_grid = H5Gcreate(group_id,gridtype.c_str(),0);
  Vector<double> grid_temp(Psi.m_grid->size());
  Vector<double> psi_temp(Psi.m_grid->size());
  for(int i=0; i<grid_temp.size(); i++)
    grid_temp[i] = Psi.m_grid->r(i);
  HDFAttribIO<Vector<double> > GridHDFOut(grid_temp);
  GridHDFOut.write(group_id_grid,"griddata");
  H5Gclose(group_id_grid);
  ofstream* osXML;
  string fnameXML = element+=".basis.xml";
  osXML = new ofstream(fnameXML.c_str());
  *osXML << "<Basis type=\"HFNG\" species=\"" << element << "\">" << endl;
  for(int orb=0; orb<Psi.NumUniqueOrb; orb++){
    char grpname[128];
    sprintf(grpname,"orbital%04d",orb);
    hid_t group_id_orb = H5Gcreate(group_id,grpname,0);
    Vector<int> NLMS;
    NLMS.resize(4);
    NLMS[0] = Psi.N[orbindex];
    NLMS[1] = Psi.L[orbindex];
    NLMS[2] = Psi.M[orbindex];
    NLMS[3] = Psi.S[orbindex];
    *osXML << "<Rnl id =\"" << grpname << "\" n=\"" << NLMS[0] 
	   << "\" l=\"" << NLMS[1] << "\" m=\"" << NLMS[2] 
	   << "\" s=\"" << NLMS[3] << "\" />" << endl;
    

    HDFAttribIO<Vector<int> > QuantumNoHDFOut(NLMS);
    QuantumNoHDFOut.write(group_id_orb,"nlms");
    
    for(int i=0; i<psi_temp.size(); i++)
      psi_temp[i] = Psi(orb,i);
    
    HDFAttribIO<Vector<double> > PsiHDFOut(psi_temp);
    PsiHDFOut.write(group_id_orb,"u_r");
    
    H5Gclose(group_id_orb);   
    orbindex += Psi.IDcount[orb];  
  }
  
  *osXML << "</Basis>" << endl;
  delete osXML; 
  H5Gclose(group_id);
  H5Fclose(afile);  
  
  /*
    ofstream* out;
    out = new ofstream("out.h5");
    Psi.get(*out);
    out->close();
  */
  return 1;

}


/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

  
  
