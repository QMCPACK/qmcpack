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
/*@author Jeongnim Kim and Jordan Vincent
 */
#ifndef OHMMS_YLMRNLSET_ONGRID_IO_H
#define OHMMS_YLMRNLSET_ONGRID_IO_H

#ifdef HAVE_LIBHDF5
#include "Numerics/HDFNumericAttrib.h"
#else
#include <string>
#include "OhmmsPETE/OhmmsMatrix.h"
#endif
#include "Numerics/Transform2GridFunctor.h"

/**@brief Print  \f$ r R_{nl}(r)\f$ */
template<class GT>
bool YlmRnlSet<GT>::get(std::ostream& os) {
  os.precision(12);
  os.setf(ios::scientific, ios::floatfield);
  for(int i=0; i < m_grid->size(); i++) {
    os << setw(20) << m_grid->r(i);
    for(int ob=0; ob < psi.size(); ob++) {
      os << setw(20) << psi[ob](i);
    }
    os << endl;
  }
  return true;
}

/**
 *@param cur the current xmlNode which contains definitions 
 *for the orbital set
 *@return true if succeeds
 *@brief Parses the xml file to add any new orbitals to the
 *set of orbitals.
 *
 *Each orbital must have the quantum numbers \f$ (n,l,m,s,c) \f$,
 *where \f$c\f$ is the occupation number.
 */

template<class GT>
bool YlmRnlSet<GT>::put(xmlNodePtr cur){

  int n;
  int l;
  int m;
  int s;
  value_type occ;

  cur = cur->xmlChildrenNode;
  while(cur != NULL) {
    if (!(xmlStrcmp(cur->name, (const xmlChar *) "orbital"))) {
      XMLReport("Found Orbital");
      n = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "n"));
      if(n < 0) {
	ERRORMSG("Invalid value for n");
	return false;
      }
      l = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "l"));
      if(l < 0) {
	ERRORMSG( "Invalid value for l");
	return false;
      }
      
      m = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "m"));
      if(m < -l && m > 1) {
	ERRORMSG( "Invalid value for m");
	return false;
      }
      
      s = atoi((const char*)xmlGetProp(cur, (const xmlChar *) "s"));
      if(s != 1 && s != -1) {
	ERRORMSG( "Invalid value for spin");
	return false;
      }
      occ = atof((const char*)xmlGetProp(cur, (const xmlChar *) "c"));
      if(static_cast<int>(occ) !=0 && static_cast<int>(occ) != 1) {
	ERRORMSG( "Invalid value for occupation");
	return false;
      }

      /*check to see if orbital with same quantum numbers has
	already been added*/
      NLMSIndex nlms(n,l,m,s);
      NLMS_Map_t::iterator it = OccNo.find(nlms);
      if(it == OccNo.end()) {
	add(n,l,m,s,occ);
	XMLReport("Adding Orbital: n=" << n << ", l=" << l <<
		  ", m=" << m << ", s=" << s << ", occ=" << occ);
	OccNo[nlms] = 1;
	//count number of spin-up and down orbitals
	if(s == 1) Nup++;
	else Ndown++;

      } else {
	ERRORMSG( "Error, orbital " << n << l << m << s 
		  << " already occupied");
	return false;
      }
    }  
    cur = cur->next;  
  }

  return true;
}  

/** 
 *@param RootName name of the element
 *@param GridType the type of grid
 *@param eigVal the eigenvalues
 *@return true if succeeds
 *@brief Prints the grid and orbital information to an HDF5
 *file named "RootName.h5".  
 */
template<class GT>
bool YlmRnlSet<GT>::print_HDF5(const std::string& RootName,
			       const std::string& GridType,
			       const std::vector<value_type>& eigVal){
#ifdef HAVE_LIBHDF5
  //print to file fname
  string HDFfname = RootName + ".h5";
  hid_t afile = H5Fcreate(HDFfname.c_str(),H5F_ACC_TRUNC,
			  H5P_DEFAULT,H5P_DEFAULT);
 
  int orbindex = 0;
  //main group  
  hid_t group_id = H5Gcreate(afile,"radial_basis_states",0);
  //output grid information
  hid_t group_id_grid = H5Gcreate(group_id,"grid",0);
  Vector<int> grid_type(1);
  Vector<double> grid_params(3);
  grid_params[0] = m_grid->rmin();
  grid_params[1] = m_grid->rmax();
  grid_params[2] = static_cast<double>(m_grid->size());

  if(GridType == "linear")
    grid_type[0] = 1;
  else if(GridType == "log")
    grid_type[0] = 2;

  HDFAttribIO<Vector<int> > GridHDFTypeOut(grid_type);
  GridHDFTypeOut.write(group_id_grid,"grid_type");
  HDFAttribIO<Vector<double> > GridHDFParamOut(grid_params);
  GridHDFParamOut.write(group_id_grid,"params");
  H5Gclose(group_id_grid);
  //output orbital information
  //output only unique radial functions
  for(int orb=0; orb<NumUniqueOrb; orb++){
    char grpname[128];
    sprintf(grpname,"orbital%04d",orb);
    hid_t group_id_orb = H5Gcreate(group_id,grpname,0);
   
    //output the eigenvalue
    Vector<double> EigVal(1);
    EigVal[0] = eigVal[orbindex];
    HDFAttribIO<Vector<double> > EigValHDFOut(EigVal);
    EigValHDFOut.write(group_id_orb,"eigenvalue");

    //for qmc code, the radial function S(r) = R(r)/r^l = u(r)/r^{l+1}
    Vector<int> Power(1);
    Power[0] = 1+L[orbindex];
    HDFAttribIO<Vector<int> > PowerHDFOut(Power);
    PowerHDFOut.write(group_id_orb,"power");

    //ouptput quantum numbers
    Vector<int> QuantumNo(3);
    QuantumNo[0] = N[orbindex];
    QuantumNo[1] = L[orbindex];
    //only single-zeta basis
    QuantumNo[2] = 1;
    HDFAttribIO<Vector<int> > QuantumNoHDFOut(QuantumNo);
    QuantumNoHDFOut.write(group_id_orb,"quantum_numbers");

    //vector to store the radial orbital
    Vector<double> rad_orb;  
    //determine the maximum radius of the orbital
    int max_rad = m_grid->size()-1;
    while(fabs(psi[orbindex](max_rad)) < 1e-12) max_rad--;   
    max_rad += 2;
    rad_orb.resize(max_rad);
 
    for(int i=0; i<rad_orb.size(); i++)
      rad_orb(i) = psi[orbindex](i);
    //output radial orbital
    HDFAttribIO<Vector<double> > RadOrbOut(rad_orb);
    RadOrbOut.write(group_id_orb,"radial_orbital");
    
    H5Gclose(group_id_orb);   
    orbindex += IDcount[orb];  

  }

  H5Gclose(group_id);
  H5Fclose(afile); 
#endif
  return true;

}

/** 
 *@param elementName name of the element
 *@param RootName name of the element
 *@param GridType the type of grid
 *@return true if succeeds
 *@brief Prints the basis information to an xml file named
 *"RootName.basis.xml" for use in qmcPlusPlus.
 */
template<class GT>
bool YlmRnlSet<GT>::print_basis(const std::string& elementName,
				const std::string& RootName,
				const std::string& GridType){
  //matrices for spin-up and down orbitals
  Matrix<double> Mup;
  Matrix<double> Mdown;
  if(Restriction == "none"){
    Mup.resize(Nup,NumUniqueOrb);
    Mdown.resize(Ndown,NumUniqueOrb);
  } else {
    if(Nup >= Ndown){
      Mup.resize(Nup,Nup);
      Mdown.resize(Ndown,Nup);
    } else {
      Mup.resize(Nup,Ndown);
      Mdown.resize(Ndown,Ndown);
    }
  }
  
  Mup = 0.0;  Mdown = 0.0;

  string fnameXML = RootName + ".basis.xml";
  string fnameHDF5 = RootName + ".h5";

  ofstream osXML(fnameXML.c_str());
  osXML << "<determinantset type=\"MolecularOrbital\">" << endl;
  osXML << "<basisset>" << endl;
  osXML << "<basis type=\"HFNG\" species=\"" << elementName 
	<< "\" file=\"" << fnameHDF5 << "\">" << endl;

  //if there is no restriction all the orbitals are basis functions
  if(Restriction == "none"){
    int nup = 0; int ndown = 0;
    for(int orb=0; orb<size(); orb++){
      if(S[orb] == 1) Mup(nup++,orb) = 1.0;
      else Mdown(ndown++,orb) = 1.0;
      char idname[128];
      sprintf(idname,"R%03d",ID[orb]);
      if(fabs(psi[orb](0)) > 0.0){ 
	osXML << "<phi id=\"" << idname << "\" n=\"" << N[orb] 
	      << "\" l=\"" << L[orb] << "\" m=\"" << M[orb] 
	      << "\" s=\"" << S[orb] << "\" zeta=\"1\"/>" << endl;
      } else {
      osXML << "<phi id=\"" << idname << "\" n=\"" << N[orb] 
	    << "\" l=\"" << L[orb] << "\" m=\"" << M[orb] 
	    << "\" s=\"" << S[orb] << "\" zeta=\"1\" imin=\"1\"/>" << endl;
      }
    }
  } else {
    /*if there is a restriction then the basis consists of one orbital 
      for each unique (n,l,m), since for a given (n,l,m) spin-up is equal 
      to spin-down, assign each basis function an up spin*/
    int nup = 0; int ndown = 0; int nbasis = -1;   
    for(int orb=0; orb<size(); orb++){
      NLMIndex nlm(N[orb],L[orb],M[orb]);
      NLM_Map_t::iterator it = NLM.find(nlm); 
      if(it == NLM.end()) {
	NLM[nlm] = nbasis;
	nbasis++;
	char idname[128];
	sprintf(idname,"R%03d",ID[orb]);
	if(fabs(psi[orb](0)) > 0.0){ 
	  osXML << "<phi id=\"" << idname << "\" n=\"" << N[orb] 
		<< "\" l=\"" << L[orb] << "\" m=\"" << M[orb] 
		<< "\" s=\"1\" zeta=\"1\"/>" << endl;
	} else {
	  osXML << "<phi id=\"" << idname << "\" n=\"" << N[orb] 
		<< "\" l=\"" << L[orb] << "\" m=\"" << M[orb] 
		<< "\" s=\"" << S[orb] << "\" zeta=\"1\" imin=\"1\"/>" << endl;
	}
	nbasis++;
      }
    }

    for(int i=0; i<Mup.rows(); i++) Mup(i,i) = 1.0;
    for(int i=0; i<Mdown.rows(); i++) Mdown(i,i) = 1.0;

  }
 
  osXML << "</basis>" << endl;
  osXML << "</basisset>" << endl;
  osXML << "<slaterdeterminant>" << endl;
  osXML << "<determinant spin=\"1cu\" orbitals=\"" << Mup.rows() 
	<< "\">" << endl;
  osXML << "<parameter id=\"1cu\" type=\"Array\">" << endl;
  osXML << Mup;
  osXML << "</parameter>" << endl;
  osXML << "</determinant>" << endl;
  osXML << "<determinant spin=\"-1\" orbitals=\"" << Mdown.rows() 
	<< "\">" << endl;
  osXML << "<parameter id=\"1cd\" type=\"Array\">" << endl;
  osXML << Mdown;
  osXML << "</parameter>" << endl;
  osXML << "</determinant>" << endl;
  osXML << "</slaterdeterminant>" << endl;
  osXML << "</determinantset>" << endl;

  osXML.close(); 

  return true;

}

#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

