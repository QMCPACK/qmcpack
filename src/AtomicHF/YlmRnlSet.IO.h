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
    
    


#ifndef OHMMS_YLMRNLSET_ONGRID_IO_H
#define OHMMS_YLMRNLSET_ONGRID_IO_H

#include "Numerics/HDFNumericAttrib.h"
#include "Numerics/Transform2GridFunctor.h"

/**@brief Print  \f$ r R_{nl}(r)\f$ */
template<class GT>
bool YlmRnlSet<GT>::get(std::ostream& os)
{
  os.precision(12);
  os.setf(std::ios::scientific, std::ios::floatfield);
  for(int i=0; i < m_grid->size(); i++)
  {
    os << std::setw(20) << m_grid->r(i);
    for(int ob=0; ob < psi.size(); ob++)
    {
      os << std::setw(20) << psi[ob](i);
    }
    os << std::endl;
  }
  return true;
}

/*! \fn template<class GT>
  bool YlmRnlSet<GT>::put(xmlNodePtr cur)
  * \param cur
  * \return true if succeeds
  * \brief Loops through the XML tree to find all the orbitals,
  checks to see if the orbital has not already been added,
  and adds any new orbitals to the orbital list.
*/

template<class GT>
bool YlmRnlSet<GT>::put(xmlNodePtr cur)
{
  int n;
  int l;
  int m;
  int s;
  value_type occ;
  cur = cur->xmlChildrenNode;
  while(cur != NULL)
  {
    if (!(xmlStrcmp(cur->name, (const xmlChar *) "Orbital")))
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
      NLMSIndex nlms(n,l,m,s);
      NLMS_Map_t::iterator it = OccNo.find(nlms);
      if(it == OccNo.end())
      {
        add(n,l,m,s,occ);
        LOGMSG("Adding Orbital: n=" << n << ", l=" << l <<
               ", m=" << m << ", s=" << s << ", occ=" << occ);
        //	OccNo[nlmn] = Num-1;
        OccNo[nlms] = 1;
      }
      else
      {
        ERRORMSG( "Error, orbital " << n << l << m << s
                  << " already occupied");
        return false;
      }
    }
    cur = cur->next;
  }
  return true;
}

/*! \fn template<class GT>
  bool YlmRnlSet<GT>::print(const std::string& elementName)
  * \param elementName name of the element
  * \return true if succeeds
  * \brief The output files produced by this function are
  for use in a Quantum Monte Carlo code and for graphical
  purposes only.  Prints the grid and orbital information to an HDF5
  file named "elementName.h5".  Also prints the basis information
  to a file named "elementName.basis.xml" and the individual
  orbitals to ASCII files named "ORB.N.L.(M).(S).elementName",
  where N,L,M,S are the quantum numbers (M,S are only printed
  if the orbitals are unrestricted).  In the ASCII files, the
  orbitals are first interpolated onto a linear grid of 500
  points for reasons of disk-space.
*/

template<class GT>
bool YlmRnlSet<GT>::print(const std::string& elementName)
{
/// print to file fname
  std::string HDFfname = elementName + ".h5";
  hid_t afile = H5Fcreate(HDFfname.c_str(),H5F_ACC_TRUNC,
                          H5P_DEFAULT,H5P_DEFAULT);
  OneDimGridBase<double>* reduced_grid;
  reduced_grid = new LinearGrid<double>;
  std::ofstream* out;
  int orbindex = 0;
/// main group
  hid_t group_id = H5Gcreate(afile,"radial_basis_states",0);
/// output grid information
  hid_t group_id_grid = H5Gcreate(group_id,"grid",0);
  Vector<int> grid_type(1);
  Vector<double> grid_params(3);
  grid_params[0] = m_grid->rmin();
  grid_params[1] = m_grid->rmax();
  grid_params[2] = static_cast<double>(m_grid->size());
// if(gridtype == "linear")
//  grid_type[0] = 1;
// else if(gridtype == "log")
  grid_type[0] = 2;
//HDFAttribIO<std::string> GridHDFTypeOut(gridtype);
//GridHDFTypeOut.write(group_id_grid,"grid_type");
  HDFAttribIO<Vector<int> > GridHDFTypeOut(grid_type);
  GridHDFTypeOut.write(group_id_grid,"grid_type");
  HDFAttribIO<Vector<double> > GridHDFParamOut(grid_params);
  GridHDFParamOut.write(group_id_grid,"params");
  H5Gclose(group_id_grid);
/// output orbital information
/// output only unique radial functions
  for(int orb=0; orb<NumUniqueOrb; orb++)
  {
    char grpname[128];
    sprintf(grpname,"orbital%04d",orb);
    hid_t group_id_orb = H5Gcreate(group_id,grpname,0);
    Vector<int> Power(1);
    Power[0] = 1+L[orbindex];
    HDFAttribIO<Vector<int> > PowerHDFOut(Power);
    PowerHDFOut.write(group_id_orb,"power");
    Vector<int> QuantumNo(3);
    QuantumNo[0] = N[orbindex];
    QuantumNo[1] = L[orbindex];
    // only single-zeta basis
    QuantumNo[2] = 1;
    HDFAttribIO<Vector<int> > QuantumNoHDFOut(QuantumNo);
    QuantumNoHDFOut.write(group_id_orb,"quantum_numbers");
    ///vector to store the radial orbital
    Vector<double> rad_orb;
    ///create the ASCII file
    char sfname[128];
    if(Restriction == "none")
      sprintf(sfname,"ORB.%01d.%01d.%01d.%01d.%01s",N[orbindex],L[orbindex],
              M[orbindex],S[orbindex],elementName.c_str());
    else
      sprintf(sfname,"ORB.%01d.%01d.%01s",N[orbindex],L[orbindex],
              elementName.c_str());
    out = new std::ofstream(sfname);
    int max_rad = m_grid->size()-1;
    while(std::abs(psi[orbindex](max_rad)) < 1e-12)
      max_rad--;
    max_rad += 2;
    rad_orb.resize(max_rad);
    reduced_grid->set(0.0,m_grid->r(max_rad),max_rad);
    RadialOrbital_t temp_orb(reduced_grid);
    for(int i=0; i<rad_orb.size(); i++)
      rad_orb(i) = psi[orbindex](i);
    ///calculate the boundary conditions
    int imin = 0;
    value_type deriv = (psi[orbindex](imin+1) - psi[orbindex](imin))/m_grid->dr(imin);
    psi[orbindex].spline(imin,deriv,max_rad,0.0);
    Transform2GridFunctor<OneDimGridFunctor<double>,OneDimGridFunctor<double> >
    transform(psi[orbindex],temp_orb);
    transform.generate();
    for(int j=0; j<reduced_grid->size(); j++)
      *out << std::setw(15) << temp_orb.r(j) << std::setw(15) << temp_orb(j) << std::endl;
    out->close();
    HDFAttribIO<Vector<double> > RadOrbOut(rad_orb);
    RadOrbOut.write(group_id_orb,"radial_orbital");
    H5Gclose(group_id_orb);
    orbindex += IDcount[orb];
  }
  H5Gclose(group_id);
  H5Fclose(afile);
/// create input for qmc
  int nUp = 0;
  int nDn = 0;
  int nBasis = 0;
  std::ofstream* osXML;
  std::string fnameXML = elementName + ".basis.xml";
  osXML = new std::ofstream(fnameXML.c_str());
  *osXML << "<DeterminantSet type=\"MolecularOrbital\">" << std::endl;
  *osXML << "<BasisSet>" << std::endl;
  *osXML << "<Basis type=\"HFNG\" species=\"" << elementName
         << "\" file=\"" << HDFfname << "\">" << std::endl;
  if(Restriction == "none")
  {
    for(int orb=0; orb<size(); orb++)
    {
      (S[orb] == 1) ? nUp++ : nDn++;
      nBasis++;
      char idname[128];
      sprintf(idname,"R%03d",ID[orb]);
      *osXML << "<psi id=\"" << idname << "\" n=\"" << N[orb]
             << "\" l=\"" << L[orb] << "\" m=\"" << M[orb]
             << "\" s=\"" << S[orb] << "\" zeta=\"1\"/>" << std::endl;
    }
  }
  else
  {
    for(int orb=0; orb<size(); orb++)
    {
      (S[orb] == 1) ? nUp++ : nDn++;
      if(S[orb] == 1)
      {
        nBasis++;
        char idname[128];
        sprintf(idname,"R%03d",ID[orb]);
        *osXML << "<psi id=\"" << idname << "\" n=\"" << N[orb]
               << "\" l=\"" << L[orb] << "\" m=\"" << M[orb]
               << "\" s=\"" << S[orb] << "\" zeta=\"1\"/>" << std::endl;
      }
    }
  }
  *osXML << "</Basis>" << std::endl;
  *osXML << "</BasisSet>" << std::endl;
  *osXML << "<SlaterDeterminant>" << std::endl;
  *osXML << "<Determinant spin=\"1\" orbitals=\"" << nUp << "\">" << std::endl;
  *osXML << "<Var type=\"Array\">" << std::endl;
  Matrix<double> M_Up(nUp,nBasis);
  if(Restriction == "none")
  {
    int row = 0;
    orbindex = 0;
    for(int j=0; j<M_Up.cols(); j++)
    {
      if(S[orbindex] == 1)
      {
        M_Up(row,j) = 1.0;
        row++;
      }
      orbindex += IDcount[j];
    }
  }
  else
  {
    int row = 0;
    int col = 0;
    while(row < M_Up.rows())
    {
      M_Up(row,col) = 1.0;
      row++;
      col++;
    }
  }
  *osXML << M_Up;
  *osXML << "</Var>" << std::endl;
  *osXML << "</Determinant>" << std::endl;
  *osXML << "<Determinant spin=\"-1\" orbitals=\"" << nDn << "\">" << std::endl;
  *osXML << "<Var type=\"Array\">" << std::endl;
  Matrix<double> M_Dn(nDn,nBasis);
  if(Restriction == "none")
  {
    orbindex = 0;
    int row = 0;
    for(int j=0; j<M_Dn.cols(); j++)
    {
      if(S[orbindex] == -1)
      {
        M_Dn(row,j) = 1.0;
        row++;
      }
      orbindex += IDcount[j];
    }
  }
  else
  {
    int row = 0;
    int col = 0;
    while(row < M_Dn.rows())
    {
      M_Dn(row,col) = 1.0;
      row++;
      col++;
    }
  }
  *osXML << M_Dn;
  *osXML << "</Var>" << std::endl;
  *osXML << "</Determinant>" << std::endl;
  *osXML << "</SlaterDeterminant>" << std::endl;
  *osXML << "</DeterminantSet>" << std::endl;
  delete osXML;
  return true;
}

#endif


