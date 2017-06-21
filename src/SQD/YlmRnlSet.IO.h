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
    
    


#ifndef OHMMS_YLMRNLSET_ONGRID_IO_H
#define OHMMS_YLMRNLSET_ONGRID_IO_H

#include <io/hdf_archive.h>
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Numerics/Transform2GridFunctor.h"
#include <fstream>

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
    if (!(xmlStrcmp(cur->name, (const xmlChar *) "orbital")))
    {
      XMLReport("Found Orbital");
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
      /*check to see if orbital with same quantum numbers has
      already been added*/
      NLMSIndex nlms(n,l,m,s);
      NLMS_Map_t::iterator it = OccNo.find(nlms);
      if(it == OccNo.end())
      {
        add(n,l,m,s,occ);
        XMLReport("Adding Orbital: n=" << n << ", l=" << l <<
                  ", m=" << m << ", s=" << s << ", occ=" << occ);
        OccNo[nlms] = 1;
        //count number of spin-up and down orbitals
        if(s == 1)
          Nup++;
        else
          Ndown++;
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
                               const std::vector<value_type>& eigVal)
{
#ifdef HAVE_LIBHDF5
  //print to file fname
  std::string HDFfname = RootName + ".h5";
  hdf_archive afile;
  afile.create(HDFfname);
  TinyVector<int,2> res_version(0,1);
  afile.write(res_version,"version");
  int orbindex = 0;
  //main group
  afile.push("radial_basis_states");
  afile.write(CuspParam,"cusp");
  //output grid information
  //Vector<double> grid_params(3);
  //grid_params[0] = m_grid->rmin();
  //grid_params[1] = m_grid->rmax();
  //grid_params[2] = static_cast<double>(m_grid->size());
  //HDFAttribIO<Vector<int> > GridHDFTypeOut(grid_type);
  //GridHDFTypeOut.write(group_id_grid,"grid_type");
  //HDFAttribIO<Vector<double> > GridHDFParamOut(grid_params);
  //GridHDFParamOut.write(group_id_grid,"params");
  //output orbital information
  //output only unique radial functions
  int rmax_safe=0;
  for(int orb=0; orb<NumUniqueOrb; orb++)
  {
    char grpname[128];
    sprintf(grpname,"orbital%04d",orb);
    afile.push(grpname);
    //double tt=0.0;
    //afile.write(tt,"cusp");
    //if(!L[orbindex]) tt=-CuspParam;
    afile.write(eigVal[orbindex],"eigenvalue");
    //for qmc code, the radial function u(r)/r^{l+1}
    //Note that we are print out u(r) not
    //int power_raised=L[orbindex]+1;
    //HDFAttribIO<int> i_out(power_raised);
    //i_out.write(group_id_orb,"power");
    //output the eigenvalue
    //Vector<double> EigVal(1);
    //EigVal[0] = eigVal[orbindex];
    //HDFAttribIO<Vector<double> > EigValHDFOut(EigVal);
    //EigValHDFOut.write(group_id_orb,"eigenvalue");
    //ouptput quantum numbers
    TinyVector<int,3> QuantumNo(N[orbindex],L[orbindex],1);
    afile.write(QuantumNo,"quantum_numbers");
    //vector to store the radial orbital
    //determine the maximum radius of the orbital
    //max_rad += 2;
    //rad_orb.resize(max_rad);
    //
    std::vector<double> uofr(m_grid->size(),0.0);
    for(int i=0; i<m_grid->size(); ++i)
      uofr[i] = (*psi[orbindex])(i);
    afile.write(uofr,"uofr");
    //test numerical derivatives to detect the uptake
    int max_rad = m_grid->size()-2;
    while(std::abs((uofr[max_rad]-uofr[max_rad-1])/((*m_grid)[max_rad]-(*m_grid)[max_rad-1])) < 1e-4)
      max_rad--;
    rmax_safe=std::max(rmax_safe,max_rad);
    std::cout << "Safe cutoff " << (*m_grid)[max_rad] << std::endl;
    int l_plus_one= 1+L[orbindex];
    for(int i=1; i<m_grid->size(); ++i)
      uofr[i] = uofr[i]/pow((*m_grid)[i],l_plus_one);
    uofr[0]=uofr[1];
    if(L[orbindex]==0)
      uofr[0] += CuspParam*((*m_grid)[1]-(*m_grid)[0]);
    afile.write(uofr,"radial_orbital");
    afile.write(l_plus_one,"power");
    afile.pop();
    orbindex += IDcount[orb];
  }
  afile.push("grid");
  {
    std::string gtype_copy=GridType+'\0';
    afile.write(gtype_copy,"type");
    double tt=m_grid->rmin();
    afile.write(tt,"ri");
    tt=m_grid->rmax();
    afile.write(tt,"rf");
    //write out the safe cutoff
    tt=(*m_grid)[rmax_safe];
    afile.write(tt,"rmax_safe");
    int npts=m_grid->size();
    afile.write(npts,"npts");
    m_grid->locate(rmax_safe);
    npts=m_grid->Loc+1;
    afile.write(npts,"npts_safe");
  }
  afile.pop();
  afile.pop();
  afile.close();
#endif
  return true;
}

/**
 *@param elementName name of the element
 *@param RootName name of the element
 *@param GridType the type of grid
 *@return true if succeeds
 *@brief Prints the basis information to an xml file named
 *"RootName.qmc.xml" for use in qmcPlusPlus.
 */
template<class GT>
bool YlmRnlSet<GT>::print_basis(const std::string& elementName,
                                const std::string& RootName,
                                const std::string& GridType)
{
    //matrices for spin-up and down orbitals
  Matrix<double> Mup;
  Matrix<double> Mdown;
  if(Restriction == "none")
  {
    Mup.resize(Nup,NumUniqueOrb);
    Mdown.resize(Ndown,NumUniqueOrb);
  }
  else
  {
    if(Nup >= Ndown)
    {
      Mup.resize(Nup,Nup);
      Mdown.resize(Ndown,Nup);
    }
    else
    {
      Mup.resize(Nup,Ndown);
      Mdown.resize(Ndown,Ndown);
    }
  }
  Mup = 0.0;
  Mdown = 0.0;
  std::string fnameXML = RootName + ".qmc.xml";
  std::string fnameHDF5 = RootName + ".h5";
  //construct xml tree
  xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
  xmlNodePtr qmc_root = xmlNewNode(NULL, BAD_CAST "qmcsystem");
  {
    //particleset for atom
    xmlNodePtr p1 = xmlNewNode(NULL,(const xmlChar*)"particleset");
    xmlNewProp(p1,(const xmlChar*)"name",(const xmlChar*)"atom");
    xmlNewProp(p1,(const xmlChar*)"size",(const xmlChar*)"1");
    xmlNodePtr p2 = xmlNewNode(NULL,(const xmlChar*)"group");
    xmlNewProp(p2,(const xmlChar*)"name",(const xmlChar*)elementName.c_str());
    std::ostringstream c;
    //c<<Nup+Ndown;
    c << Charge; //use charge
    xmlNodePtr p3 =xmlNewTextChild(p2,NULL,
                                   (const xmlChar*)"parameter", (const xmlChar*)c.str().c_str());
    xmlNewProp(p3,(const xmlChar*)"name",(const xmlChar*)"charge");
    xmlAddChild(p1,p2);
    xmlAddChild(qmc_root,p1);
  }
  {
    //particleset for electrons
    xmlNodePtr p1 = xmlNewNode(NULL,(const xmlChar*)"particleset");
    xmlNewProp(p1,(const xmlChar*)"name",(const xmlChar*)"e");
    if(Nup>0)
    {
      xmlNodePtr p2 = xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(p2,(const xmlChar*)"name",(const xmlChar*)"u");
      std::ostringstream s;
      s << Nup;
      xmlNewProp(p2,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
      xmlNodePtr p3 =xmlNewTextChild(p2,NULL, (const xmlChar*)"parameter", (const xmlChar*)"-1");
      xmlNewProp(p3,(const xmlChar*)"name",(const xmlChar*)"charge");
      xmlAddChild(p2,p3);
      xmlAddChild(p1,p2);
    }
    if(Ndown>0)
    {
      //down electrons
      xmlNodePtr p2 = xmlNewNode(NULL,(const xmlChar*)"group");
      xmlNewProp(p2,(const xmlChar*)"name",(const xmlChar*)"d");
      std::ostringstream s;
      s << Ndown;
      xmlNewProp(p2,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
      xmlNodePtr p3 =xmlNewTextChild(p2,NULL, (const xmlChar*)"parameter", (const xmlChar*)"-1");
      xmlNewProp(p3,(const xmlChar*)"name",(const xmlChar*)"charge");
      xmlAddChild(p2,p3);
      xmlAddChild(p1,p2);
    }
    xmlAddChild(qmc_root,p1);
  }
  {
    //wavefunction
    xmlNodePtr p1 = xmlNewNode(NULL,(const xmlChar*)"wavefunction");
    xmlNewProp(p1,(const xmlChar*)"name",(const xmlChar*)"psi0");
    xmlNewProp(p1,(const xmlChar*)"target",(const xmlChar*)"e");
    //determinantset
    xmlNodePtr p2 = xmlNewNode(NULL, (const xmlChar*) "determinantset");
    xmlNewProp(p2,(const xmlChar*)"type",(const xmlChar*)"MO");
    xmlNewProp(p2,(const xmlChar*)"transform",(const xmlChar*)"yes");
    xmlNewProp(p2,(const xmlChar*)"source",(const xmlChar*)"atom");
    {
      //basiset
      xmlNodePtr p3 = xmlNewNode(NULL, (const xmlChar*) "basisset");
      xmlNodePtr p4 = xmlNewNode(NULL, (const xmlChar*) "atomicBasisSet");
      xmlNewProp(p4,(const xmlChar*)"type",(const xmlChar*)"Numerical");
      xmlNewProp(p4,(const xmlChar*)"elementType",(const xmlChar*)elementName.c_str());
      xmlNewProp(p4,(const xmlChar*)"expandYlm",(const xmlChar*)"no");
      xmlNewProp(p4,(const xmlChar*)"href",(const xmlChar*)fnameHDF5.c_str());
      //{//grid
      //  std::ostringstream s1,s2,s3;
      //  s1<<m_grid->rmin(); s2<<m_grid->rmax(); s3<< m_grid->size();
      //  xmlNodePtr p5 = xmlNewNode(NULL, (const xmlChar*) "grid");
      //  xmlNewProp(p5,(const xmlChar*)"type",(const xmlChar*)GridType.c_str());
      //  xmlNewProp(p5,(const xmlChar*)"ri",(const xmlChar*)s1.str().c_str());
      //  xmlNewProp(p5,(const xmlChar*)"rf",(const xmlChar*)s2.str().c_str());
      //  xmlNewProp(p5,(const xmlChar*)"npts",(const xmlChar*)s3.str().c_str());
      //  xmlAddChild(p4,p5);
      //}
      //basisGroups
      int nup = 0;
      int ndown = 0;
      int nbasis = -1;
      for(int orb=0; orb<size(); orb++)
      {
        NLMIndex nlm(N[orb],L[orb],M[orb]);
        NLM_Map_t::iterator it = NLM.find(nlm);
        if(it == NLM.end())
        {
          NLM[nlm] = ++nbasis;
          char idname[128];
          char dsname[128];
          sprintf(idname,"R%03d",ID[orb]);
          sprintf(dsname,"orbital%04d",ID[orb]);
          xmlNodePtr p5 = xmlNewNode(NULL, (const xmlChar*) "basisGroup");
          xmlNewProp(p5,(const xmlChar*)"rid",(const xmlChar*)idname);
          xmlNewProp(p5,(const xmlChar*)"ds",(const xmlChar*)dsname);
          std::ostringstream s1,s2,s3;
          s1 << N[orb];
          s2<< L[orb];
          s3 << M[orb];
          xmlNewProp(p5,(const xmlChar*)"n",(const xmlChar*)s1.str().c_str());
          xmlNewProp(p5,(const xmlChar*)"l",(const xmlChar*)s2.str().c_str());
          xmlNewProp(p5,(const xmlChar*)"m",(const xmlChar*)s3.str().c_str());
          xmlNewProp(p5,(const xmlChar*)"zeta",(const xmlChar*)"1");
          if((*m_grid)(0)<std::numeric_limits<value_type>::epsilon())
            xmlNewProp(p5,(const xmlChar*)"imin",(const xmlChar*)"1");
          xmlAddChild(p4,p5);
        }
      }
      xmlAddChild(p3,p4);//add atomicBasisSet to basiset
      xmlAddChild(p2,p3);//basisset to determinantset
    }
    {
      xmlNodePtr p3 = xmlNewNode(NULL, (const xmlChar*) "slaterdeterminant");
      if(Nup>0)
      {
        std::ostringstream s;
        s << Nup;
        xmlNodePtr p4 = xmlNewNode(NULL,(const xmlChar*)"determinant");
        xmlNewProp(p4,(const xmlChar*)"spin",(const xmlChar*)"1");
        xmlNewProp(p4,(const xmlChar*)"id",(const xmlChar*)"detu");
        xmlNewProp(p4,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
        xmlAddChild(p3,p4);
      }
      if(Ndown>0)
      {
        std::ostringstream s;
        s << Ndown;
        xmlNodePtr p4 = xmlNewNode(NULL,(const xmlChar*)"determinant");
        xmlNewProp(p4,(const xmlChar*)"spin",(const xmlChar*)"-1");
        xmlNewProp(p4,(const xmlChar*)"id",(const xmlChar*)"detd");
        xmlNewProp(p4,(const xmlChar*)"size",(const xmlChar*)s.str().c_str());
        xmlAddChild(p3,p4);
      }
      xmlAddChild(p2,p3);//slaterdeterminant to determinantset
    }
    xmlAddChild(p1,p2);//determinantset to wavefunction
    xmlAddChild(qmc_root,p1);//wavefunction to qmcsystem
  }
  xmlDocSetRootElement(doc, qmc_root);
  xmlSaveFormatFile(fnameXML.c_str(),doc,1);
  xmlFreeDoc(doc);
  return true;
}

#endif


