//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_ATOMICORBITALBUILDER_H
#define QMCPLUSPLUS_ATOMICORBITALBUILDER_H


#include "Message/MPIObjectBase.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/lcao/RadialOrbitalSetBuilder.h"
#include "io/hdf_archive.h"

namespace qmcplusplus
{
/** atomic basisset builder
   * @tparam COT, CenteredOrbitalType = SoaAtomicBasisSet<RF,SH>
   *
   * Reimplement AtomiSPOSetBuilder.h
   */
template<typename COT>
class AOBasisBuilder : public MPIObjectBase
{
public:
  enum
  {
    DONOT_EXPAND    = 0,
    GAUSSIAN_EXPAND = 1,
    NATURAL_EXPAND,
    CARTESIAN_EXPAND,
    MOD_NATURAL_EXPAND
  };

private:
  //builder for a set of radial functors for this atom
  RadialOrbitalSetBuilder<COT> radFuncBuilder;

  bool addsignforM;
  int expandlm;
  std::string Morder;
  std::string sph;
  std::string basisType;
  std::string elementType;

  ///map for the radial orbitals
  std::map<std::string, int> RnlID;

  ///map for (n,l,m,s) to its quantum number index
  std::map<std::string, int> nlms_id;

public:
  AOBasisBuilder(const std::string& eName, Communicate* comm);

  bool put(xmlNodePtr cur);
  bool putH5(hdf_archive& hin);

  SPOSet* createSPOSetFromXML(xmlNodePtr cur) { return 0; }

  COT* createAOSet(xmlNodePtr cur);
  COT* createAOSetH5(hdf_archive& hin);

  int expandYlm(COT* aos, std::vector<int>& all_nl, int expandlm = DONOT_EXPAND);
};

template<typename COT>
AOBasisBuilder<COT>::AOBasisBuilder(const std::string& eName, Communicate* comm)
    : MPIObjectBase(comm),
      radFuncBuilder(comm),
      addsignforM(false),
      expandlm(GAUSSIAN_EXPAND),
      Morder("gaussian"),
      sph("default"),
      basisType("Numerical"),
      elementType(eName)
{
  // mmorales: for "Cartesian Gaussian", m is an integer that maps
  //           the component to Gamess notation, see Numerics/CartesianTensor.h
  nlms_id["n"] = q_n;
  nlms_id["l"] = q_l;
  nlms_id["m"] = q_m;
  nlms_id["s"] = q_s;
}

template<class COT>
bool AOBasisBuilder<COT>::put(xmlNodePtr cur)
{
  ReportEngine PRE("AtomicBasisBuilder", "put(xmlNodePtr)");
  std::string Normalized("yes");
  //Register valid attributes attributes
  OhmmsAttributeSet aAttrib;
  aAttrib.add(basisType, "type");
  aAttrib.add(sph, "angular");
  aAttrib.add(addsignforM, "expM");
  aAttrib.add(Morder, "expandYlm");
  aAttrib.add(Normalized, "normalized");
  aAttrib.put(cur);
  PRE.echo(cur);
  if (sph == "spherical")
    addsignforM = 1; //include (-1)^m
  if (Morder == "gaussian")
  {
    expandlm = GAUSSIAN_EXPAND;
  }
  else if (Morder == "natural")
  {
    expandlm = NATURAL_EXPAND;
  }
  else if (Morder == "no")
  {
    expandlm = DONOT_EXPAND;
  }
  else if (Morder == "pyscf")
  {
    expandlm    = MOD_NATURAL_EXPAND;
    addsignforM = 1;
    if (sph != "spherical")
    {
      APP_ABORT(" Error: expandYlm='pyscf' only compatible with angular='spherical'. Aborting.\n");
    }
  }

  if (sph == "cartesian" || Morder == "Gamess")
  {
    expandlm    = CARTESIAN_EXPAND;
    addsignforM = 0;
  }

  radFuncBuilder.Normalized = (Normalized == "yes");

  // Numerical basis is a special case
  if (basisType == "Numerical")
    radFuncBuilder.openNumericalBasisH5(cur);
  return true;
}

template<class COT>
bool AOBasisBuilder<COT>::putH5(hdf_archive& hin)
{
  ReportEngine PRE("AtomicBasisBuilder", "putH5(hin)");
  std::string CenterID, basisName, Normalized("yes");

  if (myComm->rank() == 0)
  {
    hin.read(sph, "angular");
    hin.read(CenterID, "elementType");
    hin.read(Normalized, "normalized");
    hin.read(Morder, "expandYlm");
    hin.read(basisName, "name");
    hin.read(basisType, "type");
    hin.read(addsignforM, "expM");
  }

  myComm->bcast(sph);
  myComm->bcast(Morder);
  myComm->bcast(CenterID);
  myComm->bcast(Normalized);
  myComm->bcast(basisName);
  myComm->bcast(basisType);
  myComm->bcast(addsignforM);

  app_log() << "<input node=\"atomicBasisSet\" name=\"" << basisName << "\" expandYlm=\"" << Morder << "\" angular=\""
            << sph << "\" elementType=\"" << CenterID << "\" normalized=\"" << Normalized << "\" type=\"" << basisType
            << "\" expM=\"" << addsignforM << "\" />" << std::endl;
  if (sph == "spherical")
    addsignforM = 1; //include (-1)^m
  if (Morder == "gaussian")
  {
    expandlm = GAUSSIAN_EXPAND;
  }
  else if (Morder == "natural")
  {
    expandlm = NATURAL_EXPAND;
  }
  else if (Morder == "no")
  {
    expandlm = DONOT_EXPAND;
  }
  else if (Morder == "pyscf")
  {
    expandlm    = MOD_NATURAL_EXPAND;
    addsignforM = 1;
    if (sph != "spherical")
    {
      APP_ABORT(" Error: expandYlm='pyscf' only compatible with angular='spherical'. Aborting.\n");
    }
  }

  if (sph == "cartesian" || Morder == "Gamess")
  {
    expandlm    = CARTESIAN_EXPAND;
    addsignforM = 0;
  }

  radFuncBuilder.Normalized = (Normalized == "yes");

  return true;
}


template<typename COT>
COT* AOBasisBuilder<COT>::createAOSet(xmlNodePtr cur)
{
  ReportEngine PRE("AtomicBasisBuilder", "createAOSet(xmlNodePtr)");
  app_log() << "  AO BasisSet for " << elementType << "\n";
  if (expandlm != CARTESIAN_EXPAND)
  {
    if (addsignforM)
      app_log() << "   Spherical Harmonics contain (-1)^m factor" << std::endl;
    else
      app_log() << "   Spherical Harmonics  DO NOT contain (-1)^m factor" << std::endl;
  }
  switch (expandlm)
  {
  case (GAUSSIAN_EXPAND):
    app_log() << "   Angular momentum m expanded according to Gaussian" << std::endl;
    break;
  case (NATURAL_EXPAND):
    app_log() << "   Angular momentum m expanded as -l, ... ,l" << std::endl;
    break;
  case (MOD_NATURAL_EXPAND):
    app_log() << "   Angular momentum m expanded as -l, ... ,l, with the exception of L=1 (1,-1,0)" << std::endl;
    break;
  case (CARTESIAN_EXPAND):
    app_log() << "   Angular momentum expanded in cartesian functions x^lx y^ly z^lz according to Gamess" << std::endl;
    break;
  default:
    app_log() << "   Angular momentum m is explicitly given." << std::endl;
  }
  QuantumNumberType nlms;
  std::string rnl;
  int Lmax(0); //maxmimum angular momentum of this center
  int num(0);  //the number of localized basis functions of this center
  //process the basic property: maximun angular momentum, the number of basis functions to be added
  std::vector<xmlNodePtr> radGroup;
  xmlNodePtr cur1 = cur->xmlChildrenNode;
  xmlNodePtr gptr = 0;
  while (cur1 != NULL)
  {
    std::string cname1((const char*)(cur1->name));
    if (cname1 == "basisGroup")
    {
      radGroup.push_back(cur1);
      const int l = std::stoi(XMLAttrString{cur1, "l"});
      Lmax  = std::max(Lmax, l);
      //expect that only Rnl is given
      if (expandlm == CARTESIAN_EXPAND)
        num += (l + 1) * (l + 2) / 2;
      else if (expandlm)
        num += 2 * l + 1;
      else
        num++;
    }
    else if (cname1 == "grid")
    {
      gptr = cur1;
    }
    cur1 = cur1->next;
  }
  //create a new set of atomic orbitals sharing a center with (Lmax, num)
  //if(addsignforM) the basis function has (-1)^m sqrt(2)Re(Ylm)
  COT* aos = new COT(Lmax, addsignforM);
  aos->LM.resize(num);
  aos->NL.resize(num);
  //Now, add distinct Radial Orbitals and (l,m) channels
  radFuncBuilder.setOrbitalSet(aos, elementType); //assign radial orbitals for the new center
  radFuncBuilder.addGrid(gptr, basisType);        //assign a radial grid for the new center
  std::vector<xmlNodePtr>::iterator it(radGroup.begin());
  std::vector<xmlNodePtr>::iterator it_end(radGroup.end());
  std::vector<int> all_nl;
  while (it != it_end)
  {
    cur1           = (*it);
    xmlAttrPtr att = cur1->properties;
    while (att != NULL)
    {
      std::string aname((const char*)(att->name));
      if (aname == "rid" || aname == "id")
      //accept id/rid
      {
        rnl = (const char*)(att->children->content);
      }
      else
      {
        std::map<std::string, int>::iterator iit = nlms_id.find(aname);
        if (iit != nlms_id.end())
        //valid for n,l,m,s
        {
          nlms[(*iit).second] = atoi((const char*)(att->children->content));
        }
      }
      att = att->next;
    }
    //add Ylm channels
    app_log() << "   R(n,l,m,s) " << nlms[0] << " " << nlms[1] << " " << nlms[2] << " " << nlms[3] << std::endl;
    std::map<std::string, int>::iterator rnl_it = RnlID.find(rnl);
    if (rnl_it == RnlID.end())
    {
      int nl = aos->RnlID.size();
      if (radFuncBuilder.addRadialOrbital(cur1, basisType, nlms))
        RnlID[rnl] = nl;
      all_nl.push_back(nl);
    }
    else
    {
      all_nl.push_back((*rnl_it).second);
    }
    ++it;
  }

  if (expandYlm(aos, all_nl, expandlm) != num)
    APP_ABORT("expandYlm doesn't match the number of basis.");
  radFuncBuilder.finalize();
  //aos->Rmax can be set small
  //aos->setRmax(0);
  aos->setBasisSetSize(-1);
  app_log() << "   Maximum Angular Momentum  = " << aos->Ylm.lmax() << std::endl
            << "   Number of Radial functors = " << aos->RnlID.size() << std::endl
            << "   Basis size                = " << aos->getBasisSetSize() << "\n\n";
  return aos;
}


template<typename COT>
COT* AOBasisBuilder<COT>::createAOSetH5(hdf_archive& hin)
{
  ReportEngine PRE("AOBasisBuilder:", "createAOSetH5(std::string)");
  app_log() << "  AO BasisSet for " << elementType << "\n";

  if (expandlm != CARTESIAN_EXPAND)
  {
    if (addsignforM)
      app_log() << "   Spherical Harmonics contain (-1)^m factor" << std::endl;
    else
      app_log() << "   Spherical Harmonics  DO NOT contain (-1)^m factor" << std::endl;
  }
  switch (expandlm)
  {
  case (GAUSSIAN_EXPAND):
    app_log() << "   Angular momentum m expanded according to Gaussian" << std::endl;
    break;
  case (NATURAL_EXPAND):
    app_log() << "   Angular momentum m expanded as -l, ... ,l" << std::endl;
    break;
  case (MOD_NATURAL_EXPAND):
    app_log() << "   Angular momentum m expanded as -l, ... ,l, with the exception of L=1 (1,-1,0)" << std::endl;
    break;
  case (CARTESIAN_EXPAND):
    app_log() << "   Angular momentum expanded in cartesian functions x^lx y^ly z^lz according to Gamess" << std::endl;
    break;
  default:
    app_log() << "   Angular momentum m is explicitly given." << std::endl;
  }

  QuantumNumberType nlms;
  std::string rnl;
  int Lmax(0); //maxmimum angular momentum of this center
  int num(0);  //the number of localized basis functions of this center

  int numbasisgroups(0);
  if (myComm->rank() == 0)
  {
    if (!hin.readEntry(numbasisgroups, "NbBasisGroups"))
      PRE.error("Could not read NbBasisGroups in H5; Probably Corrupt H5 file", true);
  }
  myComm->bcast(numbasisgroups);

  for (int i = 0; i < numbasisgroups; i++)
  {
    std::string basisGroupID = "basisGroup" + std::to_string(i);
    int l(0);
    if (myComm->rank() == 0)
    {
      hin.push(basisGroupID);
      hin.read(l, "l");
      hin.pop();
    }
    myComm->bcast(l);

    Lmax = std::max(Lmax, l);
    //expect that only Rnl is given
    if (expandlm == CARTESIAN_EXPAND)
      num += (l + 1) * (l + 2) / 2;
    else if (expandlm)
      num += 2 * l + 1;
    else
      num++;
  }

  COT* aos = new COT(Lmax, addsignforM);
  aos->LM.resize(num);
  aos->NL.resize(num);
  //Now, add distinct Radial Orbitals and (l,m) channels
  radFuncBuilder.setOrbitalSet(aos, elementType); //assign radial orbitals for the new center
  radFuncBuilder.addGridH5(hin);                  //assign a radial grid for the new center
  std::vector<int> all_nl;
  for (int i = 0; i < numbasisgroups; i++)
  {
    std::string basisGroupID = "basisGroup" + std::to_string(i);
    if (myComm->rank() == 0)
    {
      hin.push(basisGroupID);
      hin.read(rnl, "rid");
      hin.read(nlms[0], "n");
      hin.read(nlms[1], "l");
    }
    myComm->bcast(rnl);
    myComm->bcast(nlms[0]);
    myComm->bcast(nlms[1]);

    //add Ylm channels
    app_log() << "   R(n,l,m,s) " << nlms[0] << " " << nlms[1] << " " << nlms[2] << " " << nlms[3] << std::endl;
    std::map<std::string, int>::iterator rnl_it = RnlID.find(rnl);
    if (rnl_it == RnlID.end())
    {
      int nl = aos->RnlID.size();
      if (radFuncBuilder.addRadialOrbitalH5(hin, basisType, nlms))
        RnlID[rnl] = nl;
      all_nl.push_back(nl);
    }
    else
    {
      all_nl.push_back((*rnl_it).second);
    }

    if (myComm->rank() == 0)
      hin.pop();
  }

  if (expandYlm(aos, all_nl, expandlm) != num)
    APP_ABORT("expandYlm doesn't match the number of basis.");
  radFuncBuilder.finalize();
  //aos->Rmax can be set small
  //aos->setRmax(0);
  aos->setBasisSetSize(-1);
  app_log() << "   Maximum Angular Momentum  = " << aos->Ylm.lmax() << std::endl
            << "   Number of Radial functors = " << aos->RnlID.size() << std::endl
            << "   Basis size                = " << aos->getBasisSetSize() << "\n\n";
  return aos;
}


template<typename COT>
int AOBasisBuilder<COT>::expandYlm(COT* aos, std::vector<int>& all_nl, int expandlm)
{
  int num = 0;
  if (expandlm == GAUSSIAN_EXPAND)
  {
    app_log() << "Expanding Ylm according to Gaussian98" << std::endl;
    for (int nl = 0; nl < aos->RnlID.size(); nl++)
    {
      int l = aos->RnlID[nl][q_l];
      app_log() << "Adding " << 2 * l + 1 << " spherical orbitals for l= " << l << std::endl;
      switch (l)
      {
      case (0):
        aos->LM[num] = aos->Ylm.index(0, 0);
        aos->NL[num] = nl;
        num++;
        break;
      case (1): //px(1),py(-1),pz(0)
        aos->LM[num] = aos->Ylm.index(1, 1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1, -1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1, 0);
        aos->NL[num] = nl;
        num++;
        break;
      default: //0,1,-1,2,-2,...,l,-l
        aos->LM[num] = aos->Ylm.index(l, 0);
        aos->NL[num] = nl;
        num++;
        for (int tm = 1; tm <= l; tm++)
        {
          aos->LM[num] = aos->Ylm.index(l, tm);
          aos->NL[num] = nl;
          num++;
          aos->LM[num] = aos->Ylm.index(l, -tm);
          aos->NL[num] = nl;
          num++;
        }
        break;
      }
    }
  }
  else if (expandlm == MOD_NATURAL_EXPAND)
  {
    app_log() << "Expanding Ylm as L=1 as (1,-1,0) and L>1 as -l,-l+1,...,l-1,l" << std::endl;
    for (int nl = 0; nl < aos->RnlID.size(); nl++)
    {
      int l = aos->RnlID[nl][q_l];
      app_log() << "   Adding " << 2 * l + 1 << " spherical orbitals" << std::endl;
      if (l == 1)
      {
        //px(1),py(-1),pz(0)
        aos->LM[num] = aos->Ylm.index(1, 1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1, -1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1, 0);
        aos->NL[num] = nl;
        num++;
      }
      else
      {
        for (int tm = -l; tm <= l; tm++, num++)
        {
          aos->LM[num] = aos->Ylm.index(l, tm);
          aos->NL[num] = nl;
        }
      }
    }
  }
  else if (expandlm == NATURAL_EXPAND)
  {
    app_log() << "Expanding Ylm as -l,-l+1,...,l-1,l" << std::endl;
    for (int nl = 0; nl < aos->RnlID.size(); nl++)
    {
      int l = aos->RnlID[nl][q_l];
      app_log() << "   Adding " << 2 * l + 1 << " spherical orbitals" << std::endl;
      for (int tm = -l; tm <= l; tm++, num++)
      {
        aos->LM[num] = aos->Ylm.index(l, tm);
        aos->NL[num] = nl;
      }
    }
  }
  else if (expandlm == CARTESIAN_EXPAND)
  {
    app_log() << "Expanding Ylm (angular function) according to Gamess using cartesian gaussians" << std::endl;
    for (int nl = 0; nl < aos->RnlID.size(); nl++)
    {
      int l = aos->RnlID[nl][q_l];
      app_log() << "Adding " << (l + 1) * (l + 2) / 2 << " cartesian gaussian orbitals for l= " << l << std::endl;
      int nbefore = 0;
      for (int i = 0; i < l; i++)
        nbefore += (i + 1) * (i + 2) / 2;
      for (int i = 0; i < (l + 1) * (l + 2) / 2; i++)
      {
        aos->LM[num] = nbefore + i;
        aos->NL[num] = nl;
        num++;
      }
    }
  }
  else
  {
    for (int ind = 0; ind < all_nl.size(); ind++)
    {
      int nl = all_nl[ind];
      int l  = aos->RnlID[nl][q_l];
      int m  = aos->RnlID[nl][q_m];
      //assign the index for real Spherical Harmonic with (l,m)
      aos->LM[num] = aos->Ylm.index(l, m);
      //assign the index for radial orbital with (n,l)
      aos->NL[num] = nl;
      //increment number of basis functions
      num++;
    }
  }
  return num;
}

} // namespace qmcplusplus
#endif
