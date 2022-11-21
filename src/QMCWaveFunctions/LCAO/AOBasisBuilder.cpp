//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "AOBasisBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "RadialOrbitalSetBuilder.h"
#include "SoaAtomicBasisSet.h"
#include "MultiQuinticSpline1D.h"
#include "MultiFunctorAdapter.h"
#include "Numerics/SoaCartesianTensor.h"
#include "Numerics/SoaSphericalTensor.h"

namespace qmcplusplus
{
template<typename COT>
AOBasisBuilder<COT>::AOBasisBuilder(const std::string& eName, Communicate* comm)
    : MPIObjectBase(comm),
      addsignforM(false),
      expandlm(GAUSSIAN_EXPAND),
      Morder("gaussian"),
      sph("default"),
      basisType("Numerical"),
      elementType(eName),
      Normalized("yes")
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
    expandlm = GAUSSIAN_EXPAND;
  else if (Morder == "natural")
    expandlm = NATURAL_EXPAND;
  else if (Morder == "no")
    expandlm = DONOT_EXPAND;
  else if (Morder == "pyscf")
  {
    expandlm    = MOD_NATURAL_EXPAND;
    addsignforM = 1;
    if (sph != "spherical")
    {
      myComm->barrier_and_abort(" Error: expandYlm='pyscf' only compatible with angular='spherical'. Aborting.\n");
    }
  }

  if (sph == "cartesian" || Morder == "Gamess")
  {
    expandlm    = CARTESIAN_EXPAND;
    addsignforM = 0;
  }

  if (Morder == "Dirac")
  {
    expandlm    = DIRAC_CARTESIAN_EXPAND;
    addsignforM = 0;
    if (sph != "cartesian")
      myComm->barrier_and_abort(" Error: expandYlm='Dirac' only compatible with angular='cartesian'. Aborting\n");
  }

  // Numerical basis is a special case
  if (basisType == "Numerical")
    myComm->barrier_and_abort("Purely numerical atomic orbitals are not supported any longer.");

  return true;
}

template<class COT>
bool AOBasisBuilder<COT>::putH5(hdf_archive& hin)
{
  ReportEngine PRE("AtomicBasisBuilder", "putH5(hin)");
  std::string CenterID, basisName;

  if (myComm->rank() == 0)
  {
    hin.read(sph, "angular");
    hin.read(CenterID, "elementType");
    hin.read(Normalized, "normalized");
    hin.read(Morder, "expandYlm");
    hin.read(basisName, "name");
  }

  myComm->bcast(sph);
  myComm->bcast(Morder);
  myComm->bcast(CenterID);
  myComm->bcast(Normalized);
  myComm->bcast(basisName);
  myComm->bcast(basisType);
  myComm->bcast(addsignforM);

  if (sph == "spherical")
    addsignforM = 1; //include (-1)^m

  if (Morder == "gaussian")
    expandlm = GAUSSIAN_EXPAND;
  else if (Morder == "natural")
    expandlm = NATURAL_EXPAND;
  else if (Morder == "no")
    expandlm = DONOT_EXPAND;
  else if (Morder == "pyscf")
  {
    expandlm    = MOD_NATURAL_EXPAND;
    addsignforM = 1;
    if (sph != "spherical")
    {
      myComm->barrier_and_abort(" Error: expandYlm='pyscf' only compatible with angular='spherical'. Aborting.\n");
    }
  }

  if (sph == "cartesian" || Morder == "Gamess")
  {
    expandlm    = CARTESIAN_EXPAND;
    addsignforM = 0;
  }

  if (Morder == "Dirac")
  {
    expandlm    = DIRAC_CARTESIAN_EXPAND;
    addsignforM = 0;
    if (sph != "cartesian")
      myComm->barrier_and_abort(" Error: expandYlm='Dirac' only compatible with angular='cartesian'. Aborting\n");
  }
  app_log() << R"(<input node="atomicBasisSet" name=")" << basisName << "\" expandYlm=\"" << Morder << "\" angular=\""
            << sph << "\" elementType=\"" << CenterID << "\" normalized=\"" << Normalized << "\" type=\"" << basisType
            << "\" expM=\"" << addsignforM << "\" />" << std::endl;

  return true;
}


template<typename COT>
std::unique_ptr<COT> AOBasisBuilder<COT>::createAOSet(xmlNodePtr cur)
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
  case (DIRAC_CARTESIAN_EXPAND):
    app_log() << "   Angular momentum expanded in cartesian functions in DIRAC ordering" << std::endl;
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
      const int l = std::stoi(getXMLAttributeValue(cur1, "l"));
      Lmax        = std::max(Lmax, l);
      //expect that only Rnl is given
      if (expandlm == CARTESIAN_EXPAND || expandlm == DIRAC_CARTESIAN_EXPAND)
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
  auto aos = std::make_unique<COT>(Lmax, addsignforM);
  aos->LM.resize(num);
  aos->NL.resize(num);

  //Now, add distinct Radial Orbitals and (l,m) channels
  RadialOrbitalSetBuilder<COT> radFuncBuilder(myComm, *aos);
  radFuncBuilder.Normalized = (Normalized == "yes");
  radFuncBuilder.addGrid(gptr, basisType); //assign a radial grid for the new center
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

  if (expandYlm(aos.get(), all_nl, expandlm) != num)
    myComm->barrier_and_abort("expandYlm doesn't match the number of basis.");
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
std::unique_ptr<COT> AOBasisBuilder<COT>::createAOSetH5(hdf_archive& hin)
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
  case (DIRAC_CARTESIAN_EXPAND):
    app_log() << "   Angular momentum expanded in cartesian functions in DIRAC ordering" << std::endl;
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
    if (expandlm == CARTESIAN_EXPAND || expandlm == DIRAC_CARTESIAN_EXPAND)
      num += (l + 1) * (l + 2) / 2;
    else if (expandlm)
      num += 2 * l + 1;
    else
      num++;
  }

  //create a new set of atomic orbitals sharing a center with (Lmax, num)
  //if(addsignforM) the basis function has (-1)^m sqrt(2)Re(Ylm)
  auto aos = std::make_unique<COT>(Lmax, addsignforM);
  aos->LM.resize(num);
  aos->NL.resize(num);

  //Now, add distinct Radial Orbitals and (l,m) channels
  RadialOrbitalSetBuilder<COT> radFuncBuilder(myComm, *aos);
  radFuncBuilder.Normalized = (Normalized == "yes");
  radFuncBuilder.addGridH5(hin); //assign a radial grid for the new center
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

  if (expandYlm(aos.get(), all_nl, expandlm) != num)
    myComm->barrier_and_abort("expandYlm doesn't match the number of basis.");
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
  else if (expandlm == DIRAC_CARTESIAN_EXPAND)
  {
    app_log() << "Expanding Ylm (angular function) according to DIRAC using cartesian gaussians" << std::endl;
    for (int nl = 0; nl < aos->RnlID.size(); nl++)
    {
      int l = aos->RnlID[nl][q_l];
      app_log() << "Adding " << (l + 1) * (l + 2) / 2 << " cartesian gaussian orbitals for l= " << l << std::endl;
      int nbefore = 0;
      for (int i = 0; i < l; i++)
        nbefore += (i + 1) * (i + 2) / 2;
      switch (l)
      {
      case (0):
        aos->LM[num] = nbefore + 0;
        aos->NL[num] = nl;
        num++;
        break;
      case (1):
        aos->LM[num] = nbefore + 0;
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1;
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2;
        aos->NL[num] = nl;
        num++;
        break;
      case (2):
        aos->LM[num] = nbefore + 0; //xx
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 3; //xy
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 4; //xz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1; //yy
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 5; //yz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2; //zz
        aos->NL[num] = nl;
        num++;
        break;
      case (3):
        aos->LM[num] = nbefore + 0; //xxx
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 3; //xxy
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 4; //xxz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 5; //xyy
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 9; //xyz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 7; //xzz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1; //yyy
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 6; //yyz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 8; //yzz
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2; //zzz
        aos->NL[num] = nl;
        num++;
        break;
      case (4):
        aos->LM[num] = nbefore + 0; //400
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 3; //310
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 4; //301
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 9; //220
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 12; //211
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 10; //202
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 5; //130
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 13; //121
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 14; //112
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 7; //103
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1; //040
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 6; //031
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 11; //022
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 8; //013
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2; //004
        aos->NL[num] = nl;
        num++;
        break;
      case (5):
        aos->LM[num] = nbefore + 0; //500
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 3; //410
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 4; //401
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 9; //320
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 15; //311
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 10; //302
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 11; //230
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 18; //221
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 19; //212
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 13; //203
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 5; //140
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 16; //131
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 20; //122
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 17; //113
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 7; //104
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1; //050
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 6; //041
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 12; //032
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 14; //023
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 8; //014
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2; //005
        aos->NL[num] = nl;
        num++;
        break;
      case (6):
        aos->LM[num] = nbefore + 0; //600
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 3; //510
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 4; //501
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 9; //420
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 15; //411
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 10; //402
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 18; //330
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 21; //321
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 22; //312
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 19; //303
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 11; //240
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 23; //231
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 27; //222
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 25; //213
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 13; //204
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 5; //150
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 16; //141
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 24; //132
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 26; //123
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 17; //114
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 7; //105
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 1; //060
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 6; //051
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 12; //042
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 20; //033
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 14; //024
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 8; //015
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = nbefore + 2; //006
        aos->NL[num] = nl;
        num++;
        break;
      default:
        myComm->barrier_and_abort("Cartesian Tensor only defined up to Lmax=6. Aborting\n");
        break;
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

template class AOBasisBuilder<
    SoaAtomicBasisSet<MultiQuinticSpline1D<QMCTraits::RealType>, SoaCartesianTensor<QMCTraits::RealType>>>;
template class AOBasisBuilder<
    SoaAtomicBasisSet<MultiQuinticSpline1D<QMCTraits::RealType>, SoaSphericalTensor<QMCTraits::RealType>>>;
template class AOBasisBuilder<SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<QMCTraits::RealType>>,
                                                SoaCartesianTensor<QMCTraits::RealType>>>;
template class AOBasisBuilder<SoaAtomicBasisSet<MultiFunctorAdapter<GaussianCombo<QMCTraits::RealType>>,
                                                SoaSphericalTensor<QMCTraits::RealType>>>;
template class AOBasisBuilder<
    SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<QMCTraits::RealType>>, SoaCartesianTensor<QMCTraits::RealType>>>;
template class AOBasisBuilder<
    SoaAtomicBasisSet<MultiFunctorAdapter<SlaterCombo<QMCTraits::RealType>>, SoaSphericalTensor<QMCTraits::RealType>>>;

} // namespace qmcplusplus
