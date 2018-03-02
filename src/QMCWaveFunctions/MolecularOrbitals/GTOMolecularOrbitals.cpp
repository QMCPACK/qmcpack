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


#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/DetSetBuilderWithBasisSet.h"
#include "QMCWaveFunctions/MolecularOrbitals/GTOMolecularOrbitals.h"

namespace qmcplusplus
{

GTOMolecularOrbitals::GTOMolecularOrbitals(ParticleSet& els, TrialWaveFunction& psi,
    ParticleSet& ions):
  OrbitalBuilderBase(els,psi), IonSys(ions), Normalized(false),BasisSet(0), d_table(0)
{
  //int d_ie = DistanceTable::add(ions,els);
  d_table = DistanceTable::add(ions,els);
  nlms_id["n"] = q_n;
  nlms_id["l"] = q_l;
  nlms_id["m"] = q_m;
  nlms_id["s"] = q_s;
}

bool GTOMolecularOrbitals::put(xmlNodePtr cur)
{
  LOGMSG("GTOMolecularOrbitals::put")
  DetSetBuilderWithBasisSet<GTOMolecularOrbitals> spobuilder(targetPtcl,targetPsi,*this);
  if(spobuilder.put(cur))
  {
    BasisSet->resize(spobuilder.NumPtcl);
    return true;
  }
  else
  {
    return false;
  }
}

GTOMolecularOrbitals::BasisSetType*
GTOMolecularOrbitals::addBasisSet(xmlNodePtr cur)
{
  if(!BasisSet)
    BasisSet = new BasisSetType(IonSys.getSpeciesSet().getTotalNum());
  QuantumNumberType nlms;
  std::string rnl;
  //current number of centers
  int ncenters = CenterID.size();
  int activeCenter;
  int gridmode = -1;
  bool addsignforM = true;
  std::string  sph("spherical"), Morder("gaussian");
  //go thru the tree
  cur = cur->xmlChildrenNode;
  while(cur!=NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == basis_tag || cname == "atomicBasisSet")
    {
      int expandlm = GAUSSIAN_EXPAND;
      std::string abasis("invalid"), norm("no");
      //Register valid attributes attributes
      OhmmsAttributeSet aAttrib;
      aAttrib.add(abasis,"elementType");
      aAttrib.add(abasis,"species");
      aAttrib.add(sph,"angular");
      aAttrib.add(addsignforM,"expM");
      aAttrib.add(Morder,"expandYlm");
      aAttrib.add(norm,"normalized");
      aAttrib.put(cur);
      if(norm == "yes")
        Normalized=true;
      else
        Normalized=false;
      if(abasis == "invalid")
        continue;
      if(sph == "spherical")
        addsignforM=true; //include (-1)^m
      if(sph == "cartesian")
        addsignforM=false;
      if(Morder == "gaussian")
      {
        expandlm = GAUSSIAN_EXPAND;
      }
      else
        if(Morder == "natural")
        {
          expandlm = NATURAL_EXPAND;
        }
        else
          if(Morder == "no")
          {
            expandlm = DONOT_EXPAND;
          }
      if(addsignforM)
        LOGMSG("Spherical Harmonics contain (-1)^m factor")
            else
              LOGMSG("Spherical Harmonics  DO NOT contain (-1)^m factor")
              std::map<std::string,int>::iterator it = CenterID.find(abasis); //search the species name
      if(it == CenterID.end())
        //add the name to the map CenterID
      {
        //CenterID[abasis] = activeCenter = ncenters++;
        CenterID[abasis]=activeCenter=IonSys.getSpeciesSet().findSpecies(abasis);
        int Lmax(0); //maxmimum angular momentum of this center
        int num(0);//the number of localized basis functions of this center
        //process the basic property: maximun angular momentum, the number of basis functions to be added
        std::vector<xmlNodePtr> radGroup;
        xmlNodePtr cur1 = cur->xmlChildrenNode;
        xmlNodePtr gptr=0;
        while(cur1 != NULL)
        {
          std::string cname1((const char*)(cur1->name));
          if(cname1 == basisfunc_tag || cname1 == "basisGroup")
          {
            radGroup.push_back(cur1);
            int l=atoi((const char*)(xmlGetProp(cur1, (const xmlChar *)"l")));
            Lmax = std::max(Lmax,l);
            if(expandlm)
              num += 2*l+1;
            else
              num++;
          }
          cur1 = cur1->next;
        }
        LOGMSG("Adding a center " << abasis << " centerid "<< CenterID[abasis])
        LOGMSG("Maximum angular momentum    = " << Lmax)
        LOGMSG("Number of centered orbitals = " << num)
        //create a new set of atomic orbitals sharing a center with (Lmax, num)
        //if(addsignforM) the basis function has (-1)^m sqrt(2)Re(Ylm)
        CenteredOrbitalType* aos = new CenteredOrbitalType(Lmax,addsignforM);
        aos->LM.resize(num);
        aos->NL.resize(num);
        //Now, add distinct Radial Orbitals and (l,m) channels
        num=0;
        std::vector<xmlNodePtr>::iterator it(radGroup.begin());
        std::vector<xmlNodePtr>::iterator it_end(radGroup.end());
        while(it != it_end)
        {
          cur1 = (*it);
          xmlAttrPtr att = cur1->properties;
          while(att != NULL)
          {
            std::string aname((const char*)(att->name));
            if(aname == "rid" || aname == "id")
              //accept id/rid
            {
              rnl = (const char*)(att->children->content);
            }
            else
            {
              std::map<std::string,int>::iterator iit = nlms_id.find(aname);
              if(iit != nlms_id.end())
                //valid for n,l,m,s
              {
                nlms[(*iit).second] = atoi((const char*)(att->children->content));
              }
            }
            att = att->next;
          }
          LOGMSG("\n(n,l,m,s) " << nlms[0] << " " << nlms[1] << " " << nlms[2] << " " << nlms[3])
          //add Ylm channels
          num = expandYlm(rnl,nlms,num,aos,cur1,expandlm);
          ++it;
        }
        LOGMSG("Checking the order of angular momentum ")
        copy(aos->LM.begin(), aos->LM.end(), std::ostream_iterator<int>(app_log()," "));
        app_log() << std::endl;
        //add the new atomic basis to the basis set
        BasisSet->add(aos,activeCenter);
      }
      else
      {
        WARNMSG("Species " << abasis << " is already initialized. Ignore the input.")
      }
    }
    cur = cur->next;
  }
  if(BasisSet)
  {
    BasisSet->setTable(d_table);
    XMLReport("The total number of basis functions " << BasisSet->TotalBasis)
    return BasisSet;
  }
  else
  {
    ERRORMSG("BasisSet is not initialized.")
    return NULL;
  }
}

int
GTOMolecularOrbitals::expandYlm(const std::string& rnl, const QuantumNumberType& nlms,
                                int num, CenteredOrbitalType* aos, xmlNodePtr cur1,
                                int expandlm)
{
  if(Normalized)
  {
    LOGMSG("The gaussian-type orbitals are normalized")
  }
  else
  {
    LOGMSG("The gaussian-type orbitals are Not normalized")
  }
  if(expandlm == GAUSSIAN_EXPAND)
  {
    std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
    if(rnl_it == RnlID.end())
    {
      int l = nlms[q_l];
      int nl = aos->Rnl.size();
      RadialOrbitalType *sto = new RadialOrbitalType(l,Normalized);
      sto->putBasisGroup(cur1);
      aos->Rnl.push_back(sto);
      aos->RnlID.push_back(nlms);
      RnlID[rnl] = nl;
      XMLReport("Adding a Radial Orbital and Expanding Ylm according to Gaussian98.")
      XMLReport("Adding " << 2*l+1 << " spherical orbitals for l= " << l)
      switch (l)
      {
      case(0):
        aos->LM[num] = aos->Ylm.index(0,0);
        aos->NL[num] = nl;
        num++;
        break;
      case(1)://px(1),py(-1),pz(0)
        aos->LM[num] = aos->Ylm.index(1,1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1,-1);
        aos->NL[num] = nl;
        num++;
        aos->LM[num] = aos->Ylm.index(1,0);
        aos->NL[num] = nl;
        num++;
        break;
      default://0,1,-1,2,-2,...,l,-l
        aos->LM[num] = aos->Ylm.index(l,0);
        aos->NL[num] = nl;
        num++;
        for(int tm=1; tm<=l; tm++)
        {
          aos->LM[num] = aos->Ylm.index(l,tm);
          aos->NL[num] = nl;
          num++;
          aos->LM[num] = aos->Ylm.index(l,-tm);
          aos->NL[num] = nl;
          num++;
        }
        break;
      }
    }
  }
  else
    if(expandlm == NATURAL_EXPAND)
    {
      std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
      if(rnl_it == RnlID.end())
        //only when rid is different
      {
        XMLReport("Expanding Ylm as -l,-l+1,...,l-1,l")
        int nl = aos->Rnl.size();
        int l = nlms[q_l];
        RadialOrbitalType *sto = new RadialOrbitalType(l,Normalized);
        sto->putBasisGroup(cur1);
        aos->Rnl.push_back(sto);
        aos->RnlID.push_back(nlms);
        RnlID[rnl] = nl;
        XMLReport("Adding " << 2*l+1 << " spherical orbitals")
        for(int tm=-l; tm<=l; tm++,num++)
        {
          aos->LM[num] = aos->Ylm.index(l,tm);
          aos->NL[num] = nl;
        }
      }
    }
    else
    {
      XMLReport("Ylm is NOT expanded.")
      aos->LM[num] = aos->Ylm.index(nlms[q_l],nlms[q_m]);
      std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
      if(rnl_it == RnlID.end())
      {
        int nl = aos->Rnl.size();
        RadialOrbitalType *sto = new RadialOrbitalType(nlms[q_l],Normalized);
        sto->putBasisGroup(cur1);
        aos->Rnl.push_back(sto);
        aos->RnlID.push_back(nlms);
        RnlID[rnl] = nl;
        aos->NL[num] = nl;
      }
      else
      {
        aos->NL[num]=(*rnl_it).second;
      }
      num++; //increment by 1
    }
  return num;
}
}
