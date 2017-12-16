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
    
    



#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/DetSetBuilderWithBasisSet.h"
#include "QMCTools/GridMolecularOrbitals.h"
#include "QMCTools/RGFBuilderBase.h"
#include "QMCTools/STO2GridBuilder.h"
#include "QMCTools/GTO2GridBuilder.h"
#include "QMCTools/Any2GridBuilder.h"
#include "QMCTools/NumericalRGFBuilder.h"

namespace qmcplusplus
{

GridMolecularOrbitals::GridMolecularOrbitals(ParticleSet& els, TrialWaveFunction& psi,
    ParticleSet& ions):
  OrbitalBuilderBase(els,psi), IonSys(ions), BasisSet(0), d_table(0), rbuilder(0)
{
  //int d_ie = DistanceTable::add(ions,els);
  d_table = DistanceTable::add(ions,els);
  nlms_id["n"] = q_n;
  nlms_id["l"] = q_l;
  nlms_id["m"] = q_m;
  nlms_id["s"] = q_s;
}

bool GridMolecularOrbitals::put(xmlNodePtr cur)
{
  LOGMSG("GridMolecularOrbitals::put")
  DetSetBuilderWithBasisSet<GridMolecularOrbitals> spobuilder(targetPtcl,targetPsi,*this);
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

GridMolecularOrbitals::BasisSetType*
GridMolecularOrbitals::addBasisSet(xmlNodePtr cur)
{
  if(!BasisSet)
    BasisSet = new BasisSetType(IonSys.getSpeciesSet().getTotalNum());
  QuantumNumberType nlms;
  std::string rnl;
  //current number of centers
  int ncenters = CenterID.size();
  int activeCenter;
  int gridmode = -1;
  bool addsignforM = false;
  std::string  sph("default"), Morder("gaussian");
  //go thru the tree
  cur = cur->xmlChildrenNode;
  std::map<std::string,RGFBuilderBase*> rbuilderlist;
  while(cur!=NULL)
  {
    std::string cname((const char*)(cur->name));
    if(cname == basis_tag || cname == "atomicBasisSet")
    {
      int expandlm = GAUSSIAN_EXPAND;
      std::string abasis("invalid"), btype("Numerical");
      //Register valid attributes attributes
      OhmmsAttributeSet aAttrib;
      aAttrib.add(abasis,"elementType");
      aAttrib.add(abasis,"species");
      aAttrib.add(btype,"type");
      aAttrib.add(sph,"angular");
      aAttrib.add(addsignforM,"expM");
      aAttrib.add(Morder,"expandYlm");
      aAttrib.put(cur);
      if(abasis == "invalid")
        continue;
      if(sph == "spherical")
        addsignforM=1; //include (-1)^m
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
              //search the species name
              std::map<std::string,int>::iterator it = CenterID.find(abasis);
      if(it == CenterID.end())
        //add the name to the map CenterID
      {
        if(btype == "Numerical" || btype == "NG" || btype == "HFNG")
        {
          rbuilder = new NumericalRGFBuilder(cur);
        }
        else
        {
          rbuilder = new Any2GridBuilder(cur);
        }
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
            //expect that only Rnl is given
            if(expandlm)
              num += 2*l+1;
            else
              num++;
          }
          else
            if(cname1 == "grid")
            {
              gptr = cur1;
            }
          cur1 = cur1->next;
        }
        XMLReport("Adding a center " << abasis << " centerid "<< CenterID[abasis])
        XMLReport("Maximum angular momentum    = " << Lmax)
        XMLReport("Number of centered orbitals = " << num)
        //create a new set of atomic orbitals sharing a center with (Lmax, num)
        //if(addsignforM) the basis function has (-1)^m sqrt(2)Re(Ylm)
        CenteredOrbitalType* aos = new CenteredOrbitalType(Lmax,addsignforM);
        aos->LM.resize(num);
        aos->NL.resize(num);
        //Now, add distinct Radial Orbitals and (l,m) channels
        num=0;
        rbuilder->setOrbitalSet(aos,abasis); //assign radial orbitals for the new center
        rbuilder->addGrid(gptr); //assign a radial grid for the new center
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
          XMLReport("\n(n,l,m,s) " << nlms[0] << " " << nlms[1] << " " << nlms[2] << " " << nlms[3])
          //add Ylm channels
          num = expandYlm(rnl,nlms,num,aos,cur1,expandlm);
          ++it;
        }
        //add the new atomic basis to the basis set
        BasisSet->add(aos,activeCenter);
#if !defined(HAVE_MPI)
        rbuilder->print(abasis,1);
#endif
        if(rbuilder)
        {
          delete rbuilder;
          rbuilder=0;
        }
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
    LOGMSG("The total number of basis functions " << BasisSet->TotalBasis)
    return BasisSet;
  }
  else
  {
    ERRORMSG("BasisSet is not initialized.")
    return NULL;
  }
}

int
GridMolecularOrbitals::expandYlm(const std::string& rnl, const QuantumNumberType& nlms,
                                 int num, CenteredOrbitalType* aos, xmlNodePtr cur1,
                                 int expandlm)
{
  if(expandlm == GAUSSIAN_EXPAND)
  {
    XMLReport("Expanding Ylm according to Gaussian98")
    std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
    if(rnl_it == RnlID.end())
    {
      int nl = aos->Rnl.size();
      if(rbuilder->addRadialOrbital(cur1,nlms))
      {
        RnlID[rnl] = nl;
        int l = nlms[q_l];
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
  }
  else
    if(expandlm == NATURAL_EXPAND)
    {
      XMLReport("Expanding Ylm as -l,-l+1,...,l-1,l")
      std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
      if(rnl_it == RnlID.end())
      {
        int nl = aos->Rnl.size();
        if(rbuilder->addRadialOrbital(cur1,nlms))
        {
          RnlID[rnl] = nl;
          int l = nlms[q_l];
          XMLReport("Adding " << 2*l+1 << " spherical orbitals")
          for(int tm=-l; tm<=l; tm++,num++)
          {
            aos->LM[num] = aos->Ylm.index(l,tm);
            aos->NL[num] = nl;
          }
        }
      }
    }
    else
    {
      XMLReport("Ylm is NOT expanded.")
      //assign the index for real Spherical Harmonic with (l,m)
      aos->LM[num] = aos->Ylm.index(nlms[q_l],nlms[q_m]);
      //radial orbitals: add only distinct orbitals
      std::map<std::string,int>::iterator rnl_it = RnlID.find(rnl);
      if(rnl_it == RnlID.end())
      {
        int nl = aos->Rnl.size();
        if(rbuilder->addRadialOrbital(cur1,nlms))
          //assign the index for radial orbital with (n,l)
        {
          aos->NL[num] = nl;
          RnlID[rnl] = nl;
        }
      }
      else
      {
        //assign the index for radial orbital with (n,l) if repeated
        XMLReport("Already added radial function for id: " << rnl)
        aos->NL[num] = (*rnl_it).second;
      }
      //increment number of basis functions
      num++;
    }
  return num;
}
}
