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
    
    
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/LocalCorePolPotential.h"
#include "Utilities/IteratorUtility.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

LocalCorePolPotential::LocalCorePolPotential(ParticleSet& ions,
    ParticleSet& els):
  FirstTime(true), eCoreCore(0.0), IonConfig(ions), d_ie(0), d_ii(0)
{
  //set the distance tables
  d_ie = DistanceTable::add(ions,els,DT_AOS);
  d_ii = DistanceTable::add(ions,DT_AOS);
  nCenters = ions.getTotalNum();
  nParticles = els.getTotalNum();
  InpCPP.resize(IonConfig.getSpeciesSet().getTotalNum(),0);
  Centers.resize(nCenters,0);
  CoreCoreDipole.resize(nCenters,0.0);
  CoreElDipole.resize(nCenters,nParticles);
  CoreElDipole = 0.0;
}

/** destructor
 *
 * Delete InpCPP.
 */
LocalCorePolPotential::~LocalCorePolPotential()
{
  for(int i=0; i<InpCPP.size(); i++)
    if(InpCPP[i])
      delete InpCPP[i];
}


void LocalCorePolPotential::resetTargetParticleSet(ParticleSet& P)
{
  d_ie = DistanceTable::add(IonConfig,P,DT_AOS);
}

/** process xml node for each element
 * @param cur xmlnode <element name="string" alpha="double" rb="double"/>
 */
bool LocalCorePolPotential::CPP_Param::put(xmlNodePtr cur)
{
  OhmmsAttributeSet att;
  att.add(alpha,"alpha");
  att.add(r_b,"rb");
  att.put(cur);
  //const xmlChar* a_ptr = xmlGetProp(cur,(const xmlChar *)"alpha");
  //const xmlChar* b_ptr = xmlGetProp(cur,(const xmlChar *)"rb");
  //if(a_ptr) alpha = atof((const char*)a_ptr);
  //if(b_ptr) r_b = atof((const char*)b_ptr);
  C = -0.5*alpha;
  one_over_rr = 1.0/r_b/r_b;
  app_log() << "\talpha = " << alpha << " rb = " << r_b << std::endl;
  return true;
}

/** process xml node for CPP
 * @param cur xmlnode containing element+
 *
 * element/@name is used to find the index of the element of the
 * IonConfig::SpeciesSet. The size of InpCPP is the number of species.
 * The size of Centers is the number of ions.
 */
bool LocalCorePolPotential::put(xmlNodePtr cur)
{
  bool success(true);
  if(cur!= NULL)//input is provided
  {
    std::string ename;
    cur= cur->children;
    while(cur != NULL)
    {
      std::string cname((const char*)cur->name);
      if(cname == "element")
      {
        std::string species_name;
        OhmmsAttributeSet att;
        att.add(species_name,"name");
        att.put(cur);
        if(species_name.size())
        {
          int itype = IonConfig.getSpeciesSet().addSpecies(species_name); //(const char*)e_ptr);
          if(InpCPP[itype]==0)
            InpCPP[itype] = new CPP_Param;
          app_log() << "CPP parameters for " << IonConfig.getSpeciesSet().speciesName[itype] << std::endl;
          success = InpCPP[itype]->put(cur);
        }
      }
      cur=cur->next;
    }
  }
  for(int iat=0; iat<nCenters; iat++)
    Centers[iat]=InpCPP[IonConfig.GroupID[iat]];
  return success;
}

LocalCorePolPotential::Return_t
LocalCorePolPotential::evaluate(ParticleSet& P)
{
  if(FirstTime)
  {
    //index for attribute charge
    SpeciesSet& Species(IonConfig.getSpeciesSet());
    int iz = Species.addAttribute("charge");
    //calculate the Core-Core Dipole matrix
    for(int iat=0; iat<nCenters; iat++)
    {
      for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++)
      {
        int jat(d_ii->J[nn]);
        RealType rinv3 = std::pow(d_ii->rinv(nn),3);//(1/R_{JI}^3) R_{JI} = R_J-R_I
        PosType dipole(rinv3*d_ii->dr(nn));//(\vec{R_{JI}}/R_{JI}^3)
        //Sign and the charge of the paired ion are taken into account here
        CoreCoreDipole[iat] -= dipole*Species(iz,IonConfig.GroupID[jat]);
        CoreCoreDipole[jat] += dipole*Species(iz,IonConfig.GroupID[iat]);
      }
    }
    RealType corecore(0.0);
    for(int iat=0; iat<nCenters; iat++)
    {
      //app_log() << "Checking CPP = " << Centers[iat] << std::endl;
      if(Centers[iat])
        corecore+= Centers[iat]->C*dot(CoreCoreDipole[iat],CoreCoreDipole[iat]);
    }
//      LOGMSG("Core-Core Dipole = " << corecore);
    FirstTime=false;
  }
  //calculate the Electron-Core Dipole matrix
  //CoreElDipole=0.0;
  RealType e = 0.0;
  for(int iat=0; iat<nCenters; iat++)
  {
    if(Centers[iat])
    {
      PosType cc(CoreCoreDipole[iat]);
      for(int nn=d_ie->M[iat]; nn<d_ie->M[iat+1]; nn++)
      {
        int eid(d_ie->J[nn]);
        RealType rinv3 = std::pow(d_ie->rinv(nn),3);//(1/r^3)
        PosType dipole = rinv3*d_ie->dr(nn);//(\vec{r}/r^3)
        //cc +=  dipole*fcpp(d_ie->r(nn)*r_binv);
        cc += dipole*((*Centers[iat])(d_ie->r(nn)));
      }
      e += Centers[iat]->C*dot(cc,cc);
    }
  }
  return Value=e;
}

QMCHamiltonianBase* LocalCorePolPotential::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  LocalCorePolPotential* myclone=new LocalCorePolPotential(IonConfig,qp);
  //copy cpp parameters
  for(int i=0; i<InpCPP.size(); ++i)
    if(InpCPP[i])
      myclone->InpCPP[i]=new CPP_Param(*InpCPP[i]);
  myclone->put(NULL);
  return myclone;
}

}
