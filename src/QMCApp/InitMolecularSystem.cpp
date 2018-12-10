//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jordan E. Vincent, University of Illinois at Urbana-Champaign
//                    Luke Shulenburger, lshulen@sandia.gov, Sandia National Laboratories
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Mark Dewing, markdewing@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/**@file InitMolecularSystem.cpp
 * @brief Implements InitMolecuarSystem operators.
 */
#include "QMCApp/InitMolecularSystem.h"
#include "QMCApp/ParticleSetPool.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{

typedef QMCTraits::RealType RealType;

InitMolecularSystem::InitMolecularSystem(ParticleSetPool* pset,
    const char* aname):
  OhmmsElementBase(aname), ptclPool(pset) { }

bool InitMolecularSystem::put(xmlNodePtr cur)
{
  std::string target("e"), source("i"), volume("no");
  OhmmsAttributeSet hAttrib;
  hAttrib.add(target,"target");
  hAttrib.add(source,"source");
  hAttrib.add(volume,"use_volume");
  hAttrib.put(cur);
  ParticleSet* els=ptclPool->getParticleSet(target);
  if(els == 0)
  {
    ERRORMSG("No target particle " << target << " exists.")
    return false;
  }
  ParticleSet* ions=ptclPool->getParticleSet(source);
  if(ions == 0)
  {
    ERRORMSG("No target particle " << target << " exists.")
    return false;
  }

  app_log() << "<init source=\"" << source << "\" target=\""<< target << "\">" << std::endl;

  if(volume=="yes")
    initWithVolume(ions,els);
  else
    initMolecule(ions,els);

  app_log() << "</init>" << std::endl;
  app_log().flush();

  return true;
}

void InitMolecularSystem::initAtom(ParticleSet* ions, ParticleSet* els)
{
  //3N-dimensional Gaussian
  ParticleSet::ParticlePos_t chi(els->getTotalNum());
  makeGaussRandom(chi);
  RealType q=std::sqrt(static_cast<RealType>(els->getTotalNum()))*0.5;
  int nel(els->getTotalNum()), items(0);
  while(nel)
  {
    els->R[items]=ions->R[0]+q*chi[items];
    --nel;
    ++items;
  }
}

struct LoneElectron
{
  int ID;
  RealType BondLength;
  inline LoneElectron(int id, RealType bl):ID(id),BondLength(bl) {}
};

void InitMolecularSystem::initMolecule(ParticleSet* ions, ParticleSet* els)
{
  if(ions->getTotalNum()==1)
    return initAtom(ions,els);

  DistanceTableData* d_ii = DistanceTable::add(*ions,DT_SOA_PREFERRED);
  //d_ii->create(1);
  d_ii->evaluate(*ions);
  const ParticleSet::ParticleIndex_t& grID(ions->GroupID);
  SpeciesSet& Species(ions->getSpeciesSet());
  int Centers = ions->getTotalNum();
  std::vector<int> Qtot(Centers), Qcore(Centers), Qval(Centers,0);
  //use charge as the core electrons first
  int icharge = Species.addAttribute("charge");
  //Assign default core charge
  for(int iat=0; iat<Centers; iat++)
    Qtot[iat] = static_cast<int>(Species(icharge,grID[iat]));
  //cutoff radius (Bohr) this a random choice
  RealType cutoff = 4.0;
  ParticleSet::ParticlePos_t chi(els->getTotalNum());
  //makeGaussRandom(chi);
  makeSphereRandom(chi);
  int numUp =els->last(0);
  int numDown = els->last(1)-els->first(0);
  int item = 0;
  int nup_tot=0, ndown_tot=numUp;
  std::vector<LoneElectron> loneQ;
  RealType rmin=cutoff;
  ParticleSet::SingleParticlePos_t cm;

  if(d_ii->DTType == DT_SOA)
  {
    for(size_t iat=0; iat<Centers; iat++)
    {
      cm += ions->R[iat];
      const RealType* restrict dist=d_ii->Distances[iat];
      for(size_t jat=iat+1; jat<Centers; ++jat)
      {
        rmin = std::min(rmin,dist[jat]);
      }
      //use 40% of the minimum bond
      RealType sep=rmin*0.4;
      int v2=Qtot[iat]/2;
      if(Qtot[iat]>v2*2)
      {
        loneQ.push_back(LoneElectron(iat,sep));
      }
      for(int k=0; k<v2; k++)
      {
        if(nup_tot<numUp)
          els->R[nup_tot++]=ions->R[iat]+sep*chi[item++];
        if(ndown_tot<numDown)
          els->R[ndown_tot++]=ions->R[iat]+sep*chi[item++];
      }
    }
  }
  else
  {
    for(int iat=0; iat<Centers; iat++)
    {
      cm += ions->R[iat];
      for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++)
      {
        rmin = std::min(rmin,d_ii->r(nn));
      }
      //use 40% of the minimum bond
      RealType sep=rmin*0.4;
      int v2=Qtot[iat]/2;
      if(Qtot[iat]>v2*2)
      {
        loneQ.push_back(LoneElectron(iat,sep));
      }
      for(int k=0; k<v2; k++)
      {
        if(nup_tot<numUp)
          els->R[nup_tot++]=ions->R[iat]+sep*chi[item++];
        if(ndown_tot<numDown)
          els->R[ndown_tot++]=ions->R[iat]+sep*chi[item++];
      }
    }
  }

  // mmorales: changed order of spin assignment to help with spin
  // imbalances in molecules at large distances. 
  // Not guaranteed to work, but should help in most cases  
  // as long as atoms in molecules are defined sequencially 
  std::vector<LoneElectron>::iterator it(loneQ.begin());
  std::vector<LoneElectron>::iterator it_end(loneQ.end());
  //while(nup_tot<numUp && it != it_end)
  //{
  //  els->R[nup_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
  //  ++it;
  //}
  //while(it != it_end)
  //{
  //  els->R[ndown_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
  //  ++it;
  //}
  while(it != it_end && nup_tot != numUp && ndown_tot != numDown)
  {
    if(nup_tot<numUp) {
      els->R[nup_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
      ++it;
    }
    if(ndown_tot<numDown && it != it_end) {
      els->R[ndown_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
      ++it;
    }
  }
  //extra electrons around the geometric center
  RealType cnorm=1.0/static_cast<RealType>(Centers);
  cm=cnorm*cm;
  if(nup_tot<numUp)
  {
    RealType sep=rmin*2;
    int iu=0;
    while(nup_tot<numUp)
    {
      els->R[nup_tot++]=cm+sep*chi[iu++];
    }
  }
  if(ndown_tot<numDown)
  {
    RealType sep=rmin*2;
    int iu=0;
    while(ndown_tot<numDown)
    {
      els->R[ndown_tot++]=cm+sep*chi[iu++];
    }
  }
  //put all the electrons in a unit box
  if(els->Lattice.SuperCellEnum != SUPERCELL_OPEN)
  {
    els->R.setUnit(PosUnit::CartesianUnit);
    els->applyBC(els->R);
    els->update(0);
  }
}

///helper function to determine the lower bound of a domain (need to move up)
template<typename T>
inline TinyVector<T,3> lower_bound(const TinyVector<T,3>& a, const TinyVector<T,3>& b)
{
  return TinyVector<T,3>(std::min(a[0],b[0]),std::min(a[1],b[1]),std::min(a[2],b[2]));
}

///helper function to determine the upper bound of a domain (need to move up)
template<typename T>
inline TinyVector<T,3> upper_bound(const TinyVector<T,3>& a, const TinyVector<T,3>& b)
{
  return TinyVector<T,3>(std::max(a[0],b[0]),std::max(a[1],b[1]),std::max(a[2],b[2]));
}

void InitMolecularSystem::initWithVolume(ParticleSet* ions, ParticleSet* els)
{
  TinyVector<RealType,OHMMS_DIM> start(1.0);
  TinyVector<RealType,OHMMS_DIM> end(0.0);

  ParticleSet::ParticlePos_t Ru(ions->getTotalNum());
  Ru.setUnit(PosUnit::LatticeUnit);
  ions->applyBC(ions->R,Ru);

  for(int iat=0; iat<Ru.size(); iat++)
  {
    start=lower_bound(Ru[iat],start);
    end=upper_bound(Ru[iat],end);
  }

  TinyVector<RealType,OHMMS_DIM> shift;
  Tensor<RealType,OHMMS_DIM> newbox(ions->Lattice.R);

  RealType buffer=2.0; //buffer 2 bohr
  for(int idim=0; idim<OHMMS_DIM; ++idim)
  {
    //if(ions->Lattice.BoxBConds[idim]) 
    //{
    //  start[idim]=0.0; 
    //  end[idim]=1.0;
    //  shift[idim]=0.0;
    //}
    //else
    {
      RealType buffer_r=buffer*ions->Lattice.OneOverLength[idim];
      start[idim]=std::max((RealType)0.0,(start[idim]-buffer_r));
      end[idim]  =std::min((RealType)1.0,(end[idim]  +buffer_r));
      shift[idim]=start[idim]*ions->Lattice.Length[idim];
      if(std::abs(end[idim]=start[idim])<buffer)
      {//handle singular case
        start[idim]= std::max(0.0,start[idim]-buffer_r/2.0);
        end[idim]  = std::min(1.0,end[idim]  +buffer_r/2.0);
      }

      newbox(idim,idim)=(end[idim]-start[idim])*ions->Lattice.Length[idim];
    }
  }

  ParticleSet::ParticleLayout_t slattice(ions->Lattice);
  slattice.set(newbox);

  app_log() << "  InitMolecularSystem::initWithVolume " << std::endl;
  app_log() << "  Effective Lattice shifted by  " << shift << std::endl;
  app_log() <<newbox<< std::endl;

  Ru.resize(els->getTotalNum());
  makeUniformRandom(Ru);
  for(int iat=0; iat<Ru.size(); ++iat)
    els->R[iat]=slattice.toCart(Ru[iat])+shift;
  els->R.setUnit(PosUnit::CartesianUnit);
}

bool InitMolecularSystem::put(std::istream& is)
{
  return true;
}

bool InitMolecularSystem::get(std::ostream& os) const
{
  return true;
}

void InitMolecularSystem::reset()
{
}
}//namespace
