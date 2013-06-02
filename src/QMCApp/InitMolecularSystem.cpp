//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file InitMolecularSystem.cpp
 * @brief Implements InitMolecuarSystem operators.
 */
#include "QMCApp/InitMolecularSystem.h"
#include "QMCApp/ParticleSetPool.h"
using namespace std;
#include "OhmmsData/AttributeSet.h"
#include "Utilities/OhmmsInfo.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "ParticleBase/RandomSeqGenerator.h"

namespace qmcplusplus
{

InitMolecularSystem::InitMolecularSystem(ParticleSetPool* pset,
    const char* aname):
  OhmmsElementBase(aname), ptclPool(pset) { }

bool InitMolecularSystem::put(xmlNodePtr cur)
{
  string target("e"), source("i"), volume("no");
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

  app_log() << "<init source=\"" << source << "\" target=\""<< target << "\">" << endl;

  if(volume=="yes")
    initWithVolume(ions,els);
  else
    initMolecule(ions,els);

  app_log() << "</init>" << endl;
  app_log().flush();

  return true;
}

void InitMolecularSystem::initAtom(ParticleSet* ions, ParticleSet* els)
{
  //3N-dimensional Gaussian
  ParticleSet::ParticlePos_t chi(els->getTotalNum());
  makeGaussRandom(chi);
  double q=std::sqrt(static_cast<double>(els->getTotalNum()))*0.5;
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
  double BondLength;
  inline LoneElectron(int id, double bl):ID(id),BondLength(bl) {}
};

void InitMolecularSystem::initMolecule(ParticleSet* ions, ParticleSet* els)
{
  if(ions->getTotalNum()==1)
    return initAtom(ions,els);

  DistanceTableData* d_ii = DistanceTable::add(*ions);
  //d_ii->create(1);
  d_ii->evaluate(*ions);
  const ParticleSet::ParticleIndex_t& grID(ions->GroupID);
  SpeciesSet& Species(ions->getSpeciesSet());
  int Centers = ions->getTotalNum();
  vector<int> Qtot(Centers), Qcore(Centers), Qval(Centers,0);
  //use charge as the core electrons first
  int icharge = Species.addAttribute("charge");
  //Assign default core charge
  for(int iat=0; iat<Centers; iat++)
    Qtot[iat] = static_cast<int>(Species(icharge,grID[iat]));
  //cutoff radius (Bohr) this a random choice
  double cutoff = 4.0;
  ParticleSet::ParticlePos_t chi(els->getTotalNum());
  //makeGaussRandom(chi);
  makeSphereRandom(chi);
  int numUp =els->last(0);
  int numDown = els->last(1)-els->first(0);
  int item = 0;
  int nup_tot=0, ndown_tot=numUp;
  vector<LoneElectron> loneQ;
  double rmin=cutoff;
  ParticleSet::SingleParticlePos_t cm;
  for(int iat=0; iat<Centers; iat++)
  {
    cm += ions->R[iat];
    for(int nn=d_ii->M[iat]; nn<d_ii->M[iat+1]; nn++)
    {
      rmin = std::min(rmin,d_ii->r(nn));
    }
    //use 40% of the minimum bond
    double sep=rmin*0.4;
    int v2=Qtot[iat]/2;
    if(Qtot[iat]>v2*2)
    {
      loneQ.push_back(LoneElectron(iat,sep));
    }
    for(int k=0; k<v2; k++)
    {
      els->R[nup_tot++]=ions->R[iat]+sep*chi[item++];
      els->R[ndown_tot++]=ions->R[iat]+sep*chi[item++];
    }
  }
  vector<LoneElectron>::iterator it(loneQ.begin());
  vector<LoneElectron>::iterator it_end(loneQ.end());
  while(nup_tot<numUp && it != it_end)
  {
    els->R[nup_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
    ++it;
  }
  while(it != it_end)
  {
    els->R[ndown_tot++]=ions->R[(*it).ID]+(*it).BondLength*chi[item++];
    ++it;
  }
  //extra electrons around the geometric center
  double cnorm=1.0/static_cast<double>(Centers);
  cm=cnorm*cm;
  if(nup_tot<numUp)
  {
    double sep=rmin*2;
    int iu=0;
    while(nup_tot<numUp)
    {
      els->R[nup_tot++]=cm+sep*chi[iu++];
    }
  }
  if(ndown_tot<numDown)
  {
    double sep=rmin*2;
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
  TinyVector<double,OHMMS_DIM> start(1.0);
  TinyVector<double,OHMMS_DIM> end(0.0);

  ParticleSet::ParticlePos_t Ru(ions->getTotalNum());
  Ru.setUnit(PosUnit::LatticeUnit);
  ions->applyBC(ions->R,Ru);

  for(int iat=0; iat<Ru.size(); iat++)
  {
    start=lower_bound(Ru[iat],start);
    end=upper_bound(Ru[iat],end);
  }

  TinyVector<double,OHMMS_DIM> shift;
  Tensor<double,OHMMS_DIM> newbox(ions->Lattice.R);

  double buffer=2.0; //buffer 2 bohr
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
      double buffer_r=buffer*ions->Lattice.OneOverLength[idim];
      start[idim]=max(0.0,(start[idim]-buffer_r));
      end[idim]  =min(1.0,(end[idim]  +buffer_r));
      shift[idim]=start[idim]*ions->Lattice.Length[idim];
      if(std::abs(end[idim]=start[idim])<buffer)
      {//handle singular case
        start[idim]=max(0.0,start[idim]-buffer_r/2.0);
        end[idim]  =min(1.0,end[idim]  +buffer_r/2.0);
      }

      newbox(idim,idim)=(end[idim]-start[idim])*ions->Lattice.Length[idim];
    }
  }

  ParticleSet::ParticleLayout_t slattice(ions->Lattice);
  slattice.set(newbox);

  app_log() << "  InitMolecularSystem::initWithVolume " << endl;
  app_log() << "  Effective Lattice shifted by  " << shift << endl;
  app_log() <<newbox<< endl;

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
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
