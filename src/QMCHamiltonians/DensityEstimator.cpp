//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    John R. Gergely,  University of Illinois at Urbana-Champaign
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCHamiltonians/DensityEstimator.h>
#include <OhmmsData/AttributeSet.h>
#include "LongRange/LRCoulombSingleton.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus
{

typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
typedef LRCoulombSingleton::GridType       GridType;
typedef LRCoulombSingleton::RadFunctorType RadFunctorType;

DensityEstimator::DensityEstimator(ParticleSet& elns) : rVs(0)
{
  UpdateMode.set(COLLECTABLE,1);
  Periodic=(elns.Lattice.SuperCellEnum != SUPERCELL_OPEN);
  for(int dim=0; dim<OHMMS_DIM; ++dim)
  {
    density_max[dim]=elns.Lattice.Length[dim];
    ScaleFactor[dim]=1.0/elns.Lattice.Length[dim];
  }
  //    InitPotential(elns);
}

void DensityEstimator::resetTargetParticleSet(ParticleSet& P)
{
}

DensityEstimator::Return_t DensityEstimator::evaluate(ParticleSet& P)
{
  RealType wgt=tWalker->Weight;
  if (Periodic)
  {
    for(int iat=0; iat<P.getTotalNum(); ++iat)
    {
      PosType ru;
      ru=P.Lattice.toUnit(P.R[iat]);
      int i=static_cast<int>(DeltaInv[0]*(ru[0]-std::floor(ru[0])));
      int j=static_cast<int>(DeltaInv[1]*(ru[1]-std::floor(ru[1])));
      int k=static_cast<int>(DeltaInv[2]*(ru[2]-std::floor(ru[2])));
      P.Collectables[getGridIndex(i,j,k)]+=wgt;//1.0;
      //	P.Collectables[getGridIndexPotential(i,j,k)]-=1.0;
      //HACK!	P.Collectables[getGridIndexPotential(i,j,k)]+=evalSR(P,iat)+evalLR(P,iat);
    }
  }
  else
  {
    for(int iat=0; iat<P.getTotalNum(); ++iat)
    {
      PosType ru;
      for (int dim=0; dim<OHMMS_DIM; dim++)
      {
        //ru[dim]=(P.R[iat][dim]-density_min[dim])/(density_max[dim]-density_min[dim]);
        ru[dim]=(P.R[iat][dim]-density_min[dim])*ScaleFactor[dim];
      }
      if (ru[0]>0.0 && ru[1]>0.0 && ru[2]>0.0 &&
          ru[0]<1.0 && ru[1]<1.0 && ru[2]<1.0)
      {
        int i=static_cast<int>(DeltaInv[0]*(ru[0]-std::floor(ru[0])));
        int j=static_cast<int>(DeltaInv[1]*(ru[1]-std::floor(ru[1])));
        int k=static_cast<int>(DeltaInv[2]*(ru[2]-std::floor(ru[2])));
        P.Collectables[getGridIndex(i,j,k)]+=wgt;//1.0;
        //	  P.Collectables[getGridIndexPotential(i,j,k)]-=1.0;
        //HACK!	  P.Collectables[getGridIndexPotential(i,j,k)]+=evalSR(P,iat)+evalLR(P,iat);
      }
    }
  }
  return 0.0;
}

void
DensityEstimator::addEnergy(MCWalkerConfiguration &W,
                            std::vector<RealType> &LocalEnergy)
{
  int nw = W.WalkerList.size();
  int N = W.getTotalNum();
  if (Periodic)
  {
    for (int iw=0; iw<nw; iw++)
    {
      Walker_t &w = *W.WalkerList[iw];
      RealType weight=w.Weight/nw;
      for (int iat=0; iat<N; iat++)
      {
        PosType ru;
        ru=W.Lattice.toUnit(w.R[iat]);
        // for (int dim=0; dim<OHMMS_DIM; dim++)
        // {
        //   ru[dim]=(w.R[iat][dim]-density_min[dim])*ScaleFactor[dim];
        // }
        int i=static_cast<int>(DeltaInv[0]*(ru[0]-std::floor(ru[0])));
        int j=static_cast<int>(DeltaInv[1]*(ru[1]-std::floor(ru[1])));
        int k=static_cast<int>(DeltaInv[2]*(ru[2]-std::floor(ru[2])));
        W.Collectables[getGridIndex(i,j,k)]+=weight;
      }
    }
  }
}

DensityEstimator::RealType
DensityEstimator::evalLR(ParticleSet &P,int iat)
{
  RealType LR=0.0;
  const StructFact& PtclRhoK(*(P.SK));
  //    for(int spec1=0; spec1<NumSpecies; spec1++) {
  //      RealType Z1 = Zspec[spec1];
  for(int spec2=0; spec2<NumSpecies; spec2++)
  {
    RealType Z2 = Zspec[spec2];
    //RealType temp=AA->evaluate(PtclRhoK.KLists.minusk, PtclRhoK.rhok[spec1], PtclRhoK.rhok[spec2]);
#if defined(USE_REAL_STRUCT_FACTOR)
    RealType temp=AA->evaluate(PtclRhoK.KLists.kshell, iat, PtclRhoK.rhok_r[spec2], PtclRhoK.rhok_i[spec2],P);
#else
    RealType temp=AA->evaluate(PtclRhoK.KLists.kshell, iat, PtclRhoK.rhok[spec2],P);
#endif
    int spec1=spec2; ///BUG: REALLY NEED TO FIGURE OUT HOW TO GET SPEC OF IAT!
    RealType Z1 = Zspec[spec1];
    if(spec2==spec1)
      LR += 0.5*Z1*Z2*temp;
    else
      LR += Z1*Z2*temp;
  } //spec2
  //    }//spec1
  //LR*=0.5;
  return LR;
}

void DensityEstimator::InitPotential(ParticleSet &P)
{
  SpeciesSet& tspecies(P.getSpeciesSet());
  int ChargeAttribIndx = tspecies.addAttribute("charge");
  int MemberAttribIndx = tspecies.addAttribute("membersize");
  NumCenters = P.getTotalNum();
  NumSpecies = tspecies.TotalNum;
  Zat.resize(NumCenters);
  Zspec.resize(NumSpecies);
  for(int spec=0; spec<NumSpecies; spec++)
  {
    Zspec[spec] = tspecies(ChargeAttribIndx,spec);
    //      NofSpecies[spec] = static_cast<int>(tspecies(MemberAttribIndx,spec));
  }
  for(int iat=0; iat<NumCenters; iat++)
  {
    //      SpeciesID[iat]=P.GroupID[iat];
    Zat[iat] = Zspec[P.GroupID[iat]];
  }
  GridType* myGrid(0);
  RealType myRcut;
  AA = LRCoulombSingleton::getHandler(P);
  //AA->initBreakup(*PtclRef);
  myRcut=AA->get_rc();
  if(rVs==0)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA,myRcut,myGrid);
  }
}


DensityEstimator::RealType
DensityEstimator::evalSR(ParticleSet& P,int ipart)
{
  const DistanceTableData *d_aa = P.DistTables[0];
  RealType SR=0.0;
  //     for(int ipart=0; ipart<NumCenters; ipart++)
  {
    RealType esum = 0.0;
    for(int nn=d_aa->M[ipart],jpart=ipart+1; nn<d_aa->M[ipart+1]; nn++,jpart++)
    {
      esum += Zat[jpart]*d_aa->rinv(nn)*rVs->splint(d_aa->r(nn));
    }
    //Accumulate pair sums...species charge for atom i.
    SR += Zat[ipart]*esum;
  }
  return SR;
}



void DensityEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
{
  //current index
  myIndex=collectables.current();
  std::vector<RealType> tmp(NumGrids[OHMMS_DIM]);
  collectables.add(tmp.begin(),tmp.end());
  //potentialIndex=collectables.current();
  //vector<RealType> tmp2(NumGrids[OHMMS_DIM]);
  //collectables.add(tmp2.begin(),tmp2.end());
}

void DensityEstimator::registerCollectables(std::vector<observable_helper*>& h5desc
    , hid_t gid) const
{
  int loc=h5desc.size();
  std::vector<int> ng(OHMMS_DIM);
  for(int i=0; i<OHMMS_DIM; ++i)
    ng[i]=NumGrids[i];
  observable_helper* h5o=new observable_helper(myName);
  h5o->set_dimensions(ng,myIndex);
  h5o->open(gid);
  h5desc.push_back(h5o);
  //h5o=new observable_helper("Potential");
  //h5o->set_dimensions(ng,potentialIndex);
  //h5o->open(gid);
  //h5desc.push_back(h5o);
}

void DensityEstimator::setObservables(PropertySetType& plist)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex);
}

void DensityEstimator::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  //std::copy(density.first_address(),density.last_address(),plist.begin()+myDebugIndex+offset);
}

/** check xml elements
 *
 * <estimator name="density" debug="no" delta="0.1 0.1 0.1"/>
 */
bool DensityEstimator::put(xmlNodePtr cur)
{
  Delta=0.1;
  std::vector<double> delta;
  std::string debug("no");
  std::string potential("no");
  OhmmsAttributeSet attrib;
  attrib.add(debug,"debug");
  attrib.add(potential,"potential");
  attrib.add(density_min[0],"x_min");
  attrib.add(density_min[1],"y_min");
  attrib.add(density_min[2],"z_min");
  attrib.add(density_max[0],"x_max");
  attrib.add(density_max[1],"y_max");
  attrib.add(density_max[2],"z_max");
  attrib.add(Delta,"delta");
  attrib.put(cur);
  if(!Periodic)
  {
    for(int dim=0; dim<OHMMS_DIM; ++dim)
      ScaleFactor[dim]=1.0/(density_max[dim]-density_min[dim]);
  }
  resize();
  return true;
}

bool DensityEstimator::get(std::ostream& os) const
{
  os << myName << " bin =" << Delta << " bohrs " << std::endl;
  return true;
}

QMCHamiltonianBase* DensityEstimator::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  //default constructor is sufficient
  return new DensityEstimator(*this);
}

void  DensityEstimator::resize()
{
  for(int i=0; i<OHMMS_DIM; ++i)
  {
    DeltaInv[i]=1.0/Delta[i];
    NumGrids[i]=static_cast<int>(DeltaInv[i]);
    if(NumGrids[i]<2)
    {
      APP_ABORT("DensityEstimator::resize invalid bin size");
    }
  }
  app_log() << " DensityEstimator bin_size= " <<NumGrids << " delta = " << Delta << std::endl;
  NumGrids[OHMMS_DIM]=NumGrids[0]*NumGrids[1]*NumGrids[2];
}

}

