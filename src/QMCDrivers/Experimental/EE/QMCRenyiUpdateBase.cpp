//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
#include "QMCDrivers/EE/QMCRenyiUpdateBase.h"
#include "QMCDrivers/DriftOperators.h"


namespace qmcplusplus
{
void add_ee_timers(std::vector<NewTimer*>& timers)
{
  timers.push_back(new NewTimer("QMCRenyiUpdateBase::advance")); //timer for the walker loop
  timers.push_back(new NewTimer("QMCRenyiUpdateBase::movePbyP")); //timer for MC, ratio etc
  timers.push_back(new NewTimer("QMCRenyiUpdateBase::updateMBO")); //timer for measurements
  timers.push_back(new NewTimer("QMCRenyiUpdateBase::energy")); //timer for measurements
  for (int i=0; i<timers.size(); ++i)
    TimerManager.addTimer(timers[i]);
}

QMCRenyiUpdateBase::QMCRenyiUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                       QMCHamiltonian& h, RandomGenerator_t& rg,int order):
  QMCUpdateBase(w,psi,h,rg)
{
  add_ee_timers(myTimers);
  RenyiOrder=order;
  W_vec.resize(order*2),
               W_vec[0]=&w;
  W_vec[0]->UseBoundBox=false;
  Psi_vec.resize(order*2),
                 Psi_vec[0]=&psi;
  for(int i(1); i<RenyiOrder*2; i++)
  {
    W_vec[i]=new MCWalkerConfiguration(w);
    W_vec[i]->UseBoundBox=false;
    Psi_vec[i]=psi.makeClone(*W_vec[i]);
  }
  NumPtcl = W.getTotalNum();
  deltaR_vec.resize(2*RenyiOrder);
  deltaR_vec[0]=&deltaR;
  for(int i(1); i<2*RenyiOrder; i++)
    deltaR_vec[i]= new ParticleSet::ParticlePos_t(NumPtcl);
  r_map.resize(2*RenyiOrder,NumPtcl);
}

QMCRenyiUpdateBase::~QMCRenyiUpdateBase()
{
}

void QMCRenyiUpdateBase::initWalkers(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=false;
  sort_regions(true);
  update_regions();
  for (; it != it_end; it++)
  {
    W_vec[0]->R = (*it)->R;
    W_vec[0]->update();
    RealType logpsi(Psi_vec[0]->evaluateLog((*W_vec[0])));
    //setScaledDriftPbyP(Tau*m_oneovermass,W.G,(*it)->Drift);
    RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,(*W_vec[0]).G,drift);
//         RealType ene = H.evaluate(W);
    (*it)->resetProperty(logpsi,Psi_vec[0]->getPhase(),0,0.0,0.0, nodecorr);
    (*it)->Weight=1;
    for(int i(1); i<RenyiOrder*2; i++)
    {
      it++;
      (*it)->R = W_vec[0]->R;
      (*it)->G = W_vec[0]->G;
      (*it)->L = W_vec[0]->L;
      (*it)->resetProperty(logpsi,Psi_vec[0]->getPhase(),0,0.0,0.0, nodecorr);
      (*it)->Weight=1;
    }
  }
}


void QMCRenyiUpdateBase::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
{
  UpdatePbyP=true;
  sort_regions_by_r();
  for (; it != it_end; it++)
  {
    Walker_t& awalker(**it);
    (*W_vec[0]).R=awalker.R;
    (*W_vec[0]).update();
    //W.loadWalker(awalker,UpdatePbyP);
    if (awalker.DataSet.size())
      awalker.DataSet.clear();
    awalker.DataSet.rewind();
    RealType logpsi=Psi_vec[0]->registerData((*W_vec[0]),awalker.DataSet);
    RealType logpsi2=Psi_vec[0]->updateBuffer((*W_vec[0]),awalker.DataSet,false);
    W_vec[0]->saveWalker(awalker);
    for(int i(1); i<RenyiOrder*2; i++)
    {
      it++;
      Walker_t& bwalker(**it);
//           bwalker.makeCopy(awalker);
      (*W_vec[i]).R=awalker.R;
      (*W_vec[i]).update();
      if (bwalker.DataSet.size())
        bwalker.DataSet.clear();
      bwalker.DataSet.rewind();
      logpsi=Psi_vec[i]->registerData((*W_vec[i]),bwalker.DataSet);
      logpsi2=Psi_vec[i]->updateBuffer((*W_vec[i]),bwalker.DataSet,false);
      W_vec[i]->saveWalker(bwalker);
    }
  }
}

void QMCRenyiUpdateBase::check_region(WalkerIter_t it, WalkerIter_t it_end, RealType v, std::string shape, ParticleSet::ParticlePos_t& ed, ParticleSet::ParticlePos_t& Center, int maxN, int minN, bool pbyp)
{
  computeEE=shape;
  vsize=v;
  C.resize(Center.size());
  for (int x=0; x<DIM; x++)
    C[0][x]=Center[0][x];
  Edge.resize(ed.size());
  for (int x=0; x<DIM; x++)
    Edge[0][x]=ed[0][x];
  Walker_t& bwalker(**it);
  int dim1(1);
  if (!pbyp )
    dim1=2*RenyiOrder;
  regions.resize(dim1,NumPtcl+4);
  tmp_regions.resize(dim1,NumPtcl+4);
  regions=0;
  for(int d(0); d<dim1; d++)
    regions(d,NumPtcl+2)=1;
  mxN=maxN+1;
  mnN=minN;
  n_region.resize(maxN-minN,0);
//     centering
  while (it!=it_end)
  {
    Walker_t& awalker(**it);
    awalker.R=bwalker.R;
    for (int g=0; g<W.groups(); ++g)
      for (int iat=W.first(g); iat<W.last(g); ++iat)
      {
        for(int x(0); x<DIM; x++)
        {
          while (awalker.R[iat][x]<0)
            awalker.R[iat][x]+=Edge[0][x];
          while (awalker.R[iat][x]>Edge[0][x])
            awalker.R[iat][x]-=Edge[0][x];
        }
      }
    it++;
  }
  for (int g=0; g<W.groups(); ++g)
    for (int iat=W.first(g); iat<W.last(g); ++iat)
    {
      //several volumes we can integrate out
      RealType z;
      if(computeEE=="square")
      {
        z=(bwalker.R[iat][0]);
        for (int x=0; x<DIM; x++)
          z=std::max(z,bwalker.R[iat][x]);
      }
      else
        if(computeEE=="sphere")
        {
          z=0;
          for (int x=0; x<DIM; x++)
            z+=std::pow(bwalker.R[iat][x]-C[0][x],2);
          z=std::sqrt(z);
        }
        else
          if(computeEE=="halfspace")
          {
            z=(bwalker.R[iat][DIM-1]);
          }
      if (z<vsize)
        regions[0][iat]=1;
      else
        regions[0][iat]=0;
      regions[0][NumPtcl+regions[0][iat]]+=1;
    }
  if((regions[0][NumPtcl+1]<mxN)and(regions[0][NumPtcl+1]>=mnN))
    n_region[regions[0][NumPtcl+1]-mnN]+=1;
  for(int d(1); d<dim1; d++)
    for(int e(0); e<NumPtcl+4; e++)
      regions(d,e)=regions(0,e);
}

int QMCRenyiUpdateBase::get_region_all(ParticleSet::ParticlePos_t& Pos, int th)
{
  int region_alpha(0);
  for (int g=0; g<W.groups(); ++g)
    for (int iat=W.first(g); iat<W.last(g); ++iat)
    {
      //several volumes we can integrate out
      RealType z;
      if(computeEE=="square")
      {
        z=(Pos[iat][0]);
        for (int x=0; x<DIM; x++)
          z=std::max(z,Pos[iat][x]);
      }
      else
        if(computeEE=="sphere")
        {
          z=0;
          for (int x=0; x<DIM; x++)
            z+=std::pow(Pos[iat][x]-C[0][x],2);
          z=std::sqrt(z);
        }
        else
          if(computeEE=="halfspace")
          {
            z=(Pos[iat][DIM-1]);
          }
      if (z<vsize)
      {
        tmp_regions[th][iat]=1;
        region_alpha+=1;
      }
      else
        tmp_regions[th][iat]=0;
    }
  return region_alpha;
}

void QMCRenyiUpdateBase::put_in_box(PosType& Pos)
{
  for(int x(0); x<DIM; x++)
  {
    while (Pos[x]<0)
      Pos[x]+=Edge[0][x];
    while (Pos[x]>Edge[0][x])
      Pos[x]-=Edge[0][x];
  }
}

int QMCRenyiUpdateBase::get_region(ParticleSet::ParticlePos_t& Pos,int iat)
{
  PosType tmpP(Pos[iat]);
  put_in_box(tmpP);
  RealType z;
  if(computeEE=="square")
  {
    z=(Pos[iat][0]);
    for (int x=0; x<DIM; x++)
      z=std::max(z,tmpP[x]);
  }
  else
    if(computeEE=="sphere")
    {
      z=0;
      for (int x=0; x<DIM; x++)
        z+=std::pow(tmpP[x]-C[0][x],2);
      z=std::sqrt(z);
    }
    else
      if(computeEE=="halfspace")
      {
        z=(tmpP[DIM-1]);
      }
  if (z<vsize)
    return 1;
  else
    return 0;
}

int QMCRenyiUpdateBase::sort_regions_by_r( )
{
  for(int i(0); i<2*RenyiOrder; i++)
  {
    std::vector<std::pair<RealType,int> > w_i_r(NumPtcl);
    ParticleSet::ParticlePos_t Pos(W_vec[i]->R);
    for (int g=0; g<W.groups(); ++g)
      for (int iat=W.first(g); iat<W.last(g); ++iat)
      {
        PosType tmpP(Pos[iat]);
        put_in_box(tmpP);
        //several volumes we can integrate out
        RealType z;
        if(computeEE=="square")
        {
          z=(tmpP[0]);
          for (int x=0; x<DIM; x++)
            z=std::max(z,tmpP[x]);
        }
        else
          if(computeEE=="sphere")
          {
            z=0;
            for (int x=0; x<DIM; x++)
              z+=std::pow(tmpP[x]-C[0][x],2);
            z=std::sqrt(z);
          }
          else
            if(computeEE=="halfspace")
            {
              z=(tmpP[DIM-1]);
            }
        w_i_r[iat].first=z;
        w_i_r[iat].second=iat;
      }
    std::sort(w_i_r.begin(),w_i_r.end());
    for (int iat(0); iat<NumPtcl; iat++)
      r_map(i,iat)=w_i_r[iat].second;
//       for (int iat(0);iat<NumPtcl;iat++)
//         std::cerr <<r_map(i,iat)<<" ";
//       std::cerr << std::endl;
  }
  return 0;
}

void QMCRenyiUpdateBase::print_all()
{
  for(int x(0); x<NumPtcl; x++)
  {
    for(int i(0); i<RenyiOrder*2; i++)
      std::cerr <<W_vec[i]->R[x]<<" ";
    for(int i(0); i<regions.size1(); i++)
      std::cerr <<regions[i][x]<<" ";
    std::cerr << std::endl;
  }
  std::cerr <<regions[0][NumPtcl+2]<< std::endl;
  std::cerr <<regions[0][NumPtcl+3]<< std::endl;
}

}

