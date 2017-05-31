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
    
    


#ifndef QMCPLUSPLUS_QMCRENYIUPDATEBASE_H
#define QMCPLUSPLUS_QMCRENYIUPDATEBASE_H
#include "QMCDrivers/QMCDriver.h"
#include "QMCDrivers/CloneManager.h"
namespace qmcplusplus
{
class QMCRenyiUpdateBase: public QMCUpdateBase
{
public:
  /// Constructor.
  QMCRenyiUpdateBase(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                     QMCHamiltonian& h, RandomGenerator_t& rg, int order=2);

  ~QMCRenyiUpdateBase();

  void initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end);
  void initWalkers(WalkerIter_t it, WalkerIter_t it_end);
  virtual void plotSwapAmplitude(WalkerIter_t it, WalkerIter_t it_end,Matrix<RealType>& averageSwaps) {};

  virtual void check_region(WalkerIter_t it, WalkerIter_t it_end, RealType v, std::string shape, ParticleSet::ParticlePos_t& ed, ParticleSet::ParticlePos_t& Center, int maxN, int minN, bool pbyp=true);
  virtual int get_region(ParticleSet::ParticlePos_t& Pos,int iat);
  virtual int get_region_all(ParticleSet::ParticlePos_t& Pos, int th=0);

  void put_in_box(PosType& Pos);

  void takestats()
  {
    regions[0][NumPtcl+3] += regions[0][NumPtcl+2];
    if((regions[0][NumPtcl+1]<mxN)and(regions[0][NumPtcl+1]>=mnN))
      n_region[regions[0][NumPtcl+1]-mnN]+=1;
    cnt++;
  }

  RealType get_stats(std::vector<RealType>& n)
  {
    RealType nrm=1.0/cnt;
    for (int i(0); i<n.size(); i++)
      n[i]+=nrm*n_region[i];
    return nrm*regions[0][NumPtcl+3];
  }

  void clear_stats()
  {
    std::fill(n_region.begin(),n_region.end(),0);
    cnt=0;
    regions[0][NumPtcl+3]=0;
  }

  void update_regions()
  {
    regions=tmp_regions;
  }

  int sort_regions(bool initial=false)
  {
    if (initial)
      for(int i(0); i<RenyiOrder; i++)
        get_region_all(W_vec[i]->R,i);
    std::vector<int> miat(RenyiOrder,0);
    int swaps(0), all_swaps(0);
//       first include all region 0
    for(int i(0); i<RenyiOrder; i++,swaps=0)
      for(int j(0); j<NumPtcl; j++)
        if(tmp_regions(i,j)==0)
        {
//             use other walkers for tmp storage
          W_vec[i+RenyiOrder]->R[miat[i]++] = W_vec[i]->R[j];
          W_vec[i+RenyiOrder]->G[miat[i]]   = W_vec[i]->G[j];
          W_vec[i+RenyiOrder]->L[miat[i]]   = W_vec[i]->L[j];
          all_swaps+=swaps;
        }
        else
          swaps++;
    for(int i(0); i<RenyiOrder; i++)
      for(int j(0); j<NumPtcl; j++)
        if(tmp_regions(i,j)==1)
        {
          W_vec[i+RenyiOrder]->R[miat[i]++] = W_vec[i]->R[j];
          W_vec[i+RenyiOrder]->G[miat[i]]   = W_vec[i]->G[j];
          W_vec[i+RenyiOrder]->L[miat[i]]   = W_vec[i]->L[j];
        }
// sorted into regions now
    for(int i(0); i<RenyiOrder; i++)
    {
      W_vec[i]->R = W_vec[i+RenyiOrder]->R;
      W_vec[i]->G = W_vec[i+RenyiOrder]->G;
      W_vec[i]->L = W_vec[i+RenyiOrder]->L;
    }
// update the tmp_regions lists
    for(int i(0); i<RenyiOrder; i++)
      get_region_all(W_vec[i]->R,i);
    return (all_swaps&1)?-1:1;
    //return std::pow(-1,all_swaps);
  }

  int sort_regions_by_r();

  int sort_regions_and_dr()
  {
    std::vector<int> miat(RenyiOrder,0);
    int swaps(0), all_swaps(0);
//       first include all region 0
    for(int i(0); i<RenyiOrder; i++,swaps=0)
      for(int j(0); j<NumPtcl; j++)
        if(tmp_regions(i,j)==0)
        {
//             use other walkers for tmp storage
          W_vec[i+RenyiOrder]->R[miat[i]] = W_vec[i]->R[j];
          W_vec[i+RenyiOrder]->G[miat[i]]   = W_vec[i]->G[j];
          W_vec[i+RenyiOrder]->L[miat[i]]   = W_vec[i]->L[j];
          (*deltaR_vec[i+RenyiOrder])[miat[i]++]=(*deltaR_vec[i])[j];
          all_swaps+=swaps;
        }
        else
          swaps++;
    for(int i(0); i<RenyiOrder; i++)
      for(int j(0); j<NumPtcl; j++)
        if(tmp_regions(i,j)==1)
        {
          W_vec[i+RenyiOrder]->R[miat[i]] = W_vec[i]->R[j];
          W_vec[i+RenyiOrder]->G[miat[i]]   = W_vec[i]->G[j];
          W_vec[i+RenyiOrder]->L[miat[i]]   = W_vec[i]->L[j];
          (*deltaR_vec[i+RenyiOrder])[miat[i]++]=(*deltaR_vec[i])[j];
        }
// sorted into regions now
    for(int i(0); i<RenyiOrder; i++)
    {
      W_vec[i]->R = W_vec[i+RenyiOrder]->R;
      W_vec[i]->G = W_vec[i+RenyiOrder]->G;
      W_vec[i]->L = W_vec[i+RenyiOrder]->L;
      *deltaR_vec[i] = *deltaR_vec[i+RenyiOrder];
    }
// update the tmp_regions lists
    for(int i(0); i<RenyiOrder; i++)
      get_region_all(W_vec[i]->R,i);
    return (all_swaps&1)?-1:1;
    //return std::pow(-1,all_swaps);
  }


  void print_all();


//     RealType advanceWalkerForEE(Walker_t& w1, std::vector<PosType>& dR, std::vector<int>& iats, std::vector<int>& rs, std::vector<RealType>& ratios);

  std::vector<NewTimer*> myTimers;
  std::vector<MCWalkerConfiguration*> W_vec;
  std::vector<TrialWaveFunction*> Psi_vec;
  std::vector<ParticleSet::ParticlePos_t*> deltaR_vec;


  int RenyiOrder;
  int mxN,mnN;
  int cnt;
  ///patch geometry (sphere,...
  std::string computeEE;
  ///patch size (sphere:r^2,...
  RealType vsize;
  ParticleSet::ParticlePos_t C,Edge;
  std::vector<RealType> n_region;
  Matrix<int> regions, tmp_regions;
  Matrix<int> r_map;

  std::vector<PosType> plottingGrid;
};

}

#endif
