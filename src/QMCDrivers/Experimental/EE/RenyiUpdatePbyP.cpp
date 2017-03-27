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
    
    


#include "QMCDrivers/EE/RenyiUpdatePbyP.h"
#include "QMCDrivers/DriftOperators.h"
#include "Message/OpenMP.h"

namespace qmcplusplus
{

RenyiUpdatePbyP::RenyiUpdatePbyP(MCWalkerConfiguration& w, TrialWaveFunction& psi,
                                 QMCHamiltonian& h, RandomGenerator_t& rg,int order):
  QMCRenyiUpdateBase(w,psi,h,rg,order)
{
}

RenyiUpdatePbyP::~RenyiUpdatePbyP()
{
}

void RenyiUpdatePbyP::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
{
  myTimers[0]->start();
  while (it != it_end)
  {
    WalkerIter_t begin_it(it);
    for (int i(0); i<2*RenyiOrder; i++)
    {
      Walker_t& thisWalker(**it);
      W_vec[i]->loadWalker(thisWalker,true);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      Psi_vec[i]->copyFromBuffer((*W_vec[i]),w_buffer);
      it++;
    }
//       this tells us what order to move particles in
    sort_regions_by_r();
    myTimers[1]->start();
    for (int iter=0; iter<nSubSteps; ++iter)
    {
      for (int iat=0; iat<W.getTotalNum(); ++iat)
      {
        int firstmove(-2);
        for (int i(0); i<RenyiOrder; i++)
          makeGaussRandomWithEngine(*deltaR_vec[i],RandomGen);
//           int starting0(get_region((*W_vec[0]).R,r_map(0,iat)));
        PosType dr = m_sqrttau* (*deltaR_vec[0])[r_map(0,iat)];
        if (W_vec[0]->makeMoveAndCheck(r_map(0,iat),dr))
          firstmove=get_region((*W_vec[0]).R,r_map(0,iat));
        bool goodmove(true);
//           only propose half the moves, all other walkers R depends on the first half
        for (int th(1); th<RenyiOrder; th++)
        {
//             //          all subsequent moves must cross or not cross as the first one
//             int startingth(get_region((*W_vec[th]).R,r_map(th,iat)));
//             if(startingth!=starting0)
//               std::cerr <<"Starting from different regions"<< std::endl;
          dr = m_sqrttau*(*deltaR_vec[th])[r_map(th,iat)];
          int nextmove(-1);
          if (W_vec[th]->makeMoveAndCheck(r_map(th,iat),dr))
            nextmove=get_region( (*W_vec[th]).R,r_map(th,iat));
          goodmove=((firstmove==nextmove) and goodmove);
        }
        if(goodmove)
        {
          for (int th(0); th<RenyiOrder; th++)
          {
            //          all subsequent moves must cross or not cross as the first one
//               psi(r_1)psi(r_2)...psi(r_0^a r_1^b)psi(r_1^a r_2^b)...)
            int indx=th+RenyiOrder+firstmove;
            indx=(indx==RenyiOrder*2?RenyiOrder:indx);
//               int startingth(get_region((*W_vec[indx]).R,r_map(indx,iat)));
//               if(startingth!=starting0)
//                 std::cerr <<"Starting Ths from different regions"<< std::endl;
            dr = (*W_vec[th]).R[r_map(th,iat)] - (*W_vec[indx]).R[r_map(indx,iat)];
            bool x=W_vec[indx]->makeMoveAndCheck(r_map(indx,iat),dr);
//               if (not x)
//                 std::cerr <<"rejected x move ?!"<< std::endl;
          }
          RealType rate(1);
          for (int th(0); th<2*RenyiOrder; th++)
            rate *= Psi_vec[th]->ratio((*W_vec[th]),r_map(th,iat));
          if (RandomGen() < std::abs(rate))
          {
//               for (int th(0);th<2*RenyiOrder;th++)
//                 put_in_box((*W_vec[th]).R[iat]);
            for (int th(0); th<2*RenyiOrder; th++)
              (*W_vec[th]).acceptMove(r_map(th,iat));
            for (int th(0); th<2*RenyiOrder; th++)
              (*Psi_vec[th]).acceptMove((*W_vec[th]),r_map(th,iat));
//               changed regions
            if (regions[0][r_map(0,iat)]!=firstmove)
            {
              regions[0][NumPtcl+firstmove] +=1;
              regions[0][NumPtcl+regions[0][r_map(0,iat)]] -=1;
              regions[0][r_map(0,iat)]=firstmove;
            }
            RealType sgn;
            sgn=(rate>0?+1:-1);
            regions[0][NumPtcl+2] *= sgn;
// //               sanity check
//               RealType sumphase(0);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase+= Psi_vec[th]->getPhase();
//               regions[NumPtcl+2] = std::cos(sumphase);
            takestats();
            ++nAccept;
          }
          else
          {
            for (int th(0); th<2*RenyiOrder; th++)
              (*W_vec[th]).rejectMove(r_map(th,iat));
            for (int th(0); th<2*RenyiOrder; th++)
              (*Psi_vec[th]).rejectMove(r_map(th,iat));
//               RealType sumphase(0);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase+= Psi_vec[th]->getPhase();
//               regions[NumPtcl+2] = std::cos(sumphase);
            takestats();
            ++nReject;
          }
        }
        else
        {
          for (int th(0); th<RenyiOrder; th++)
            (*W_vec[th]).rejectMove(r_map(th,iat));
//             RealType sumphase(0);
//             for (int th(0);th<2*RenyiOrder;th++)
//               sumphase+= Psi_vec[th]->getPhase();
//             regions[NumPtcl+2] = std::cos(sumphase);
          takestats();
          ++nReject;
        }
      }
    }
    myTimers[1]->stop();
//
    myTimers[2]->start();
    it=begin_it;
    for (int i(0); i<2*RenyiOrder; i++,it++)
    {
      Walker_t::Buffer_t& w_buffer((**it).DataSet);
      RealType logpsi = Psi_vec[i]->updateBuffer(*W_vec[i],w_buffer,false);
//         RealType p0=Psi_vec[i]->getPhase();
//         Psi_vec[i]->evaluateLog(*W_vec[i]);
//         RealType p1=Psi_vec[i]->getPhase();
//         Psum+=p1;
      W_vec[i]->saveWalker((**it));
      (**it).resetProperty(logpsi,Psi_vec[i]->getPhase(),0);
    }
//       double_check_region(begin_it,it_end);
    myTimers[2]->stop();
  }
  myTimers[0]->stop();
}


void RenyiUpdatePbyP::plotSwapAmplitude(WalkerIter_t it, WalkerIter_t it_end, Matrix<RealType>& averageSwaps)
{
  int ctr=(averageSwaps.size1())/2;
  int N_A(0);
  while (it != it_end)
  {
    for (int i(0); i<2*RenyiOrder; i++)
    {
      Walker_t& thisWalker(**it);
      W_vec[i]->loadWalker(thisWalker,true);
      Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
      Psi_vec[i]->copyFromBuffer((*W_vec[i]),w_buffer);
      it++;
    }
    sort_regions_by_r();
    std::vector<RealType> rate_i(RenyiOrder,0);
    for (int i(0); i<RenyiOrder; i++)
    {
      int indx=i+RenyiOrder+1;
      indx=(indx==RenyiOrder*2?RenyiOrder:indx);
      rate_i[i]=std::cos(Psi_vec[indx]->getPhase()+Psi_vec[i]->getPhase());
    }
    for (int th(0); th<RenyiOrder; th++)
      for (int iat=0; iat<W.getTotalNum(); ++iat)
      {
//           all in region 1
        if(get_region(W_vec[th]->R,r_map(th,iat))==1)
        {
          N_A++;
          for (int i(0); i<averageSwaps.size1(); i++)
            for (int j(0); j<averageSwaps.size2(); j++)
            {
              if((i-ctr)*(i-ctr)+(j-ctr)*(j-ctr)<=(ctr-1)*(ctr-1))
              {
                PosType X;
                X[0]=(i-ctr)*vsize+C[0][0];
                X[1]=(j-ctr)*vsize+C[0][1];
                PosType dr = X - (*W_vec[th]).R[r_map(th,iat)];
                W_vec[th]->makeMoveAndCheck(r_map(th,iat),dr);
                double r=Psi_vec[th]->ratio((*W_vec[th]),r_map(th,iat));
                int indx=th+RenyiOrder+1;
                indx=(indx==RenyiOrder*2?RenyiOrder:indx);
                dr = X - (*W_vec[indx]).R[r_map(indx,iat)];
                bool x=W_vec[indx]->makeMoveAndCheck(r_map(indx,iat),dr);
                r*=Psi_vec[indx]->ratio((*W_vec[indx]),r_map(indx,iat));
//                   averageSwaps(i,j)+=(rate_i[th]*r >0?1:-1);
                averageSwaps(i,j)+=regions[0][NumPtcl+2]*r;
//                   averageSwaps(i,j)+=(avgtmp*rate_i[th]>0?1:-1);
                (*W_vec[th]).rejectMove(r_map(th,iat));
                (*W_vec[indx]).rejectMove(r_map(indx,iat));
                (*Psi_vec[th]).rejectMove(r_map(th,iat));
                (*Psi_vec[indx]).rejectMove(r_map(indx,iat));
              }
            }
        }
      }
  }
  averageSwaps*=1.0/N_A;
}


//   RenyiUpdateAll::RenyiUpdateAll(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//       QMCHamiltonian& h, RandomGenerator_t& rg,int order):
//     QMCRenyiUpdateBase(w,psi,h,rg,order)
//   {
//   }
//
//   RenyiUpdateAll::~RenyiUpdateAll()
//   {
//   }
//
//   void RenyiUpdateAll::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     myTimers[0]->start();
//     while (it != it_end)
//     {
//       std::vector<Walker_t* > walkers(2*RenyiOrder);
//       for (int i(0);i<2*RenyiOrder;i++,it++)
//       {
//         walkers[i] = *it;
//         (*W_vec[i]).R=(**it).R;
//         (*W_vec[i]).update();
//       }
//
//       myTimers[1]->start();
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//           int firstmove(-2);
//           for (int i(0);i<RenyiOrder;i++)
//             makeGaussRandomWithEngine(*deltaR_vec[i],RandomGen);
//
//           if (W_vec[0]->makeMove(*walkers[0] ,*deltaR_vec[0], m_sqrttau))
//             firstmove=get_region_all((*W_vec[0]).R,0);
//
//           bool goodmove(true);
// //           only propose half the moves, all other walkers R depends on the first half
//           for (int th(1);th<RenyiOrder;th++)
//           {
//             int nextmove(-1);
//             if (W_vec[th]->makeMove(*walkers[th] ,*deltaR_vec[th], m_sqrttau))
//               nextmove=get_region_all( (*W_vec[th]).R,th);
//             goodmove=((firstmove==nextmove) and goodmove);
//           }
//
//           if(goodmove)
//           {
//             sort_regions_by_r();
//             for (int th(0);th<RenyiOrder;th++)
//             {
//               for( int iat(0);iat<NumPtcl;iat++)
//               {
//                 int indx=th+tmp_regions[th][iat];
//                 indx=(indx==RenyiOrder?0:indx);
//                 (*W_vec[indx+RenyiOrder]).R[iat]= (*W_vec[th]).R[iat];
//               }
//             }
//             for (int th(0);th<2*RenyiOrder;th++)
//               (*W_vec[th]).update();
//
//             RealType rate(0);
//             for (int th(0);th<2*RenyiOrder;th++)
//               rate += Psi_vec[th]->evaluateLog((*W_vec[th])) - walkers[th]->Properties(LOGPSI);
//             rate=std::exp(rate);
//             if (RandomGen() < std::abs(rate))
//             {
//               for (int th(0);th<2*RenyiOrder;th++)
//               {
//                 walkers[th]->R =  W_vec[th]->R;
//                 walkers[th]->resetProperty(Psi_vec[th]->getLogPsi(),Psi_vec[th]->getPhase(),0);
//               }
//               update_regions();
//
//               regions[0][NumPtcl+1]=firstmove;
//               regions[0][NumPtcl]=NumPtcl-firstmove;
//
//               RealType sumphase(1);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase *= std::cos(Psi_vec[th]->getPhase());
//               regions[0][NumPtcl+2] = sumphase;
//
//               takestats();++nAccept;
//             }
//             else
//             {
//               takestats();++nReject;
//             }
//           }
//           else
//           {
//             takestats();++nReject;
//           }
//       }
//       myTimers[1]->stop();
//     }
//     myTimers[0]->stop();
//   }
//
//   /// Constructor.
//   RenyiUpdateAllWithDrift::RenyiUpdateAllWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//       QMCHamiltonian& h, RandomGenerator_t& rg, int order): QMCRenyiUpdateBase(w,psi,h,rg,order)
//   {
//     drifts.resize(2*RenyiOrder);
//     for(int i(0);i<2*RenyiOrder;i++)
//       drifts[i]=new ParticleSet::ParticlePos_t(NumPtcl);
//   }
//
//   RenyiUpdateAllWithDrift::~RenyiUpdateAllWithDrift()
//   {
//   }
//
//   void RenyiUpdateAllWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     myTimers[0]->start();
//     while (it != it_end)
//     {
//       std::vector<Walker_t* > walkers(2*RenyiOrder);
//       for (int i(0);i<2*RenyiOrder;i++,it++)
//       {
//         walkers[i] = * it;
//         (*W_vec[i]).loadWalker((*walkers[i]),false);
//       }
//
//       myTimers[1]->start();
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//           for (int i(0);i<2*RenyiOrder;i++)
//             RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,(*W_vec[i]).G,*drifts[i]);
//
//           for (int th(0);th<RenyiOrder;th++)
//           {
//             for( int iat(0);iat<NumPtcl;iat++)
//             {
//               int indx=th+regions[th][iat]+RenyiOrder;
//               indx=(indx==2*RenyiOrder?RenyiOrder:indx);
//               (*drifts[th])[iat] += (*drifts[indx])[iat];
//               (*drifts[th])[iat] *= 0.5;
// //               (*drifts[indx])[iat] = (*drifts[th])[iat];
//             }
//           }
//
//           for (int i(0);i<RenyiOrder;i++)
//             makeGaussRandomWithEngine(*deltaR_vec[i],RandomGen);
//
//           int firstmove(-2);
//           if (W_vec[0]->makeMoveWithDrift( *walkers[0], *drifts[0], *deltaR_vec[0], m_sqrttau))
//             firstmove=get_region_all((*W_vec[0]).R,0);
//
//           bool goodmove(true);
// //           only propose half the moves, all other walkers R depends on the first half
//           for (int th(1);th<RenyiOrder;th++)
//           {
//             int nextmove(-1);
//             if (W_vec[th]->makeMoveWithDrift( *walkers[th], *drifts[th], *deltaR_vec[th], m_sqrttau))
//               nextmove=get_region_all( (*W_vec[th]).R,th);
//             goodmove=((firstmove==nextmove) and goodmove);
//           }
//
//
//           if(goodmove)
//           {
//             RealType logGf(0);
//             for (int th(0);th<RenyiOrder;th++)
//               logGf -= 0.5*Dot(*deltaR_vec[th],*deltaR_vec[th]);
//
//             for (int th(0);th<RenyiOrder;th++)
//               *deltaR_vec[th] = walkers[th]->R;
//
//             sort_regions_and_dr();
//             for (int th(0);th<RenyiOrder;th++)
//             {
//               for( int iat(0);iat<NumPtcl;iat++)
//               {
//                 int indx=th+tmp_regions[th][iat];
//                 indx=(indx==RenyiOrder?0:indx);
//                 (*W_vec[indx+RenyiOrder]).R [iat] = (*W_vec[th]).R[iat];
//               }
//             }
//             for (int th(0);th<2*RenyiOrder;th++)
//               (*W_vec[th]).update();
//
//             RealType rate(0);
//             for (int th(0);th<2*RenyiOrder;th++)
//               rate += Psi_vec[th]->evaluateLog((*W_vec[th])) - walkers[th]->Properties(LOGPSI);
//
//             for (int i(0);i<2*RenyiOrder;i++)
//               RealType nodecorr=setScaledDriftPbyPandNodeCorr(m_tauovermass,(*W_vec[i]).G,*drifts[i]);
//
//             for (int th(0);th<RenyiOrder;th++)
//             {
//               for( int iat(0);iat<NumPtcl;iat++)
//               {
//                 int indx=th+tmp_regions[th][iat]+RenyiOrder;
//                 indx=(indx==2*RenyiOrder?RenyiOrder:indx);
//                 (*drifts[th])[iat] += (*drifts[indx])[iat];
//
// //                 RealType d12=dot((*W_vec[th]).R[iat] - (*W_vec[indx]).R[iat],(*W_vec[th]).R[iat] - (*W_vec[indx]).R[iat]);
// //                 if( d12>0.001 )
// //                   std::cerr <<"broken"<< std::endl;
//
//                 (*drifts[th])[iat] *= 0.5;
// //                 (*drifts[indx])[iat] = (*drifts[th])[iat];
//               }
//             }
//
//             RealType logGb(0);
//             for (int th(0);th<RenyiOrder;th++)
//             {
//               (*deltaR_vec[RenyiOrder]) = *deltaR_vec[th] - (*W_vec[th]).R - *drifts[th];
//               logGb -= m_oneover2tau*Dot(*deltaR_vec[RenyiOrder],*deltaR_vec[RenyiOrder]);
//             }
//
//             rate=std::exp(rate + logGb-logGf);
//             if (RandomGen() < std::abs(rate))
//             {
//               for (int th(0);th<2*RenyiOrder;th++)
//               {
//                 W_vec[th]->saveWalker(*walkers[th]);
//                 walkers[th]->resetProperty(Psi_vec[th]->getLogPsi(),Psi_vec[th]->getPhase(),0);
//               }
//               update_regions();
//
//               regions[0][NumPtcl+1]=firstmove;
//               regions[0][NumPtcl]=NumPtcl-firstmove;
// //               update sign
//               RealType sumphase(1);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase *= std::cos(Psi_vec[th]->getPhase());
//               regions[0][NumPtcl+2] = sumphase;
//
//               takestats();++nAccept;
//             }
//             else
//             {
//               takestats();++nReject;
//             }
//           }
//           else
//           {
//             takestats();++nReject;
//           }
//       }
//       myTimers[1]->stop();
//     }
//     myTimers[0]->stop();
//   }
//
}

