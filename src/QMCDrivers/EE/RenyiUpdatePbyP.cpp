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
      for (int i(0);i<2*RenyiOrder;i++)
      {
        Walker_t& thisWalker(**it);
        W_vec[i]->loadWalker(thisWalker,true);
        Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
        Psi_vec[i]->copyFromBuffer((*W_vec[i]),w_buffer);
        it++;
      }
      
      myTimers[1]->start();
      for (int iter=0; iter<nSubSteps; ++iter)
      {
        for (int iat=0; iat<W.getTotalNum(); ++iat)
        {
//           cerr<<regions[NumPtcl+2]<<endl;
          int firstmove(-2);
          for (int i(0);i<RenyiOrder;i++)
            makeGaussRandomWithEngine(*deltaR_vec[i],RandomGen);
          
// //           propose first move. it may cross regions
//           int initialregion=get_region((*W_vec[0]).R,iat);
//           if(initialregion!=regions[iat])
//             cerr<<"wrong region"<<endl;
          PosType dr = m_sqrttau* (*deltaR_vec[0])[iat];
          if (W_vec[0]->makeMoveAndCheck(iat,dr))
            firstmove=get_region((*W_vec[0]).R,iat);
          
          bool goodmove(true);
//           only propose half the moves, all other walkers R depends on the first half
          for (int th(1);th<RenyiOrder;th++)
          {
//             int initialregion2=get_region((*W_vec[th]).R,iat);
//             if(initialregion2!=regions[iat])
//               cerr<<"WHOA! wrong region"<<endl;
//             //          all subsequent moves must cross or not cross as the first one
            dr = m_sqrttau*(*deltaR_vec[th])[iat];
            int nextmove(-1);
            if (W_vec[th]->makeMoveAndCheck(iat,dr))
              nextmove=get_region( (*W_vec[th]).R,iat);
            goodmove=((firstmove==nextmove) and goodmove);
          }
          
          if(goodmove)
          {
            for (int th(0);th<RenyiOrder;th++)
            {
  //          all subsequent moves must cross or not cross as the first one
//               psi(r_1)psi(r_2)...psi(r_0^a r_1^b)psi(r_1^a r_2^b)...)

              int indx=th+RenyiOrder+firstmove;
              indx=(indx==RenyiOrder*2?RenyiOrder:indx);
              dr = (*W_vec[th]).R[iat] - (*W_vec[indx]).R[iat];
              bool x=W_vec[indx]->makeMoveAndCheck(iat,dr);
//               if (not x)
//                 cerr<<"rejected x move ?!"<<endl;
            }
            
            RealType rate(1);
//             std::vector<RealType> ratios(2*RenyiOrder,0);
            for (int th(0);th<2*RenyiOrder;th++)
              rate *= Psi_vec[th]->ratio((*W_vec[th]),iat);
//             for (int th(0);th<2*RenyiOrder;th++) 
//               cerr<<ratios[th]<<" "; 
//             cerr<<rate<<endl;
            
            if (RandomGen() < std::abs(rate))
            {
//               for (int th(0);th<2*RenyiOrder;th++)
//                 put_in_box((*W_vec[th]).R[iat]);
              
              for (int th(0);th<2*RenyiOrder;th++)
                (*W_vec[th]).acceptMove(iat);
              for (int th(0);th<2*RenyiOrder;th++)
                (*Psi_vec[th]).acceptMove((*W_vec[th]),iat);
              
//               changed regions
              if (regions[iat]!=firstmove)
              {
                regions[NumPtcl+firstmove] +=1;
                regions[NumPtcl+regions[iat]] -=1;
                regions[iat]=firstmove;
              }
              
              
              RealType sgn; sgn=(rate>0?+1:-1);
              regions[NumPtcl+2] *= sgn;
// //               sanity check
//               RealType sumphase(0);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase+= Psi_vec[th]->getPhase();
//               regions[NumPtcl+2] = std::cos(sumphase);
              
              regions[NumPtcl+3] += regions[NumPtcl+2];
              if((regions[NumPtcl+1]<mxN)and(regions[NumPtcl+1]>=mnN))
                n_region[regions[NumPtcl+1]-mnN]+=1;
              cnt++; ++nAccept;

            }
            else
            {
              for (int th(0);th<2*RenyiOrder;th++)
                (*W_vec[th]).rejectMove(iat);
              for (int th(0);th<2*RenyiOrder;th++)
                (*Psi_vec[th]).rejectMove(iat);
              
//               RealType sumphase(0);
//               for (int th(0);th<2*RenyiOrder;th++)
//                 sumphase+= Psi_vec[th]->getPhase();
//               regions[NumPtcl+2] = std::cos(sumphase);
              
              regions[NumPtcl+3] += regions[NumPtcl+2];
              if((regions[NumPtcl+1]<mxN)and(regions[NumPtcl+1]>=mnN))
                n_region[regions[NumPtcl+1]-mnN]+=1;
              cnt++;++nReject;
            }
          }
          else
          {
            for (int th(0);th<RenyiOrder;th++)
              (*W_vec[th]).rejectMove(iat);
            
//             RealType sumphase(0);
//             for (int th(0);th<2*RenyiOrder;th++)
//               sumphase+= Psi_vec[th]->getPhase();
//             regions[NumPtcl+2] = std::cos(sumphase);
            
            regions[NumPtcl+3] += regions[NumPtcl+2];
            if((regions[NumPtcl+1]<mxN)and(regions[NumPtcl+1]>=mnN))
              n_region[regions[NumPtcl+1]-mnN]+=1;
            cnt++;++nReject;
          }
        }
      }
      myTimers[1]->stop();
// 
      myTimers[2]->start();

      it=begin_it;
//       RealType Psum(0);
      for (int i(0);i<2*RenyiOrder;i++,it++)
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

//   RenyiUpdatePbyP::RealType RenyiUpdatePbyP::advanceWalkerForEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios) 
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     
//     vector<RealType> runningratio(3,1.0); int sh(0);
//     int nshells(ratios.size()/3);
//     vector<RealType> orb_ratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
//       runningratio[0]*=ratio;
//       runningratio[1]*=orb_ratios[0];
//       runningratio[2]*=ratio/orb_ratios[0];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//       {
//         ratios[sh]=runningratio[0];
//         ratios[nshells+sh]=runningratio[1];
//         ratios[nshells*2+sh++]=runningratio[2];
//       }
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
//       runningratio[0]*=ratio;
//       runningratio[1]*=orb_ratios[0];
//       runningratio[2]*=ratio/orb_ratios[0];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(nshells*2+sh < ratios.size())      
//       {
//         ratios[sh]=runningratio[0];
//         ratios[nshells+sh]=runningratio[1];
//         ratios[nshells*2+sh++]=runningratio[2];
//       }
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//     
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//        
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//     
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//       
//     return runningratio[0];
//   }


//   VMCUpdatePbyP::RealType VMCUpdatePbyP::advanceWalkerForCSEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios, vector<RealType>& weights, vector<RealType>& logs) 
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     Psi.getLogs(logs);
//     
//     RealType runningratio(std::exp(2.0*logs[1]));
//     int sh(0);
//     vector<RealType> lzratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=lzratios[1];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//         ratios[sh++]=runningratio;
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=lzratios[1];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(sh<ratios.size())
//         ratios[sh++]=runningratio;
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//     
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//        
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//     
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//       
//     return runningratio;
//   }
// 
//   void VMCUpdatePbyP::advanceCSWalkers(vector<TrialWaveFunction*>& pclone
//       , vector<MCWalkerConfiguration*>& wclone
//       , vector<QMCHamiltonian*>& hclone
//       , vector<RandomGenerator_t*>& rng, vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
// 
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ipx=0; ipx<NumThreads; ++ipx)
//       psi2_i_new[ipx] = 2.0*W[ipx]->getPropertyBase()[LOGPSI] + c_i[ipx];
// 
// #pragma omp parallel 
//     {
//       RealType psi2_new=1.0;
//       for (int ipx=1; ipx<NumThreads; ++ipx)
//         psi2_new += std::exp(psi2_i_new[ipx]-psi2_i_new[0]);
// 
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
//       
//       //copy the new to now
//       RealType psi2_now=psi2_new;
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
// 
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
// 
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           PosType dr=m_sqrttau*deltaR[iat];
// 
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
// #pragma omp barrier
//           if (!movePtcl) continue;
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// // #pragma omp barrier
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
//           psi2_new=1.0;
//           for(int jp=1; jp<NumThreads; ++jp)
//             psi2_new+=std::exp(psi2_i_new[jp]-psi2_i_new[0]);
//           RealType prob = std::exp(psi2_i_new[0]-psi2_0_now)*(psi2_new/psi2_now);  
//           if (prob > rng_loc())
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
// 
//             psi2_i_now=psi2_i_new[ip];
//             psi2_0_now=psi2_i_new[0];
//             psi2_now=psi2_new;
//           }
//           else
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }//loop over iat
//       }
// 
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       wclone[ip]->saveWalker(*W[ip]);
//       RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//       thisWalker.Weight=1.0;
//       hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
// 
//     myTimers[0]->stop();
//   }
// 
//   void VMCUpdatePbyP::estimateNormWalkers(vector<TrialWaveFunction*>& pclone
//       , vector<MCWalkerConfiguration*>& wclone
//       , vector<QMCHamiltonian*>& hclone
//       , vector<RandomGenerator_t*>& rng
//       , vector<RealType>& ratio_i_0)
//   {
//     int NumThreads(pclone.size());
// 
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
// 
//     RealType nn = -std::log(1.0*nSubSteps*W.getTotalNum());
// #pragma omp parallel 
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
// 
//       //copy the new to now
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
// 
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
// 
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           PosType dr=m_sqrttau*deltaR[iat];
// 
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
// 
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*logl(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
// #pragma omp critical
//           {
//             ratio_i_0[ip] += expl( psi2_i_new[ip]-psi2_i_new[0] + nn);
//           }
//           wclone[ip]->rejectMove(iat);
//           pclone[ip]->rejectMove(iat);
//         }
//       }
// 
//       //     Walker_t& thisWalker(*W[ip]);
//       //     Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       //     RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       //     wclone[ip]->saveWalker(*W[ip]);
//       //     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //     //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       //     thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
// 
//       //     hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       //     hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
// 
//     myTimers[0]->stop();
//   }
// 
//   /// Constructor.
//   VMCUpdatePbyPWithDrift::VMCUpdatePbyPWithDrift(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//       QMCHamiltonian& h, RandomGenerator_t& rg):
//     QMCUpdateBase(w,psi,h,rg)
//   {
//     add_vmc_timers(myTimers);
//   }
// 
//   VMCUpdatePbyPWithDrift::~VMCUpdatePbyPWithDrift()
//   {
//   }
// 
//   void VMCUpdatePbyPWithDrift::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     myTimers[0]->start();
//     for (; it != it_end; ++it)
//     {
//       Walker_t& thisWalker(**it);
// 
//       W.loadWalker(thisWalker,true);
// 
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       Psi.copyFromBuffer(W,thisWalker.DataSet);
//       //create a 3N-Dimensional Gaussian with variance=1
//       bool moved = false;
//       myTimers[1]->start();
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
// 
//         makeGaussRandomWithEngine(deltaR,RandomGen);
// 
//         moved = false;
// 
//         for (int iat=0; iat<W.getTotalNum(); ++iat)
//         {
//           PosType dr;
//           ///dr = m_sqrttau*deltaR[iat]+thisWalker.Drift[iat];
//           //RealType sc=getDriftScale(m_tauovermass,W.G[iat]);
//           //PosType dr(m_sqrttau*deltaR[iat]+sc*real(W.G[iat]));
//           getScaledDrift(m_tauovermass,W.G[iat],dr);
//           dr += m_sqrttau*deltaR[iat];
// 
//           //reject illegal moves
//           if (!W.makeMoveAndCheck(iat,dr))
//           {
//             ++nReject;
//             continue;
//           }
//           //PosType newpos=W.R[iat];
//           //PosType newpos = W.makeMove(iat,dr);
// 
//           RealType ratio = Psi.ratio(W,iat,dG,dL);
//           RealType prob = ratio*ratio;
// 
//           //zero is always rejected
//           if (prob<numeric_limits<RealType>::epsilon())
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             continue;
//           }
// 
//           G = W.G+dG;
// 
//           //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
//           //dr = (*it)->R[iat]-newpos-Tau*G[iat];
//           //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
//           RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
// 
//           //RealType scale=getDriftScale(m_tauovermass,G[iat]);
//           //dr = thisWalker.R[iat]-W.R[iat]-scale*real(G[iat]);
//           getScaledDrift(m_tauovermass,G[iat],dr);
//           dr = thisWalker.R[iat]-W.R[iat]-dr;
// 
//           RealType logGb = -m_oneover2tau*dot(dr,dr);
// 
//           //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
//           if (RandomGen() < prob*std::exp(logGb-logGf))
//           {
//             moved = true;
//             ++nAccept;
//             W.acceptMove(iat);
//             Psi.acceptMove(W,iat);
//             W.G = G;
//             W.L += dL;
// 
//             //do not need to update Drift
//             //assignDrift(scale,G,thisWalker.Drift);
// 
//           }
//           else
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//           }
//         }
//       }
//       myTimers[1]->stop();
// 
//       if (moved)
//       {
//         myTimers[2]->start();
//         //thisWalker.R = W.R;
//         //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
//         //w_buffer.rewind();
//         //W.copyToBuffer(w_buffer);
// 
//         RealType logpsi = Psi.evaluateLog(W,w_buffer);
//         W.saveWalker(thisWalker);
// 
//         myTimers[2]->stop();
// 
//         myTimers[3]->start();
//         RealType eloc=H.evaluate(W);
//         myTimers[3]->stop();
//         //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
//         thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
//         H.auxHevaluate(W,thisWalker);
//         H.saveProperty(thisWalker.getPropertyBase());
//       }
//       else
//       {
//         ++nAllRejected;
//         H.rejectedMove(W,thisWalker);
//       }
//     }
//     myTimers[0]->stop();
//   }
// 
//   /// Constructor.
//   VMCUpdatePbyPWithDriftFast::VMCUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//       QMCHamiltonian& h, RandomGenerator_t& rg) :
//     QMCUpdateBase(w,psi,h,rg)
//   {
//     add_vmc_timers(myTimers);
//   }
// 
//   VMCUpdatePbyPWithDriftFast::~VMCUpdatePbyPWithDriftFast()
//   {
//   }
// 
//   VMCUpdatePbyPWithDriftFast::RealType VMCUpdatePbyPWithDriftFast::advanceWalkerForEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios) 
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     std::vector<RealType> logs;
//     Psi.getLogs(logs);
//     
//     
//     vector<RealType> runningratio(4,1.0); int sh(0);
// //     runningratio[3]=std::exp(2.0*(logs[0]-csoffset));
//     int nshells(rs.size());
//     
//     vector<RealType> orb_ratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
//       runningratio[0]*=ratio;
//       runningratio[1]*=orb_ratios[0];
//       runningratio[2]*=orb_ratios[1];
//       runningratio[3]*=ratio/orb_ratios[0];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//       {
//         ratios[sh]=runningratio[0];
//         ratios[nshells+sh]=runningratio[1];
//         ratios[nshells*2+sh]=runningratio[2];
//         ratios[nshells*3+sh++]=runningratio[3];
//       }
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,orb_ratios);
//       runningratio[0]*=ratio;
//       runningratio[1]*=orb_ratios[0];
//       runningratio[2]*=orb_ratios[1];
//       runningratio[3]*=ratio/orb_ratios[0];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(nshells*3+sh < ratios.size())      
//       {
//         ratios[sh]=runningratio[0];
//         ratios[nshells+sh]=runningratio[1];
//         ratios[nshells*2+sh]=runningratio[2];
//         ratios[nshells*3+sh++]=runningratio[3];
//       }
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//     
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//        
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//     
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//       
//     return std::exp(2.0*(logs[0]-csoffset));
//   }
// 
//   VMCUpdatePbyPWithDriftFast::RealType VMCUpdatePbyPWithDriftFast::advanceWalkerForCSEE(Walker_t& w1, vector<PosType>& dR, vector<int>& iats, vector<int>& rs, vector<RealType>& ratios, vector<RealType>& weights, vector<RealType>& logs ) 
//   {
//     Walker_t::Buffer_t& w_buffer(w1.DataSet);
//     W.loadWalker(w1,true);
//     Psi.copyFromBuffer(W,w_buffer);
//     Psi.getLogs(logs);
//     RealType logpsi = Psi.getLogPsi();
//     RealType logpsiCS = -(logpsi-logs[0]-csoffset);
//     
//     RealType runningratio(1);
//     RealType runningweight(1);
// //     RealType runningweight(std::exp(2.0*logpsiCS));
//     logs[0]=runningweight;
//     int sh(0);
//     vector<RealType> lzratios;
// //                 accumulate ratios on the way there
//     for(int itz(0); itz<(iats.size()-1); itz++)
//     {
//       int iat=iats[itz];
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=ratio;
//       runningweight*=lzratios[0];
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//       while(itz+1==rs[sh])
//       {
//         ratios[sh]=runningratio;
//         weights[sh++]=runningweight;
//       }
//     }
//     //we get to reject the last move
//     {
//       int iat=iats[iats.size()-1];
//       W.makeMove(iat,dR[iats.size()-1]);
//       RealType ratio = Psi.ratioVector(W,iat,lzratios);
//       runningratio*=ratio;
//       runningweight*=lzratios[0];
//       W.rejectMove(iat);
//       Psi.rejectMove(iat);
//       while(sh<ratios.size())      
//       {
//         ratios[sh]=runningratio;
//         weights[sh++]=runningweight;
//       }
//     }
// //                 put the particles back
//     for(int itz(0); itz<iats.size(); itz++)
//       dR[itz]*=-1;
//     
//     for(int itz(iats.size()-2); itz>=0; itz--)
//     {
//       int iat=iats[itz];
//        
//       W.makeMove(iat,dR[itz]);
//       RealType ratio = Psi.ratio(W,iat);
//       W.acceptMove(iat);
//       Psi.acceptMove(W,iat);
//     }
//     
// //     Psi.updateBuffer(W,w_buffer,false);
// //     W.saveWalker(thisWalker);
//       
//     return runningratio;
//   }
// 
//   void VMCUpdatePbyPWithDriftFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     myTimers[0]->start();
//     for (; it != it_end; ++it)
//     {
//       Walker_t& thisWalker(**it);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
// 
//       W.loadWalker(thisWalker,true);
//       //W.R = thisWalker.R;
//       //w_buffer.rewind();
//       //W.copyFromBuffer(w_buffer);
//       Psi.copyFromBuffer(W,w_buffer);
// 
// 
//       myTimers[1]->start();
//       bool moved = false;
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//         moved = false;
//         for (int iat=0; iat<W.getTotalNum(); ++iat)
//         {
// 
//           GradType grad_now=Psi.evalGrad(W,iat), grad_new;
//           PosType dr;
//           //    //RealType sc=getDriftScale(m_tauovermass,grad_now);
//           //    //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_now));
//           getScaledDrift(m_tauovermass,grad_now,dr);
//           dr += m_sqrttau*deltaR[iat];
//           //    PosType dr(m_sqrttau*deltaR[iat]+m_tauovermass*grad_now);
//           if (!W.makeMoveAndCheck(iat,dr))
//           {
//             ++nReject;
//             continue;
//           }
// 
//           //PosType newpos = W.makeMove(iat,dr);
//           RealType ratio = Psi.ratioGrad(W,iat,grad_new);
//           RealType prob = ratio*ratio;
// 
//           //zero is always rejected
//           if (prob<numeric_limits<RealType>::epsilon())
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             continue;
//           }
// 
//           //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
//           //dr = (*it)->R[iat]-newpos-Tau*G[iat];
//           //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
//           RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
// 
//           //    //sc=getDriftScale(m_tauovermass,grad_new);
//           //    //dr = thisWalker.R[iat]-W.R[iat]-sc*real(grad_new);
//           getScaledDrift(m_tauovermass,grad_new,dr);
//           dr = thisWalker.R[iat]-W.R[iat]-dr;
//           //    dr = thisWalker.R[iat]-W.R[iat]-m_tauovermass*grad_new;
//           RealType logGb = -m_oneover2tau*dot(dr,dr);
// 
//           //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
//           if (RandomGen() < prob*std::exp(logGb-logGf))
//           {
//             moved = true;
//             ++nAccept;
//             W.acceptMove(iat);
//             Psi.acceptMove(W,iat);
//           }
//           else
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//           }
//         }
//         //for subSteps must update thiswalker
//         thisWalker.R=W.R;
//         thisWalker.G=W.G;
//         thisWalker.L=W.L;
//       }
//       myTimers[1]->stop();
// 
//       //Always compute the energy
//       //if(moved)
//       {
//         myTimers[2]->start();
//         //thisWalker.R = W.R;
//         //w_buffer.rewind();
//         //W.updateBuffer(w_buffer);
//         RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
//         W.saveWalker(thisWalker);
//         myTimers[2]->stop();
// 
//         myTimers[3]->start();
//         RealType eloc=H.evaluate(W);
//         myTimers[3]->stop();
//         //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
//         thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
//         H.auxHevaluate(W,thisWalker);
//         H.saveProperty(thisWalker.getPropertyBase());
//       }
// 
//       if(!moved) ++nAllRejected;
//       //else
//       //{
//       //  ++nAllRejected;
//       //  H.rejectedMove(W,thisWalker);
//       //}
// 
//     }
//     myTimers[0]->stop();
//   }
// 
// 
//   /// Constructor.
//   VMCUpdateRenyiWithDriftFast::VMCUpdateRenyiWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//       QMCHamiltonian& h, RandomGenerator_t& rg) :
//     QMCUpdateBase(w,psi,h,rg)
//   {
//     add_vmc_timers(myTimers);
//   }
// 
//   VMCUpdateRenyiWithDriftFast::~VMCUpdateRenyiWithDriftFast()
//   {
//   }
// 
// 
//   void VMCUpdateRenyiWithDriftFast::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
//     myTimers[0]->start();
//     WalkerIter_t begin(it);
//     
//     for (; it != it_end; ++it)
//     {
//       Walker_t& thisWalker(**it);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
// 
//       W.loadWalker(thisWalker,true);
//       Psi.copyFromBuffer(W,w_buffer);
// 
// 
//       myTimers[1]->start();
//       bool moved = false;
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//         moved = false;
//         for (int iat=0; iat<W.getTotalNum(); ++iat)
//         {
// 
//           GradType grad_now=Psi.evalGrad(W,iat), grad_new;
//           PosType dr;
//           //    //RealType sc=getDriftScale(m_tauovermass,grad_now);
//           //    //PosType dr(m_sqrttau*deltaR[iat]+sc*real(grad_now));
//           getScaledDrift(m_tauovermass,grad_now,dr);
//           dr += m_sqrttau*deltaR[iat];
//           //    PosType dr(m_sqrttau*deltaR[iat]+m_tauovermass*grad_now);
//           if (!W.makeMoveAndCheck(iat,dr))
//           {
//             ++nReject;
//             continue;
//           }
// 
//           //PosType newpos = W.makeMove(iat,dr);
//           RealType ratio = Psi.ratioGrad(W,iat,grad_new);
//           RealType prob = ratio*ratio;
// 
//           //zero is always rejected
//           if (prob<numeric_limits<RealType>::epsilon())
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             continue;
//           }
// 
//           //RealType forwardGF = exp(-0.5*dot(deltaR[iat],deltaR[iat]));
//           //dr = (*it)->R[iat]-newpos-Tau*G[iat];
//           //RealType backwardGF = exp(-oneover2tau*dot(dr,dr));
//           RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
// 
//           //    //sc=getDriftScale(m_tauovermass,grad_new);
//           //    //dr = thisWalker.R[iat]-W.R[iat]-sc*real(grad_new);
//           getScaledDrift(m_tauovermass,grad_new,dr);
//           dr = thisWalker.R[iat]-W.R[iat]-dr;
//           //    dr = thisWalker.R[iat]-W.R[iat]-m_tauovermass*grad_new;
//           RealType logGb = -m_oneover2tau*dot(dr,dr);
// 
//           //RealType prob = std::min(1.0e0,ratio*ratio*std::exp(logGb-logGf));
//           if (RandomGen() < prob*std::exp(logGb-logGf))
//           {
//             moved = true;
//             ++nAccept;
//             W.acceptMove(iat);
//             Psi.acceptMove(W,iat);
//           }
//           else
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//           }
//         }
//         //for subSteps must update thiswalker
//         thisWalker.R=W.R;
//         thisWalker.G=W.G;
//         thisWalker.L=W.L;
//       }
//       myTimers[1]->stop();
// 
//       //Always compute the energy
//       //if(moved)
//       {
//         myTimers[2]->start();
//         //thisWalker.R = W.R;
//         //w_buffer.rewind();
//         //W.updateBuffer(w_buffer);
//         RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
//         W.saveWalker(thisWalker);
//         myTimers[2]->stop();
// 
//         myTimers[3]->start();
//         RealType eloc=H.evaluate(W);
//         myTimers[3]->stop();
//         //thisWalker.resetProperty(std::log(abs(psi)), psi,eloc);
//         thisWalker.resetProperty(logpsi,Psi.getPhase(), eloc);
//         H.auxHevaluate(W,thisWalker);
//         H.saveProperty(thisWalker.getPropertyBase());
//       }
// 
//       if(!moved) ++nAllRejected;
//       //else
//       //{
//       //  ++nAllRejected;
//       //  H.rejectedMove(W,thisWalker);
//       //}
// 
//     }
//     myTimers[0]->stop();
//   }
// 
// 
//   /// Constructor.
//   // VMCCSUpdatePbyPWithDriftFast::VMCCSUpdatePbyPWithDriftFast(MCWalkerConfiguration& w, TrialWaveFunction& psi,
//   //     QMCHamiltonian& h, RandomGenerator_t& rg) :
//   //     QMCUpdateBase(w,psi,h,rg)
//   // {
//   //   add_vmc_timers(myTimers);
//   // }
//   // 
//   // VMCCSUpdatePbyPWithDriftFast::~VMCCSUpdatePbyPWithDriftFast()
//   // {
//   // }
// 
//   void VMCUpdatePbyPWithDriftFast::advanceCSWalkers(vector<TrialWaveFunction*>& pclone, vector<MCWalkerConfiguration*>& wclone, vector<QMCHamiltonian*>& hclone, vector<RandomGenerator_t*>& rng, vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
//     myTimers[0]->start();
//     myTimers[1]->start();
//     bool moved = false;
// 
//     std::vector<RealType> psi2_i_now(NumThreads);
//     RealType psi2_now=0;
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       psi2_i_now[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
//       psi2_now += std::exp(psi2_i_now[ip]-psi2_i_now[0]);
//     }
// 
//     for (int iter=0; iter<nSubSteps; ++iter)
//     {
//       //create a 3N-Dimensional Gaussian with variance=1
//       makeGaussRandomWithEngine(deltaR, RandomGen);
//       moved = false;
// 
//       for (int iat=0; iat<W.getTotalNum(); ++iat)
//       {
//         std::vector<GradType> grad_i_now(NumThreads);
// 
//         // #pragma omp parallel for  
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_i_now[ip] = pclone[ip]->evalGrad(*wclone[ip],iat);
// 
// 
//         GradType grad_now;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_now += std::exp(psi2_i_now[ip]-psi2_i_now[0])/psi2_now*grad_i_now[ip];
// 
//         PosType dr;
//         getScaledDrift(m_tauovermass,grad_now,dr);
//         dr += m_sqrttau*deltaR[iat];
// 
// 
//         for (int ip=1; ip<NumThreads; ++ip)
//           wclone[ip]->makeMoveAndCheck(iat,dr);
// 
// 
//         if (!wclone[0]->makeMoveAndCheck(iat,dr))
//         {
//           ++nReject;
//           continue;
//         }
// 
//         RealType psi2_new(0);
//         std::vector<GradType> grad_i_new(NumThreads);
//         std::vector<RealType> psi2_i_new(NumThreads);
//         // #pragma omp parallel for  
//         for (int ip=0; ip<NumThreads; ++ip)
//         {
//           RealType ratio = pclone[ip]->ratioGrad(*wclone[ip],iat,grad_i_new[ip]);
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now[ip];
//         }
//         for (int ip=0; ip<NumThreads; ++ip) psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
// 
//         RealType prob = std::exp(psi2_i_new[0]-psi2_i_now[0])*(psi2_new/psi2_now);
// 
//         //zero is always rejected
//         if (prob < numeric_limits<RealType>::epsilon())
//         {
//           ++nReject;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//           continue;
//         }
// 
//         RealType logGf = -0.5e0*dot(deltaR[iat],deltaR[iat]);
// 
//         GradType grad_new;
//         for (int ip=0; ip<NumThreads; ++ip)
//           grad_new += std::exp(psi2_i_new[ip]-psi2_i_new[0])/psi2_new*grad_i_new[ip];
// 
//         getScaledDrift(m_tauovermass,grad_new,dr);
//         dr += W[0]->R[iat] - wclone[0]->R[iat] -dr;
// 
//         RealType logGb = -m_oneover2tau*dot(dr,dr);
// 
//         if ( RandomGen() < prob*std::exp(logGb-logGf))
//         {
//           moved = true;
//           ++nAccept;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
//             psi2_i_now[ip]=psi2_i_new[ip];
//           }
//           psi2_now = psi2_new;
//         } 
//         else
//         {
//           ++nReject;
//           for (int ip=0; ip<NumThreads; ++ip)
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }
//       }
//       //for subSteps must update walkers
//       for (int ip=0; ip<NumThreads; ++ip) wclone[ip]->saveWalker(*W[ip]);
//     }
// 
//     myTimers[1]->stop();
// 
//     // #pragma omp parallel for  
//     for (int ip=0; ip<NumThreads; ++ip)
//     {
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       if (moved)
//       {
//         RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//         wclone[ip]->saveWalker(*W[ip]);
//         RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//         //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//         thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//         hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//         hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//       }
//       else
//       {
//         ++nAllRejected;
//         hclone[ip]->rejectedMove(*wclone[ip],thisWalker);
//       }
//     }
// 
//     myTimers[0]->stop();
//   }
// 
//   VMCUpdatePbyPSampleRN::VMCUpdatePbyPSampleRN(MCWalkerConfiguration& w, TrialWaveFunction& psi, TrialWaveFunction& guide,
//       QMCHamiltonian& h, RandomGenerator_t& rg):
//     QMCUpdateBase(w,psi,guide,h,rg), logEpsilon(0.0)
//   {
//     add_vmc_timers(myTimers);
//     //   dummy= new MCWalkerConfiguration(w);
//   }
// 
//   VMCUpdatePbyPSampleRN::~VMCUpdatePbyPSampleRN()
//   {
//   }
// 
//   void VMCUpdatePbyPSampleRN::initWalkersForPbyP(WalkerIter_t it, WalkerIter_t it_end)
//   {
//     UpdatePbyP=true;
//     //   Guide.resizeTempP(W);
// 
//     for (;it != it_end; ++it)
//     {
//       Walker_t& thisWalker(**it);
//       W.loadWalker(thisWalker,UpdatePbyP);
// 
//       Walker_t::Buffer_t tbuffer;
//       RealType logguide=Guide.registerData(W,tbuffer)+logEpsilon;
//       RealType logpsi=Psi.registerData(W,tbuffer);
//       thisWalker.DataSet=tbuffer;
//       RealType ene = H.evaluate(W);
// 
//       thisWalker.resetProperty(logpsi,Psi.getPhase(),ene);
//       H.saveProperty(thisWalker.getPropertyBase());
//       thisWalker.ReleasedNodeAge=0;
//       thisWalker.ReleasedNodeWeight=0;
//       //       thisWalker.Weight=1.0;
//       thisWalker.Weight=std::exp(2.0*(logpsi-logguide));
//     }
//   }
// 
//   void VMCUpdatePbyPSampleRN::advanceWalkers(WalkerIter_t it, WalkerIter_t it_end, bool measure)
//   {
// 
//     myTimers[0]->start();
//     for (; it != it_end; ++it)
//     {
// 
//       Walker_t& thisWalker(**it);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
// 
//       W.loadWalker(thisWalker,true);
//       //       dummy.loadWalker(thisWalker,true);
//       //       RealType logguide2_now=2.0*(logEpsilon+Guide.evaluateLog(W));
//       Guide.copyFromBuffer(W,w_buffer);
//       RealType logguide2_now = 2.0*(logEpsilon+Guide.getLogPsi());
// 
//       Psi.copyFromBuffer(W,w_buffer);
//       RealType logpsi2_now = 2.0*Psi.getLogPsi();
// 
//       myTimers[1]->start();
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         makeGaussRandomWithEngine(deltaR,RandomGen);
//         bool stucked=true;
//         for (int iat=0; iat<W.getTotalNum(); ++iat)
//         {
//           PosType dr = m_sqrttau*deltaR[iat];
//           //ignore illegal moves
//           if (!W.makeMoveAndCheck(iat,dr))
//           {
//             ++nReject;
//             continue;
//           }
//           //               dummy.makeMoveAndCheck(iat,dr);
//           //PosType newpos = W.makeMove(iat,dr);
//           RealType psi_ratio = Psi.ratio(W,iat);
//           RealType logpsi2_new = 2.0*std::log(std::abs(psi_ratio))+logpsi2_now;
//           RealType logguide2_new = 2.0*std::log(Guide.ratio(W,iat))+ logguide2_now;
// 
//           long double prob = psi_ratio*psi_ratio*(1.0 + expl(logguide2_new-logpsi2_new))/(1.0+ expl(logguide2_now-logpsi2_now));
//           //               app_log()<<prob<<endl;
//           //RealType prob = std::min(1.0e0,ratio*ratio);
//           if (RandomGen() < prob)
//           {
//             stucked=false;
//             ++nAccept;
//             W.acceptMove(iat);
//             Guide.acceptMove(W,iat);
//             Psi.acceptMove(W,iat);
//             logpsi2_now=logpsi2_new;
//             logguide2_now=logguide2_new;
//           }
//           else
//           {
//             ++nReject;
//             W.rejectMove(iat);
//             Psi.rejectMove(iat);
//             Guide.rejectMove(iat);
//           }
//         }
//         if (stucked)
//         {
//           ++nAllRejected;
//           H.rejectedMove(W,thisWalker);
//         }
//       }
//       myTimers[1]->stop();
// 
//       myTimers[2]->start();
//       //thisWalker.R = W.R;
//       //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
//       //w_buffer.rewind();
//       //W.updateBuffer(w_buffer);
//       //       RealType logguide = logguide2_now;
//       RealType logguide = Guide.updateBuffer(W,w_buffer,true);
//       RealType logpsi = Psi.updateBuffer(W,w_buffer,true);
//       W.saveWalker(thisWalker);
// 
//       myTimers[2]->stop();
// 
//       myTimers[3]->start();
//       RealType eloc=H.evaluate(W);
//       myTimers[3]->stop();
// 
//       thisWalker.resetProperty(logpsi,Psi.getPhase(),eloc);
//       //       thisWalker.resetProperty(0.5*logpsi2_now,Psi.getPhase(),eloc);
//       H.auxHevaluate(W,thisWalker);
//       //       thisWalker.Weight = 1.0;
//       //       thisWalker.Weight = expl(logpsi2_now-logguide2_now);
//       thisWalker.Weight = 1.0/(1.0+expl(2.0*(logguide+logEpsilon-logpsi)));
//       H.saveProperty(thisWalker.getPropertyBase());
//     }
//     myTimers[0]->stop();
//   }
// 
//   void VMCUpdatePbyPSampleRN::advanceCSWalkers(vector<TrialWaveFunction*>& pclone
//       , vector<MCWalkerConfiguration*>& wclone
//       , vector<QMCHamiltonian*>& hclone
//       , vector<RandomGenerator_t*>& rng, vector<RealType>& c_i)
//   {
//     int NumThreads(pclone.size());
// 
//     //this can be modified for cache etc
//     long double psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI] + c_i[ip];
// 
// #pragma omp parallel 
//     {
//       long double psi2_new=1.0;
//       for (int ip=1; ip<NumThreads; ++ip)
//         psi2_new += std::exp(psi2_i_new[ip]-psi2_i_new[0]);
// 
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
// 
//       //copy the new to now
//       long double psi2_now=psi2_new;
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
// 
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
// 
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           PosType dr=m_sqrttau*deltaR[iat];
// 
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
// 
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*std::log(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
//           psi2_new=1.0;
//           for(int jp=1; jp<NumThreads; ++jp)
//             psi2_new+= expl(psi2_i_new[jp]-psi2_i_new[0]);
//           long double prob = expl(psi2_i_new[0]-psi2_0_now)*(psi2_new/psi2_now);  
//           if (prob > rng_loc())
//           {
//             wclone[ip]->acceptMove(iat);
//             pclone[ip]->acceptMove(*wclone[ip],iat);
// 
//             psi2_i_now=psi2_i_new[ip];
//             psi2_0_now=psi2_i_new[0];
//             psi2_now=psi2_new;
//           }
//           else
//           {
//             wclone[ip]->rejectMove(iat);
//             pclone[ip]->rejectMove(iat);
//           }
//         }//loop over iat
//       }
// 
//       Walker_t& thisWalker(*W[ip]);
//       Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       wclone[ip]->saveWalker(*W[ip]);
//       RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
//       thisWalker.Weight=1.0;
//       hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
// 
//     myTimers[0]->stop();
//   }
// 
//   void VMCUpdatePbyPSampleRN::estimateNormWalkers(vector<TrialWaveFunction*>& pclone
//       , vector<MCWalkerConfiguration*>& wclone
//       , vector<QMCHamiltonian*>& hclone
//       , vector<RandomGenerator_t*>& rng
//       , vector<RealType>& ratio_i_0)
//   {
//     int NumThreads(pclone.size());
// 
//     //this can be modified for cache etc
//     RealType psi2_i_new[64];
//     for (int ip=0; ip<NumThreads; ++ip)
//       psi2_i_new[ip] = 2.0*W[ip]->getPropertyBase()[LOGPSI];
// 
// //     RealType nn = -std::log(1.0*nSubSteps*W.getTotalNum());
// #pragma omp parallel 
//     {
//       int nptcl=W.getTotalNum();
//       int ip=omp_get_thread_num();
//       RandomGenerator_t& rng_loc(*rng[ip]);
// 
//       //copy the new to now
//       RealType psi2_i_now=psi2_i_new[ip];
//       RealType psi2_0_now=psi2_i_new[0];
// 
//       for (int iter=0; iter<nSubSteps; ++iter)
//       {
//         //create a 3N-Dimensional Gaussian with variance=1
//         makeGaussRandomWithEngine(deltaR,rng_loc);
// 
//         for (int iat=0; iat<nptcl; ++iat)
//         {
//           PosType dr=m_sqrttau*deltaR[iat];
// 
//           bool movePtcl = wclone[ip]->makeMoveAndCheck(iat,dr);
//           //everyone should skip this; could be a problem with compilers
//           if (!movePtcl) continue;
// 
//           RealType ratio = pclone[ip]->ratio(*wclone[ip],iat);
// #pragma omp barrier
//           psi2_i_new[ip] = 2.0*logl(std::abs(ratio)) + psi2_i_now;
// #pragma omp barrier
// #pragma omp critical
//           {
//             ratio_i_0[ip] += expl( psi2_i_new[ip]-psi2_i_new[0]);
//           }
//           wclone[ip]->rejectMove(iat);
//           pclone[ip]->rejectMove(iat);
//         }
//       }
// 
//       //     Walker_t& thisWalker(*W[ip]);
//       //     Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
//       //     RealType logpsi = pclone[ip]->updateBuffer(*wclone[ip],w_buffer,false);
//       //     wclone[ip]->saveWalker(*W[ip]);
//       //     RealType eloc=hclone[ip]->evaluate(*wclone[ip]);
//       //     //           thisWalker.resetProperty(0.5*psi2_i_now[ip],pclone[ip]->getPhase(), eloc);
//       //     thisWalker.resetProperty(logpsi,pclone[ip]->getPhase(), eloc);
// 
//       //     hclone[ip]->auxHevaluate(*wclone[ip],thisWalker);
//       //     hclone[ip]->saveProperty(thisWalker.getPropertyBase());
//     }
//     RealType nn = 1.0/ratio_i_0[0];
//     for (int ip=0; ip<NumThreads; ++ip) ratio_i_0[ip]*=nn;
//     myTimers[0]->stop();
//   }

}

/***************************************************************************
 * $RCSfile: VMCUpdatePbyP.cpp,v $   $Author: jnkim $
 * $Revision: 1.25 $   $Date: 2006/10/18 17:03:05 $
 * $Id: VMCUpdatePbyP.cpp,v 1.25 2006/10/18 17:03:05 jnkim Exp $
 ***************************************************************************/
