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
    
    


#include "QMCDrivers/EE/VMCMultiRenyiOMP.h"
#include "QMCDrivers/VMC/VMCUpdatePbyP.h"
#include "QMCDrivers/VMC/VMCUpdateAll.h"
#include "OhmmsApp/RandomNumberControl.h"
#include "Message/OpenMP.h"
#include "Message/CommOperators.h"
#include "tau/profiler.h"

namespace qmcplusplus
{

/// Constructor.
VMCMultiRenyiOMP::VMCMultiRenyiOMP(MCWalkerConfiguration& w, TrialWaveFunction& psi, QMCHamiltonian& h,
                                   HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCDriver(w,psi,h,ppool),  CloneManager(hpool),
  myWarmupSteps(0),UseDrift("yes"), logoffset(2.0), logepsilon(0), computeEE("no"),vsize(0),EEN(4),EENR(1),EEdR(1),csEE("no"),csoffset(0)
{
  RootName = "vmc";
  QMCType ="VMCMultiRenyiOMP";
  QMCDriverMode.set(QMC_UPDATE_MODE,1);
  QMCDriverMode.set(QMC_WARMUP,0);
  m_param.add(UseDrift,"useDrift","string");
  m_param.add(UseDrift,"usedrift","string");
  m_param.add(UseDrift,"use_drift","string");
  m_param.add(logepsilon,"logepsilon","double");
  m_param.add(logoffset,"logoffset","double");
  m_param.add(computeEE,"computeEE","string");
  m_param.add(csEE,"csEE","string");
  m_param.add(vsize,"EEsize","double");
  m_param.add(EEN,"EEN","int");
  m_param.add(EENR,"EErsize","int");
  m_param.add(EEdR,"EEdr","double");
  m_param.add(myRNWarmupSteps,"rnwarmupsteps","int");
  m_param.add(myWarmupSteps,"warmupSteps","int");
  m_param.add(myWarmupSteps,"warmupsteps","int");
  m_param.add(myWarmupSteps,"warmup_steps","int");
  m_param.add(nTargetSamples,"targetWalkers","int");
  m_param.add(nTargetSamples,"targetwalkers","int");
  m_param.add(nTargetSamples,"target_walkers","int");
}

bool VMCMultiRenyiOMP::run()
{
  resetRun();
  if (EEN>W.getActiveWalkers()-1)
    EEN=W.getActiveWalkers()-1;
  if(myComm->rank()==0)
  {
    N_dat.str("");
    s_dat.str("");
    ee_dat.str("");
    N_dat<<RootName.c_str()<<".N.dat";
    s_dat<<RootName.c_str()<<".scaling.dat";
    ee_dat<<RootName.c_str()<<".ee.dat";
    file_out.open(ee_dat.str().c_str(),fstream::out | fstream::trunc);
    file_out<<"# indx";
    file_out<<" JV";
    for(int j(0); j<EENR; j++)
      for(int k(0); k<4; k++)
      {
        for(int i(0); i<EEN; i++)
          file_out<<"  EE"<<i<<"_A"<<k<<"_r"<<j;
        for(int i(0); i<EEN; i++)
          file_out<<"  EE"<<i<<"_V"<<k<<"_r"<<j;
        for(int i(0); i<EEN; i++)
          file_out<<"  EE"<<i<<"_S"<<k<<"_r"<<j;
        for(int i(0); i<EEN; i++)
          file_out<<"  EE"<<i<<"_N"<<k<<"_r"<<j;
      }
    file_out<< std::endl;
    file_out.close();
    file_out.open(s_dat.str().c_str(),fstream::out | fstream::trunc);
    file_out<<"# indx";
    for(int j(0); j<EENR; j++)
      for(int k(0); k<4; k++)
      {
        for(int i(1); i<EEN; i++)
          file_out<<"  EER"<<i<<"_A"<<k<<"_r"<<j;
        for(int i(1); i<EEN; i++)
          file_out<<"  EER"<<i<<"_V"<<k<<"_r"<<j;
        for(int i(1); i<EEN; i++)
          file_out<<"  EER"<<i<<"_S"<<k<<"_r"<<j;
        for(int i(1); i<EEN; i++)
          file_out<<"  EER"<<i<<"_N"<<k<<"_r"<<j;
      }
    file_out<< std::endl;
    file_out.close();
    file_out.open(N_dat.str().c_str(),fstream::out| fstream::trunc);
    file_out<<"# indx";
    for(int j(0); j<EENR; j++)
    {
      if(W.groups()>1)
        for(int g(0); g<W.groups()+1; g++)
          for(int i(0); i<10; i++)
            file_out<<"  N"<<g<<"_"<<i<<"_"<<j;
      else
        for(int i(0); i<10; i++)
          file_out<<"   N"<<"_"<<i<<"_"<<j;
    }
    file_out<< std::endl;
    file_out.close();
  }
//       if(csEE!="no")
//       {
//         get a decent normalization for a crazy ;arge log value for the trial function
  std::vector<RealType> Jvalue(NumThreads,0);
  #pragma omp parallel
  {
    int ip=omp_get_thread_num();
    int i=wPerNode[ip];
    while(i<wPerNode[ip+1])
    {
      std::vector<PosType> dR;
      std::vector<int> iat(0);
      std::vector<RealType> lratio(0);
      std::vector<RealType> logs;
      std::vector<RealType> lweights(0);
      {
        Movers[ip]->getLogs(logs);
        Jvalue[ip] += logs[0];
      }
      i++;
    }
  }
  #pragma omp barrier
  RealType nrm=1.0/(myComm->size()*W.getActiveWalkers());
//         if(csEE!="no")
//         {
  for(int i(0); i<Jvalue.size(); i++)
    csoffset += Jvalue[i];
  myComm->allreduce(csoffset);
  csoffset *= nrm;
//         }
  for (int ip=0; ip<NumThreads; ++ip)
    Movers[ip]->csoffset=csoffset;
  app_log()<<"Using csoffset:"<<csoffset<< std::endl;
//       }
  //start the main estimator
  Estimators->start(nBlocks);
  for (int ip=0; ip<NumThreads; ++ip)
    Movers[ip]->startRun(nBlocks,false);
  const bool has_collectables=W.Collectables.size();
  hpmStart(QMC_VMC_0_EVENT,"vmc::main");
  for (int block=0; block<nBlocks; ++block)
  {
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      //IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:(nBlocks+1)*nSteps;
      IndexType updatePeriod=(QMCDriverMode[QMC_UPDATE_MODE])?Period4CheckProperties:0;
      //assign the iterators and resuse them
      MCWalkerConfiguration::iterator wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);

      Movers[ip]->startBlock(nSteps);
      int now_loc=CurrentStep;

      RealType cnorm=1.0/static_cast<RealType>(wPerNode[ip+1]-wPerNode[ip]);

      for (int step=0; step<nSteps; ++step)
      {
        //collectables are reset, it is accumulated while advancing walkers
        wClones[ip]->resetCollectables();
        Movers[ip]->advanceWalkers(wit,wit_end,false);
        if(has_collectables)
          wClones[ip]->Collectables *= cnorm;
        Movers[ip]->accumulate(wit,wit_end);
        ++now_loc;
        //if (updatePeriod&& now_loc%updatePeriod==0) Movers[ip]->updateWalkers(wit,wit_end);
        if (Period4WalkerDump&& now_loc%myPeriod4WalkerDump==0)
          wClones[ip]->saveEnsemble(wit,wit_end);
        if(storeConfigs && (now_loc%storeConfigs == 0))
          ForwardWalkingHistory.storeConfigsForForwardWalking(*wClones[ip]);
      }

      Movers[ip]->stopBlock(false);
    }//end-of-parallel for
    int nnode(W.getActiveWalkers());
    //         set up ratio helpers
    std::vector<qmcplusplus::EEholder> EE_indices(EENR);
    for (int i=0; i<EENR; i++)
      EE_indices[i].resize(nnode,W.groups()+1);
    ParticleSet::ParticlePos_t L(1);
    L.setUnit(PosUnit::LatticeUnit);
    for (int i=0; i<DIM; i++)
      L[0][i]=0;
    if(W.Lattice.SuperCellEnum!=SUPERCELL_OPEN)
    {
      for (int i=0; i<DIM; i++)
        L[0][i]=1;
      W.convert2Cart(L);
//         put all in boxes so we can integrate the right areas
      #pragma omp parallel
      {
        int ip=omp_get_thread_num();
        MCWalkerConfiguration::iterator  wit(W.begin()+wPerNode[ip]), wit_end(W.begin()+wPerNode[ip+1]);
        while(wit!=wit_end)
        {
          for (int iat(0); iat<(*wit)->R.size(); iat++)
            for(int x(0); x<DIM; x++)
            {
              while ((*wit)->R[iat][x]<0)
                (*wit)->R[iat][x]+=L[0][x];
              while ((*wit)->R[iat][x]> L[0][x])
                (*wit)->R[iat][x]-=L[0][x];
            }
          wit++;
        }
      }
      if(computeEE=="sphere")
      {
//           set it to point at the middle of the simulation cell
        L.setUnit(PosUnit::LatticeUnit);
        for (int i=0; i<DIM; i++)
          L[0][i]=0.5;
        W.convert2Cart(L);
//           app_log()<<" Using EE_R: ";
//           for(int i(0);i<EENR;i++) app_log()<<std::sqrt(vsize+i*EEdR)<<" ";
//           app_log()<< std::endl;
      }
    }
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      int i=wPerNode[ip];
      while(i<wPerNode[ip+1]) //       for(int i(0); i<nnode; i++)
      {
        for (int g=0; g<W.groups(); ++g)
          for (int iat=W.first(g); iat<W.last(g); ++iat)
          {
            //               several volumes we can integrate out
            RealType z;
            if(computeEE=="square")
            {
              z=(W[i]->R[iat][0]);
              for (int x=0; x<DIM; x++)
                z=std::max(z,W[i]->R[iat][x]);
            }
            else
              if(computeEE=="sphere")
              {
                z=0;
                for (int x=0; x<DIM; x++)
                  z+=std::pow( W[i]->R[iat][x]-L[0][x],2);
                z=std::pow(z,0.5);
              }
              else
                if(computeEE=="halfspace")
                {
                  z=(W[i]->R[iat][DIM-1]);
                }
            for(int j(0); j<EENR; j++)
              if (z<(vsize+j*EEdR)) //only record smallest shell that it is in
              {
                EE_indices[j].g_w_pindex[g][i].push_back(iat);
                EE_indices[j].g_w_num[g][i]++;
                EE_indices[j].g_w_num[W.groups()][i]++;
                break;
              }
          }
        i++;
      }
    }
    #pragma omp barrier
    for(int j(1); j<EENR; j++)
      for (int g=0; g<W.groups()+1; g++)
        for(int i(0); i<nnode; i++)
          EE_indices[j].g_w_num[g][i] += EE_indices[j-1].g_w_num[g][i];
//         for(int j(0);j<EENR;j++) for (int g=0; g<W.groups()+1;g++)
//         app_log()<<" "<<EE_indices[j].g_w_num[g][0];
//         app_log()<< std::endl;
    int EENR2(2*EENR);
    int EENR3(3*EENR);
    std::vector<Matrix<RealType> > ratios(4*EENR,Matrix<RealType>(nnode,nnode));
    std::vector<RealType> Jvalue(nnode+1,0);
    #pragma omp parallel
    {
      int ip=omp_get_thread_num();
      int i=wPerNode[ip];
      while(i<wPerNode[ip+1]) //           for(int i(0); i<nnode; i++)
      {
        for(int j(0); j<nnode; j++)
        {
          if(i==j)
          {
            for(int k(0); k<EENR; k++)
            {
              ratios[EENR3+k](i,j)=ratios[EENR+k](i,j)=ratios[EENR2+k](i,j)=ratios[k](i,j)=0;
            }
          }
          else
          {
            bool one_same(false);
            std::vector<int> same_number;
            std::vector<int> number_k;
            for(int k(0); k<EENR; k++)
            {
              bool goodk(true);
              for (int g=0; g<W.groups(); g++)
                if (EE_indices[k].g_w_num[g][i] != EE_indices[k].g_w_num[g][j])
                  goodk=false;
              if (goodk)
              {
                if (EE_indices[k].g_w_num[W.groups()][i]==0)
                  ratios[EENR3+k](i,j)=ratios[EENR+k](i,j)=ratios[EENR2+k](i,j)=ratios[k](i,j)=1;
                else
                {
//                       ratios[k](i,j)=1;
                  same_number.push_back(k);
                  number_k.push_back(EE_indices[k].g_w_num[W.groups()][i]);
                  one_same=true;
                  //will get filled in later
                  ratios[EENR3+k](i,j)=ratios[EENR+k](i,j)=ratios[EENR2+k](i,j)=ratios[k](i,j)=1;
                }
              }
              else
              {
                ratios[EENR3+k](i,j)=ratios[EENR+k](i,j)=ratios[EENR2+k](i,j)=ratios[k](i,j)=0;
              }
            }
            if (one_same)
            {
              int nshells(same_number.size());
              int maxk(same_number[nshells-1]);
              int max_Sub(number_k[maxk]);
              std::vector<PosType> dR;
              std::vector<int> iat;
              for (int k=0; k<=maxk; k++)
                for (int g=0; g<W.groups(); g++)
                  for (int p=0; p<EE_indices[k].g_w_pindex[g][j].size(); p++)
                  {
                    int iat_j=EE_indices[k].g_w_pindex[g][j][p];
                    dR.push_back(W[j]->R[iat_j]);
                  }
              int iti(0);
              for (int k=0; k<=maxk; k++)
                for (int g=0; g<W.groups(); g++)
                  for (int p=0; p<EE_indices[k].g_w_pindex[g][i].size(); p++)
                  {
                    int iat_i=EE_indices[k].g_w_pindex[g][i][p];
                    dR[iti] = dR[iti] - W[i]->R[iat_i];
                    iat.push_back(iat_i);
                    iti++;
                  }
              std::vector<RealType> lratio(4*nshells);
//                   std::vector<RealType> logs;
//                   if(csEE!="no")
//                   {
//                     Movers[ip]->advanceWalkerForCSEE((*W[i]), dR, iat,number_k,lratio,lweights,logs);
// #pragma omp critical
//                     {
//                       Jvalue[i] += logs[0];
//                       Jvalue[nnode]+=1;
//                     }
//                   }
//                   else
              Jvalue[i] += Movers[ip]->advanceWalkerForEE((*W[i]), dR, iat,number_k,lratio);
              Jvalue[nnode]+=1;
//                   Movers[0]->advanceWalkerForEE((*W[i]), dR, iat,number_k,lratio);
//                   for(int a(0);a<same_number.size();a++) app_log()<<a<<" "<<i<<" "<<j<<" "<<same_number[a]<<" "<<number_k[a]<<" "<<lratio[a]<< std::endl;
              #pragma omp critical
              for (int k=0; k<nshells; k++)
              {
                ratios[same_number[k]](i,j)=lratio[k];
                ratios[same_number[k]+EENR](i,j)=lratio[nshells+k];
                ratios[same_number[k]+EENR2](i,j)=lratio[nshells*2+k];
                ratios[same_number[k]+EENR3](i,j)=lratio[nshells*3+k];
//                     weights[same_number[k]](i,j)=lweights[k];
              }
            }
          }
        }
        i++;
      }
    }
    #pragma omp barrier
    RealType csJval(0);
//         if(csEE!="no")
//         {
    for(int i(0); i<nnode; i++)
      csJval += Jvalue[i];
    myComm->allreduce(csJval);
    myComm->allreduce(Jvalue[nnode]);
    csJval *= 1.0/Jvalue[nnode];
//         }
//         matrix is computed, now we power up and take trace
    if(myComm->rank()==0)
    {
      file_out.open(ee_dat.str().c_str(),fstream::out | fstream::app);
      file_out<<block;
      file_out<<"  "<<csJval;
      file_out.close();
      file_out.open(s_dat.str().c_str(),fstream::out | fstream::app);
      file_out<<block;
      file_out.close();
      file_out.open(N_dat.str().c_str(),fstream::out | fstream::app);
      file_out<<block;
      file_out.close();
    }
    std::vector<RealType> prev_estimate(EEN,0), prev_nrm(EEN,0);
    int EEN2(EEN*2);
    int EEN3(EEN*3);
    for(int j(0); j<EENR; j++)
    {
      std::vector<RealType> estimateA(4*EEN,0);
      std::vector<RealType> estimateV(4*EEN,0);
      std::vector<RealType> estimateS(4*EEN,0);
      std::vector<RealType> estimateN(4*EEN,0);
      std::vector<RealType> tv(4*EEN,0);
      std::vector<int> itr(EEN+1,0);
      std::vector<Matrix<RealType> > ratios3(4);
      ratios3[0]=ratios[j];
      ratios3[1]=ratios[EENR+j];
      ratios3[2]=ratios[EENR2+j];
      ratios3[3]=ratios[EENR3+j];
      for(itr[ 0 ]=0; itr[ 0 ]<nnode; itr[ 0 ]++)
        for(itr[ 1 ]=itr[0]+1; itr[ 1 ]<nnode; itr[ 1 ]++)
        {
          tv[ 0 ] = ratios3[0](itr[ 0 ],itr[ 1 ]);
          tv[ EEN+0 ] = ratios3[1](itr[ 0 ],itr[ 1 ]);
          tv[ EEN2+0 ] = ratios3[2](itr[ 0 ],itr[ 1 ]);
          tv[ EEN3+0 ] = ratios3[3](itr[ 0 ],itr[ 1 ]);
          if(tv[0]!=0)
          {
            RealType x=tv[ 0 ]*ratios3[0](itr[ 1 ],itr[0]);
            RealType y=tv[ EEN+0 ]*ratios3[1](itr[ 1 ],itr[0]);
            RealType z=tv[ EEN2+0 ]*ratios3[2](itr[ 1 ],itr[0]);
            RealType a=tv[ EEN3+0 ]*ratios3[3](itr[ 1 ],itr[0]);
            estimateA[0]+=x;
            estimateV[ 0 ]+=std::abs(x);
            estimateN[0]+=1;
            estimateA[EEN+0]+=y;
            estimateV[EEN+ 0 ]+=std::abs(y);
            estimateN[EEN+0]+=1;
            estimateA[EEN2+0]+=z;
            estimateV[ EEN2+0 ]+=std::abs(z);
            estimateN[EEN2+0]+=1;
            estimateA[EEN3+0]+=a;
            estimateV[ EEN3+0 ]+=std::abs(a);
            estimateN[EEN3+0]+=1;
//                estimateA[0]+=x; estimateV[ 0 ]+=std::abs(x); estimateS[0]+=(x>0?1:-1); estimateN[0]+=1;
//                estimateA[EEN+0]+=y; estimateV[EEN+ 0 ]+=std::abs(y); estimateS[EEN+0]+=(y>0?1:-1); estimateN[EEN+0]+=1;
//                estimateA[EEN2+0]+=z; estimateV[ EEN2+0 ]+=std::abs(z); estimateS[EEN2+0]+=(z>0?1:-1); estimateN[EEN2+0]+=1;
//                estimateA[EEN3+0]+=a; estimateV[ EEN3+0 ]+=std::abs(a); estimateS[EEN3+0]+=(a>0?1:-1); estimateN[EEN3+0]+=1;
            if(EEN > 1)
              godeeper(1,ratios3,estimateA,estimateV,estimateS,estimateN,tv,itr);
          }
        }
      myComm->allreduce(estimateA);
      myComm->allreduce(estimateV);
      myComm->allreduce(estimateS);
      myComm->allreduce(estimateN);
//           csJval = 1.0/csJval;
//           for(int i(0);i<EEN;i++)
//           {
//             estimateA[EEN*3+i]*=csJval;
//             estimateV[EEN*3+i]*=csJval;
//           }
      RealType nrm=1.0/(myComm->size()*nnode);
      for(int i(0); i<EEN; i++)
      {
//combinatorial factor for each increasing power of S_N
        nrm*=(i+2.0)/(nnode-i-1.0);
        for(int j(0); j<4; j++)
        {
          estimateA[EEN*j+i]*=nrm;
          estimateV[EEN*j+i]*=nrm;
          estimateN[EEN*j+i]*=nrm;
          estimateS[EEN*j+i]*=nrm;
          if (estimateN[EEN*j+i] !=0)
          {
            estimateS[EEN*j+i]=estimateA[EEN*j+i]/estimateV[EEN*j+i];
          }
        }
      }
      Matrix<RealType> SpinN(W.groups()+1,10);
      for(int g(0); g<W.groups(); g++)
        for(int i(0); i<nnode; i++)
          for(int n(1); n<11; n++)
            SpinN(g,n-1) += pow(EE_indices[j].g_w_num[g][i],n);
      if(W.groups()>1)
      {
        for(int i(0); i<nnode; i++)
          for(int n(1); n<11; n++)
            SpinN(W.groups(),n-1) += pow(EE_indices[j].g_w_num[W.groups()][i],n);
      }
      myComm->allreduce(SpinN);
      for(int g(0); g<W.groups()+1; g++)
        for(int n(0); n<10; n++)
          SpinN(g,n) *= 1.0/(myComm->size()*nnode);
      if(myComm->rank()==0)
      {
        file_out.open(ee_dat.str().c_str(),fstream::out | fstream::app);
        for(int j(0); j<4; j++)
        {
          for(int i(0); i<EEN; i++)
            file_out<<"  "<<estimateA[EEN*j+i];
          for(int i(0); i<EEN; i++)
            file_out<<"  "<<estimateV[EEN*j+i];
          for(int i(0); i<EEN; i++)
            file_out<<"  "<<estimateS[EEN*j+i];
          for(int i(0); i<EEN; i++)
            file_out<<"  "<<estimateN[EEN*j+i];
        }
        file_out.close();
        file_out.open(s_dat.str().c_str(),fstream::out | fstream::app);
        for(int j(0); j<4; j++)
        {
          for(int i(1); i<EEN; i++)
            if (estimateA[EEN*j+i-1]!=0)
              file_out<<"  "<<(estimateA[EEN*j+i]/estimateA[EEN*j+i-1]);
            else
              file_out<<"  "<<0.0;
          for(int i(1); i<EEN; i++)
            if (estimateV[EEN*j+i-1]!=0)
              file_out<<"  "<<(estimateV[EEN*j+i]/estimateV[EEN*j+i-1]);
            else
              file_out<<"  "<<0.0;
          for(int i(1); i<EEN; i++)
            if (estimateS[EEN*j+i-1]!=0)
              file_out<<"  "<<(estimateS[EEN*j+i]/estimateS[EEN*j+i-1]);
            else
              file_out<<"  "<<0.0;
          for(int i(1); i<EEN; i++)
            if (estimateN[EEN*j+i-1]!=0)
              file_out<<"  "<<(estimateN[EEN*j+i]/estimateN[EEN*j+i-1]);
            else
              file_out<<"  "<<0.0;
        }
        file_out.close();
        file_out.open(N_dat.str().c_str(),fstream::out | fstream::app);
        if(W.groups()>1)
          for(int g(0); g<W.groups()+1; g++)
            for(int i(0); i<10; i++)
              file_out<<"  "<<SpinN(g,i);
        else
          for(int i(0); i<10; i++)
            file_out<<"  "<<SpinN(0,i);
        file_out.close();
      }
    }
    if(myComm->rank()==0)
    {
      file_out.open(ee_dat.str().c_str(),fstream::out | fstream::app);
      file_out<< std::endl;
      file_out.close();
      file_out.open(N_dat.str().c_str(),fstream::out | fstream::app);
      file_out<< std::endl;
      file_out.close();
      file_out.open(s_dat.str().c_str(),fstream::out | fstream::app);
      file_out<< std::endl;
      file_out.close();
    }
//         app_log()<<"new ratios block"<< std::endl;
//         for(int i(0);i<nnode;i++)
//         {
//           for(int j(0);j<nnode;j++) app_log()<<ratios(i,j)<<" ";
//           app_log()<< std::endl;
//         }
    //Estimators->accumulateCollectables(wClones,nSteps);
    CurrentStep+=nSteps;
    Estimators->stopBlock(estimatorClones);
    //why was this commented out? Are checkpoints stored some other way?
    if(storeConfigs)
      recordBlock(block);
  }//block
  hpmStop(QMC_VMC_0_EVENT);
  Estimators->stop(estimatorClones);
  file_out.close();
  //copy back the random states
  for (int ip=0; ip<NumThreads; ++ip)
    *(RandomNumberControl::Children[ip])=*(Rng[ip]);
  //finalize a qmc section
  return finalize(nBlocks);
}

void VMCMultiRenyiOMP::resetRun()
{
  makeClones(W,Psi,H);
  std::vector<IndexType> samples_th(omp_get_max_threads(),0);
  myPeriod4WalkerDump=(Period4WalkerDump>0)?Period4WalkerDump:(nBlocks+1)*nSteps;
  int samples_this_node = nTargetSamples/myComm->size();
  if (nTargetSamples%myComm->size() > myComm->rank())
    samples_this_node+=1;
  int samples_each_thread = samples_this_node/omp_get_max_threads();
  for (int ip=0; ip<omp_get_max_threads(); ++ip)
    samples_th[ip]=samples_each_thread;
  if(samples_this_node%omp_get_max_threads())
    for (int ip=0; ip < samples_this_node%omp_get_max_threads(); ++ip)
      samples_th[ip] +=1;
  app_log() << "  Samples are dumped every " << myPeriod4WalkerDump << " steps " << std::endl;
  app_log() << "  Total Sample Size =" << nTargetSamples << std::endl;
  app_log() << "  Nodes Sample Size =" << samples_this_node << std::endl;
  for (int ip=0; ip<NumThreads; ++ip)
    app_log()  << "    Sample size for thread " <<ip<<" = " << samples_th[ip] << std::endl;
  app_log() << "  Warmup Steps " << myWarmupSteps << std::endl;
//     if (UseDrift == "rn") makeClones( *(psiPool.getWaveFunction("guide")) );
  if (Movers.empty())
  {
    Movers.resize(NumThreads,0);
    branchClones.resize(NumThreads,0);
    estimatorClones.resize(NumThreads,0);
    Rng.resize(NumThreads,0);
    int nwtot=(W.getActiveWalkers()/NumThreads)*NumThreads;
    if(nwtot!=W.getActiveWalkers())
    {
      app_log() << "  walker number not integer multiple of threads. Removing "<< W.getActiveWalkers()-nwtot <<" walkers. ";
      W.destroyWalkers(W.getActiveWalkers()-nwtot);
      std::vector<int> nwoff(myComm->size()+1,0);
      for(int ip=0; ip<myComm->size(); ip++)
        nwoff[ip+1]=nwoff[ip]+nwtot;
      W.setGlobalNumWalkers(nwoff[myComm->size()]);
      W.setWalkerOffsets(nwoff);
    }
    FairDivideLow(nwtot,NumThreads,wPerNode);
    app_log() << "  Initial partition of walkers ";
    copy(wPerNode.begin(),wPerNode.end(),std::ostream_iterator<int>(app_log()," "));
    app_log() << std::endl;
#if !defined(BGP_BUG)
    #pragma omp parallel for
#endif
    for(int ip=0; ip<NumThreads; ++ip)
    {
      std::ostringstream os;
      estimatorClones[ip]= new EstimatorManager(*Estimators);//,*hClones[ip]);
      estimatorClones[ip]->resetTargetParticleSet(*wClones[ip]);
      estimatorClones[ip]->setCollectionMode(false);
      Rng[ip]=new RandomGenerator_t(*(RandomNumberControl::Children[ip]));
      hClones[ip]->setRandomGenerator(Rng[ip]);
      branchClones[ip] = new BranchEngineType(*branchEngine);
      if (QMCDriverMode[QMC_UPDATE_MODE])
      {
        if (UseDrift == "yes")
        {
          os <<"  PbyP moves with drift, using VMCUpdatePbyPWithDriftFast"<< std::endl;
          Movers[ip]=new VMCUpdatePbyPWithDriftFast(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
          // Movers[ip]=new VMCUpdatePbyPWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        else
        {
          os <<"  PbyP moves with |psi^2|, using VMCUpdatePbyP"<< std::endl;
          Movers[ip]=new VMCUpdatePbyP(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      else
      {
        if (UseDrift == "yes")
        {
          os <<"  walker moves with drift, using VMCUpdateAllWithDrift"<< std::endl;
          Movers[ip]=new VMCUpdateAllWithDrift(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        else
        {
          os <<"  walker moves with |psi|^2, using VMCUpdateAll"<< std::endl;
          Movers[ip]=new VMCUpdateAll(*wClones[ip],*psiClones[ip],*hClones[ip],*Rng[ip]);
        }
        //Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
      }
      if(ip==0)
        app_log() << os.str() << std::endl;
    }
  }
#if !defined(BGP_BUG)
  #pragma omp parallel for
#endif
  for(int ip=0; ip<NumThreads; ++ip)
  {
    //int ip=omp_get_thread_num();
    Movers[ip]->put(qmcNode);
    Movers[ip]->resetRun(branchClones[ip],estimatorClones[ip]);
    if (QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->initWalkersForPbyP(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    else
      Movers[ip]->initWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
    for (int prestep=0; prestep<myWarmupSteps; ++prestep)
      Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
    if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
      Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
  }
//     //JNKIM: THIS IS BAD AND WRONG
//     if (UseDrift == "rn")
//     {
//       RealType avg_w(0);
//       RealType n_w(0);
// #pragma omp parallel
//       {
//         int ip=omp_get_thread_num();
//         for (int step=0; step<myWarmupSteps; ++step)
//         {
//           avg_w=0;
//           n_w=0;
//           for (int prestep=0; prestep<myRNWarmupSteps; ++prestep)
//           {
//             Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//             #pragma omp single
//             {
//               MCWalkerConfiguration::iterator wit(W.begin()), wit_end(W.end());
//               while (wit!=wit_end)
//               {
//                 avg_w += (*wit)->Weight;
//                 n_w +=1;
//                 wit++;
//               }
//             }
//             #pragma omp barrier
//            }
//            #pragma omp single
//            {
//              avg_w *= 1.0/n_w;
//              RealType w_m = avg_w/(1.0-avg_w);
//              w_m = std::log(0.5+0.5*w_m);
//              if (std::abs(w_m)>0.01)
//                logepsilon += w_m;
//            }
//            #pragma omp barrier
//            Movers[ip]->setLogEpsilon(logepsilon);
//           }
//
//         for (int prestep=0; prestep<myWarmupSteps; ++prestep)
//           Movers[ip]->advanceWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],true);
//
//         if (myWarmupSteps && QMCDriverMode[QMC_UPDATE_MODE])
//           Movers[ip]->updateWalkers(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1]);
//       }
//     }
  for(int ip=0; ip<NumThreads; ++ip)
  {
    wClones[ip]->clearEnsemble();
    wClones[ip]->setNumSamples(samples_th[ip]);
  }
  myWarmupSteps=0;
  //Used to debug and benchmark opnemp
  //#pragma omp parallel for
  //    for(int ip=0; ip<NumThreads; ip++)
  //    {
  //      Movers[ip]->benchMark(W.begin()+wPerNode[ip],W.begin()+wPerNode[ip+1],ip);
  //    }
}

//  void VMCMultiRenyiOMP::godeeper(int lvl,int EEN,Matrix<RealType>& ratios,std::vector<RealType>& estimate,std::vector<RealType>& nrm,std::vector<RealType>& tv,std::vector<int>& itr)
//{
//  int itrl(lvl+1);
//  for(itr[ itrl ]=0;itr[ itrl ]<W.getActiveWalkers();itr[ itrl ]++)
//  {
//    bool dupe(false);
//    for(int i(0);i<itrl;i++)
//      if (itr[itrl]==itr[i])
//      {
//        dupe=true; break;
//      }
//    if (!dupe)
//    {
//      tv[ lvl ] = tv[ lvl-1 ]*ratios(itr[ lvl],itr[ itrl ]);
//      estimate[ lvl]+=tv[ lvl ]*ratios(itr[ itrl ],itr[0]); nrm[lvl]++;
//      if(EEN > lvl)
//        godeeper(lvl+1,EEN,ratios,estimate,nrm,tv,itr);
//    }
//  }
//}

void VMCMultiRenyiOMP::godeeper(int lvl,std::vector<Matrix<RealType> >& ratios, std::vector<RealType>& estimateA,std::vector<RealType>& estimateV,std::vector<RealType>& estimateS,std::vector<RealType>& estimateN,std::vector<RealType>& tv,std::vector<int>& itr)
{
  int itrl(lvl+1);
  for(itr[ itrl ]=itr[lvl]+1; itr[ itrl ]<W.getActiveWalkers(); itr[ itrl ]++)
  {
    tv[ lvl ] = tv[ lvl-1 ]*ratios[0](itr[ lvl],itr[ itrl ]);
    tv[ EEN+lvl ] = tv[EENR+ lvl-1 ]*ratios[1](itr[ lvl],itr[ itrl ]);
    tv[ 2*EEN+lvl ] = tv[2*EENR+ lvl-1 ]*ratios[2](itr[ lvl],itr[ itrl ]);
    tv[ 3*EEN+lvl ] = tv[3*EENR+ lvl-1 ]*ratios[3](itr[ lvl],itr[ itrl ]);
    if(tv[ lvl ]!=0)
    {
      RealType x=tv[ lvl ]*ratios[0](itr[ itrl ],itr[0]);
      RealType y=tv[ EEN+lvl ]*ratios[1](itr[ itrl ],itr[0]);
      RealType z=tv[ 2*EEN+lvl ]*ratios[2](itr[ itrl ],itr[0]);
      RealType a=tv[ 3*EEN+lvl ]*ratios[3](itr[ itrl ],itr[0]);
      estimateA[lvl]+=x;
      estimateV[lvl]+=std::abs(x);
      estimateN[ lvl]+=1;
      estimateA[EEN+lvl]+=y;
      estimateV[EEN+lvl]+=std::abs(y);
      estimateN[EEN+ lvl]+=1;
      estimateA[2*EEN+lvl]+=z;
      estimateV[2*EEN+lvl]+=std::abs(z);
      estimateN[2*EEN+ lvl]+=1;
      estimateA[3*EEN+lvl]+=a;
      estimateV[3*EEN+lvl]+=std::abs(a);
      estimateN[3*EEN+ lvl]+=1;
      if(EEN > lvl)
        godeeper(itrl,ratios,estimateA,estimateV,estimateS,estimateN,tv,itr);
    }
  }
}

bool
VMCMultiRenyiOMP::put(xmlNodePtr q)
{
  //nothing to add
  return true;
}
}

