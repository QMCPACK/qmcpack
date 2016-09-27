//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by:  Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////



#ifndef QMCPLUSPLUS_REPTILECORRELATION_H
#define QMCPLUSPLUS_REPTILECORRELATION_H

#include "ReptileEstimator.h"
#include <cstdlib>

namespace qmcplusplus
{

struct ReptileCorrelationEstimator: public ReptileEstimator
{

  ReptileCorrelationEstimator(MultiChain* polymer)
  {
    headO=-10;
    bodyO=-10;
    sumO=true;
  }

  ~ReptileCorrelationEstimator() {}

  void put(xmlNodePtr cur, int Rlength, int SFH)
  {
    FirstHamiltonian=SFH;
    Rsize = Rlength;
    headObservable="LocalEnergy";
    bodyObservable="PotentialEnergy";
    sumObservable="true";
    Npoints=2;
    OhmmsAttributeSet Tattrib;
    Tattrib.add(headObservable,"headObservable");
    Tattrib.add(bodyObservable,"bodyObservable");
    Tattrib.add(sumObservable,"sum");
    Tattrib.add(Npoints,"intervals");
    Tattrib.put(cur);
    if ((sumObservable=="false")||(sumObservable=="no"))
      sumO=false;
    if (Npoints>Rsize)
    {
      app_log()<<"ERROR! More points("<<Npoints<<") than beads("<<Rsize<<") in reptile"<< std::endl;
      exit(-1);
    }
    if ((Rsize-1) % Npoints !=0 )
    {
      app_log()<<"ERROR!,Number of links(Reptile Length -1="<<Rsize-1<<") must be a multiple of intervals("<<Npoints<<") "<< std::endl;
      exit(-1);
    }
    app_log()<<" Calculating Correlations at:";
    if (Npoints % 2 == 1)
    {
      //Odd
      RealType interval= (Rsize) /Npoints;
      for(RealType i=0; i<Rsize; i+=interval)
      {
        int pnt;
        if (i<=0.5*Rsize)
          pnt =static_cast<int>(i+0.5);
        else
          pnt =static_cast<int>(i);
        outputPoints.push_back(pnt);
        app_log()<<"  "<<pnt;
      }
      if (outputPoints.back()!=Rsize-1)
      {
        outputPoints.push_back(Rsize-1);
        app_log()<<"  "<<Rsize-1;
      }
    }
    else
    {
      RealType interval= (Rsize) /Npoints;
      for(RealType i=0; i<Rsize; i+=interval)
      {
        int pnt = static_cast<int>(i+0.5);
        outputPoints.push_back(pnt);
        app_log()<<"  "<<pnt;
      }
    }
    app_log()<< std::endl;
    for(int i=1; i<outputPoints.size(); i++)
    {
      outputStride.push_back( outputPoints[i]-outputPoints[i-1] );
    }
//       app_log()<<"outputStride "<<outputStride.size()<< std::endl;
  }

  //     void put(xmlNodePtr cur, MCWalkerConfiguration& refWalker) {put(cur);}

  void addNames(std::vector<std::string>& names)
  {
    for(int j=0; j<names.size(); j++)
    {
      if (headObservable==names[j])
        headO=j-1+FirstHamiltonian ;
      if (bodyObservable==names[j])
        bodyO=j-1+FirstHamiltonian ;
    }
    if(headObservable=="LocalPotential")
    {
      headO=LOCALPOTENTIAL;
      headObservable="LocPot";
    }
    else
      if(headObservable=="LogPsi")
        headO=LOGPSI;
      else
        if(headObservable=="LocalEnergy")
        {
          headObservable="LE";
          headO=LOCALENERGY;
        }
    if(bodyObservable=="LocalPotential")
    {
      bodyO=LOCALPOTENTIAL;
      bodyObservable="LocPot";
    }
    else
      if(bodyObservable=="LogPsi")
        bodyO=LOGPSI;
      else
        if(bodyObservable=="LocalEnergy")
        {
          bodyObservable="LE";
          bodyO=LOCALENERGY;
        }
    if ((headO==-10)||(bodyO==-10))
    {
      app_log()<<  "ERROR! the hamiltonian elements do not exist. Choose from the following,"<< std::endl;
      for(int j=0; j<names.size(); j++)
      {
        app_log()<<names[j]<<"  ";
      }
      app_log()<< std::endl;
      exit(-1);
    }
    app_log()<<"  Head Index, Body Index:"<<headO<<", "<<bodyO<< std::endl;
    myIndex=names.size();
    std::stringstream sstr;
    sstr<<headObservable<<"_"<<Rsize;
    names.push_back("H_"+headObservable);
    names.push_back("T_"+headObservable);
    for(int i=0; i<outputPoints.size(); i++)
    {
      std::stringstream sstrb,sstrc;
      sstrb<<headObservable<<"_"<<bodyObservable<<"_"<<outputPoints[i];
      sstrc<<bodyObservable<<"_"<<outputPoints[i];
      names.push_back(sstrb.str());
      names.push_back(sstrc.str());
    }
    if (sumO)
    {
      names.push_back(bodyObservable+"_sum");
      names.push_back(headObservable+"_"+bodyObservable+"_sum");
    }
    Values.resize(names.size()-myIndex);
  }

  void evaluate(MultiChain::iterator first, MultiChain::iterator last, int ipsi)
  {
    MultiChain::iterator end(last), begin(first);
    std::vector<RealType>::iterator Vit(Values.begin());
    RealType Ohead = (*first)->getPropertyBase(ipsi)[headO];
    (*Vit) = Ohead;
    Vit++;
    RealType Otail = (*(last-1))->getPropertyBase(ipsi)[headO];
    (*Vit) = Otail;
    Vit++;
    last--;
    (*Vit) = 0.5*(Ohead* (*first)->getPropertyBase(ipsi)[bodyO] + Otail*(*last)->getPropertyBase(ipsi)[bodyO]);
    Vit++;
    *(Vit) = (*first)->getPropertyBase(ipsi)[bodyO];
    Vit++;
    for(std::vector<int>::iterator strit(outputStride.begin()); strit!=outputStride.end() ; strit++,Vit+=2)
    {
      last-=(*strit);
      first+=(*strit);
      (*Vit)     = 0.5*(Ohead* (*first)->getPropertyBase(ipsi)[bodyO] + Otail*(*last)->getPropertyBase(ipsi)[bodyO]);
      (*(Vit+1)) = (*first)->getPropertyBase(ipsi)[bodyO];
    }
    if (sumO)
    {
      (*Vit)=0.0;
      while(begin!=end)
      {
        (*Vit)+= (*begin)->getPropertyBase(ipsi)[bodyO];
        begin++;
      }
      (*(Vit+1)) = (*Vit)*(Otail+Ohead)*0.5;
    }
  }
  //

private:
  std::string headObservable,bodyObservable,sumObservable;
  int headO ,bodyO;
  bool sumO;
  int Npoints,Rsize,FirstHamiltonian;
  std::vector<int> outputPoints,outputStride;
};
}
#endif
