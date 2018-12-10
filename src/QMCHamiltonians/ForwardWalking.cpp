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
    
    
#include "QMCHamiltonians/ForwardWalking.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "ParticleBase/ParticleAttribOps.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"
#include <cstdlib>
#include <set>
#include <string>

namespace qmcplusplus
{

bool ForwardWalking::putSpecial(xmlNodePtr cur, QMCHamiltonian& h, ParticleSet& P)
{
  FirstHamiltonian = h.startIndex();
  nObservables=0;
  nValues=0;
  blockT=1;
  //       OhmmsAttributeSet attrib;
  //       attrib.add(blockT,"blockSize");
  //       attrib.put(cur);
  xmlNodePtr tcur = cur->children;
  while(tcur != NULL)
  {
    std::string cname((const char*)tcur->name);
    if(cname == "Observable")
    {
      std::string tagName("none");
      int Hindex(-100);
      int blockSeries(0);
      int blockFreq(0);
      OhmmsAttributeSet Tattrib;
      Tattrib.add(tagName,"name");
      Tattrib.add(blockSeries,"max");
      Tattrib.add(blockFreq,"frequency");
      Tattrib.put(tcur);
      if (tagName.find("*")==std::string::npos)
      {
        //Single Observable case
        int numProps = P.PropertyList.Names.size();
        Hindex = h.getObservable(tagName)+NUMPROPERTIES;
        if(tagName=="LocalPotential")
        {
          Hindex=LOCALPOTENTIAL ;
          tagName="LocPot";
        }
        else
          if(tagName=="LocalEnergy")
          {
            Hindex=LOCALENERGY ;
            tagName="LocEn";
          }
          else
            if (Hindex==(NUMPROPERTIES-1))
            {
              app_log()<<"Not a valid H element("<<Hindex<<") Valid names are:";
              for (int jk=0; jk<h.sizeOfObservables(); jk++)
                app_log()<<"  "<<h.getObservableName(jk);
              app_log()<< std::endl;
              exit(-1);
            }
        Names.push_back(tagName);
        Hindices.push_back( Hindex);
        app_log()<<" Hamiltonian Element "<<tagName<<" was found at "<< Hindex<< std::endl;
        int numT=blockSeries/blockFreq ;
        nObservables+=1;
        nValues+=numT;
        app_log()<<"   "<<numT<<" values will be calculated every "<<blockFreq<<"*tau H^-1"<< std::endl;
        std::vector<int> pms(3);
        pms[0]=blockFreq;
        pms[1]=numT;
        pms[2]=blockSeries+2;
        walkerLengths.push_back(pms);
        int maxWsize=blockSeries+2;
        int pindx = P.addPropertyHistory(maxWsize);
        // summed values.
        //         P.addPropertyHistory(numT);
        Pindices.push_back(pindx);
      }
      else
      {
        bool FOUNDH(false);
        // 	    Multiple observables for this tag
        int found=tagName.rfind("*");
        tagName.replace (found,1,"");
        int numProps = P.PropertyList.Names.size();
        for(int j=0; j<h.sizeOfObservables(); j++)
        {
          std::string Hname = h.getObservableName(j);
          if (Hname.find(tagName) != std::string::npos)
          {
            //               std::vector<int> Parameters;
            //               if(blockSeries==0)
            //                 putContent(Parameters,tcur);
            //               else
            //                 for( int pl=blockFreq;pl<=blockSeries;pl+=blockFreq) Parameters.push_back(pl);
            FOUNDH=true;
            app_log()<<" Hamiltonian Element "<<Hname<<" was found at "<< j<< std::endl;
            Names.push_back(Hname);
            Hindex = j+NUMPROPERTIES;
            Hindices.push_back( Hindex);
            int numT=blockSeries/blockFreq ;
            nObservables+=1;
            nValues+=numT;
            app_log()<<"   "<<numT<<" values will be calculated every "<<blockFreq<<"*tau H^-1"<< std::endl;
            std::vector<int> pms(3);
            pms[0]=blockFreq;
            pms[1]=numT;
            pms[2]=blockSeries+2;
            walkerLengths.push_back(pms);
            int maxWsize=blockSeries+2;
            int pindx = P.addPropertyHistory(maxWsize);
            Pindices.push_back(pindx);
          }
        }
        //handle FOUNDH
        if (FOUNDH)
        {
          nObservables+=1;
        }
        else
        {
          app_log()<<"Not a valid H element("<<Hindex<<") Valid names are:";
          for (int jk=0; jk<h.sizeOfObservables(); jk++)
            app_log()<<"  "<<h.getObservableName(jk);
          app_log()<< std::endl;
          APP_ABORT("ForwardWalking::put");
        }
      }
    }
    tcur = tcur->next;
  }
  app_log()<<"Total number of observables calculated:"<<nObservables<< std::endl;
  app_log()<<"Total number of values calculated:"<<nValues<< std::endl;
  Values.resize(nValues);
  return true;
}

QMCHamiltonianBase* ForwardWalking::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  //nothing to worry, default copy constructor will do
  return new ForwardWalking(*this);
}

void ForwardWalking::addObservables(PropertySetType& plist)
{
  //not used
}

void ForwardWalking::addObservables(PropertySetType& plist, BufferType& collectables)
{
  myIndex=plist.size();
  int nc=0;
  for(int i=0; i<nObservables; ++i)
    for(int j=0; j<walkerLengths[i][1]; ++j,++nc)
    {
      std::stringstream sstr;
      sstr << "FWE_" << Names[i] << "_" << j*walkerLengths[i][0];
      int id=plist.add(sstr.str());
      //         myIndex=std::min(myIndex,id);
      //app_log() <<" Observables named "<<sstr.str() << " at " << id << std::endl;
    }
  app_log()<<"ForwardWalking::Observables [" << myIndex << ", " << myIndex+nc << ")" << std::endl;
}
}

