//////////////////////////////////////////////////////////////////
// (c) Copyright 2005- by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_WALKER_CONTROL_BASE_FW_STRUCTURE_H
#define QMCPLUSPLUS_WALKER_CONTROL_BASE_FW_STRUCTURE_H

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
// #include "QMCDrivers/WalkerControlBase.h"
 

namespace qmcplusplus {
  
  struct ForwardWalkingData
    { 
      typedef MCWalkerConfiguration::Walker_t Walker_t;
      typedef TinyVector<float,OHMMS_DIM>       StoredPosType;
      typedef ParticleAttrib<StoredPosType>       StoredPosVector;
      long ID;
      long ParentID;
      StoredPosVector Pos;
      
      inline ForwardWalkingData()
      {
      }
      
      inline ForwardWalkingData(const Walker_t& a)
      {
        Pos.resize(a.R.size()); 
        Pos = a.R;
        ID = a.ID;
        ParentID = a.ParentID; 
      }
       
      inline ForwardWalkingData(const ForwardWalkingData& a)
      {
        Pos.resize(a.Pos.size()); 
        Pos = a.Pos;
        ID = a.ID;
        ParentID = a.ParentID; 
      }
      
      inline ForwardWalkingData(const int a)
      {
        Pos.resize(a); 
      }
      
      inline ForwardWalkingData& operator=(const Walker_t& a) {
        Pos.resize(a.R.size()); 
        Pos = a.R;
        ID = a.ID;
        ParentID = a.ParentID;
        return *this;
      } 
       
      inline ForwardWalkingData& operator=(const ForwardWalkingData& a)
      {
        Pos.resize(a.Pos.size()); 
        Pos = a.Pos;
        ID = a.ID;
        ParentID = a.ParentID; 
        return *this;
      }
      
      inline void resize(const int a)
      {
        Pos.resize(a); 
      }
      
      inline int SizeOf()
      {
        return sizeof(long)*2 + Pos.size()*OHMMS_DIM*sizeof(float);
      }
      
      inline void toFloat(vector<float>& pout)
      {
        pout.resize(Pos.size()*OHMMS_DIM);
        vector<float>::iterator pit=pout.begin() ;
        for (int i=0;i<Pos.size();i++) for (int j=0;j<OHMMS_DIM;j++,pit++) (*pit) = Pos[i][j];
      }
      
      inline void fromFloat(vector<float>& pout)
      {
//         Pos.resize(pout.size()/OHMMS_DIM);
        assert(pout.size()/OHMMS_DIM == Pos.size() );
        vector<float>::iterator pit=pout.begin() ;
        for (int i=0;i<Pos.size();i++) for (int j=0;j<OHMMS_DIM;j++,pit++) Pos[i][j] = (*pit);
      }
      
    };
    
    class ForwardWalkingHistoryObject{
      typedef MCWalkerConfiguration::Walker_t Walker_t;
      typedef TinyVector<float,OHMMS_DIM>       StoredPosType;
      typedef ParticleAttrib<StoredPosType>       StoredPosVector;
      typedef vector<ForwardWalkingData> ForwardWalkingConfiguration;
      
      public:
        vector<ForwardWalkingConfiguration> ForwardWalkingHistory;
        
        inline ForwardWalkingHistoryObject() {}
        
        ~ForwardWalkingHistoryObject()
        {
          this->clearConfigsForForwardWalking();
        }
        
        inline void storeConfigsForForwardWalking(MCWalkerConfiguration& W)
        { 
          vector<ForwardWalkingData> ForwardWalkingHere;
          
          for(std::vector<Walker_t*>::iterator Wit(W.begin()); Wit != W.end(); Wit++ )
          {
            ForwardWalkingData fwhere( *(*Wit) );
            ForwardWalkingHere.push_back(fwhere);
          }
          
          ForwardWalkingHistory.push_back(ForwardWalkingHere);
        }
        
        inline void clearConfigsForForwardWalking()
        {
          ForwardWalkingHistory.clear();
        }
        
        inline int sizeOfConfigsForForwardWalking()
        {
          int szeFW(0);
          int singleSize = ForwardWalkingHistory[0][0].SizeOf();
          for(int i=0;i<ForwardWalkingHistory.size();i++) szeFW += ForwardWalkingHistory[i].size() * singleSize;
          return szeFW;
        }
        
        inline void layoutOfConfigsForForwardWalking(vector<int>& returnVal)
        {
          returnVal.resize(ForwardWalkingHistory.size(),0);
          for(int i=0;i<ForwardWalkingHistory.size();i++) returnVal[i]=ForwardWalkingHistory[i].size();
        }
    };
    
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: mcminis2 $
 * $Revision: 3737 $   $Date: 2009-04-06 11:40:52 -0500 (Mon, 06 Apr 2009) $
 * $Id: WalkerControlBase.h 3737 2009-04-06 16:40:52Z mcminis2 $ 
 ***************************************************************************/

