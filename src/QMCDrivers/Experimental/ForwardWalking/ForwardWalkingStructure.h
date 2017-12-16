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


#ifndef QMCPLUSPLUS_WALKER_CONTROL_BASE_FW_STRUCTURE_H
#define QMCPLUSPLUS_WALKER_CONTROL_BASE_FW_STRUCTURE_H

#include <Configuration.h>
#include <Particle/MCWalkerConfiguration.h>
#include <type_traits/scalar_traits.h>

namespace qmcplusplus
{

struct ForwardWalkingData
{
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef TinyVector<float,OHMMS_DIM>   StoredPosType;
  typedef ParticleAttrib<StoredPosType> StoredPosVector;
  long ID;
  long ParentID;
  StoredPosVector Pos;

  inline ForwardWalkingData() { }

  inline ForwardWalkingData(int a):ID(0),ParentID(0)
  {
    Pos.resize(a);
  }

  /** copy walker to the position */
  inline void copy(const Walker_t& a)
  {
    ID = a.ID;
    ParentID = a.ParentID;
    Pos.resize(a.R.size());
    for(int i=0; i<Pos.size(); ++i)
      convert(a.R[i],Pos[i]);
  }
  inline void resize(const int a)
  {
    Pos.resize(a);
  }

  inline int SizeOf()
  {
    return sizeof(long)*2 + Pos.size()*OHMMS_DIM*sizeof(float);
  }

  //template<typename T>
  //inline void toFloat(std::vector<T>& pout)
  //{
  //  pout.resize(Pos.size()*OHMMS_DIM);
  //  int i=0;
  //  for(int k=0; k<Pos.size(); ++k)
  //    for(int dim=0;dim<OHMMS_DIM; ++dim)
  //      pout[i++]=Pos[k][dim];
  //}

  template<typename T>
  inline void fromFloat(std::vector<T>& pin)
  {
    //Pos.resize(pout.size()/OHMMS_DIM);
    assert(pin.size()/OHMMS_DIM == Pos.size() );
    int i=0;
    for(int k=0; k<Pos.size(); ++k)
      for(int dim=0; dim<OHMMS_DIM; ++dim)
        Pos[k][dim]=pin[i++];
  }

  template<typename T>
  inline void append(std::vector<T>& pout) const
  {
    for(int k=0; k<Pos.size(); ++k)
      for(int dim=0; dim<OHMMS_DIM; ++dim)
        pout.push_back(Pos[k][dim]);
  }
};

/** Container for the forward walking history object
 *
 * Using std::vector<std::vector<<ForwardWalkingData>*> to limit the allocation and copy
 */
struct ForwardWalkingHistoryObject
{
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef std::vector<ForwardWalkingData> ForwardWalkingConfiguration;

  ///accumulated number of walkers
  size_t number_of_walkers;
  std::vector<ForwardWalkingConfiguration*> ForwardWalkingHistory;

  inline ForwardWalkingHistoryObject():number_of_walkers(0) {}

  ~ForwardWalkingHistoryObject()
  {
    clearConfigsForForwardWalking();
  }

  inline const ForwardWalkingData& operator()(int i, int j) const
  {
    return (*ForwardWalkingHistory[i])[j];
  }

  inline ForwardWalkingData& operator()(int i, int j)
  {
    return (*ForwardWalkingHistory[i])[j];
  }

  inline void storeConfigsForForwardWalking(MCWalkerConfiguration& W)
  {
    const int nw=W.getActiveWalkers();
    //add the vector first with the size of nw
    ForwardWalkingHistory.push_back(new ForwardWalkingConfiguration(nw));
    //use the reference
    ForwardWalkingConfiguration& now(*ForwardWalkingHistory.back());
    for(int iw=0; iw<nw; ++iw)
      now[iw].copy(*W[iw]);
    number_of_walkers+=nw;
  }

  inline void clearConfigsForForwardWalking()
  {
    for(int i=0; i<ForwardWalkingHistory.size(); ++i)
      delete ForwardWalkingHistory[i];
    ForwardWalkingHistory.clear();
    number_of_walkers=0;
  }

  inline int sizeOfConfigsForForwardWalking()
  {
    int szeFW(0);
    int singleSize = (*ForwardWalkingHistory[0])[0].SizeOf();
    for(int i=0; i<ForwardWalkingHistory.size(); i++)
      szeFW += ForwardWalkingHistory[i]->size() * singleSize;
    return szeFW;
  }

  inline void layoutOfConfigsForForwardWalking(std::vector<int>& returnVal)
  {
    returnVal.resize(ForwardWalkingHistory.size(),0);
    for(int i=0; i<ForwardWalkingHistory.size(); i++)
      returnVal[i]=ForwardWalkingHistory[i]->size();
  }
};

}
#endif

