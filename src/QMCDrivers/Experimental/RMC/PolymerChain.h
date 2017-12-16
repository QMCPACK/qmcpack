//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_POLYMERCHAIN_H
#define QMCPLUSPLUS_POLYMERCHAIN_H
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "Utilities/IteratorUtility.h"
#include <deque>
namespace qmcplusplus
{

struct PolymerChain: public std::deque<MCWalkerConfiguration::Walker_t*>
{

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::ParticleGradient_t ParticleGradient_t;
  typedef MCWalkerConfiguration::RealType RealType;

  ///Boolean for the direction of a move
  bool MoveHead;

  ///The number of beads to be cut
  int  NumCuts;

  ///The index of the middle
  int Middle;

  ///The index of the Last
  int Last;

  ///The number of H/Psi pairs
  int nPsi;

  RealType NumMoves;
  RealType AcceptedMoves;

  RealType SumRatioPrimary;
  Vector<int> RefSign;
  Vector<int> TotalSign;
  Vector<RealType> UmbrellaWeight;
  Vector<RealType> LogRatioActionIJ;

  ///The walkers for the proposed move(s).
  std::vector<Walker_t*> heads;

  ///The walkers that belong to the tails.
  std::vector<Walker_t*> tails;

  ///The walker repository to reuse Walkers.
  std::vector<Walker_t*> repository;

  Matrix<ParticleGradient_t*> Gradients;
  Matrix<ParticleGradient_t*> HeadGradients;

  /** constructor
   *@param awalker the walker which is cloned to form a chain
   *@param len the number of chains in a polymer
   *@param movables the number of heads and tails that are moved
   */
  PolymerChain(Walker_t* awalker,int len, int movables):
    MoveHead(true), NumCuts(movables), nPsi(0), NumMoves(0.0),AcceptedMoves(0.0)
  {
    //always add number of beads
    if(len%2 == 0)
      len++;
    Middle = len/2;
    Last = len-1;
    for(int i=0; i<len; i++)
    {
      Walker_t* acopy=new Walker_t(*awalker);
      acopy->ID = i;
      push_back(acopy);
    }
    tails.resize(movables+1,0);
    heads.resize(movables,0);
    for(int i=0; i<2*movables; i++)
    {
      Walker_t* acopy=new Walker_t(*awalker);
      acopy->ID = i+len;
      repository.push_back(acopy);
    }
  }

  /** destructor
   *
   * Need to clean up the walkers in the repository and the polymer chain
   */
  ~PolymerChain()
  {
    delete_iter(Gradients.begin(),Gradients.end());
    delete_iter(HeadGradients.begin(),HeadGradients.end());
    delete_iter(repository.begin(),repository.end());
    delete_iter(this->begin(),this->end());
  }

  inline size_t getID() const
  {
    return (*this)[Middle]->ID;
  }
  inline void subCuts()
  {
    if(NumCuts>1)
      NumCuts--;
  }

  inline void addCuts()
  {
    NumCuts++;
    if(heads.size()<NumCuts)
    {
      heads.resize(2*NumCuts);
      tails.resize(2*NumCuts);
      for(int i=heads.size(); i<2*NumCuts; i++)
        repository.push_back(new Walker_t(*(repository[0])));
    }
  }

  //resize containers for gradients
  inline void resizeArrays(int npsi=1)
  {
    int nptcl=this->front()->R.size();
    if(Gradients.size1() != this->size() || Gradients.size2() != npsi)
    {
      delete_iter(Gradients.begin(),Gradients.end());
      Gradients.resize(this->size(),npsi);
      for(int i=0; i<Gradients.size(); i++)
        Gradients(i) = new ParticleGradient_t(nptcl);
    }
    if(HeadGradients.size1()< NumCuts+1 || HeadGradients.size2() != npsi)
    {
      delete_iter(HeadGradients.begin(),HeadGradients.end());
      HeadGradients.resize(NumCuts+1,npsi);
      for(int i=0; i<HeadGradients.size(); i++)
        HeadGradients(i) = new ParticleGradient_t(nptcl);
    }
    RefSign.resize(npsi);
    TotalSign.resize(npsi);
    UmbrellaWeight.resize(npsi);
    LogRatioActionIJ.resize(npsi*(npsi-1)/2);
    //set the default to one
    UmbrellaWeight=1.0;
    nPsi=npsi;
  }

  /** make tails and heads to make moves
   *@return the pointer to the anchor Walker
   *
   *\if MoveHead == true,
   *tails are built from the end of the chain. The anchor is the first walker.
   *\else
   *tails are built from the start of the chain. The anchor is the last walker.
   *\endif
   *The heads are copied from the repository and the heads walkers will contain
   *the data with the new configuration by the drift-and-diffusion.
   */
  inline Walker_t* makeEnds()
  {
    NumMoves+=1.0;
    for(int i=0; i<NumCuts; i++)
      heads[i]=repository[i];
    Walker_t* anchor = 0;
    if(MoveHead)
    {
      //anchor=(*this)[0];
      anchor=this->front();
      for(int i=0, j=Last; i<NumCuts+1; i++,j--)
      {
        tails[i]=(*this)[j];
      }
    }
    else
    {
      //anchor=(*this)[Last];
      anchor=this->back();
      for(int i=0; i<NumCuts+1; i++)
      {
        tails[i]=(*this)[i];
      }
    }
    //ID has to be copied
    for(int i=0; i<NumCuts; i++)
      heads[i]->ID = tails[i]->ID;
    if(nPsi)
      //Copy the gradients of the anchor
    {
      int curID=anchor->ID;
      for(int ipsi=0; ipsi<nPsi; ipsi++)
        *HeadGradients(0,ipsi)=*Gradients(curID,ipsi);
    }
    //return the anchor
    return anchor;
  }

  inline void updateEnds()
  {
    AcceptedMoves+=1.0;
    if(nPsi)
    {
      for(int i=0, iplus=1; i<NumCuts; i++,iplus++)
      {
        int targetID=heads[i]->ID;
        for(int ipsi=0; ipsi<nPsi; ipsi++)
          *Gradients(targetID,ipsi) = *HeadGradients(iplus,ipsi);
      }
    }
    if(MoveHead)
    {
      for(int i=0; i<NumCuts; i++)
      {
        push_front(heads[i]);
        pop_back();
      }
    }
    else
    {
      for(int i=0; i<NumCuts; i++)
      {
        push_back(heads[i]);
        pop_front();
      }
    }
    //copy NumCuts of the tails to the repository for next step
    for(int i=0; i<NumCuts; i++)
    {
      repository[i] = tails[i];
    }
  }

  /** rather stupid average **/
  inline void average(Walker_t& center)
  {
    center.ID = (*this)[Middle]->ID;
    center.R = (*this)[Middle]->R;
    center.Properties = ((*this)[0]->Properties + (*this)[Last]->Properties)*0.5;
    center.Properties(LOCALPOTENTIAL) = (*this)[Middle]->Properties(LOCALPOTENTIAL);
  }

  inline void flip()
  {
    MoveHead = !MoveHead; //flip the direction
  }
};
}
#endif
