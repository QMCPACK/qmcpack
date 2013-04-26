/*
 * Reptile.h
 *
 * This class is a wrapper/utility class for groups of walkers within
 * MCWalkerConfiguration for Reptation Monte Carlo.
 *
 * Walkers are stored as usual in MCWalkerConfiguration.  However, this interface
 * represents the reptile as a circular queue within a segment of the walker list.
 */

#ifndef QMCPLUSPLUS_REPTILE_H
#define QMCPLUSPLUS_REPTILE_H

#include "QMCDrivers/DriftOperators.h"
#include "Configuration.h"

namespace qmcplusplus
{
class MCWalkerConfiguration;

class Reptile: public QMCTraits
{


public:
  typedef MCWalkerConfiguration::Walker_t Walker_t;
  //typedef Walker_t::Buffer_t              Buffer_t;
  //	typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::iterator WalkerIter_t;

  std::vector<IndexType> Action;
  std::vector<IndexType> TransProb;

  RealType forwardprob;
  RealType backwardprob;
  RealType forwardaction;
  RealType backwardaction;

  RealType r2prop;
  RealType r2accept;
  IndexType r2samp;
  IndexType maxSamp;
  RealType tauscale;
  RealType tau;

  RealType erun;
  RealType erun2;
  RealType eest;
  RealType evar;
  IndexType esamp;

  IndexType direction, headindex, nbeads;
  MCWalkerConfiguration& w;
  Walker_t* prophead;
  WalkerIter_t repstart, repend;

  inline Reptile(MCWalkerConfiguration& W, WalkerIter_t start, WalkerIter_t end):
    w(W),repstart(start),repend(end),direction(1),headindex(0),prophead(0), r2prop(0.0), r2accept(0.0),tau(0.0)
  {
    //w=W;
    //repstart(start);
    //repend(end);
    //direction(1);
    //headindex(0);
    maxSamp=10000;
    erun=0.0;
    erun2=0.0;
    eest=0.0;
    evar=100000;
    esamp=0;
    r2samp=0;
    r2accept=1;
    r2prop=1;
    Action.resize(3);
    Action[0]=w.addProperty("ActionBackward");
    Action[1]=w.addProperty("ActionForward");
    Action[2]=w.addProperty("ActionLocal");
    TransProb.resize(2);
    TransProb[0]=w.addProperty("TransProbBackward");
    TransProb[1]=w.addProperty("TransProbForward");
    //forwardaction=0.0;
    //backwardaction=0.0;
    //forwardprob=0.0;
    //backwardprob=0.0;
    nbeads=repend-repstart;
  }

  ~Reptile() {}

  inline IndexType wrapIndex(IndexType repindex)
  {
    return (repindex%nbeads + nbeads)%nbeads ;
  }

  inline Walker_t& getWalker(IndexType i)
  {
    WalkerIter_t bead = repstart + wrapIndex(i);
    return **bead;
  }

  inline IndexType getBeadIndex(IndexType i)
  {
    return wrapIndex(headindex + direction*i) ;
  }
  inline Walker_t& getBead(IndexType i)
  {
    return getWalker(getBeadIndex(i)) ;
  }
  inline Walker_t& getHead()
  {
    return getWalker(getBeadIndex(0)) ;
  }
  inline Walker_t& getTail()
  {
    return getWalker(getBeadIndex(nbeads-1)) ;
  }
  inline Walker_t& getNext()
  {
    return getWalker(getBeadIndex(nbeads-2)) ;
  }
  inline Walker_t& getCenter()
  {
    return getWalker(getBeadIndex( (nbeads-1)/2 ));
  }
  //inline void setProposedHead(){

  inline void flip()
  {
    // direction*=-1;
    // headindex = getBeadIndex(nbeads-1);
    headindex = wrapIndex(headindex - direction);
    direction*=-1;
  }

  inline void setDirection(IndexType dir)
  {
    direction=dir;
  }

  inline void setBead(Walker_t& walker, IndexType i)
  {
    IndexType index = getBeadIndex(i);
    Walker_t& newbead(getWalker(index));
    newbead=walker;								//This should be a hard copy
  }

  inline void setHead(Walker_t& overwrite)
  {
    //overwrite last element.
    headindex = getBeadIndex(nbeads-1);  //sets to position of tail.
    Walker_t& newhead(getBead(0));
    newhead=overwrite;
  }

  inline Walker_t&  getNewHead()
  {
    //overwrite last element.
    headindex = getBeadIndex(nbeads-1);  //sets to position of tail.
    return getWalker(headindex);
  }

  inline void printState()
  {
    app_log()<<"********PRINT REPTILE STATE*********\n";
    app_log()<<"Direction="<<direction<<"  Headindex="<<headindex<<"  tail="<<getBeadIndex(nbeads-1)<<"\n  next="<<getBeadIndex(nbeads-2)<<"  nbeads="<<nbeads<<endl;
    app_log()<<"BeadIndex\tWrapIndex\tEnergy\tAction[0]\tAction[1]\tAction[2]\t\n";
    for( int i=0; i<nbeads; i++)
    {
      app_log()<<i<<"\t"<<getBeadIndex(i)<<"\t"<<getBead(i).Properties(LOCALENERGY)<<"\t"<<getBead(i).Properties(Action[0])<<"\t"<<getBead(i).Properties(Action[1])<<"\t"<<getBead(i).Properties(Action[2])<<"\n";
    }
    app_log()<<"POSITIONS===============:\n";
    for( int i=0; i<nbeads; i++)
    {
      app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).R[0]<<"\n";
      app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).R[1]<<"\n";
    }
    app_log()<<"GVECS===============:\n";
    for( int i=0; i<nbeads; i++)
    {
      app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).G[0]<<"\n";
      app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).G[1]<<"\n";
    }
    app_log()<<"************************************\n";
  }
  //from i to j
  inline void resetR2Avg()
  {
    r2samp=0.0;
    r2prop=0.0;
    r2accept=0.0;
    //app_log()<<"Dumping running averages...\n";
  }

  inline void resetEAvg()
  {
    erun=0.0;
    erun2=0.0;
    esamp=0;
  }

  inline void calcTauScaling()
  {
    if (r2prop != 0.0 && r2samp > maxSamp)
    {
      tauscale=r2accept/r2prop;
      resetR2Avg();
    }
  }

  inline void accumulateE(RealType eloc)
  {
    erun+=eloc;
    erun2+=eloc*eloc;
    esamp++;
  }

  inline void calcERun()
  {
    if (esamp > maxSamp)
    {
      eest=erun/RealType(esamp);
      evar = erun2/RealType(esamp) - eest*eest;
      //app_log()<<"eest="<<eest<<"  evar="<<evar<<endl;
      //app_log()<<"tauscale = "<<tauscale<<endl;
      resetEAvg();
    }
  }


  //	inline RealType calcActionDiff()
  //	{
  //			Walker_t& head = getHead();
  //			Walker_t& newhead=*prophead;
  //			Walker_t& tail=getTail();
  ///		Walker_t& nexttail=getNext();

  //	ParticlePos_t delRfromhead;
  ///		newhead.R - head.R -

  //		return 0.0;

  //0	}








};


}
#endif
