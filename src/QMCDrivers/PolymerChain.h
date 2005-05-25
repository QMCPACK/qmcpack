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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "Particle/MCWalkerConfiguration.h"
#include <deque>
namespace ohmmsqmc {

  struct PolymerChain: public std::deque<MCWalkerConfiguration::Walker_t*> {

    typedef MCWalkerConfiguration::Walker_t Walker_t;

    ///Boolean for the direction of a move
    bool MoveHead;

    ///The number of beads to be cut
    int  NumCuts;

    double NumMoves;
    double AcceptedMoves;

    ///The index of the middle
    int Middle;    
    int Last;

    ///The walkers for the proposed move(s).
    std::vector<Walker_t*> heads; 

    ///The walkers that belong to the tails.
    std::vector<Walker_t*> tails;

    ///The walker repository to reuse Walkers.
    std::vector<Walker_t*> repository;

    /** constructor    
     *@param awalker the walker which is cloned to form a chain
     *@param len the number of chains in a polymer
     *@param movables the number of heads and tails that are moved
     */
    PolymerChain(Walker_t* awalker,int len, int movables): 
      MoveHead(true), NumCuts(movables),NumMoves(0.0),AcceptedMoves(0.0) {

      //always add number of beads
      if(len%2 == 0) len++;
      Middle = len/2;

      Last = len-1;
      for(int i=0; i<len; i++) {
	Walker_t* acopy=new Walker_t(*awalker);
	acopy->ID = i;
	push_back(acopy);
      }
     
      tails.resize(movables+1,0);
      heads.resize(movables,0);
      for(int i=0; i<2*movables; i++) {
	Walker_t* acopy=new Walker_t(*awalker);
        acopy->ID = i+len;
	repository.push_back(acopy);
      }
    }
    
    /** destructor
     *
     * Need to clean up the walkers in the repository and the polymer chain
     */
    ~PolymerChain() {
      for(int i=0; i<repository.size(); i++) delete repository[i];
      for(int i=0; i<size(); i++) delete (*this)[i];
    }

    inline size_t getID() const { return (*this)[Middle]->ID;}
    inline void subCuts() {
      if(NumCuts>1) NumCuts--;
    }

    inline void addCuts() {
      NumCuts++;
      if(heads.size()<NumCuts) {
	heads.resize(2*NumCuts);
	tails.resize(2*NumCuts);
	for(int i=heads.size(); i<2*NumCuts; i++) 
	  repository.push_back(new Walker_t(*(repository[0])));
      }
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
    inline Walker_t* makeEnds() {
      NumMoves+=1.0;
      for(int i=0; i<NumCuts; i++) heads[i]=repository[i];
      Walker_t* anchor = 0;
      if(MoveHead) {
	anchor=(*this)[0];
	for(int i=0, j=Last; i<NumCuts+1; i++,j--) {
	  tails[i]=(*this)[j];
	}
      } else {
	anchor=(*this)[Last];
	for(int i=0; i<NumCuts+1; i++) {
	  tails[i]=(*this)[i];
	}
      }
      return anchor;
    }
    
    inline void updateEnds() {
      AcceptedMoves+=1.0;
      if(MoveHead){
	for(int i=0; i<NumCuts; i++) {
	  push_front(heads[i]);
	  pop_back();
	}
      }else {
	for(int i=0; i<NumCuts; i++) {
	  push_back(heads[i]);
	  pop_front();
	}
      }
      //copy NumCuts of the tails to the repository for next step
      for(int i=0; i<NumCuts; i++) {repository[i] = tails[i];}
    }
    
    /** rather stupid average **/
    inline void average(Walker_t& center) {
      center.ID = (*this)[Middle]->ID;
      center.R = (*this)[Middle]->R;
      center.Properties = ((*this)[0]->Properties + (*this)[Last]->Properties)*0.5;
      center.Properties(LOCALPOTENTIAL) = (*this)[Middle]->Properties(LOCALPOTENTIAL);
    }

    inline void flip() {
      MoveHead = !MoveHead; //flip the direction
    }
  };

  class LocalPotentialEstimator {

    int Middle;
    int Counter;
    ofstream* fout;
    std::vector<PolymerChain*>* Polymers;
    std::vector<double> PEavg;
    std::vector<double> PE2;

  public:
  
    LocalPotentialEstimator(std::vector<PolymerChain*>* polymers):
      Middle(0), Counter(0), fout(0),Polymers(polymers) {  
      Middle=(*polymers)[0]->Middle;
      PEavg.resize(Middle+1,0.0);
      PE2.resize(Middle+1,0.0);
    }

    ~LocalPotentialEstimator() {
      if(fout) {delete fout;}
    }

    inline void resetReportSettings(const string& aname) {
      if(fout) {
	delete fout;
      }
      std::string filename(aname);
      filename.append(".pe.dat");
      fout = new ofstream(filename.c_str());
    }

    inline void reset() { 
      Counter = 0;
      for(int i=0; i<PEavg.size(); i++) PEavg[i]=0.0;
      for(int i=0; i<PE2.size(); i++) PE2[i]=0.0;
    }

    inline void report(int iter) {
      (*fout) << iter;
      double wgtinv = 1.0/static_cast<double>(Counter);
      for(int i=0; i<PEavg.size(); i++) {
	double eavg = PEavg[i]*wgtinv;
	(*fout)  <<  setw(15) << eavg;
      }
      (*fout) << endl;
      reset();
    }

    inline void accumulate() {
      Counter+=Polymers->size();
      register double e;
      for(int i=0; i<Polymers->size(); i++) {
	PolymerChain& curLink = *((*Polymers)[i]);
	for(int k=0,j=curLink.size()-1; k<Middle; k++,j--) {
	  e =0.5*(curLink[k]->Properties(LOCALPOTENTIAL)+curLink[j]->Properties(LOCALPOTENTIAL));
	  PEavg[k] += e;
	  PE2[k] += e*e;
	}
        e= curLink[Middle]->Properties(LOCALPOTENTIAL);
	PEavg[Middle] += e;
	PE2[Middle] += e*e;
      }
    }
  };

}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
