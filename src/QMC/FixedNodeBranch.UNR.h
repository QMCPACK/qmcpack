//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Jordan Vincent
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
#ifndef OHMMS_QMC_FIXEDNODE_BRANCHER_H
#define OHMMS_QMC_FIXEDNODE_BRANCHER_H
#include "OhmmsData/libxmldefs.h"
#include <deque>
#include <algorithm>
#include <libxml++/libxml++.h>

namespace ohmmsqmc {

  /**
     @brief Class to perform branching.
     *
     *For use in Fixed-Node Diffusion Monte Carlo.
   */
  template<class T>
  class FixedNodeBranch {
    
  public:
    /// the timestep
    T Tau;
    /// feedback parameter to control the population
    T feed;
    /// energy offset to control branching
    T E_T;
    /// ideal population
    int Nideal;
    /// maximum population
    int Nmax;
    /// minumum population
    int Nmin;
    /// determines when to branch
    int Stride;
    /// container to store old values of e_ref
    std::deque<T> Eg;
    /// size of container Eg
    int size;

    int Counter;
    /// Constructor
    FixedNodeBranch(T tau, int nideal): 
      Tau(tau), E_T(0.0), Nideal(nideal), 
      Stride(1), size(100), Counter(0) {
      feed = 1.0/(50.0*tau); 
      Nmax = 2*nideal;
      Nmin = static_cast<int>(0.5*nideal);
      Eg.resize(size);
      for(int i=0; i<Eg.size(); i++) Eg[i] = 0.0; 
    }
    
    ///return true if the nodal surface is crossed
    inline bool operator()(T psi0, T psi1) const { return psi0*psi1 < 0;}
    
    ///returns weight \f$\exp(-\tau \left[E_L(R)+E_L(R')-2E_T\right/2])\f$
    inline T ratio(T tau, T elold, T elnew, T reject) const { 
      return exp(-tau*(elnew+elold-2*E_T)/2.0);
    }
    
    inline void setEguess(T eg){
      E_T = eg;
      for(int i=0; i<Eg.size(); i++) Eg[i] = eg;
    } 

    inline int branch(int iter, MCWalkerConfiguration& w) {

      MCWalkerConfiguration::iterator it = w.begin();
      int iw=0, nw = w.getActiveWalkers();

      typedef MCWalkerConfiguration::Walker_t Walker_t;
      WalkerList_t good,bad,ok;

      while(it != w.end()) {
	//reweight walkers that are older than 4 time steps
	if ((*it)->Properties(Age) > 3 ) {
	  (*it)->Properties(Weight) = min(0.5,(*it)->Properties(Weight));
	}
	T wgt=(*it)->Properties(Multiplicity);
	if(wget<0.1) 
	  bad.push_back(*it);
	else if(wget>1.0) 
	  good.push_back(*it);
	else 
	  ok.push_back(*it);
	it++;
      }

      
      else if(*
	//limit maximun number of copies to 10
	int ncopy = min(static_cast<int>((*it)->Properties(Weight) + Random()),10);
	
	//(*it)->Properties(Weight) = 1.0;
	if(ncopy == 0) {
	  it = w.destroyWalker(it);
	} else {
	  if(ncopy>1) {
	    //add ncopy-1 walkers to the end of the list
	    (*it)->Properties(Weight) /= static_cast<RealType>(ncopy);
	    w.copyWalker(it,ncopy-1);
	  }
	  it++;
	}
	iw++;
      }

//       typedef MCWalkerConfiguration::Walker_t Walker_t;
//       list<Walker_t*> keep, discard;
//       while(it != w.end()) {
// 	//reweight walkers that are older than 4 time steps
// 	if ((*it)->Properties(Age) > 3 ) {
// 	  (*it)->Properties(Weight) = min(0.5,(*it)->Properties(Weight));
// 	  //another option is to kill the walker
// 	  //(*it)->Properties(Weight) = 0.0;
// 	}

// 	//limit maximun number of copies to 10
// 	int ncopy = min(static_cast<int>((*it)->Properties(Weight) + Random()),10);
	
// 	//(*it)->Properties(Weight) = 1.0;
	
// 	if(ncopy == 0) {
//           discard.push_back(*it);
// 	} else {
// 	  //add ncopy-1 walkers to the end of the list
// 	  (*it)->Properties(Weight) /= static_cast<RealType>(ncopy);
// 	  keep.push_back(*it);
// 	  for(int icopy=1; icopy<ncopy; icopy++) 
// 	    keep.push_back(new Walker_t(**it));
// 	}
// 	iw++; it++;
//       }

//       it = discard.begin();
//       while(it != discard.end()) delete *it;
//       w.clear();
//       w.copy(keep.begin(), keep.end());

      int nwalkers = w.getActiveWalkers();
      if (nwalkers > Nmax){
	//if too many walkers, kill until the population is 90%
	//of Nmax
	DEBUGMSG("Too many walkers at step " << iter)
	  int nsubtract =  nwalkers-static_cast<int>(0.9*Nmax);
	MCWalkerConfiguration::iterator itend = w.begin();
	for(int i=0; i < nsubtract; i++) itend++;
	w.destroyWalker(w.begin(), itend);
	//JVINCENT attempting to randomly destroy walkers 
	//int nsubtract = nwalkers-static_cast<int>(4*Nmax/5);
	//int i=0;	 
	// while(i<nsubtract) {
	//  advance(itstart,static_cast<int>(Random()*w.getActiveWalkers()));
	//  w.destroyWalker(itstart);
	//  i++;
	// }
	//JVINCENT
	  
      } else if(nwalkers < Nmin) {
	//if too few walkers, copy until the population is 10%
	//more than Nmin
	it = w.begin();
	int nadd = static_cast<int>(Nmin*1.1)-nwalkers;
	if(nadd < nwalkers){
	  int i=0;
	  while(i<nadd){
	    (*it)->Properties(Weight) *= 0.5;
	    w.copyWalker(it,1);
	    it++; i++;
	  }
	} else {
	  cerr << "Too few walkers to copy!" << endl;
	  exit(-1);
	}
      }

      return w.getActiveWalkers();

    }

    inline T update(int pop, T eavg) {
      Counter++;
      //pop off the last value of Eg
      Eg.pop_back();
      //insert a new value at the beggining of the deque
      Eg.push_front(eavg);
      T Esum = 0.0;
      //calculate the average
      //average over the last half of the simulation
      int limit = min(Counter/2+1,size);
      for(int i=0; i<limit; i++) Esum += Eg[i];
      T egavg = Esum/static_cast<T>(limit);
      E_T = egavg - feed*log(static_cast<T>(pop)/static_cast<T>(Nideal));
      //cout << "Input E_t " << E_T << " new average " << egavg << endl;
      //return egavg;
      return E_T;
    }


    bool put(xmlpp::Node* q){
      using namespace xmlpp;
      //  NodeSet bset = q->find("./Branching");
      //       if(bset.empty()){
      // 	WARNMSG("No branching information!")
      // 	  } else {
      // 	    Element* cur = dynamic_cast<Element*>(bset[0]);
      // 	    Element::AttributeList atts = cur->get_attributes();
      // 	    Element::AttributeList::iterator it = atts.begin();
      // 	    while(it != atts.end()) {  
      // 	      const string& aname = (*it)->get_name();
      // 	      if(aname == "e_ref"){
      // 		T eg =  atof(((*it)->get_value()).c_str());
      // 		XMLReport("Branching:Referece Energy = " << eg)
      // 		  for(int i=0; i<Eg.size(); i++)
      // 		    Eg[i] = eg;
      // 	      } else if(aname == "stride"){
      // 		int s = atoi(((*it)->get_value()).c_str());	
      // 		XMLReport("Branching:Stride = " << Stride)
      // 		  Stride = s;		 	
      // 	      } else if(aname == "n_gen"){
      // 		T n = atof(((*it)->get_value()).c_str());
      // 		XMLReport("Branching:Number of Generations = " << n)		 
      // 		  feed = 1.0/(n*Tau);
      // 	      }
      // 	      it++;
      // 	    }
      // 	  }
      NodeSet pset = q->find("./parameter");
      if(pset.empty()){
	WARNMSG("No branching information!")
 	  } else {
	    for(int j=0; j<pset.size(); j++){
	      Element* cur = dynamic_cast<Element*>(pset[j]);
	      Element::AttributeList atts = cur->get_attributes();
	      Attribute* att = cur->get_attribute("name");
	      string pname = att->get_value();
	      xmlNode* curc = cur->cobj();
	      if(pname == "en_ref") {
		T eg;
		putContent(eg,curc);
		for(int i=0; i<Eg.size(); i++) Eg[i] = eg;
		E_T = eg;
	      } else if(pname == "branch") {
		putContent(Stride,curc);
	      } else if(pname == "num_gen") {
		int n;
		putContent(n,curc);
		feed = 1.0/(static_cast<T>(n)*Tau);
	      }
	    }
	  }
      XMLReport("Branching: Referece energy = " << Eg[0])
	XMLReport("Branching: Branch every = " << Stride << " steps.")
	XMLReport("Branching: Feedback parameter = " << feed)

	return true;
    }

  private:
    ///Default Constructor (disabled)
    FixedNodeBranch(){}
  };
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

