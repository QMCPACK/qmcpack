//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim and Jordan Vincent
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

#include <deque>
#include <algorithm>
namespace ohmmsqmc {

  /** For use in Fixed-Node Diffusion Monte Carlo. 
   *
   *@brief Provides routines for the fixed-node diffusion
   Monte Carlo algorithm.  
   *
   Calculates the Branching Green's function
   \f[
   G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])
   \f]
   which is used for the weight and mulitplicity of each walker
   \f[ Weight =  G_{branch} \f]
   \f[ Mulitplicity =  G_{branch} + \nu \f]
   and to update the energy offset
   \f[ E_T = <E_G> - feed \log \(\frac{P(t)}{P_0}\) \f]
   where \f$P(t)\f$ is the current population, \f$P_0\f$ is the 
   ideal population, and \f$<E_G>\f$ is an estimate of the
   local energy. 
  */
  template<class T>
  class MolecuFixedNodeBranch {

  public:
    ///the timestep
    T Tau;
    ///feedback parameter to control the population
    T feed;
    ///energy offset to control branching
    T E_T;
    ///ideal population
    int Nideal;
    ///maximum population
    int Nmax;
    ///minumum population
    int Nmin;
    //determines when to branch
    //int Stride;
    ///container to store old values of e_ref
    std::deque<T> Eg;
    ///size of container Eg
    int size;
    ///counts the number of times update has been called
    int Counter;

    ///Constructor
    MolecuFixedNodeBranch(T tau, int nideal): 
      Tau(tau), E_T(0.0), Nideal(nideal), 
      size(100), Counter(0) {
      feed = 1.0/(50.0*tau); 
      Nmax = 2*nideal;
      Nmin = static_cast<int>(0.5*nideal);
      Eg.resize(size);
      for(int i=0; i<Eg.size(); i++) Eg[i] = 0.0; 

    }
    
    ///return true if the nodal surface is crossed
    inline bool operator()(T psi0, T psi1) const { return psi0*psi1 < 0;}
    
    /**
     *@param tau effective time step
     *@param emixed mixed energy \f$(E_L(R)+E_L(R'))/2\f$
     *@param reject rejection probability
     *@return \f$G_{branch}\f$
     *@brief Calculates the Branching Green's function
     \f[G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])\f]
     *@note Use the rejection probability \f$q\f$ to limit \f$G_{branch}\f$
     \f[ G_{branch} = \min\(\frac{1}{2q},G_{branch}\). \f]
   */
    inline T branchGF(T tau, T emixed, T reject) const { 
      return exp(-tau*(emixed-E_T));
      //return min(0.5/(reject+1e-12),exp(-tau*(emix-E_T)));
    }
    
    ///set \f$ <E_G> = eg \f$
    inline void setEguess(T eg){
      E_T = eg;
      for(int i=0; i<Eg.size(); i++) Eg[i] = eg;
    } 


    /**
     *@param iter the iteration
     *@param w the walker ensemble
     *@return the number of walkers after branching
     *@brief For each walker place /f$ncopy = /min(Multiplicity)/f$ 
     *copies back into the ensemble for the next iteration. 
     *If /f$ncopy=0,/f$ then remove the walker from the ensemble.
     *@note If the size of the population is greater than \f$Nmax,\f$
     *delete the extra walkers.  Copy walkers if the population is
     *less than \f$Nmin.\f$
    */
    inline int branch(int iter, MCWalkerConfiguration& w) {

      MCWalkerConfiguration::iterator it = w.begin();
      int iw=0, nw = w.getActiveWalkers();

      while(iw < nw && it != w.end()) {

	//limit maximun number of copies to 10
	//all copies are added at the end of the list
	int ncopy = min(static_cast<int>((*it)->Properties(Multiplicity)),10);
	(*it)->Properties(Weight) = 1.0;
	(*it)->Properties(Multiplicity) = 1.0;
	if(ncopy == 0) {
	  it = w.destroyWalker(it);
	} else {
	  if(ncopy>1) {
	    w.copyWalker(it,ncopy-1);
	  }
	  it++;
	}
	iw++;
      }

      int nwalkers = w.getActiveWalkers();
      if (nwalkers > Nmax){
	/*if too many walkers, kill until the population is 90%
	  of Nmax*/
	//ERRORMSG("Too many walkers at step " << iter)
	int nsubtract =  nwalkers-static_cast<int>(0.9*Nmax);
	MCWalkerConfiguration::iterator itend = w.begin();
	for(int i=0; i < nsubtract; i++) itend++;
	w.destroyWalker(w.begin(), itend);
      } else if(nwalkers < Nmin) {
	/*if too few walkers, copy until the population is 10%
	  more than Nmin*/
	it = w.begin();
	int nadd = static_cast<int>(Nmin*1.1)-nwalkers;
	if(nadd < nwalkers){
	  int i=0;
	  while(i<nadd){
	    //ERRORMSG("Too few walkers at step " << iter)
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

    /**
     *@param pop the current population \f$ P(t) \f$
     *@param eavg the average value of the local energy
     *@return the energy offset \f$E_T\f$
     *@brief Update the energy offset
     \f[ E_T = <E_G> - feed \log \(\frac{P(t)}{P_0}\) \f]
    */
    inline T update(T pop, T eavg) {
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
      // return egavg;
      return E_T;
    }


    /**
     *@param q the current xmlNode 
     *@brief Parse the xml file for parameters
     <ul>
     <li> en_ref: a reference energy
     <li> num_gen: number of generations \f$N_G\f$ to reach 
     equilibrium, used in the feedback parameter
     \f$ feed = \frac{1}{N_G \tau} \f$ 
     <\ul>
    */
    bool put(xmlNodePtr cur){
      cur=cur->children;
      while(cur != NULL) {
	string cname((const char*)(cur->name));
	if(cname == "parameter") {
	  xmlChar* att= xmlGetProp(cur,(const xmlChar*)"name");
	  if(att) {
	    string pname((const char*)att);
	    if(pname == "en_ref") {
	      T eg;
	      putContent(eg,cur);
	      for(int i=0; i<Eg.size(); i++) Eg[i] = eg;
	      E_T = eg;
	    } else if(pname == "num_gen") {
	      int n;
	      putContent(n,cur);
	      feed = 1.0/(static_cast<T>(n)*Tau);
	    }
	  }
	}
	cur=cur->next;
      }
      XMLReport("Branching: Referece energy = " << Eg[0])
      // XMLReport("Branching: Branch every = " << Stride << " steps.")
      XMLReport("Branching: Feedback parameter = " << feed)
      return true;
    }

  private:
    ///default constructor (disabled)
    MolecuFixedNodeBranch(){}
  };
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

