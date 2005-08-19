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
#ifndef OHMMS_QMC_SIMPLE_FIXEDNODE_BRANCHER_H
#define OHMMS_QMC_SIMPLE_FIXEDNODE_BRANCHER_H

#include <deque>
#include <algorithm>
#include <numeric>
#include "OhmmsData/ParameterSet.h"

namespace ohmmsqmc {

  /** Implements branching algorithm for the fixed-node Diffusion Monte Carlo
   *
   * Calculates the Branching Green's function
   * \f[
   * G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])
   * \f]
   * which is used for the weight and mulitplicity of each walker
   * \f[ Weight =  G_{branch} \f]
   * \f[ Mulitplicity =  G_{branch} + \nu \f]
   * and to update the energy offset
   * \f[ E_T = <E_G> - feed \log \left( \frac{P(t)}{P_0} \right) \f]
   * where \f$P(t)\f$ is the current population, \f$P_0\f$ is the 
   * ideal population, and \f$<E_G>\f$ is an estimate of the
   * local energy. 
   */
  template<class T>
  class SimpleFixedNodeBranch {

   public:
    ///counts the number of times update has been called
    int Counter;
    ///ideal population
    int Nideal;
    ///maximum population
    int Nmax;
    ///minumum population
    int Nmin;
    ///maximum copies of a walker
    int MaxCopy;
    ///control population fluctutaions
    int NumGeneration;

    ///the timestep
    T Tau;
    ///feedback parameter to control the population
    T Feed;
    ///energy offset to control branching
    T E_T;
    ///Feed*log(N)
    T logN;
    ///Accumulation of the energy
    T EavgSum;
    ///Accumulation of the weight
    T WgtSum;

    ///Constructor
    SimpleFixedNodeBranch(T tau, int nideal): Counter(0), Nideal(nideal), NumGeneration(50), MaxCopy(10),
      Tau(tau), E_T(0.0), EavgSum(0.0), WgtSum(0.0) {
    }
    
    ///return true if the nodal surface is crossed
    inline bool operator()(T psi0, T psi1) const { return psi0*psi1 < 0;}
    
    /**  Calculates the Branching Green's function
     *@param tau effective time step
     *@param emixed mixed energy \f$(E_L(R)+E_L(R'))/2\f$
     *@param reject rejection probability
     *@return \f$G_{branch}\f$
     *
     \f[G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])\f]
     *@note Use the rejection probability \f$q\f$ to limit \f$G_{branch}\f$
     \f[ G_{branch} = \min \left(\frac{1}{2q},G_{branch}\right). \f]
     */
    inline T branchGF(T tau, T emixed, T reject) const { 
      return exp(-tau*(emixed-E_T));
      //return min(0.5/(reject+1e-12),exp(-tau*(emix-E_T)));
    }
    
    ///set \f$ <E_G> = eg \f$
    inline void setEguess(T eg){
      E_T = eg;
      LOGMSG("Current Counter = " << Counter << " Trial Energy = " << E_T)
    } 


    /** call MCWalkerConfiguration::branch
     *@param iter the iteration
     *@param w the walker ensemble
     *@return the number of walkers after branching
     */
    inline int branch(int iter, MCWalkerConfiguration& w) {
      return w.branch(10,Nmax,Nmin);
    }

    inline void flush(int counter) {
      if(counter== 0) { EavgSum=0.0; WgtSum=0.0; }
      Counter=counter;
    }

    inline void accumulate(T eloc, T wgt) {
      Counter++;
      EavgSum += eloc*wgt; WgtSum += wgt;
    }

    /** Update the energy offset
     *@param pop_now current population \f$ P(t) \f$
     *@param eloc local energy
     *@return the energy offset \f$E_T\f$
     *
     * The trial energy is set according to
     *\f[ E_T = <E_G> - feed \log \left( \frac{P(t)}{P_0} \right) \f]
     *<E_G> is a running average over multiple runs.
    */
    inline T update(T pop_now) {
      return E_T = EavgSum/WgtSum-Feed*log(static_cast<T>(pop_now))+logN;
    }

    /**  Parse the xml file for parameters
     *@param cur current xmlNode 
     *@param LogOut ostream to which the run-time report is sent
     *
     * Few important parameters are:
     * <ul>
     * <li> en_ref: a reference energy
     * <li> num_gen: number of generations \f$N_G\f$ to reach  equilibrium, used in the feedback parameter
     * \f$ feed = \frac{1}{N_G \tau} \f$ 
     * </ul>
    */
    bool put(xmlNodePtr cur, OhmmsInform *LogOut){
      ParameterSet m_param;
      m_param.add(E_T,"en_ref","AU");
      m_param.add(NumGeneration,"num_gen","int");
      m_param.add(MaxCopy,"max_copy","int");
      m_param.add(Nideal,"target_walkers","int");
      m_param.add(EavgSum,"energy_sum","AU");
      m_param.add(WgtSum,"weight_sum","none");
      m_param.put(cur);
      reset();
      LogOut->getStream() << "reference energy = " << E_T << endl;
      LogOut->getStream() << "number of generations = " << NumGeneration << endl;
      LogOut->getStream() << "feedback = " << Feed << endl;
      return true;
    }

    void reset() {
      Nmax = 2*Nideal;
      Nmin = static_cast<int>(Nideal/2);
      Feed = 1.0/(static_cast<T>(NumGeneration)*Tau);
      logN = Feed*log(static_cast<T>(Nideal));
    }

  private:
    ///default constructor (disabled)
    SimpleFixedNodeBranch(){}
  };
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

