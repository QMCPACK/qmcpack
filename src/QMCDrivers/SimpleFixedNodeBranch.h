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
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H
#define QMCPLUSPLUS_SIMPLE_FIXEDNODE_BRANCHER_H

#include "Configuration.h"
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsData/HDFAttribIO.h"
#include "QMCDrivers/WalkerControlBase.h"

namespace qmcplusplus {

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
  class SimpleFixedNodeBranch: public QMCTraits {

   public:

     typedef SimpleFixedNodeBranch ThisType;

     ///use reconfiguration method for DMC
     bool FixedNumWalkers;
     ///boolean to swap walkers among processors
     int SwapMode;
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
     RealType Tau;
     ///feedback parameter to control the population
     RealType Feed;
     ///energy offset to control branching
     RealType E_T;
     ///Feed*log(N)
     RealType logN;
     ///Accumulation of the energy
     RealType EavgSum;
     ///Accumulation of the weight
     RealType WgtSum;
     ///tolerance to trigger swapping to control populations per node
     RealType PopControl;
     ///LogJacob
     RealType LogJacobRef;

    
     WalkerControlBase* WalkerController;

     ///Constructor
     SimpleFixedNodeBranch(RealType tau, int nideal);

     ///copy constructor
     SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch);

     ///return true if the nodal surface is crossed
     inline bool operator()(RealType psi0, RealType psi1) const { return psi0*psi1 < 0;}

     //inline bool operator()(complex<RealType>& psi0, complex<RealType>& psi1) const { 
     //  return true;
     //}


     /** initialize  the WalkerController 
      * @param fixW true, if reconfiguration with the fixed number of walkers is used
      */
     void initWalkerController(RealType tau, bool fixW=false);

     //void setWeights(MCWalkerConfiguration::iterator it, MCWalkerConfiguration::iterator it_end);

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
     inline RealType branchGF(RealType tau, RealType emixed, RealType reject) const { 
       //return min(0.5/(reject+1e-12),exp(-tau*(emix-E_T)));
       return exp(-tau*(emixed-E_T));
     }

     /** set the trial energy \f$ <E_G> = eg \f$
      * @param eg input trial energy
      */
     inline void setEguess(RealType eg){ E_T = eg; } 

     /** call MCWalkerConfiguration::branch
      *@param iter the iteration
      *@param w the walker ensemble
      *@return the number of walkers after branching
      */
     inline int branch(int iter, MCWalkerConfiguration& w) {
       return WalkerController->branch(iter,w,PopControl);
     }

     int branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones);

     /** restart averaging
      * @param counter Counter to determine the cummulative average will be reset.
      */
     inline void flush(int counter) {
       if(counter == 0) {
         EavgSum=0.0;
         WgtSum=0.0; 
       }
       Counter=counter;
     }

     inline void accumulate(RealType eloc, RealType wgt) {
       Counter++;
       EavgSum += eloc*wgt; WgtSum += wgt;
     }

    /** Update the energy offset
     *@param pop_now current population \f$ P(t) \f$
     *@param ecur local energy
     *@return the energy offset \f$E_T\f$
     *
     * The trial energy is set according to
     *\f[ E_T = <E_G> - feed \log \left( \frac{P(t)}{P_0} \right) \f]
     *<E_G> is a running average over multiple runs.
    */
    inline RealType update(int pop_now, RealType ecur) {
      return E_T = EavgSum/WgtSum-Feed*log(static_cast<RealType>(pop_now))+logN;
    }

    inline RealType CollectAndUpdate(int pop_now, RealType ecur) {
      return E_T = WalkerController->average(EavgSum,WgtSum)-Feed*log(static_cast<RealType>(pop_now))+logN;
    }

    /** reset the internal parameters */
    void reset();

    /**  Parse the xml file for parameters
     *@param cur current xmlNode 
     *
     * Few important parameters are:
     * <ul>
     * <li> en_ref: a reference energy
     * <li> num_gen: number of generations \f$N_G\f$ to reach  equilibrium, used in the feedback parameter
     * \f$ feed = \frac{1}{N_G \tau} \f$ 
     * </ul>
    */
    bool put(xmlNodePtr cur);

    void write(hid_t grp, bool append=false);

    void read(hid_t grp);

    void read(const string& fname);

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

