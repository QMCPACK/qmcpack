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
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/HDFAttribIO.h"
#include "Particle/MCWalkerConfiguration.h"

namespace qmcplusplus {

  class WalkerControlBase;
  class ScalarEstimatorManager;

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
      ///number of qmc sections executed
      int QMCCounter;
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
      ///control population fluctutaions
      int NumGeneration;
      ///index of the trial energy
      int ETrialIndex;
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
      ///save xml element
      xmlNodePtr myNode;
      ///LogNorm
      vector<RealType> LogNorm;
      ///WalkerController
      WalkerControlBase* WalkerController;
      ///EstimatorManager
      ScalarEstimatorManager*  MyEstimator;
      ///root name
      string RootName;
      ///(yes|no) string to determine SwapMode
      string SwapWalkers;
      ///set of parameters
      ParameterSet m_param;
      ///Constructor
      SimpleFixedNodeBranch(RealType tau, int nideal);

      ///copy constructor
      SimpleFixedNodeBranch(const SimpleFixedNodeBranch& abranch);

      inline bool phaseChanged(RealType psi0, RealType psi1) const { 
        return abs(psi0-psi1) > numeric_limits<RealType>::epsilon();
      }

      /** increment QMCCounter
       *
       * QMCCounter is the number of times any QMC section is processed.
       */
      inline void advanceQMCCounter() { QMCCounter++;}

      /** get the ScalarEstimatorManager */
      ScalarEstimatorManager* getEstimatorManager()
      {
        return MyEstimator;
      }

      /** set the ScalarEstimatorManager
       * @param est estimator created by the first QMCDriver
       * */
      void setEstimatorManager(ScalarEstimatorManager* est)
      {
        MyEstimator=est;
      }

      /** initialize  the WalkerController 
       * @param fixW true, if reconfiguration with the fixed number of walkers is used
       */
      void initWalkerController(RealType tau, bool fixW=false);

      /**  Calculates the Branching Green's function
       * @param tau effective time step
       * @param emixed mixed energy \f$(E_L(R)+E_L(R'))/2\f$
       * @param reject rejection probability
       * @return \f$G_{branch}\f$
       *
       * \f[G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])\f]
       * @note Use the rejection probability \f$q\f$ to limit \f$G_{branch}\f$
       * \f[ G_{branch} = \min \left(\frac{1}{2q},G_{branch}\right). \f]
       */
      inline RealType branchGF(RealType tau, RealType emixed, RealType reject) const { 
        //return min(0.5/(reject+1e-12),exp(-tau*(emix-E_T)));
        return std::exp(-tau*(emixed-E_T));
      }

      /** set the trial energy \f$ <E_G> = eg \f$
       * @param eg input trial energy
       */
      inline void setEguess(RealType eg){ E_T = eg; } 

      inline void setTrialEnergy(RealType etot, RealType wtot) {
        EavgSum=etot;
        WgtSum=wtot;
        E_T=etot/wtot;
      }

      /** perform branching
       * @param iter current step
       * @param w Walker configuration
       */
      void branch(int iter, MCWalkerConfiguration& w);

      /** perform branching
       * @param iter the iteration
       * @param w the walker ensemble
       * @param clones of the branch engine for OpenMP threads
       */
      void branch(int iter, MCWalkerConfiguration& w, vector<ThisType*>& clones);

      /** restart averaging
       * @param counter Counter to determine the cummulative average will be reset.
       */
      void flush(int counter);

      /** reset the internal parameters */
      void reset();

      bool put(xmlNodePtr cur);

      void write(hid_t grp, bool append=false);

      void read(hid_t grp);

      void read(const string& fname);

      /** create map between the parameter name and variables */
      void registerParameters();

      ///start a run
      void start(const string& froot, bool append);
      ///finalize the simulation
      void finalize();

    private:
      ///default constructor (disabled)
      SimpleFixedNodeBranch(){}
  };

}
#endif
/***************************************************************************
 * $RCSfile: SimpleFixedNodeBranch.h,v $   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

