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
#ifndef QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#define QMCPLUSPLUS_QMCDRIVER_FACTORY_H
#include "OhmmsData/OhmmsElementBase.h"
#include <bitset>
namespace qmcplusplus {

  //forward declaration
  class MCWalkerConfiguration;
  class QMCDriver;
  class ParticleSetPool;
  class WaveFunctionPool;
  class HamiltonianPool;

  struct QMCDriverFactory {

    /*! enum for QMC Run Type */
    enum QMCRunType {DUMMY_RUN, /*!< dummy */
      VMC_RUN, /**< VMC type: vmc, vmc-ptcl, vmc-multiple, vmc-ptcl-multiple */
      DMC_RUN, /**< DMC type: dmc, dmc-ptcl*/
      RMC_RUN, /**< RMC type: rmc, rmc-ptcl */
      OPTIMIZE_RUN /*!< Optimization */
    };

    /*! enum to set the bit to determine the QMC mode */
    enum QMCModeEnum {
      UPDATE_MODE,  /**< bit for move: walker or pbyp */
      MULTIPLE_MODE, /**< bit for multple configuration */
      SPACEWARP_MODE /**< bit for space-warping */
    };

    ///current QMC mode determined by curQmcModeBits
    unsigned long curQmcMode;

    ///3-bit (SPACEWARP_MODE, MULTIPLE_MODE, UPDATE_MODE)
    std::bitset<3> curQmcModeBits;

    ///type of qmcdriver
    QMCRunType curRunType;

    ///name of the current QMCriver
    std::string curMethod;

    /** current MCWalkerConfiguration
     */
    MCWalkerConfiguration *qmcSystem;

    /** current QMCDriver
     */
    QMCDriver *qmcDriver;

    /** ParticleSet Pool
     */
    ParticleSetPool* ptclPool;

    /** TrialWaveFunction Pool
     */
    WaveFunctionPool* psiPool;

    /** QMCHamiltonian Pool
     */
    HamiltonianPool* hamPool;

    /** default constructor **/
    QMCDriverFactory();

    /** set the active qmcDriver */
    bool setQMCDriver(int curSeries, xmlNodePtr cur);

    /** virtual destructor **/
    virtual ~QMCDriverFactory();

    /** create a new QMCDriver 
     */
    void createQMCDriver(xmlNodePtr cur);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

