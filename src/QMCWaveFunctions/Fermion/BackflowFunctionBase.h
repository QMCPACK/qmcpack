//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef QMCPLUSPLUS_BACKFLOW_FUNCTIONBASE_H
#define QMCPLUSPLUS_BACKFLOW_FUNCTIONBASE_H
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Configuration.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

  /**  Base class for backflow transformations. 
   *  FT is an optimizable functor class that implements the radial function
   *  Any class used for Jastrow functions should work
   */
  class BackflowFunctionBase: public OrbitalSetTraits<QMCTraits::ValueType> 
  {

    public:
   
      ///recasting enum of DistanceTableData to maintain consistency
      enum {SourceIndex  = DistanceTableData::SourceIndex,
            VisitorIndex = DistanceTableData::VisitorIndex,
            WalkerIndex  = DistanceTableData::WalkerIndex
           };

    DistanceTableData* myTable; 

    ///Reference to the center
    const ParticleSet& CenterSys;
    ///number of centers, e.g., ions
    int NumCenters;
    ///number of quantum particles
    int NumTargets;

    BackflowFunctionBase(ParticleSet& ions, ParticleSet& els):
     CenterSys(ions), myTable(0) {
      NumCenters=CenterSys.getTotalNum(); // in case
      NumTargets=els.getTotalNum();
    }

    BackflowFunctionBase(BackflowFunctionBase &fn):
     CenterSys(fn.CenterSys), myTable(fn.myTable),NumTargets(fn.NumTargets),NumCenters(fn.NumCenters)
    {}

    virtual
    BackflowFunctionBase* makeClone()=0;

    ~BackflowFunctionBase() {}; 
 
    /** reset the distance table with a new target P
     */
    void resetTargetParticleSet(ParticleSet& P)
    {
      myTable = DistanceTable::add(CenterSys,P);
    }

    virtual void resetParameters(const opt_variables_type& active)=0;

    virtual void checkInVariables(opt_variables_type& active)=0;

    virtual void checkOutVariables(const opt_variables_type& active)=0;

    
    /** calculate quasi-particle coordinates only
     */
    virtual inline void 
    evaluate(const ParticleSet& P, ParticleSet& QP)=0;

    /** calculate quasi-particle coordinates, Bmat and Amat 
     */
    virtual inline void
    evaluate(const ParticleSet& P, ParticleSet& QP, GradMatrix_t& Bmat, HessMatrix_t& Amat)=0;

  };

}
  
#endif
