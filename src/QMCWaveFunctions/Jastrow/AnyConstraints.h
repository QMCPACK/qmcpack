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
#ifndef QMCPLUSPLUS_ANY_CONSTRAINTS_H
#define QMCPLUSPLUS_ANY_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/Jastrow/WMFunctor.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/ComboOrbital.h"

namespace qmcplusplus {

  struct AnyConstraints: public OrbitalConstraintsBase {
    ///analytic functor
    typedef OptimizableFunctorBase<RealType> InFuncType;
    //typedef ComboFunctor<RealType> InFuncType;
    ///numerical functor
    typedef CubicBsplineSingle<RealType> OutFuncType;
    /** class to define a basisGroup to represent a radial function
     */
    struct BasisGroupType
    {
      InFuncType* In_;
      OutFuncType* Out_;
      int NumGridPoints;
      RealType Rcut;
      BasisGroupType():In_(0),Out_(0),Rcut(10.0),NumGridPoints(101){}
      void setGrid(RealType rc, int npts)
      {
        Rcut=rc;
        NumGridPoints=npts;
      }
    };
    ///typedef for BasisGroupContainer 
    typedef map<string,BasisGroupType*> BGContainerType;
    ///flag to tunr on/off spin-dependent term, always off
    bool IgnoreSpin;
    ///has all the basis groups
    BGContainerType BasisGroups;

    AnyConstraints(ParticleSet& p, TrialWaveFunction& psi);
    ~AnyConstraints();
    void apply();
    void addOptimizables(VarRegistry<RealType>& outVars);
    OrbitalBase* createTwoBody();
    OrbitalBase* createOneBody(ParticleSet& source);

    inline void addTwoBodyPart(ComboOrbital* jcombo) {
      //OrbitalBase* twobody = createTwoBody();
      //if (twobody) jcombo->Psi.push_back(twobody);
    }

    bool put(xmlNodePtr cur);

    void addSingleBasisPerSpecies(xmlNodePtr cur);
    BasisGroupType* createBasisGroup(xmlNodePtr cur, const string& elementType);
    void add2BasisGroup(BasisGroupType* curBG, xmlNodePtr cur);
  };

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1693 $   $Date: 2007-02-02 12:13:11 -0600 (Fri, 02 Feb 2007) $
 * $Id: AnyConstratints.h 1693 2007-02-02 18:13:11Z jnkim $ 
 ***************************************************************************/
