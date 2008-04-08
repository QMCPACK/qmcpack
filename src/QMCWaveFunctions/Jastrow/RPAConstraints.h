//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#define QMCPLUSPLUS_RPA_COMBO_CONSTRAINTS_H
#include "QMCWaveFunctions/OrbitalConstraintsBase.h"
#include "QMCWaveFunctions/ComboOrbital.h"
#include "LongRange/LRHandlerBase.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"


namespace qmcplusplus {

  class LRTwoBodyJastrow;

  class RPAPBCConstraints: public OrbitalConstraintsBase {

    public:
      //for optimize need to pull these out
      typedef LRHandlerBase HandlerType;    
      typedef CubicBspline<RealType,LINEAR_1DGRID,FIRSTDERIV_CONSTRAINTS> SplineEngineType;
      typedef CubicSplineSingle<RealType,SplineEngineType> FuncType;
      typedef LinearGrid<RealType> GridType;
      enum {USE_BREAKUP=0, USE_RPA, USE_NONE};

      RPAPBCConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);
      ~RPAPBCConstraints();
      void addOptimizables(OptimizableSetType& outVars);
      void resetParameters(OptimizableSetType& optVariables);
      OrbitalBase* createTwoBody();
      OrbitalBase* createOneBody(ParticleSet& source);
      void addExtra2ComboOrbital(ComboOrbital* jcombo);
      bool put(xmlNodePtr cur);

    private:
      ///@{options to enable/disable features and to set parameters
      bool IgnoreSpin;
      bool DropLongRange;
      bool DropShortRange;
      bool OwnHandler;
      int LongRangeForm;
      RealType Rs;
      RealType Kc;
      RealType Rcut;
      string ID_Rs;
      string MyName;
      ///@}

      /** main handler
      */
      HandlerType* myHandler;
      /** handler which performs a breakup in real and k-space
       *
       * ShortRangeRPA uses realHanlder
       */
      HandlerType* realHandler;
      /** hanlder which can overwrite the long-range part of a RPA jastro
      */
      HandlerType* kspaceHandler;

      ///object to handle the long-range part
      LRTwoBodyJastrow* LongRangeRPA;

      ///@{objects to handle the short-range part
      ///two-body Jastrow function
      OrbitalBase* ShortRangeRPA;
      ///numerical function owned by ShortRangeRPA
      FuncType* nfunc;
      ///adaptor function to initialize nfunc
      ShortRangePartAdapter<RealType>* SRA;
      ///@}
      
      /** container to keep track unique instantiations of HandlerType 
       *
       * The key uses jastrow/@name. This object is used to perform the breakup only once
       * for the handler with the same parameters (name). 
       */
      static map<string,HandlerType*> handlerSet;
  };

 // struct RPAConstraints: public OrbitalConstraintsBase {
 //   typedef RPAJastrow<RealType> FuncType;
 //   bool IgnoreSpin;
 //   RealType Rs;
 //   string ID_Rs;
 //   vector<FuncType*> FuncList;

 //   ~RPAConstraints();

 //   RPAConstraints(ParticleSet& p, TrialWaveFunction& psi, bool nospin=true);

 //   void addOptimizables(OptimizableSetType& outVars); 
 //   /** update the optimizable variables
 //    */
 //   void resetParameters(OptimizableSetType& optVariables);


 //   OrbitalBase* createTwoBody();
 //   OrbitalBase* createOneBody(ParticleSet& source);
 //   void addExtra2ComboOrbital(ComboOrbital* jcombo) {}
 //   bool put(xmlNodePtr cur);
 // };
     
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
