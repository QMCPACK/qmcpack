//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
/** @file OrbitalConstraintsBase.h 
 * @brief Declaration of the base class for the constraints on the wavefunctions
 */
#ifndef QMCPLUSPLUS_ORBITALCONSTRAINTSBASE_H
#define QMCPLUSPLUS_ORBITALCONSTRAINTSBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCFactory/OneDimGridFactory.h"
#include "Optimize/VarList.h"

namespace qmcplusplus {
  
  class ComboOrbital;
  
  struct OrbitalConstraintsBase: public QMCTraits {

    typedef OneDimGridFactory::GridType RadialGridType;

    /** xmlNode on which the constraints operate
     */
    xmlNodePtr myNode;
    /** a common radial grid 
     */
    RadialGridType* myGrid;
    ///print the numerical table 
    bool PrintTables;
    /** cutoff radius to truncated any radial grid functor*/
    RealType Rcut;
    /** inVars[variable name] = (variable id, value)
     */
    std::map<std::string,std::pair<std::string,RealType> > inVars;
    ///default contructor
    OrbitalConstraintsBase();
    ///virtual destructor
    virtual ~OrbitalConstraintsBase() {}
    ///Apply contraints
    virtual void apply()=0;
    /** Add named values to outVars for optimization
     * @param outVars optimizable variables
     */
    virtual void addOptimizables(VarRegistry<RealType>& outVars)=0;
    /** Create an OrbitalBase using two-body relation
     * @param target Quantum Particle Set on which an Orbital depend
     * @return A OrbitalBase*, typically ComboOrbital*
     */

    virtual void addTwoBodyPart(ParticleSet& target, ComboOrbital* j)=0;
    /** Add the appropriate orbital (or orbitals in the case of 
     *  a jastrow with a short and a long range part) to the ComboOrbital
     */
     
//    virtual OrbitalBase* createTwoBody(ParticleSet& target)=0;
    /** Create an OrbitalBase using one-body relation
     * @param target Quantum Particle Set on which an Orbital depend
     * @param source Quantum/Classical ParticleSet 
     * @return A OrbitalBase*, typically ComboOrbital*
     */
    virtual OrbitalBase* createOneBody(ParticleSet& target,ParticleSet& source)=0;
    /** Process an xml element
     */
    virtual bool put(xmlNodePtr cur)=0;

    void getParam(xmlNodePtr cur);
    bool getVariables(xmlNodePtr cur);
    void setRadialGrid(ParticleSet& target);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
