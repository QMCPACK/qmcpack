//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file OrbitalConstraintsBase.h
 * @brief Declaration of the base class for the constraints on the wavefunctions
 */
#ifndef QMCPLUSPLUS_ORBITALCONSTRAINTSBASE_H
#define QMCPLUSPLUS_ORBITALCONSTRAINTSBASE_H
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCFactory/OneDimGridFactory.h"

namespace qmcplusplus
{

class ProductOrbital;

/**  Base class to build Jastrow functions with constraints.
 *
 * OrbitalConstraintsBase is derived from OrbitalBuilderBase  and its derived
 * classes implement a builder for Jastrow functions with constraints, such
 * as cusp conditions.
 */
struct OrbitalConstraintsBase: public OrbitalBuilderBase
{

  typedef OneDimGridFactory::GridType RadialGridType;

  enum {MULTIPLE=0, LONGRANGE, ONEBODY, TWOBODY, THREEBODY, FOURBODY};
  ///bitset to set the Jastrow component
  std::bitset<8> JComponent;
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
  OrbitalConstraintsBase(ParticleSet& p, TrialWaveFunction& psi);
  ///virtual destructor
  virtual ~OrbitalConstraintsBase() {}
  /** reset parameters
   * @param optVariables optimizable variables
   */
  virtual void resetParameters(const opt_variables_type& optVaraibles)=0;

  /** Add the appropriate orbital (or orbitals in the case of
   *  a jastrow with a short and a long range part) to the ProductOrbital
   */
  virtual OrbitalBase* createTwoBody()=0;

  /** Create an OrbitalBase using two-body relation
   * @param target Quantum Particle Set on which an Orbital depend
   * @return A OrbitalBase*, typically ProductOrbital*
   */
  virtual void addExtra2ProductOrbital(ProductOrbital* j)=0;

  /** Create an OrbitalBase using one-body relation
   * @param target Quantum Particle Set on which an Orbital depend
   * @param source Quantum/Classical ParticleSet
   * @return A OrbitalBase*, typically ProductOrbital*
   */
  virtual OrbitalBase* createOneBody(ParticleSet& source)=0;


  void getParam(xmlNodePtr cur);
  bool getVariables(xmlNodePtr cur);
  void setRadialGrid();
};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
