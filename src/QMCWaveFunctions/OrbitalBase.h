//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_ORBITALBASE_H
#define QMCPLUSPLUS_ORBITALBASE_H
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/RecordProperty.h"
#include "QMCWaveFunctions/OrbitalTraits.h"
#include "Optimize/VarList.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#if defined(ENABLE_SMARTPOINTER)
#include <boost/shared_ptr.hpp>
#endif

/**@file OrbitalBase.h
 *@brief Declaration of OrbitalBase
 */
namespace qmcplusplus {

  ///forward declaration of OrbitalBase
  class OrbitalBase;
  ///forward declaration of DiffOrbitalBase
  class DiffOrbitalBase;

#if defined(ENABLE_SMARTPOINTER)
  typedef boost::shared_ptr<OrbitalBase>     OrbitalBasePtr;
  typedef boost::shared_ptr<DiffOrbitalBase> DiffOrbitalBasePtr;
#else
  typedef OrbitalBase*                       OrbitalBasePtr;
  typedef DiffOrbitalBase*                   DiffOrbitalBasePtr;
#endif

  /**@defgroup OrbitalComponent Orbital group
   * @brief Classes which constitute a many-body trial wave function
   *
   * A many-body trial wave function is 
   * \f[
   * \Psi(\{ {\bf R}\}) = \prod_i \psi_{i}(\{ {\bf R}\}),
   * \f]
   * where \f$\psi\f$s are represented by 
   * the derived classes from OrbtialBase.
   */
  /** @ingroup OrbitalComponent
   * @brief An abstract class for a component of a many-body trial wave function
   */
  struct OrbitalBase: public QMCTraits 
  {

    ///recasting enum of DistanceTableData to maintain consistency
    enum {SourceIndex  = DistanceTableData::SourceIndex, 
	  VisitorIndex = DistanceTableData::VisitorIndex, 
	  WalkerIndex  = DistanceTableData::WalkerIndex};

    typedef ParticleAttrib<ValueType> ValueVectorType;
    typedef ParticleAttrib<GradType>  GradVectorType;
    typedef PooledData<RealType>      BufferType;

    /** boolean to set optimization
     *
     * If true, this object is actively modified during optimization
     */
    bool Optimizable;
    /** boolean to turn on/off the buffer
     */
    bool UseBuffer;
    ///** integer to keep track its usage
    // */
    //int Counter;
    /** current \f$\log\phi \f$
     */
    ValueType LogValue;
    /** current phase 
     */
    RealType PhaseValue;
    /** Pointer to the differential orbital of this object
     *
     * If dPsi=0, this orbital is constant with respect to the optimizable variables
     */
    DiffOrbitalBasePtr dPsi;
    /** A vector for \f$ \frac{\partial \nabla \log\phi}{\partial \alpha} \f$
     */
    GradVectorType dLogPsi;
    /** A vector for \f$ \frac{\partial \nabla^2 \log\phi}{\partial \alpha} \f$
     */
    ValueVectorType d2LogPsi;
    /** Name of this orbital
     */
    string OrbitalName;
    ///list of variables this orbital handles
    opt_variables_type myVars;

    /// default constructor
    OrbitalBase();
    //OrbitalBase(const OrbitalBase& old);

    ///default destructor
    virtual ~OrbitalBase() { }

    inline void setOptimizable(bool optimizeit) { Optimizable = optimizeit;}

    ///assign a differential orbital
    void setDiffOrbital(DiffOrbitalBasePtr d);

    /** check in optimizable parameters
     * @param active a super set of optimizable variables
     *
     * Add the paramemters this orbital manage to active.
     */
    virtual void checkInVariables(opt_variables_type& active)=0;

    /** check out optimizable variables
     *
     * Update myVars index map
     */
    virtual void checkOutVariables(const opt_variables_type& active)=0;

    /** reset the parameters during optimizations
     */
    virtual void resetParameters(const opt_variables_type& active)=0;

    /** print the state, e.g., optimizables */
    virtual void reportStatus(ostream& os)=0;

    /** reset properties, e.g., distance tables, for a new target ParticleSet
     * @param P ParticleSet
     */
    virtual void resetTargetParticleSet(ParticleSet& P)=0;

    /** evaluate the value of the orbital for a configuration P.R
     *@param P  active ParticleSet
     *@param G  Gradients
     *@param L  Laplacians
     *@return the value
     *
     *Mainly for walker-by-walker move. The initial stage of particle-by-particle
     *move also uses this.
     */
    virtual ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L) = 0;

    /** evaluate the value of the orbital
     * @param P active ParticleSet
     * @param G Gradients, \f$\nabla\ln\Psi\f$
     * @param L Laplacians, \f$\nabla^2\ln\Psi\f$
     *
     */
    virtual ValueType
    evaluateLog(ParticleSet& P, 
        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L) = 0;


    /** evaluate the ratio of the new to old orbital value
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *@param dG the differential gradient
     *@param dL the differential laplacian
     *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'} \}) \f$
     *
     *Paired with acceptMove(ParticleSet& P, int iat).
     */
    virtual ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL) = 0;

    /** evaluate the log ratio of the new to old orbital value
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *@param dG the differential gradient
     *@param dL the differential laplacian
     *@return \f$ \log\[\psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'} \})\] \f$
     *
     *Paired with update(ParticleSet& P, int iat).
     */
    virtual ValueType logRatio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL) = 0;

    /** a move for iat-th particle is accepted. Update the content for the next moves
     * @param P target ParticleSet
     * @param iat index of the particle whose new position was proposed
     */
    virtual void acceptMove(ParticleSet& P, int iat) =0;

    /** a move for iat-th particle is reject. Restore to the content.
     * @param iat index of the particle whose new position was proposed
     */
    virtual void restore(int iat) = 0;

    /** evalaute the ratio of the new to old orbital value
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'}\})\f$
     *
     *Specialized for particle-by-particle move.
     */
    virtual ValueType ratio(ParticleSet& P, int iat) =0;

    /** update the gradient and laplacian values by accepting a move
     *@param P the active ParticleSet
     *@param dG the differential gradients
     *@param dL the differential laplacians
     *@param iat the index of a particle
     *
     *Specialized for particle-by-particle move. Each Hamiltonian 
     *updates its data for next update and evaluates differential gradients
     *and laplacians.
     */
    virtual void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& dG, 
			ParticleSet::ParticleLaplacian_t& dL,
			int iat) =0;


    /** equivalent to evaluate(P,G,L) with write-back function */
    virtual ValueType evaluate(ParticleSet& P,BufferType& buf)=0;

    /** add temporary data reserved for particle-by-particle move.
     *
     * Return the log|psi|  like evalaute evaluateLog
     */
    virtual ValueType registerData(ParticleSet& P, BufferType& buf) =0;

    /** re-evaluate the content and buffer data
     * @param P particle set
     * @param buf Anonymous storage
     *
     * This function is introduced to update the data periodically for particle-by-particle move.
     */
    virtual ValueType updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch=false) =0;

    /** copy the internal data saved for particle-by-particle move.*/
    virtual void copyFromBuffer(ParticleSet& P, BufferType& buf)=0;

    /** dump the internal data to buf for optimizations 
     *
     * Implments the default function that does nothing
     */
    virtual void dumpToBuffer(ParticleSet& P, BufferType& buf) {}

    /** copy the internal data from buf for optimizations
     *
     * Implments the default function that does nothing
     */
    virtual void dumpFromBuffer(ParticleSet& P, BufferType& buf){}

    /** return a proxy orbital of itself
     */
    OrbitalBasePtr makeProxy(ParticleSet& tqp);
    /** make clone 
     * @param tqp target Quantum ParticleSet
     * @param deepcopy if true, make a decopy
     *
     * If not true, return a proxy class
     */
    virtual OrbitalBasePtr makeClone(ParticleSet& tqp) const;
    ///** copy data members from old
    // * @param old existing OrbitalBase from which all the data members are copied.
    // *
    // * It is up to the derived classes to determine to use deep, shallow and mixed copy methods.
    // */
    //virtual void copyFrom(const OrbitalBase& old);
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

