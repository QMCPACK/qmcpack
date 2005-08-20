// -*- C++ -*-
#ifndef OHMMS_QMC_ORBITALBASE_H
#define OHMMS_QMC_ORBITALBASE_H
#include "Configuration.h"
#include "Particle/ParticleSet.h"
#include "Particle/DistanceTableData.h"
#include "OhmmsData/RecordProperty.h"

/**@file OrbitalBase.h
 *@brief Declaration of OrbitalBase
 */
namespace ohmmsqmc {

  class WalkerSetRef;

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
  struct OrbitalBase: public QMCTraits {

    ///recasting enum of DistanceTableData to maintain consistency
    enum {SourceIndex  = DistanceTableData::SourceIndex, 
	  VisitorIndex = DistanceTableData::VisitorIndex, 
	  WalkerIndex  = DistanceTableData::WalkerIndex};

    typedef ParticleAttrib<ValueType>    ValueVectorType;
    typedef ParticleAttrib<GradType>     GradVectorType;
    typedef PooledData<RealType>         BufferType;

    bool Optimizable;
    ValueType LogValue;
    ValueType SignValue;

    ///default constructor
    inline OrbitalBase(): Optimizable(true),LogValue(1.0),SignValue(1.0){ }

    ///default destructor
    virtual ~OrbitalBase() { }

    inline void setOptimizable(bool optimizeit) { Optimizable = optimizeit;}

    /**reset the Orbital
     *
     *During wavefunction optimizations, the TrailWaveFunction calls
     *the reset function with new functional variables.
     */
    virtual void reset() = 0;

    /**resize the internal storage with multiple walkers
     *@param nwalkers the number of walkers
     *
     *Specialized for all-walker move.
     */
    virtual void resizeByWalkers(int nwalkers) {}

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

    /** evaluate the orbital values of for all the walkers
     *@param W a collection of walkers
     *@param psi an array of orbital values
     *@param G the Gradients
     *@param L the Laplacians
     *@return the value
     *
     *Specialized for all-walker move. G and L are two-dimensional arrays.
     */
    virtual void 
    evaluate(WalkerSetRef& W, 
	     ValueVectorType& psi,
	     WalkerSetRef::WalkerGradient_t& G,
	     WalkerSetRef::WalkerLaplacian_t& L) = 0;


    /** evaluate the ratio of the new to old orbital value
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *@param dG the differential gradient
     *@param dL the differential laplacian
     *@return \f$ \psi( \{ {\bf R}^{'} \} )/ \psi( \{ {\bf R}^{'} \}) \f$
     *
     *Paired with update(ParticleSet& P, int iat).
     */
    virtual ValueType ratio(ParticleSet& P, int iat,
			    ParticleSet::ParticleGradient_t& dG,
			    ParticleSet::ParticleLaplacian_t& dL) = 0;

    virtual void update(ParticleSet& P, int iat) =0;

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
    virtual ValueType updateBuffer(ParticleSet& P, BufferType& buf) =0;

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
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

