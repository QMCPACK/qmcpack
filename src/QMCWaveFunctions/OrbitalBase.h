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

  /** An abstract class for many-body orbitals
   *
   * TrialWaveFunction is a product of  OrbitalBase objects.
   */
  struct OrbitalBase: public QMCTraits {

    ///recasting enum of DistanceTableData to maintain consistency
    enum {SourceIndex  = DistanceTableData::SourceIndex, 
	  VisitorIndex = DistanceTableData::VisitorIndex, 
	  WalkerIndex  = DistanceTableData::WalkerIndex};

    typedef ParticleAttrib<ValueType>    ValueVectorType;
    typedef ParticleAttrib<GradType>     GradVectorType;

    ///default constructor
    inline OrbitalBase(){ }

    ///default destructor
    virtual ~OrbitalBase() { }

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
     *@param P the active ParticleSet
     *@param G the Gradients
     *@param L the Laplacians
     *@return the value
     *
     *Mainly for walker-by-walker move. The initial stage of particle-by-particle
     *move also uses this.
     */
    virtual ValueType
    evaluate(ParticleSet& P, 
	     ParticleSet::ParticleGradient_t& G, 
	     ParticleSet::ParticleLaplacian_t& L) = 0;

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


    /** evalaute the ratio of the new to old orbital value
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *@return \f$\phi(\{\bf R}^{'}\})/\phi(\{\bf R}^{'}\})\f$.
     *
     *Specialized for particle-by-particle move.
     */
    virtual ValueType ratio(ParticleSet& P, int iat) =0;

    /** update the gradient and laplacian values by accepting a move
     *@param P the active ParticleSet
     *@param iat the index of a particle
     *
     *Specialized for particle-by-particle move.
     */
    virtual void update(ParticleSet& P, 
			ParticleSet::ParticleGradient_t& G, 
			ParticleSet::ParticleLaplacian_t& L,
			int iat) =0;

    /** equivalent to evaluate(P,G,L) with write-back function */
    virtual ValueType evaluate(ParticleSet& P,PooledData<RealType>& buf)=0;

    /** add temporary data reserved for particle-by-particle move.*/
    virtual void registerData(ParticleSet& P, PooledData<RealType>& buf) =0;

    /** copy the internal data saved for particle-by-particle move.*/
    virtual void putData(ParticleSet& P, PooledData<RealType>& buf)=0;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

