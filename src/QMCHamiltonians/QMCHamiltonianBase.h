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
/**@file QMCHamiltonianBase.h
 *@brief Declaration of QMCHamiltonianBase
 */
#ifndef QMCPLUSPLUS_HAMILTONIANBASE_H
#define QMCPLUSPLUS_HAMILTONIANBASE_H

#include <Particle/ParticleSet.h>
#include <OhmmsData/RecordProperty.h>
#include <Utilities/RandomGenerator.h>
#include <QMCHamiltonians/observable_helper.h>
#include <bitset>

namespace qmcplusplus {

  /**@defgroup hamiltonian Hamiltonian group
   * @brief QMCHamiltonian and its component, QMCHamiltonianBase
   */

  class DistanceTableData;
  class TrialWaveFunction;
  class QMCHamiltonian;

  struct NonLocalData: public QMCTraits {
    IndexType PID;
    RealType Weight;
    PosType Delta;
    inline NonLocalData():PID(-1),Weight(1.0){}
    inline NonLocalData(IndexType id, RealType w, const PosType& d):PID(id),Weight(w),Delta(d) {}
  };

  /** @ingroup hamiltonian
   * @brief An abstract class for Local Energy operators 
   *
   * Return_t is defined as RealTye. 
   * The types should be checked when using complex wave functions.
   */ 
  struct QMCHamiltonianBase: public QMCTraits {
    
    /** type of return value of evaluate
     */
    typedef RealType Return_t;
    typedef ParticleSet::Walker_t::Buffer_t  BufferType;

    enum {PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL};
    ///set the current update mode
    bitset<4> UpdateMode;
    ///starting index of this object
    int myIndex;
    ///number of dependents: to be removed
    int Dependants;
    ///current value
    RealType Value;
    ///a new value for a proposed move
    RealType NewValue;
    /// This is used to store the value for force on the source
    /// ParticleSet.  It is accumulated if setComputeForces(true).
    ParticleSet::ParticlePos_t IonForce;
    ///what is this????
    Walker<Return_t, ParticleSet::ParticleGradient_t>* tWalker;
    ///name of this object
    string myName;
    ///name of dependent object: to be removed
    string depName;
   
    ///constructor
    QMCHamiltonianBase():myIndex(-1),Value(0.0),Dependants(0)
    {
      UpdateMode.set(PRIMARY,1);
    }

    ///virtual destructor
    virtual ~QMCHamiltonianBase() { }

    /** return the mode i
     * @param i index among PRIMARY, OPTIMIZABLE, RATIOUPDATE, PHYSICAL
     */
    inline bool getMode(int i) { return UpdateMode[i];}

    /** default implementation to add named values to a list
     * @param plist RecordNameProperty
     */
    virtual void addObservables(PropertySetType& plist)
    {
      myIndex=plist.add(myName.c_str());
    }

    /*** add to observable descriptor for hdf5 
     * @param h5list contains a set of hdf5 descriptors for observables
     *
     * Default implementation is for scalar observables.
     */
    virtual void registerObservables(vector<observable_helper*>& h5list
        , hid_t gid) const ;

    /** set the values evaluated by this object to plist
     * @param plist RecordNameProperty
     *
     * Default implementation is to assign Value which is updated
     * by evaluate  function using myIndex.
     */
    virtual void setObservables(PropertySetType& plist)
    {
      plist[myIndex]=Value;
    }
    
    virtual void setParticlePropertyList(PropertySetType& plist
        , int offset)
    {
      plist[myIndex+offset]=Value;
    }

    virtual void setHistories(
        Walker<Return_t, ParticleSet::ParticleGradient_t>& ThisWalker)
    {
       tWalker = &(ThisWalker);
    }
    
    /** reset the data with the target ParticleSet
     * @param P new target ParticleSet
     */
    virtual void resetTargetParticleSet(ParticleSet& P) = 0;

    /** Evaluate the local energies of an N-particle configuration
     *@param P input configuration containing N particles
     *@return the value of the Hamiltonian
     */
    virtual Return_t evaluate(ParticleSet& P) = 0; 
    virtual Return_t rejectedMove(ParticleSet& P){ return 0; }
    virtual Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) = 0; 

    /*@{
     * @brief Functions to handle particle-by-particle move
     *
     * Default implementations use evaluate.
     */
    virtual Return_t registerData(ParticleSet& P, BufferType& buffer)
    {
      return evaluate(P);
    }
    virtual Return_t updateBuffer(ParticleSet& P, BufferType& buf)
    {
      return evaluate(P);
    }
    virtual void copyFromBuffer(ParticleSet& P, BufferType& buf)
    {
      Value=evaluate(P);
    }
    virtual void copyToBuffer(ParticleSet& P, BufferType& buf)
    {
    }
    virtual Return_t evaluatePbyP(ParticleSet& P, int active)
    {
      APP_ABORT(myName + " missing evaluatePbyP");
      return NewValue;
    }
    virtual void acceptMove(int active)
    { 
      Value=NewValue;
    }
    virtual void rejectMove(int active)
    {
    }
    /*@}*/

    /** return an average value by collective operation
     */ 
    virtual Return_t getEnsembleAverage() { return 0.0;}

    /** read the input parameter
     * @param cur xml node for a QMCHamiltonianBase object
     */
    virtual bool put(xmlNodePtr cur)=0;

    /** write about the class */
    virtual bool get(std::ostream& os) const =0;
    
    virtual QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)=0;

    virtual void setRandomGenerator(RandomGenerator_t* rng)
    {
      //empty
    }
    
    virtual void add2Hamiltonian(ParticleSet& qp, TrialWaveFunction& psi
        , QMCHamiltonian& targetH);
    //virtual QMCHamiltonianBase* makeDependants(ParticleSet& qp )
    //{
    //  return 0;
    //}

    virtual void setComputeForces(bool compute) 
    {
      // empty
    }

  };
}
#endif

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

