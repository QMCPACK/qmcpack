//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
/**@file ParticleSetPool.h
 * @brief Declaration of ParticleSetPool
 */
#ifndef QMCPLUSPLUS_PARTICLESETPOOL_H
#define QMCPLUSPLUS_PARTICLESETPOOL_H

#include "OhmmsData/OhmmsElementBase.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Message/MPIObjectBase.h"

namespace qmcplusplus {

  /** @ingroup qmcapp
   * @brief Manage a collection of ParticleSet objects
   *
   * This object handles \<particleset\> elements and
   * functions as a builder class for ParticleSet objects.
   */
  class ParticleSetPool : public MPIObjectBase
  {

  public:

    typedef map<string,ParticleSet*> PoolType;

    /** constructor
     * @param aname xml tag
     */
    ParticleSetPool(Communicate* c, const char* aname = "particleset");

    bool put(xmlNodePtr cur);
    bool get(ostream& os) const;
    void reset();

    ///assign TileMatrix
    bool putTileMatrix(xmlNodePtr cur);

    /** initialize the supercell shared by all the particle sets 
     */
    bool putLattice(xmlNodePtr cur);
    ///return true, if the pool is empty
    inline bool empty() const { return myPool.empty();}

    ///add a ParticleSet* to the pool
    void addParticleSet(ParticleSet* p);
    /** get a named ParticleSet
     * @param pname name of the ParticleSet
     * @return a MCWalkerConfiguration object with pname
     *
     * When the named ParticleSet is not in this object, return 0.
     */
    ParticleSet* getParticleSet(const string& pname);

    /** get a named MCWalkerConfiguration
     * @param pname name of the MCWalkerConfiguration
     * @return a MCWalkerConfiguration object with pname
     *
     * When the named MCWalkerConfiguration is not in this object, return 0.
     */
    MCWalkerConfiguration* getWalkerSet(const string& pname);

    /** get the Pool object
     */
    inline PoolType& getPool() { return myPool;}

    /** create a target particleset and other associated particlesets 
     * @param cur xml node
     * @return A ParticleSet
     *
     * Introduced to avoid conflicting definitions of the particlesets
     */
    ParticleSet* createESParticleSet(xmlNodePtr cur, const string& target);

    /** randomize a particleset particleset/@random='yes' && particleset@random_source exists
     */
    void randomize();

  private:
    /** global SimulationCell
     *
     * SimulationCell cannot not modified once it is initialized by
     * - <simulationcell> element
     * - the first particleset created with ES-HDF
     */
    ParticleSet::ParticleLayout_t* SimulationCell;
    /** tiling matrix 
     */
    Tensor<int,OHMMS_DIM> TileMatrix;
    /** List of ParticleSet 
     *
     * Each ParticleSet has to have a unique name which is used as a key for the map
     */
    map<string,ParticleSet*> myPool;
    /** xml node for random initialization. 
     *
     * randomize() process initializations just before starting qmc sections
     */
    vector<xmlNodePtr> randomize_nodes;
  };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
