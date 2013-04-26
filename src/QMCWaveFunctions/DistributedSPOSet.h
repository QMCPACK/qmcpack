//////////////////////////////////////////////////////////////////
// (c) Copyright 2007-  by Jeongnim Kim
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
#ifndef QMCPLUSPLUS_DISTRIBUTEDORBITALSET_H
#define QMCPLUSPLUS_DISTRIBUTEDORBITALSET_H
/**@file DistributedSPOSet.h
 * @brief Declaration of distributed single-particle orbital set class
 */
#include "QMCWaveFunctions/SPOSetBase.h"
#include "Utilities/PooledData.h"

namespace qmcplusplus
{

/** class specialized for small-memory distributed machines
 *
 * Limitations: Only valid for the same number of walkers per node.
 */
struct DistributedSPOSet: public SPOSetBase
{

  ///use orthogonal and non-truncated orbitals
  typedef TricubicBsplineSPOSet<ValueType,true,false> BspilneSetType;
  typedef PooledBuffer<RealType> BufferType;

  BspineSetType* Phi;
  ///wavefunction
  //SPOSetBase* Phi;
  ///index map
  vector<int> OrbitalIndex;
  ///communicator for a mpi group sharing the wavefunctions
  Communicate* myComm;
  ///node id for the remote nodes
  vector<int> RemoteNodes;
  ///current position
  vector<PosType> Rnow;
  ///orbital indices
  vector<int> OrbitalCount;
  ///orbital indices
  vector<int> OrbitalOffset;
  ///send buffer : maybe better using local variables
  vector<BufferType*> SendBuffer;
  ///recv buffer : maybe better using local variables
  vector<BufferType*> RecvBuffer;

  ///constructor
  DistributedSPOSet(int norbs=0);
  ///destructor
  ~DistributedSPOSet();

  ///set the communicate
  void setCommunicator(Communicate* c);

  void setOrbitalSetSize(int norbs);

  void resetParameters(VarRegistry<RealType>& optVariables);

  void resetTargetParticleSet(ParticleSet& P);

  void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

  void evaluate(const ParticleSet& P, int iat,
                ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi);

  void evaluate_notranspose(const ParticleSet& P, int first, int last,
                            ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet);

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    APP_ABORT("Need specialization of DistributedSPOSet::evaluate_notranspose() for grad_grad_logdet. \n");
  }

  void evaluate_notranspose(const ParticleSet& P, int first, int last
                            , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
    APP_ABORT("Need specialization of DistributedSPOSet::evaluate_notranspose() for grad_grad_grad_logdet. \n");
  }

};
}
#endif
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 1772 $   $Date: 2007-02-17 17:47:37 -0600 (Sat, 17 Feb 2007) $
 * $Id: DistributedSPOSet.h 1772 2007-02-17 23:47:37Z jnkim $
 ***************************************************************************/
