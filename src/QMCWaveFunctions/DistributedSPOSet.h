//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
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
  std::vector<int> OrbitalIndex;
  ///communicator for a mpi group sharing the wavefunctions
  Communicate* myComm;
  ///node id for the remote nodes
  std::vector<int> RemoteNodes;
  ///current position
  std::vector<PosType> Rnow;
  ///orbital indices
  std::vector<int> OrbitalCount;
  ///orbital indices
  std::vector<int> OrbitalOffset;
  ///send buffer : maybe better using local variables
  std::vector<BufferType*> SendBuffer;
  ///recv buffer : maybe better using local variables
  std::vector<BufferType*> RecvBuffer;

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
