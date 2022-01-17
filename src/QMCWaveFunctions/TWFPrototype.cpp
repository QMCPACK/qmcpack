//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//
// File created by:   Raymond Clay III, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/TWFPrototype.h"
#include "Numerics/DeterminantOperators.h"
#include <iostream>
namespace qmcplusplus
{

TWFPrototype::IndexType TWFPrototype::getTWFGroupIndex(const IndexType gid)
{
  IndexType return_group_index(-1);
  for (IndexType i = 0; i < groups_.size(); i++)
    if (gid == groups_[i])
      return_group_index=i;
  
  assert(return_group_index != -1);

  return return_group_index;
}

void TWFPrototype::addGroup(const ParticleSet& P, const IndexType gid, SPOSet* spo)
{
  if (std::find(groups_.begin(), groups_.end(), gid) == groups_.end())
  {
    groups_.push_back(gid);
    spos_.push_back(spo);
    IndexType first = P.first(gid);
    IndexType last  = P.last(gid);
    IndexType norbs = spo->getOrbitalSetSize();
  }
}

void TWFPrototype::getM(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec)
{
  IndexType ndets  = spos_.size();
  IndexType norbs  = 0;
  IndexType nptcls = 0;
  IndexType gid    = 0;
  IndexType first  = 0;
  IndexType last   = 0;
  for (IndexType i = 0; i < ndets; i++)
  {
    gid     = groups_[i];
    first   = P.first(i);
    last    = P.last(i);
    nptcls  = last - first;
    norbs   = spos_[i]->getOrbitalSetSize();
    mvec[i] = 0;
    GradMatrix_t tmpgmat;
    ValueMatrix_t tmplmat;
    tmpgmat.resize(nptcls, norbs);
    tmplmat.resize(nptcls, norbs);
    spos_[i]->evaluate_notranspose(P, first, last, mvec[i], tmpgmat, tmplmat);
  }
}

void TWFPrototype::getEGradELaplM(const ParticleSet& P,
                                     std::vector<ValueMatrix_t>& mvec,
                                     std::vector<GradMatrix_t>& gmat,
                                     std::vector<ValueMatrix_t>& lmat)
{
  IndexType ndets  = mvec.size();
  IndexType norbs  = 0;
  IndexType nptcls = 0;
  IndexType gid    = 0;
  IndexType first  = 0;
  IndexType last   = 0;
  for (IndexType i = 0; i < ndets; i++)
  {
    gid     = groups_[i];
    first   = P.first(i);
    last    = P.last(i);
    nptcls  = last - first;
    norbs   = spos_[i]->getOrbitalSetSize();
    mvec[i] = 0;
    gmat[i] = 0;
    lmat[i] = 0;
    spos_[i]->evaluate_notranspose(P, first, last, mvec[i], gmat[i], lmat[i]);
  }
}

void TWFPrototype::getIonGradM(const ParticleSet& P,
                               const ParticleSet& source,
                               const int iat,
                               std::vector<std::vector<ValueMatrix_t>>& dmvec)
{
  IndexType ndets  = dmvec[0].size();
  IndexType norbs  = 0;
  IndexType nptcls = 0;
  IndexType gid    = 0;
  IndexType first  = 0;
  IndexType last   = 0;
  for (IndexType i = 0; i < ndets; i++)
  {
    gid    = groups_[i];
    first  = P.first(i);
    last   = P.last(i);
    nptcls = last - first;
    norbs  = spos_[i]->getOrbitalSetSize();

    GradMatrix_t grad_phi;

    grad_phi.resize(nptcls, norbs);

    spos_[i]->evaluateGradSource(P, first, last, source, iat, grad_phi);

    for (IndexType idim = 0; idim < OHMMS_DIM; idim++)
      for (IndexType iptcl = 0; iptcl < nptcls; iptcl++)
        for (IndexType iorb = 0; iorb < norbs; iorb++)
        {
          dmvec[idim][i][iptcl][iorb] += grad_phi[iptcl][iorb][idim];
        }
  }
}

void TWFPrototype::getIonGradIonGradELaplM(const ParticleSet& P,
                                          const ParticleSet& source,
                                          int iat,
                                          std::vector<std::vector<ValueMatrix_t>>& dmvec,
                                          std::vector<std::vector<ValueMatrix_t>>& dlmat)
{
  IndexType ndets  = dmvec[0].size();
  IndexType norbs  = 0;
  IndexType nptcls = 0;
  IndexType gid    = 0;
  IndexType first  = 0;
  IndexType last   = 0;
  for (IndexType i = 0; i < ndets; i++)
  {
    gid    = groups_[i];
    first  = P.first(i);
    last   = P.last(i);
    nptcls = last - first;
    norbs  = spos_[i]->getOrbitalSetSize();

    GradMatrix_t grad_phi;
    HessMatrix_t grad_grad_phi;
    GradMatrix_t grad_lapl_phi;

    grad_phi.resize(nptcls, norbs);
    grad_grad_phi.resize(nptcls, norbs);
    grad_lapl_phi.resize(nptcls, norbs);

    spos_[i]->evaluateGradSource(P, first, last, source, iat, grad_phi, grad_grad_phi, grad_lapl_phi);

    for (IndexType idim = 0; idim < OHMMS_DIM; idim++)
      for (IndexType iptcl = 0; iptcl < nptcls; iptcl++)
        for (IndexType iorb = 0; iorb < norbs; iorb++)
        {
          dmvec[idim][i][iptcl][iorb] += grad_phi[iptcl][iorb][idim];
          dlmat[idim][i][iptcl][iorb] += grad_lapl_phi[iptcl][iorb][idim];
        }
  }
}

TWFPrototype::ValueType TWFPrototype::computeGSDerivative(const std::vector<ValueMatrix_t>& Minv,
                                                            const std::vector<ValueMatrix_t>& X,
                                                            const std::vector<ValueMatrix_t>& dM,
                                                            const std::vector<ValueMatrix_t>& dB)
{
  IndexType nspecies = Minv.size();
  ValueType dval     = 0.0;
  for (int id = 0; id < nspecies; id++)
  {
    int ptclnum       = Minv[id].rows();
    ValueType dval_id = 0.0;
    for (int i = 0; i < ptclnum; i++)
      for (int j = 0; j < ptclnum; j++)
      {
        //Tr[M^{-1} dB - X * dM ]
        dval_id += Minv[id][i][j] * dB[id][j][i] - X[id][i][j] * dM[id][j][i];
      }
    dval += dval_id;
  }
  return dval;
}

void TWFPrototype::invertMatrix(const std::vector<ValueMatrix_t>& M, std::vector<ValueMatrix_t>& Minv)
{
  IndexType nspecies = M.size();
  for (IndexType id = 0; id < nspecies; id++)
  {
    assert(M[id].cols() == M[id].rows());
    Minv[id] = M[id];
    invert_matrix(Minv[id]);
  }
}

void TWFPrototype::buildX(const std::vector<ValueMatrix_t>& Minv,
                           const std::vector<ValueMatrix_t>& B,
                           std::vector<ValueMatrix_t>& X)
{
  IndexType nspecies = Minv.size();

  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum = Minv[id].rows();
    assert(Minv[id].rows()==Minv[id].cols());
    ValueMatrix_t tmpmat;
    X[id].resize(ptclnum, ptclnum);
    tmpmat.resize(ptclnum, ptclnum);
    //(B*A^-1)
    for (int i = 0; i < ptclnum; i++)
      for (int j = 0; j < ptclnum; j++)
        for (int k = 0; k < ptclnum; k++)
        {
          tmpmat[i][j] += B[id][i][k] * Minv[id][k][j];
        }
    //A^{-1}*B*A^{-1}
    for (int i = 0; i < ptclnum; i++)
      for (int j = 0; j < ptclnum; j++)
        for (int k = 0; k < ptclnum; k++)
        {
          X[id][i][j] += Minv[id][i][k] * tmpmat[k][j];
        }
  }
}

void TWFPrototype::wipeMatrix(std::vector<ValueMatrix_t>& A)
{

  for (IndexType id = 0; id < A.size(); id++)
  {
    A[id] = 0.0;
  }
}

TWFPrototype::ValueType TWFPrototype::trAB(const std::vector<ValueMatrix_t>& A, const std::vector<ValueMatrix_t>& B)
{
  IndexType nspecies = A.size();
  assert(A.size() == B.size());
  ValueType val      = 0.0;
  //Now to compute the kinetic energy
  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum      = A[id].rows();
    ValueType val_id = 0.0;
    assert(A[id].cols() == B[id].rows() && A[id].rows() == B[id].cols());
    for (int i = 0; i < A[id].rows(); i++)
      for (int j = 0; j < A[id].cols(); j++)
      {
        val_id += A[id][i][j] * B[id][j][i];
      }
    val += val_id;
  }

  return val;
}

void TWFPrototype::getGSMatrix(const std::vector<ValueMatrix_t>& A, std::vector<ValueMatrix_t>& Aslice)
{
  IndexType nspecies = A.size();
  Aslice.resize(nspecies);
  for (IndexType id = 0; id < nspecies; id++)
  {
    IndexType ptclnum = A[id].rows();
    Aslice[id].resize(ptclnum, ptclnum);
    for (IndexType i = 0; i < ptclnum; i++)
      for (IndexType j = 0; j < ptclnum; j++)
        Aslice[id][i][j] = A[id][i][j];
  }
}


TWFPrototype::IndexType TWFPrototype::getRowM(const ParticleSet& P, const IndexType iel, ValueVector_t& val)
{
  IndexType gid      = P.getGroupID(iel);
  IndexType sid = getTWFGroupIndex(gid);

  GradVector_t tempg;
  ValueVector_t templ;

  IndexType norbs = spos_[sid]->getOrbitalSetSize();

  tempg.resize(norbs);
  templ.resize(norbs);

  spos_[sid]->evaluateVGL(P, iel, val, tempg, templ);

  return sid;
}


} // namespace qmcplusplus
