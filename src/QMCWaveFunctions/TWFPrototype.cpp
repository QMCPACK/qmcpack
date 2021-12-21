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
TWFPrototype::TWFPrototype() : initialized(false) { std::cout << "TWFPrototype Constructed\n"; }

TWFPrototype::IndexType TWFPrototype::get_group_index(const IndexType gid)
{
  IndexType ndets = groups.size();
  for (IndexType i = 0; i < ndets; i++)
    if (gid == groups[i])
      return i;

  return IndexType(-1);
}
void TWFPrototype::add_determinant(const ParticleSet& P, const IndexType gid, SPOSet* spo)
{
  if (std::find(groups.begin(), groups.end(), gid) == groups.end())
  {
    groups.push_back(gid);
    spos.push_back(spo);
    IndexType first = P.first(gid);
    IndexType last  = P.last(gid);
    IndexType norbs = spo->getOrbitalSetSize();
    num_orbs.push_back(norbs);
    num_ptcls.push_back(last - first);
  }
  initialized = true;
}

TWFPrototype::IndexType TWFPrototype::get_det_id(const IndexType species_id)
{
  IndexType detIndex = -1;
  for (IndexType i = 0; i < groups.size(); i++)
  {
    if (species_id == groups[i])
      detIndex = i;
  }
  assert(detIndex != -1);

  return detIndex;
}

void TWFPrototype::get_M(const ParticleSet& P, std::vector<ValueMatrix_t>& mvec)
{
  IndexType ndets  = spos.size();
  IndexType norbs  = 0;
  IndexType nptcls = 0;
  IndexType gid    = 0;
  IndexType first  = 0;
  IndexType last   = 0;
  for (IndexType i = 0; i < ndets; i++)
  {
    gid     = groups[i];
    first   = P.first(i);
    last    = P.last(i);
    nptcls  = last - first;
    norbs   = spos[i]->getOrbitalSetSize();
    mvec[i] = 0;
    GradMatrix_t tmpgmat;
    ValueMatrix_t tmplmat;
    tmpgmat.resize(nptcls, norbs);
    tmplmat.resize(nptcls, norbs);
    spos[i]->evaluate_notranspose(P, first, last, mvec[i], tmpgmat, tmplmat);
  }
}

void TWFPrototype::get_egrad_elapl_M(const ParticleSet& P,
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
    gid     = groups[i];
    first   = P.first(i);
    last    = P.last(i);
    nptcls  = last - first;
    norbs   = spos[i]->getOrbitalSetSize();
    mvec[i] = 0;
    gmat[i] = 0;
    lmat[i] = 0;
    spos[i]->evaluate_notranspose(P, first, last, mvec[i], gmat[i], lmat[i]);
  }
}

void TWFPrototype::get_igrad_M(const ParticleSet& P,
                               const ParticleSet& source,
                               int iat,
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
    gid    = groups[i];
    first  = P.first(i);
    last   = P.last(i);
    nptcls = last - first;
    norbs  = spos[i]->getOrbitalSetSize();

    GradMatrix_t grad_phi;

    grad_phi.resize(nptcls, norbs);

    spos[i]->evaluateGradSource(P, first, last, source, iat, grad_phi);

    for (IndexType idim = 0; idim < OHMMS_DIM; idim++)
      for (IndexType iptcl = 0; iptcl < nptcls; iptcl++)
        for (IndexType iorb = 0; iorb < norbs; iorb++)
        {
          dmvec[idim][i][iptcl][iorb] += grad_phi[iptcl][iorb][idim];
        }
  }
}

void TWFPrototype::get_igrad_igradelapl_M(const ParticleSet& P,
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
    gid    = groups[i];
    first  = P.first(i);
    last   = P.last(i);
    nptcls = last - first;
    norbs  = spos[i]->getOrbitalSetSize();

    GradMatrix_t grad_phi;
    HessMatrix_t grad_grad_phi;
    GradMatrix_t grad_lapl_phi;

    grad_phi.resize(nptcls, norbs);
    grad_grad_phi.resize(nptcls, norbs);
    grad_lapl_phi.resize(nptcls, norbs);

    spos[i]->evaluateGradSource(P, first, last, source, iat, grad_phi, grad_grad_phi, grad_lapl_phi);

    for (IndexType idim = 0; idim < OHMMS_DIM; idim++)
      for (IndexType iptcl = 0; iptcl < nptcls; iptcl++)
        for (IndexType iorb = 0; iorb < norbs; iorb++)
        {
          dmvec[idim][i][iptcl][iorb] += grad_phi[iptcl][iorb][idim];
          dlmat[idim][i][iptcl][iorb] += grad_lapl_phi[iptcl][iorb][idim];
        }
  }
}

TWFPrototype::ValueType TWFPrototype::compute_gs_derivative(const std::vector<ValueMatrix_t>& Minv,
                                                            const std::vector<ValueMatrix_t>& X,
                                                            const std::vector<ValueMatrix_t>& dM,
                                                            const std::vector<ValueMatrix_t>& dB)
{
  IndexType nspecies = num_species();
  ValueType dval     = 0.0;
  for (int id = 0; id < nspecies; id++)
  {
    int ptclnum       = num_particles(id);
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

void TWFPrototype::invert_M(const std::vector<ValueMatrix_t>& M, std::vector<ValueMatrix_t>& Minv)
{
  IndexType nspecies = num_species();
  for (IndexType id = 0; id < nspecies; id++)
  {
    assert(M[id].cols() == M[id].rows());
    Minv[id] = M[id];
    invert_matrix(Minv[id]);
  }
}

void TWFPrototype::build_X(const std::vector<ValueMatrix_t>& Minv,
                           const std::vector<ValueMatrix_t>& B,
                           std::vector<ValueMatrix_t>& X)
{
  IndexType nspecies = num_species();

  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum = num_particles(id);
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

void TWFPrototype::wipe_matrix(std::vector<ValueMatrix_t>& A)
{
  IndexType nspecies = num_species();

  for (IndexType id = 0; id < nspecies; id++)
  {
    A[id] = 0.0;
  }
}

TWFPrototype::ValueType TWFPrototype::trAB(const std::vector<ValueMatrix_t>& A, const std::vector<ValueMatrix_t>& B)
{
  IndexType nspecies = num_species();
  ValueType val      = 0.0;
  //Now to compute the kinetic energy
  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum      = num_particles(id);
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

void TWFPrototype::get_gs_matrix(const std::vector<ValueMatrix_t>& A, std::vector<ValueMatrix_t>& Aslice)
{
  IndexType nspecies = num_species();
  Aslice.resize(nspecies);
  for (IndexType id = 0; id < nspecies; id++)
  {
    IndexType ptclnum = num_particles(id);
    Aslice[id].resize(ptclnum, ptclnum);
    for (IndexType i = 0; i < ptclnum; i++)
      for (IndexType j = 0; j < ptclnum; j++)
        Aslice[id][i][j] = A[id][i][j];
  }
}

TWFPrototype::IndexType TWFPrototype::get_igrad_row(const ParticleSet& P,
                                                    const ParticleSet& source,
                                                    IndexType iel,
                                                    IndexType iat_source,
                                                    std::vector<ValueVector_t>& dval)
{
  return -1;
}

TWFPrototype::IndexType TWFPrototype::get_M_row(const ParticleSet& P, IndexType iel, ValueVector_t& val)
{
  IndexType gid      = P.getGroupID(iel);
  IndexType detIndex = get_det_id(gid);

  GradVector_t tempg;
  ValueVector_t templ;

  IndexType norbs = spos[detIndex]->getOrbitalSetSize();

  tempg.resize(norbs);
  templ.resize(norbs);

  spos[detIndex]->evaluateVGL(P, iel, val, tempg, templ);

  return detIndex;
}


TWFPrototype::RealType TWFPrototype::evaluateLog(ParticleSet& P) { return 0; }

} // namespace qmcplusplus
