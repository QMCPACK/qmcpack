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

#include "QMCWaveFunctions/TWFFastDerivWrapper.h"
#include "Numerics/DeterminantOperators.h"
#include "type_traits/ConvertToReal.h"
#include <iostream>
namespace qmcplusplus
{
TWFFastDerivWrapper::IndexType TWFFastDerivWrapper::getTWFGroupIndex(const IndexType gid) const
{
  IndexType return_group_index(-1);
  for (IndexType i = 0; i < groups_.size(); i++)
    if (gid == groups_[i])
      return_group_index = i;

  assert(return_group_index != -1);

  return return_group_index;
}

void TWFFastDerivWrapper::addGroup(const ParticleSet& P, const IndexType gid, SPOSet* spo)
{
  if (std::find(groups_.begin(), groups_.end(), gid) == groups_.end())
  {
    groups_.push_back(gid);
    spos_.push_back(spo);
  }
}

void TWFFastDerivWrapper::getM(const ParticleSet& P, std::vector<ValueMatrix>& mvec) const
{
  IndexType ngroups = spos_.size();
  for (IndexType i = 0; i < ngroups; i++)
  {
    const IndexType gid    = groups_[i];
    const IndexType first  = P.first(i);
    const IndexType last   = P.last(i);
    const IndexType nptcls = last - first;
    const IndexType norbs  = spos_[i]->getOrbitalSetSize();
    GradMatrix tmpgmat;
    ValueMatrix tmplmat;
    tmpgmat.resize(nptcls, norbs);
    tmplmat.resize(nptcls, norbs);
    spos_[i]->evaluate_notranspose(P, first, last, mvec[i], tmpgmat, tmplmat);
  }
}

TWFFastDerivWrapper::RealType TWFFastDerivWrapper::evaluateJastrowVGL(const ParticleSet& P,
                                                                      ParticleSet::ParticleGradient& G,
                                                                      ParticleSet::ParticleLaplacian& L) const
{
  WaveFunctionComponent::LogValue logpsi = 0.0;
  G                                      = 0.0;
  L                                      = 0.0;
  for (int i = 0; i < jastrow_list_.size(); ++i)
  {
    logpsi += jastrow_list_[i]->evaluateLog(P, G, L);
  }
  RealType rval = std::real(logpsi);
  return rval;
}

TWFFastDerivWrapper::RealType TWFFastDerivWrapper::evaluateJastrowRatio(ParticleSet& P, const int iel) const
{
  //legacy calls are hit and miss with const.  Remove const for index.
  int iel_(iel);
  WaveFunctionComponent::PsiValue r(1.0);
  for (int i = 0; i < jastrow_list_.size(); ++i)
  {
    r *= jastrow_list_[i]->ratio(P, iel_);
  }

  RealType ratio_return(1.0);
  convertToReal(r, ratio_return);
  return ratio_return;
}

TWFFastDerivWrapper::RealType TWFFastDerivWrapper::calcJastrowRatioGrad(ParticleSet& P,
                                                                        const int iel,
                                                                        GradType& grad) const
{
  int iel_(iel);
  WaveFunctionComponent::PsiValue r(1.0);
  for (int i = 0; i < jastrow_list_.size(); ++i)
  {
    r *= jastrow_list_[i]->ratioGrad(P, iel_, grad);
  }
  RealType ratio_return(1.0);
  convertToReal(r, ratio_return);
  return ratio_return;
}

TWFFastDerivWrapper::GradType TWFFastDerivWrapper::evaluateJastrowGradSource(ParticleSet& P,
                                                                             ParticleSet& source,
                                                                             const int iat) const
{
  GradType grad_iat = GradType();
  for (int i = 0; i < jastrow_list_.size(); ++i)
    grad_iat += jastrow_list_[i]->evalGradSource(P, source, iat);
  return grad_iat;
}

TWFFastDerivWrapper::GradType TWFFastDerivWrapper::evaluateJastrowGradSource(
    ParticleSet& P,
    ParticleSet& source,
    const int iat,
    TinyVector<ParticleSet::ParticleGradient, OHMMS_DIM>& grad_grad,
    TinyVector<ParticleSet::ParticleLaplacian, OHMMS_DIM>& lapl_grad) const
{
  GradType grad_iat = GradType();
  for (int dim = 0; dim < OHMMS_DIM; dim++)
    for (int i = 0; i < grad_grad[0].size(); i++)
    {
      grad_grad[dim][i] = GradType();
      lapl_grad[dim][i] = 0.0;
    }
  for (int i = 0; i < jastrow_list_.size(); ++i)
    grad_iat += jastrow_list_[i]->evalGradSource(P, source, iat, grad_grad, lapl_grad);
  return grad_iat;
}

void TWFFastDerivWrapper::getEGradELaplM(const ParticleSet& P,
                                         std::vector<ValueMatrix>& mvec,
                                         std::vector<GradMatrix>& gmat,
                                         std::vector<ValueMatrix>& lmat) const
{
  IndexType ngroups = mvec.size();
  for (IndexType i = 0; i < ngroups; i++)
  {
    const IndexType gid    = groups_[i];
    const IndexType first  = P.first(i);
    const IndexType last   = P.last(i);
    const IndexType nptcls = last - first;
    const IndexType norbs  = spos_[i]->getOrbitalSetSize();
    spos_[i]->evaluate_notranspose(P, first, last, mvec[i], gmat[i], lmat[i]);
  }
}

void TWFFastDerivWrapper::getIonGradM(const ParticleSet& P,
                                      const ParticleSet& source,
                                      const int iat,
                                      std::vector<std::vector<ValueMatrix>>& dmvec) const
{
  IndexType ngroups = dmvec[0].size();
  for (IndexType i = 0; i < ngroups; i++)
  {
    const IndexType gid    = groups_[i];
    const IndexType first  = P.first(i);
    const IndexType last   = P.last(i);
    const IndexType nptcls = last - first;
    const IndexType norbs  = spos_[i]->getOrbitalSetSize();

    GradMatrix grad_phi;

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

void TWFFastDerivWrapper::getIonGradIonGradELaplM(const ParticleSet& P,
                                                  const ParticleSet& source,
                                                  int iat,
                                                  std::vector<std::vector<ValueMatrix>>& dmvec,
                                                  std::vector<std::vector<GradMatrix>>& dgmat,
                                                  std::vector<std::vector<ValueMatrix>>& dlmat) const
{
  IndexType ngroups = dmvec[0].size();
  for (IndexType i = 0; i < ngroups; i++)
  {
    const IndexType gid    = groups_[i];
    const IndexType first  = P.first(i);
    const IndexType last   = P.last(i);
    const IndexType nptcls = last - first;
    const IndexType norbs  = spos_[i]->getOrbitalSetSize();

    GradMatrix grad_phi;
    HessMatrix grad_grad_phi;
    GradMatrix grad_lapl_phi;

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
          for (IndexType ielec = 0; ielec < OHMMS_DIM; ielec++)
            dgmat[idim][i][iptcl][iorb][ielec] += grad_grad_phi[iptcl][iorb](idim, ielec);
        }
  }
}

TWFFastDerivWrapper::ValueType TWFFastDerivWrapper::computeGSDerivative(const std::vector<ValueMatrix>& Minv,
                                                                        const std::vector<ValueMatrix>& X,
                                                                        const std::vector<ValueMatrix>& dM,
                                                                        const std::vector<ValueMatrix>& dB) const
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

// only doing single species here
void TWFFastDerivWrapper::computeMDDerivative(const std::vector<ValueMatrix>& Minv,
                                              const std::vector<ValueMatrix>& X,
                                              const std::vector<ValueMatrix>& dM,
                                              const std::vector<ValueMatrix>& dB,
                                              const std::vector<ValueMatrix>& B,
                                              const std::vector<ValueMatrix>& Mvirt,
                                              IndexType id;
                                              const std::vector<int>& excdata,
                                              const std::vector<ValueType>& excsign,
                                              std::vector<ValueType>& dvals) const
{
  IndexType ndet = excsign.size();

  for (int idet = 0; idet < ndet; idet++)
  {
    dvals[idet] = 0.0;
  }

  // input mats:
  // Minv   [occ, elec]
  // X      [occ, elec]
  // dM     [elec, occ+virt]
  // dB     [elec, occ+virt]
  // B      [elec, virt]
  // Mvirt  [elec, virt]

  // O is all occ orbs; o is occ orbs that appear as holes in exc. list; e is all elecs; v is virt orbs that appear as particles in exc. list

  // dvals[idet] = dval_ref + tr(
  //     - inv({Minv[o,e].M[e,v]}).{Minv[o,e].dM[e,v] - Minv[o,e].dM[e,O].Minv[O,e].M[e,v]}.inv({Minv[o,e].M[e,v]}).{S} ...
  //     + inv({Minv[o,e].M[e,v]}).{Minv[o,e].dB[e,v] - X[o,e].dM[e,v] - (Minv[o,e].dB[e.O] - X[o,e].dM[e,O]).Minv[O,e].M[e,v] - Minv[o,e].dM[e,O].S[O,v]}
  //   )
  
  // S = (Minv[O,e].B[e,v] - X[O,e].M[e,v]) ("M" from paper)
  
  // X = Minv[O,e].B[e,O].Minv[O,e]

  // a = Minv[O,e].M[e,v] (slice this for alpha from paper)
  
  // mat1b  = Minv[o,e].dM[e,O].Minv[O,e].M[e,v]
  // mat1   = Minv[o,e].dM[e,v] - mat1b[o,v]
  
  // mat2ba = Minv[o,e].dB[e.O] - X[o,e].dM[e,O]
  // mat2b  = mat2ba[o,O].a[O,v]
  // mat2c  = Minv[o,e].dM[e,O].S[O,v]
  // mat2   = Minv[o,e].dB[e,v] - X[o,e].dM[e,v] - mat2b[o,v] - mat2c[o,v]

  //  dvals[idet] = dval_ref + tr(-inv({a}).{mat1}.inv({a}).{S} + inv({a}).{mat2})

  // dval_ref = tr(Minv[O,e].dB[e,O] - X[O,e].dM[e,O])


  /// FIXME: only keep relevant row/col subsets

  // store this and update as single-elec moves are made
  // Minv[occ,elec].M[elec,virt]
  ValueMatrix Minv_Mvirt_ov(nocc, nvirt);

  // Minv,Mvirt are row-major, but gemm assumes col-major layout, so reorder args
  // Minv_Mvirt_ov.T = (Mvirt.T) . (Minv.T)
  BLAS::gemm('n', 'n', nvirt, nocc, nelec, 1.0, Mvirt[id].data(), Mvirt[id].cols(), Minv[id].data(), Minv[id].cols(),
             0.0, Minv_Mvirt_ov.data(), Minv_Mvirt_ov.cols());

  // combine these two, just offset correctly later
  // Minv[occ,elec].dM[elec,occ]
  // ValueMatrix Minv_dM_oo(nocc, nocc);
  // Minv[occ,elec].dM[elec,virt]
  // ValueMatrix Minv_dM_ov(nocc, nvirt);

  ValueMatrix Minv_dM_on(nocc, nocc + nvirt);
  BLAS::gemm('n', 'n', nocc + nvirt, nocc, nelec, 1.0, dM[id].data(), dM[id].cols(), Minv[id].data(), Minv[id].cols(),
             0.0, Minv_dM_on.data(), Minv_dM_on.cols());


  // Minv[occ,elec].B[elec,virt]
  ValueMatrix Minv_B_ov(nocc, nvirt);
  BLAS::gemm('n', 'n', nvirt, nocc, nelec, 1.0, B[id].data(), B[id].cols(), Minv[id].data(), Minv[id].cols(), 0.0,
             Minv_B_ov.data(), Minv_B_ov.cols());

  // X[occ,elec].M[elec,virt]
  ValueMatrix X_Mvirt_ov(nocc, nvirt);
  BLAS::gemm('n', 'n', nvirt, nocc, nelec, 1.0, Mvirt[id].data(), Mvirt[id].cols(), X[id].data(), X[id].cols(), 0.0,
             X_Mvirt_ov.data(), X_Mvirt_ov.cols());


  // combine these two
  // X[occ,elec].dM[elec,occ]
  // ValueMatrix X_dM_oo(nocc, nocc);
  // X[occ,elec].dM[elec,virt]
  // ValueMatrix X_dM_ov(nocc, nvirt);

  ValueMatrix X_dM_on(nocc, nocc + nvirt);
  BLAS::gemm('n', 'n', nocc + nvirt, nocc, nelec, 1.0, dM[id].data(), dM[id].cols(), X[id].data(), X[id].cols(), 0.0,
             X_dM_on.data(), X_dM_on.cols());

  // combine these two
  // Minv[occ,elec].dB[elec,occ]
  // ValueMatrix Minv_dB_oo(nocc, nocc);
  // Minv[occ,elec].dB[elec,virt]
  // ValueMatrix Minv_dB_ov(nocc, nvirt);

  ValueMatrix Minv_dB_on(nocc, nocc + nvirt);
  BLAS::gemm('n', 'n', nocc + nvirt, nocc, nelec, 1.0, dB[id].data(), dB[id].cols(), Minv[id].data(), Minv[id].cols(),
             0.0, Minv_dB_on.data(), Minv_dB_on.cols());


  // S = (Minv_B[v] - X_Mvirt[v]) ("M" from paper)
  ValueMatrix S(nocc, nvirt);
  for (size_t i = 0; i < S.rows(); i++)
    for (size_t j = 0; j < S.cols(); j++)
      S(i, j) = Minv_B_ov(i, j) - X_Mvirt_ov(i, j);

  // Minv_dM[o,o].Minv_Mv[o,v]
  ValueMatrix mat1b(nocc, nvirt);
  // mat1b.T = Minv_Mv.T @ Minv_dM.T
  BLAS::gemm('n', 'n', nvirt, nocc, nocc, 1.0, Minv_Mvirt_ov.data(), Minv_Mvirt_ov.cols(), Minv_dM_on.data(),
             Minv_dM_on.cols(), 0.0, mat1b.data(), mat1b.cols());

  // Minv_dM[o,v] - Minv_dM[o,o].Minv_Mv[o,v]
  ValueMatrix mat1(nocc, nvirt);
  for (size_t i = 0; i < mat1.rows(); i++)
    for (size_t j = 0; j < mat1.cols(); j++)
      mat1(i, j) = Minv_dM_on(i, j) - mat1b(i, j);


  // Minv_dB[o,o]-X_dM[o,o]
  ValueMatrix mat2ba(nocc, nocc);
  for (size_t i = 0; i < mat2ba.rows(); i++)
    for (size_t j = 0; j < mat2ba.cols(); j++)
      mat2ba(i, j) = Minv_dB_on(i, j) - X_dM_on(i, j);

  // mat2ba[o,o].Minv_Mv[o,v]
  ValueMatrix mat2b(nocc, nvirt);
  BLAS::gemm('n', 'n', nvirt, nocc, nocc, 1.0, Minv_Mvirt_ov.data(), Minv_Mvirt_ov.cols(), mat2ba.data(), mat2ba.cols(),
             0.0, mat2b.data(), mat2b.cols());


  // Minv_dM[o,o].S[o,v]
  ValueMatrix mat2c(nocc, nvirt);
  BLAS::gemm('n', 'n', nvirt, nocc, nocc, 1.0, S.data(), S.cols(), Minv_dM_on.data(), Minv_dM_on.cols(), 0.0,
             mat2c.data(), mat2c.cols());


  // Minv_dB - X_dM - mat2b - mat2c
  ValueMatrix mat2(nocc, nvirt);
  for (size_t i = 0; i < mat2.rows(); i++)
    for (size_t j = 0; j < mat2.cols(); j++)
      mat2(i, j) = Minv_dB_on(i, j) - X_dM_on(i, j) - mat2b(i, j) - mat2c(i, j);


  // tr(Minv_dB - X_dM)
  ValueType dval_refdet = 0.0;
  for (size_t i = 0; i < mat2ba.rows(); i++)
  {
    dval_refdet += mat2ba(i, i);
  }


  // TODO: Eq. 43
  // dvals[idet] = dval_refdet + tr(-inv(Minv_M[P,Q]).mat1[P,Q].inv(Minv_M[P,Q]).S[P,Q] + inv(Minv_M[P,Q]).mat2[P,Q])
  // build full mats for all required pairs, then select submatrices for each excited det
  // use exc index data to map into full arrays and create [k,k] tables


  size_t det_offset  = 1;
  size_t data_offset = 1;

  auto update_offsets = [&](size_t ext_level) {
    det_offset += (*ndets_per_excitation_level_)[ext_level];
    data_offset += (*ndets_per_excitation_level_)[ext_level] * (3 * ext_level + 1);
  };

  for (size_t iexc_level = 1; iexc_level <= max_ext_level; iexc_level++)
  {
    // calc vals
    
    // see MultiDiracDeterminant::mw_updateRatios, MultiDiracDeterminant::mw_buildTableMatrix_calculateRatios_impl
    // for examples of excited det idx layout

    // get hole/particle idx
    // take corresponding rows cols of Minv_M, S, mat1, mat2
    // a = inv(Minv_M[P,Q])
    // dvals[idet] = dval_refdet + tr(-a.mat1[P,Q].a.S[P,Q] + a.S[P,Q])
    // dvals[idet] = dval_refdet + tr(a.(I - mat1[P,Q].a).S[P,Q])
    // dvals[idet] = dval_refdet + tr((I - a.mat1[P,Q]).a.S[P,Q])

    // a.S[P,Q] also needed for Eq. 29


    // update offsets
    update_offsets(iexc_level);
  }

  return;
}

void TWFFastDerivWrapper::invertMatrices(const std::vector<ValueMatrix>& M, std::vector<ValueMatrix>& Minv)
{
  IndexType nspecies = M.size();
  for (IndexType id = 0; id < nspecies; id++)
  {
    assert(M[id].cols() == M[id].rows());
    Minv[id] = M[id];
    invert_matrix(Minv[id]);
  }
}

void TWFFastDerivWrapper::buildX(const std::vector<ValueMatrix>& Minv,
                                 const std::vector<ValueMatrix>& B,
                                 std::vector<ValueMatrix>& X)
{
  IndexType nspecies = Minv.size();

  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum = Minv[id].rows();
    assert(Minv[id].rows() == Minv[id].cols());
    ValueMatrix tmpmat;
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

void TWFFastDerivWrapper::wipeMatrices(std::vector<ValueMatrix>& A)
{
  for (IndexType id = 0; id < A.size(); id++)
  {
    A[id] = 0.0;
  }
}

TWFFastDerivWrapper::ValueType TWFFastDerivWrapper::trAB(const std::vector<ValueMatrix>& A,
                                                         const std::vector<ValueMatrix>& B)
{
  IndexType nspecies = A.size();
  assert(A.size() == B.size());
  ValueType val = 0.0;
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

void TWFFastDerivWrapper::getGSMatrices(const std::vector<ValueMatrix>& A, std::vector<ValueMatrix>& Aslice) const
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


TWFFastDerivWrapper::IndexType TWFFastDerivWrapper::getRowM(const ParticleSet& P,
                                                            const IndexType iel,
                                                            ValueVector& val) const
{
  IndexType gid = P.getGroupID(iel);
  IndexType sid = getTWFGroupIndex(gid);

  GradVector tempg;
  ValueVector templ;

  IndexType norbs = spos_[sid]->getOrbitalSetSize();

  tempg.resize(norbs);
  templ.resize(norbs);

  spos_[sid]->evaluateVGL(P, iel, val, tempg, templ);

  return sid;
}

} // namespace qmcplusplus
