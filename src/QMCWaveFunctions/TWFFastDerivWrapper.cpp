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
void TWFFastDerivWrapper::addMultiSlaterDet(const ParticleSet& P, MultiSlaterDetTableMethod* msd)
{
  // add msd to slaterdet list
  slaterdets_.push_back(msd);
  std::vector<IndexType> groupids(msd->Dets.size(), 0);

  /// NOTE: we could call `addGroup` for the msd->Dets here. The singledet version does that by
  /// registering the diracdets when registerTWFFastDerivWrapper is called, so this is consistent with that behavior
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

void TWFFastDerivWrapper::computeMDDerivative(const std::vector<ValueMatrix>& Minv,
                                              const std::vector<ValueMatrix>& X,
                                              const std::vector<ValueMatrix>& dM,
                                              const std::vector<ValueMatrix>& dB,
                                              const std::vector<ValueMatrix>& B,
                                              const std::vector<ValueMatrix>& M,
                                              const std::vector<IndexType>& mdd_spo_ids,
                                              const std::vector<MultiDiracDeterminant*>& mdds,
                                              std::vector<ValueType>& dvals) const
{
  // mdd_id is multidiracdet id
  // sid is sposet id (index into first dim of M, X, B, etc.)
  for (size_t mdd_id = 0; mdd_id < mdd_spo_ids.size(); mdd_id++)
  {
    IndexType sid = mdd_spo_ids[mdd_id];

    // MultiDiracDeterminant* multidiracdet_i = mdds[mdd_id];
    const MultiDiracDeterminant& multidiracdet_i = *mdds[mdd_id];

    // number of DiracDets in this MultiDiracDet
    // IndexType ndet = multidiracdet_i->getNumDets();
    IndexType ndet = multidiracdet_i.getNumDets();

    /// TODO: just do this before calling
    for (int idet = 0; idet < ndet; idet++)
    {
      dvals[mdd_id][idet] = 0.0;
    }

    // input mats:
    // Minv   [occ, elec]
    // X      [occ, elec]
    // dM     [elec, occ+virt]
    // dB     [elec, occ+virt]
    // B      [elec, virt]
    // M      [elec, occ+virt]
    // Mvirt  [elec, virt] (subset of M)
    // Mvirt.data() = M.data() + nOcc (ldim is occ+virt)

    // O is all occ orbs; o is occ orbs that appear as holes in exc. list; e is all elecs; v is virt orbs that appear as particles in exc. list

    // dvals[idet] = dval_ref + tr(
    //     - inv({Minv[o,e].M[e,v]}).{Minv[o,e].dM[e,v] - Minv[o,e].dM[e,O].Minv[O,e].M[e,v]}.inv({Minv[o,e].M[e,v]}).{S} ...
    //     + inv({Minv[o,e].M[e,v]}).{Minv[o,e].dB[e,v] - X[o,e].dM[e,v] - (Minv[o,e].dB[e.O] - X[o,e].dM[e,O]).Minv[O,e].M[e,v] - Minv[o,e].dM[e,O].S[O,v]}
    //   )


    // X = Minv[O,e].B[e,O].Minv[O,e]

    // a = Minv[O,e].M[e,v] (slice this for alpha from paper)
    // b = Minv[o,e].dM[e,Ov]
    // c = Minv[o,e].dB[e,Ov] (first dim O for ref det)
    // d = X[o,e].dM[e,Ov]    (first dim O for ref det)
    // f = Minv[O,e].B[e,v]
    // g = X[O,e].M[e,v]
    // h = c[o,Ov] - d[o,Ov]  (refdet is tr(h[O,O]) )

    // S = (Minv[O,e].B[e,v] - X[O,e].M[e,v]) ("M" from paper)
    //   = Minv[O,e].(B[e,v] - B[e,O].Minv[O,e].M[e,v])
    //   = f[O,v] - g[O,v]

    // mat1b  = Minv[o,e].dM[e,O].Minv[O,e].M[e,v]
    //        = Minv[o,e].dM[e,O].a[O,v]
    //        = b[o,O].a[O,v]

    // mat1   = Minv[o,e].dM[e,v] - mat1b[o,v]
    //        = b[o,v] - mat1b[o,v]

    /// NOTE: also need c[o,v] - d[o,v] for mat2, so can construct (c-d)[o,Ov] and slice later
    // mat2ba = Minv[o,e].dB[e,O] - X[o,e].dM[e,O]
    //        = c[o,O] - d[o,O]
    //        = h[o,O]

    // mat2b  = (Minv[o,e].dB[e,O] - X[o,e].dM[e,O]).a[O,v]
    //        = h[o,O].a[O,v]

    // mat2c  = Minv[o,e].dM[e,O].S[O,v]
    //        = b[o,O].S[O,v]

    // mat2   = Minv[o,e].dB[e,v] - X[o,e].dM[e,v] - mat2b[o,v] - mat2c[o,v]
    //        = c[o,v] - d[o,v] - mat2b[o,v] - mat2c[o,v]
    //        = (c-d)[o,v] - (c-d)[o,O].a[O,v] - mat2c[o,v]
    //        = h[o,v] - h[o,O].a[O,v] - mat2c[o,v]

    //  dvals[idet] = dval_ref + tr(-inv({a}).{mat1}.inv({a}).{S} + inv({a}).{mat2})

    // dval_ref = tr(Minv[O,e].dB[e,O] - X[O,e].dM[e,O])
    //          = tr(c[O,O] - d[O,O])


    /// FIXME: set these correctly; decide how to index into nocc if not contiguous within nOcc
    size_t nOcc;               // total occ orbs in refdet (should be same as n_elec)
    size_t nocc;               // number of orbs that appear as holes in exc. list
    size_t nvirt;              // should be number of virtuals that appear as particles in exc. list
    size_t virt_offset = nOcc; // first col of virt orbs in M

    /// NOTE: gemms below are reversed from what is written in comments
    // a_Ov = Minv.Mvirt
    // a_Ov.T = (Mvirt.T) . (Minv.T)
    // Minv,M are row-major, but gemm assumes col-major layout, so reorder args
    // if we want output to be col-major, we can reverse inputs and use 't','t'

    /// FIXME: only keep relevant row/col subsets

    /// TODO: store this and update as single-elec moves are made
    // a[O,v] = Minv[O,e].M[e,v]
    ValueMatrix a_Ov(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nelec, 1.0, M[sid].data() + virt_offset, M[sid].cols(), Minv[sid].data(),
               Minv[sid].cols(), 0.0, a_Ov.data(), a_Ov.cols());


    /// FIXME: we only need o for first dim
    // b[O,Ov] = Minv[O,e].dM[e,Ov]
    ValueMatrix b_On(nOcc, nOcc + nvirt);
    BLAS::gemm('n', 'n', nOcc + nvirt, nOcc, nelec, 1.0, dM[sid].data(), dM[sid].cols(), Minv[sid].data(),
               Minv[sid].cols(), 0.0, b_On.data(), b_On.cols());


    /// FIXME: for MD term, we only need o for first dim (for refdet, we need full O first dim)
    // c[O,Ov] = Minv[O,e].dB[elec,Ov]
    ValueMatrix c_On(nOcc, nOcc + nvirt);
    BLAS::gemm('n', 'n', nOcc + nvirt, nOcc, nelec, 1.0, dB[sid].data(), dB[sid].cols(), Minv[sid].data(),
               Minv[sid].cols(), 0.0, c_On.data(), c_On.cols());

    /// FIXME: for MD term, we only need o for first dim (for refdet, we need full O first dim)
    // d[O,Ov] = X[O,e].dM[e,Ov]
    ValueMatrix d_On(nOcc, nOcc + nvirt);
    BLAS::gemm('n', 'n', nOcc + nvirt, nOcc, nelec, 1.0, dM[sid].data(), dM[sid].cols(), X[sid].data(), X[sid].cols(),
               0.0, d_On.data(), d_On.cols());

    // f[O,v] = Minv[O,e].B[e,v]
    ValueMatrix f_Ov(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nelec, 1.0, B[sid].data() + virt_offset, B[sid].cols(), Minv[sid].data(),
               Minv[sid].cols(), 0.0, f_Ov.data(), f_Ov.cols());

    // g[O,v] = X[O,e].M[e,v]
    ValueMatrix g_Ov(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nelec, 1.0, M[sid].data() + virt_offset, M[sid].cols(), X[sid].data(),
               X[sid].cols(), 0.0, g_Ov.data(), g_Ov.cols());

    /// FIXME: for MD term, we only need o for first dim (for refdet, we need full O first dim)
    // h[O,Ov] = c[O,Ov] - d[O,Ov]
    ValueMatrix h_On(nOcc, nOcc + nvirt);
    for (size_t i = 0; i < h_On.rows(); i++)
      for (size_t j = 0; j < h_On.cols(); j++)
        h_On(i, j) = c_On(i, j) - d_On(i, j);

    // S = (Minv[O,e].B[e,v] - X[O,e].M[e,v]) ("M" from paper)
    //   = f[O,v] - g[O,v]
    ValueMatrix S_Ov(nOcc, nvirt);
    for (size_t i = 0; i < S_Ov.rows(); i++)
      for (size_t j = 0; j < S_Ov.cols(); j++)
        S_Ov(i, j) = f_Ov(i, j) - g_Ov(i, j);


    /// FIXME: only need o for first dim
    /// NOTE: second dim of b is Ov
    // mat1b[O,v] = b[O,O].a[O,v]
    ValueMatrix mat1b(nocc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nOcc, 1.0, a_Ov.data(), a_Ov.cols(), b_On.data(), b_On.cols(), 0.0, mat1b.data(),
               mat1b.cols());

    /// FIXME: only need o for first dim
    /// NOTE: second dim of b is Ov
    //  mat1[O,v] = b[O,v] - mat1b[O,v]
    ValueMatrix mat1(nOcc, nvirt);
    for (size_t i = 0; i < mat1.rows(); i++)
      for (size_t j = 0; j < mat1.cols(); j++)
        mat1(i, j) = b_On(i, j + virt_offset) - mat1b(i, j);


    /// FIXME: only need o for first dim
    // mat2b[O,v] = h[O,O].a[O,v]
    ValueMatrix mat2b(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nOcc, 1.0, a_Ov.data(), a_Ov.cols(), h_On.data(), h_On.cols(), 0.0, mat2b.data(),
               mat2b.cols());


    /// FIXME: only need o for first dim
    // mat2c[O,v] = b[O,O].S[O,v]
    ValueMatrix mat2c(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nOcc, 1.0, S_Ov.data(), S_Ov.cols(), b_On.data(), b_On.cols(), 0.0, mat2c.data(),
               mat2c.cols());

    /// FIXME: only need o for first dim
    // mat2[O,v] = h[O,v] - mat2b[O,v] - mat2c[O,v]
    ValueMatrix mat2(nOcc, nvirt);
    for (size_t i = 0; i < mat2.rows(); i++)
      for (size_t j = 0; j < mat2.cols(); j++)
        mat2(i, j) = h_On(i, j + virt_offset) - mat2b(i, j) - mat2c(i, j);

    // tr(Minv_dB - X_dM)
    ValueType dval_refdet = 0.0;
    for (size_t i = 0; i < h_On.rows(); i++)
    {
      dval_refdet += h_On(i, i);
    }

    dvals[mdd_id][0] = dval_refdet;


    // TODO: Eq. 43
    // {} signify submatrix corresponding to holes/particles rows/cols
    //
    // dvals[idet] = dval_refdet + tr(-inv({a}).{mat1}.inv({a}).{S} + inv({a}).{mat2})
    // build full mats for all required pairs, then select submatrices for each excited det
    // use exc index data to map into full arrays and create [k,k] tables


    size_t det_offset  = 1;
    size_t data_offset = 1;

    const std::vector<int>& ndets_per_exc_lvl = *multidiracdet_i->ndets_per_excitation_level_;
    const OffloadVector<int>& excdata         = *multidiracdet_i->detData;
    const OffloadVector<RealType>& det_signs  = *multidiracdet_i->DetSigns;

    auto update_offsets = [&](size_t ext_level) {
      det_offset += ndets_per_exc_lvl[ext_level];
      data_offset += ndets_per_exc_lvl[ext_level] * (3 * ext_level + 1);
    };

    // dets ordered as ref, all singles, all doubles, all triples, etc.
    //      shift by ndet_per_exc_level after each exc_level
    // data is ordered in same way, but each det has several contiguous elements
    //      each det of exc_level k has (3*k+1) elements: [k, {pos}, {ocp}, {uno}]
    //      k: exc level
    //      pos : positions of holes in refdet
    //      ocp : MO idx of holes
    //      uno : MO idx of particles

    for (size_t k = 1; k <= max_ext_level; k++)
    {
      size_t ndet_exc = ndets_per_exc_lvl[k];
      // calc vals

      // see MultiDiracDeterminant::mw_updateRatios, MultiDiracDeterminant::mw_buildTableMatrix_calculateRatios_impl
      // for examples of excited det idx layout

      // get hole/particle idx
      // take corresponding rows/cols of a, S, mat1, mat2

      // a' = inv({a})
      // dvals[idet] = dval_refdet + tr(-a'.{mat1}.a'.{S} + a'.{mat2})
      // dvals[idet] = dval_refdet + tr(-a'.({mat1}.a'.{S} - {mat2}))

      // a'.{S} also needed for Eq. 29

      for (size_t idet = 0; idet < ndet_exc; idet++)
      {
        std::vector<IndexType> hlist(k, 0);
        std::vector<IndexType> plist(k, 0);
        ValueMatrix a(k, k);
        ValueMatrix m1(k, k);
        ValueMatrix m2(k, k);
        ValueMatrix s(k, k);
        size_t excdata_offset = data_offset + idet * (3 * k + 1);
        for (size_t i = 0; i < k; i++)
        {
          hlist[i] = excdata[excdata_offset + k + 1 + i];
          plist[i] = excdata[excdata_offset + 2 * k + 1 + i];
        }

        // construct submatrices
        for (size_t i = 0; i < k; i++)
          for (size_t j = 0; j < k; j++)
          {
            a(i, j)  = a_Ov(hlist[i], plist[j]);
            m1(i, j) = mat1(hlist[i], plist[j]);
            m2(i, j) = mat2(hlist[i], plist[j]);
            s(i, j)  = S_Ov(hlist[i], plist[j]);
          }

        ValueType dval;
        /// TODO: compute dval = tr(-ainv.m1.ainv.s + ainv.m2)

        dvals[mdd_id][det_offset + idet] = dval_refdet + dval;
      }


      // update offsets
      update_offsets(k);
    }
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
