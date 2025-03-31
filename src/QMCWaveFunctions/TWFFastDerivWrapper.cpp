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
#include "QMCWaveFunctions/Fermion/MultiSlaterDetTableMethod.h"
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

void TWFFastDerivWrapper::addMultiSlaterDet(const ParticleSet& P, const WaveFunctionComponent* msd)
{
  // add msd to slaterdet list
  slaterdets_.push_back(msd);
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

void TWFFastDerivWrapper::computeMDDerivatives_Obs(const std::vector<ValueMatrix>& Minv_Mv,
                                                   const std::vector<ValueMatrix>& Minv_B,
                                                   const std::vector<IndexType>& mdd_spo_ids,
                                                   const std::vector<const WaveFunctionComponent*>& mdds,
                                                   std::vector<ValueVector>& dvals_O) const
{
  // mdd_id is multidiracdet id
  // sid is sposet id (index into first dim of M, X, B, etc.)
  for (size_t mdd_id = 0; mdd_id < mdd_spo_ids.size(); mdd_id++)
  {
    IndexType sid = mdd_spo_ids[mdd_id];

    const auto& multidiracdet_i = static_cast<const MultiDiracDeterminant&>(*mdds[mdd_id]);

    // number of DiracDets in this MultiDiracDet
    // IndexType ndet = multidiracdet_i->getNumDets();
    IndexType ndet = multidiracdet_i.getNumDets();

    /// FIXME: set these correctly; decide how to index into nocc if not contiguous within nOcc
    size_t nelec = multidiracdet_i.getNumPtcls(); // total occ orbs in refdet (should be same as n_elec)

    /// should be same as spos.getOrbitalSetSize
    size_t nOcc = nelec; // total occ orbs in refdet (should be same as n_elec)

    /// TODO: can change this to only handle holes (i.e. orbs that we excite from in excited dets)
    ///       would need to be careful about ordering
    size_t nocc = nelec; // number of orbs that appear as holes in exc. list

    size_t norb_tot = multidiracdet_i.getNumOrbitals();

    size_t nvirt       = norb_tot - nOcc; // should be number of virtuals that appear as particles in exc. list
    size_t virt_offset = nOcc;            // first col of virt orbs in M

    // S = Minv[o,e].B[e,v] - Minv[o,e].B[e,o].Minv[o,e].Mv[e,v] ("M" from paper)


    // input mats:
    // Minv_Mv   [occ, virt]
    // Minv_B    [occ, occ]

    // O is all occ orbs; o is occ orbs that appear as holes in exc. list; e is all elecs; v is virt orbs that appear as particles in exc. list


    /// NOTE: gemms below are reversed from what is written in comments
    // a_Ov = Minv.Mvirt
    // a_Ov.T = (Mvirt.T) . (Minv.T)
    // Minv,M are row-major, but gemm assumes col-major layout, so reorder args
    // if we want output to be col-major, we can reverse inputs and use 't','t'

    /// FIXME: only keep relevant row/col subsets

    /// TODO: store this and update as single-elec moves are made
    // s[O,v] = -Minv_B[O,O].Minv_Mv[O,v]
    ValueMatrix s_Ov(nOcc, nvirt);
    BLAS::gemm('n', 'n', nvirt, nOcc, nOcc, -1.0, Minv_Mv[sid].data(), Minv_Mv[sid].cols(), Minv_B[sid].data(),
               Minv_B[sid].cols(), 0.0, s_Ov.data(), s_Ov.cols());

    // a1[O,v] = Minv_B[O,v] - Minv_B[O,O].Minv_Mv[O,v]
    for (size_t i = 0; i < s_Ov.rows(); i++)
      for (size_t j = 0; j < s_Ov.cols(); j++)
        s_Ov(i, j) += Minv_B[sid](i, j + virt_offset);


    // compute difference from ref here and add that term back later
    dvals_O[mdd_id][0] = 0.0;

    // TODO: Eq. 43
    // {} signify submatrix corresponding to holes/particles rows/cols
    //
    // dvals_O[idet] = tr(inv({Minv_Mv}).{a1})
    // build full mats for all required pairs, then select submatrices for each excited det
    // use exc index data to map into full arrays and create [k,k] tables


    size_t det_offset  = 1;
    size_t data_offset = 1;
    int max_exc_level  = multidiracdet_i.getMaxExcLevel();

    // const std::vector<int>& ndets_per_exc_lvl = *multidiracdet_i.ndets_per_excitation_level_;
    const OffloadVector<int>& excdata = multidiracdet_i.getDetData();

    auto update_offsets = [&](size_t ext_level) {
      det_offset += multidiracdet_i.getNdetPerExcLevel(ext_level);
      data_offset += multidiracdet_i.getNdetPerExcLevel(ext_level) * (3 * ext_level + 1);
    };

    // dets ordered as ref, all singles, all doubles, all triples, etc.
    //      shift by ndet_per_exc_level after each exc_level
    // data is ordered in same way, but each det has several contiguous elements
    //      each det of exc_level k has (3*k+1) elements: [k, {pos}, {uno}, {ocp}]
    //      k: exc level
    //      pos : positions of holes in refdet
    //      uno : MO idx of particles
    //      ocp : MO idx of holes

    for (size_t k = 1; k <= max_exc_level; k++)
    {
      size_t ndet_exc = multidiracdet_i.getNdetPerExcLevel(k);
      // calc vals
      // get hole/particle idx
      // take corresponding rows/cols

      /// TODO: do we need to zero these out at each iter?
      ///       should be assigning to all elements for each exc_det, so should be fine
      std::vector<IndexType> hlist(k, 0);
      std::vector<IndexType> plist(k, 0);

      // also ainv (invert in place)
      ValueMatrix a(k, k); // Minv_Mv
      ValueMatrix s(k, k); // s_Ov

      for (size_t idet = 0; idet < ndet_exc; idet++)
      {
        size_t excdata_offset = data_offset + idet * (3 * k + 1);
        assert(excdata[excdata_offset] == k);
        for (size_t i = 0; i < k; i++)
        {
          plist[i] = excdata[excdata_offset + k + 1 + i];
          hlist[i] = excdata[excdata_offset + 2 * k + 1 + i];
        }

        // construct submatrices
        for (size_t i = 0; i < k; i++)
          for (size_t j = 0; j < k; j++)
          {
            a(i, j) = Minv_Mv[sid](hlist[i], plist[j] - virt_offset);
            s(i, j) = s_Ov(hlist[i], plist[j] - virt_offset);
          }

        // invert a in place
        invert_matrix(a, false);

        // dval_O = dval_O_refdet + tr(ainv.s)

        ValueType dval_O = 0.0;

        // tr(ainv.s)
        for (size_t i = 0; i < k; i++)
          for (size_t j = 0; j < k; j++)
            dval_O += a(i, j) * s(j, i);

        dvals_O[mdd_id][det_offset + idet] = dval_O;
      }

      // update offsets
      update_offsets(k);
    }
  }

  return;
}

void TWFFastDerivWrapper::transform_Av_AoBv(const ValueMatrix& A, const ValueMatrix& B, ValueMatrix& X) const
{
  // A [nocc,nocc+nvirt]
  // B [nocc,nvirt]

  // return A[o,v] - A[o,o].B[o,v]

  const int nocc  = A.rows();
  const int nvirt = B.cols();

  for (size_t i = 0; i < nocc; i++)
    for (size_t j = 0; j < nvirt; j++)
      X(i, j) = A(i, j + nocc);

  // X = X - A[o,o].B[o,v]
  BLAS::gemm('n', 'n', nvirt, nocc, nocc, -1.0, A.data(), A.cols(), B.data(), B.cols(), 1.0, X.data(), X.cols());
  return;
}

void TWFFastDerivWrapper::computeMDDerivatives_dmu(const std::vector<ValueMatrix>& Minv_Mv,
                                                   const std::vector<ValueMatrix>& Minv_B,
                                                   const std::vector<ValueMatrix>& Minv_dM,
                                                   const std::vector<ValueMatrix>& Minv_dB,
                                                   const std::vector<IndexType>& mdd_spo_ids,
                                                   const std::vector<const WaveFunctionComponent*>& mdds,
                                                   std::vector<ValueVector>& dvals_dmu_O,
                                                   std::vector<ValueVector>& dvals_dmu) const
{
  // mdd_id is multidiracdet id
  // sid is sposet id (index into first dim of M, X, B, etc.)
  for (size_t mdd_id = 0; mdd_id < mdd_spo_ids.size(); mdd_id++)
  {
    IndexType sid = mdd_spo_ids[mdd_id];

    const auto& multidiracdet_i = static_cast<const MultiDiracDeterminant&>(*mdds[mdd_id]);

    // number of DiracDets in this MultiDiracDet
    // IndexType ndet = multidiracdet_i->getNumDets();
    IndexType ndet = multidiracdet_i.getNumDets();


    /// FIXME: set these correctly; decide how to index into nocc if not contiguous within nOcc
    size_t nelec = multidiracdet_i.getNumPtcls(); // total occ orbs in refdet (should be same as n_elec)

    /// should be same as spos.getOrbitalSetSize
    size_t nOcc = nelec; // total occ orbs in refdet (should be same as n_elec)

    /// TODO: can change this to only handle holes (i.e. orbs that we excite from in excited dets)
    ///       would need to be careful about ordering
    size_t nocc = nelec; // number of orbs that appear as holes in exc. list

    size_t norb_tot = multidiracdet_i.getNumOrbitals();

    size_t nvirt       = norb_tot - nOcc; // should be number of virtuals that appear as particles in exc. list
    size_t virt_offset = nOcc;            // first col of virt orbs in M

    // input mats:
    // a = Minv_Mv  [occ, virt]
    // X2 = Minv_dM  [occ, occ+virt]
    // X3 = Minv_B   [occ, occ+virt]
    // X4 = Minv_dB  [occ, occ+virt]

    // O is all occ orbs; o is occ orbs that appear as holes in exc. list; e is all elecs; v is virt orbs that appear as particles in exc. list

    // dvals_dmu_O[idet] = tr(
    //     - inv({a}).{X2[o,v] - X2[o,o].a[o,v]}.inv({a}).{X3[o,v] - X3[o,o].a[o,v]} ...
    //     + inv({a}).{(X4[o,v] - X3[o,o].X2[o,v]) - (X4[o,o] - X3[o,o].X2[o,o]).a[o,v] - X2[o,o].(X3[o,v] - X3[o,o].a)}
    //   )

    // X32[o,o+v] = X3[o,o].X2[o,o+v]
    // X32[o,o+v] = Minv_B[o,o].Minv_dM[o,o+v]
    ValueMatrix X32_On(nOcc, nOcc + nvirt);
    BLAS::gemm('n', 'n', nOcc + nvirt, nOcc, nOcc, 1.0, Minv_dM[sid].data(), Minv_dM[sid].cols(), Minv_B[sid].data(),
               Minv_B[sid].cols(), 0.0, X32_On.data(), X32_On.cols());

    // Minv_dB - X32
    ValueMatrix X432_On(nOcc, nOcc + nvirt);
    X432_On = Minv_dB[sid] - X32_On;

    ValueMatrix X432b_Ov(nOcc, nvirt);
    transform_Av_AoBv(X432_On, Minv_Mv[sid], X432b_Ov);

    ValueMatrix X2b_Ov(nOcc, nvirt);
    transform_Av_AoBv(Minv_dM[sid], Minv_Mv[sid], X2b_Ov);

    ValueMatrix X3b_Ov(nOcc, nvirt);
    transform_Av_AoBv(Minv_B[sid], Minv_Mv[sid], X3b_Ov);

    // X432b = X432b - Minv_dM[o,o].X3b[o,v]
    BLAS::gemm('n', 'n', nvirt, nOcc, nOcc, -1.0, X3b_Ov.data(), X3b_Ov.cols(), Minv_dM[sid].data(),
               Minv_dM[sid].cols(), 1.0, X432b_Ov.data(), X432b_Ov.cols());


    // compute difference from ref here and add that term back later
    dvals_dmu_O[mdd_id][0] = 0.0;
    dvals_dmu[mdd_id][0]   = 0.0;


    // TODO: Eq. 43
    // {} signify submatrix corresponding to holes/particles rows/cols
    //
    // dvals_dmu_O[idet] = tr(-inv({a}).{X2b}.inv({a}).{X3b} + inv({a}).{X432b})
    // build full mats for all required pairs, then select submatrices for each excited det
    // use exc index data to map into full arrays and create [k,k] tables


    size_t det_offset  = 1;
    size_t data_offset = 1;
    int max_exc_level  = multidiracdet_i.getMaxExcLevel();

    // const std::vector<int>& ndets_per_exc_lvl = *multidiracdet_i.ndets_per_excitation_level_;
    const OffloadVector<int>& excdata = multidiracdet_i.getDetData();

    auto update_offsets = [&](size_t ext_level) {
      det_offset += multidiracdet_i.getNdetPerExcLevel(ext_level);
      data_offset += multidiracdet_i.getNdetPerExcLevel(ext_level) * (3 * ext_level + 1);
    };

    // dets ordered as ref, all singles, all doubles, all triples, etc.
    //      shift by ndet_per_exc_level after each exc_level
    // data is ordered in same way, but each det has several contiguous elements
    //      each det of exc_level k has (3*k+1) elements: [k, {pos}, {uno}, {ocp}]
    //      k: exc level
    //      pos : positions of holes in refdet
    //      uno : MO idx of particles
    //      ocp : MO idx of holes

    for (size_t k = 1; k <= max_exc_level; k++)
    {
      size_t ndet_exc = multidiracdet_i.getNdetPerExcLevel(k);

      // get hole/particle idx
      // take corresponding rows/cols

      // a' = inv({a})
      // dvals_dmu_O[idet] = tr(-a'.{X2b}.a'.{X3b} + a'.{X432b})
      // dvals_dmu[idet]   = tr(a'.{X2b})

      // a'.{S} also needed for Eq. 29

      /// TODO: do we need to zero these out at each iter?
      ///       should be assigning to all elements for each exc_det, so should be fine
      std::vector<IndexType> hlist(k, 0);
      std::vector<IndexType> plist(k, 0);

      // also ainv (invert in place)
      ValueMatrix a(k, k);
      ValueMatrix x2(k, k);
      // also (m2 - m1.ainv.s) (update in place)
      ValueMatrix x3(k, k);
      ValueMatrix x4(k, k);

      ValueMatrix ainv_x2(k, k); // tr -> d(logD)
      ValueMatrix ainv_x3(k, k);
      ValueMatrix ainv_x4(k, k);


      for (size_t idet = 0; idet < ndet_exc; idet++)
      {
        size_t excdata_offset = data_offset + idet * (3 * k + 1);
        assert(excdata[excdata_offset] == k);
        for (size_t i = 0; i < k; i++)
        {
          plist[i] = excdata[excdata_offset + k + 1 + i];
          hlist[i] = excdata[excdata_offset + 2 * k + 1 + i];
        }

        // construct submatrices
        for (size_t i = 0; i < k; i++)
          for (size_t j = 0; j < k; j++)
          {
            a(i, j)  = Minv_Mv[sid](hlist[i], plist[j] - virt_offset);
            x2(i, j) = X2b_Ov(hlist[i], plist[j] - virt_offset);
            x3(i, j) = X3b_Ov(hlist[i], plist[j] - virt_offset);
            x4(i, j) = X432b_Ov(hlist[i], plist[j] - virt_offset);
          }

        // invert a in place
        invert_matrix(a, false);

        // ainv_x4 = ainv.x4
        BLAS::gemm('n', 'n', k, k, k, 1.0, x4.data(), x4.cols(), a.data(), a.cols(), 0.0, ainv_x4.data(),
                   ainv_x4.cols());

        // ainv_x3 = ainv.x3
        BLAS::gemm('n', 'n', k, k, k, 1.0, x3.data(), x3.cols(), a.data(), a.cols(), 0.0, ainv_x3.data(),
                   ainv_x3.cols());

        // ainv_x2 = ainv.x2
        BLAS::gemm('n', 'n', k, k, k, 1.0, x2.data(), x2.cols(), a.data(), a.cols(), 0.0, ainv_x2.data(),
                   ainv_x2.cols());

        // ainv_x4 = ainv_x4 - ainv_x2.ainv_x3
        BLAS::gemm('n', 'n', k, k, k, -1.0, ainv_x3.data(), ainv_x3.cols(), ainv_x2.data(), ainv_x2.cols(), 1.0,
                   ainv_x4.data(), ainv_x4.cols());


        // dval_dmu_O = tr(ainv_x4 - ainv_x2.ainv_x3)
        // dval_dmu = tr(ainv_x2)

        ValueType dval_dmu_O = 0.0;
        ValueType dval_dmu   = 0.0;


        for (size_t i = 0; i < k; i++)
        {
          dval_dmu += ainv_x2(i, i);
          dval_dmu_O += ainv_x4(i, i);
        }

        dvals_dmu_O[mdd_id][det_offset + idet] = dval_dmu_O;
        dvals_dmu[mdd_id][det_offset + idet]   = dval_dmu;
      }

      // update offsets
      update_offsets(k);
    }
  }

  return;
}

void TWFFastDerivWrapper::computeMDDerivatives_ExcDets(const std::vector<ValueMatrix>& Minv,
                                                       const std::vector<ValueMatrix>& X,
                                                       const std::vector<ValueMatrix>& dM,
                                                       const std::vector<ValueMatrix>& dB,
                                                       const std::vector<ValueMatrix>& B,
                                                       const std::vector<ValueMatrix>& M,
                                                       const std::vector<IndexType>& mdd_spo_ids,
                                                       const std::vector<const WaveFunctionComponent*>& mdds,
                                                       std::vector<ValueVector>& dvals_dmu_O,
                                                       std::vector<ValueVector>& dvals_O,
                                                       std::vector<ValueVector>& dvals_dmu) const
{
  // mdd_id is multidiracdet id
  // sid is sposet id (index into first dim of M, X, B, etc.)
  for (size_t mdd_id = 0; mdd_id < mdd_spo_ids.size(); mdd_id++)
  {
    IndexType sid = mdd_spo_ids[mdd_id];

    const auto& multidiracdet_i = static_cast<const MultiDiracDeterminant&>(*mdds[mdd_id]);

    // number of DiracDets in this MultiDiracDet
    // IndexType ndet = multidiracdet_i->getNumDets();
    IndexType ndet = multidiracdet_i.getNumDets();

    /// TODO: just do this before calling
    for (int idet = 0; idet < ndet; idet++)
    {
      dvals_dmu_O[mdd_id][idet] = 0.0;
      dvals_O[mdd_id][idet]     = 0.0;
      dvals_dmu[mdd_id][idet]   = 0.0;
    }


    /// FIXME: set these correctly; decide how to index into nocc if not contiguous within nOcc
    size_t nelec = multidiracdet_i.getNumPtcls(); // total occ orbs in refdet (should be same as n_elec)

    /// should be same as spos.getOrbitalSetSize
    size_t nOcc = nelec; // total occ orbs in refdet (should be same as n_elec)

    /// TODO: can change this to only handle holes (i.e. orbs that we excite from in excited dets)
    ///       would need to be careful about ordering
    size_t nocc = nelec; // number of orbs that appear as holes in exc. list

    size_t norb_tot = multidiracdet_i.getNumOrbitals();

    size_t nvirt       = norb_tot - nOcc; // should be number of virtuals that appear as particles in exc. list
    size_t virt_offset = nOcc;            // first col of virt orbs in M

    // input mats:
    // Minv   [occ, elec]
    // X      [occ, elec] (Minv.B(Occ).Minv)
    // dM     [elec, occ+virt]
    // dB     [elec, occ+virt]
    // B      [elec, (occ+?)virt]
    // M      [elec, occ+virt]
    // Mvirt  [elec, virt] (subset of M)
    // Mvirt.data() = M.data() + nOcc (ldim is occ+virt)

    // O is all occ orbs; o is occ orbs that appear as holes in exc. list; e is all elecs; v is virt orbs that appear as particles in exc. list

    // dvals_dmu_O[idet] = dval_ref + tr(
    //     - inv({Minv[o,e].M[e,v]}).{Minv[o,e].dM[e,v] - Minv[o,e].dM[e,O].Minv[O,e].M[e,v]}.inv({Minv[o,e].M[e,v]}).{S} ...
    //     + inv({Minv[o,e].M[e,v]}).{Minv[o,e].dB[e,v] - X[o,e].dM[e,v] - (Minv[o,e].dB[e.O] - X[o,e].dM[e,O]).Minv[O,e].M[e,v] - Minv[o,e].dM[e,O].S[O,v]}
    //   )


    // X = Minv[O,e].B[e,O].Minv[O,e]

    // a = Minv[O,e].M[e,v] (slice this for alpha from paper)
    // b = Minv[o,e].dM[e,Ov]
    // c = Minv[o,e].dB[e,Ov] (first dim O for ref det)
    // d = X[o,e].dM[e,Ov]    (first dim O for ref det)
    // f = Minv[O,e].B[e,v/O] (only needed v for exc dmu; need O for refdet Od)
    // g = X[O,e].M[e,v]
    // h = c[o,Ov] - d[o,Ov]  (refdet is tr(h[O,O]) )

    // S = (Minv[O,e].B[e,v] - X[O,e].M[e,v]) ("M" from paper)
    //   = Minv[O,e].(B[e,v] - B[e,O].Minv[O,e].M[e,v])
    //   = f[O,v] - g[O,v]

    // mat1b  = Minv[o,e].dM[e,O].Minv[O,e].M[e,v]
    //        = Minv[o,e].dM[e,O].a[O,v]
    //        = b[o,O].a[O,v]

    /// NOTE: use for dmu_log
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

    //  dvals_dmu_O[idet] = dval_dmu_ref + tr(-inv({a}).{mat1}.inv({a}).{S} + inv({a}).{mat2})

    // dval_ref = tr(Minv[O,e].dB[e,O] - X[O,e].dM[e,O])
    //          = tr(c[O,O] - d[O,O])


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

    // for refdet Od/d
    // f[O,O] = Minv[O,e].B[e,O]
    ValueMatrix f_OO(nOcc, nOcc);
    BLAS::gemm('n', 'n', nOcc, nOcc, nelec, 1.0, B[sid].data(), B[sid].cols(), Minv[sid].data(), Minv[sid].cols(), 0.0,
               f_OO.data(), f_OO.cols());

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
    ValueType dval_dmu_O_refdet = 0.0;
    ValueType dval_O_refdet     = 0.0;
    ValueType dval_dmu_refdet   = 0.0;

    for (size_t i = 0; i < h_On.rows(); i++)
    {
      dval_dmu_O_refdet += h_On(i, i);
      dval_O_refdet += f_OO(i, i);
      dval_dmu_refdet += b_On(i, i);
    }

    // dvals_dmu_O[mdd_id][0] = dval_dmu_O_refdet;
    // dvals_O[mdd_id][0]     = dval_O_refdet;
    // dvals_dmu[mdd_id][0]   = dval_dmu_refdet;
    // debugging: compute difference from ref here and add that term back later
    dvals_dmu_O[mdd_id][0] = 0.0;
    dvals_O[mdd_id][0]     = 0.0;
    dvals_dmu[mdd_id][0]   = 0.0;


    // TODO: Eq. 43
    // {} signify submatrix corresponding to holes/particles rows/cols
    //
    // dvals_dmu_O[idet] = dval_refdet + tr(-inv({a}).{mat1}.inv({a}).{S} + inv({a}).{mat2})
    // build full mats for all required pairs, then select submatrices for each excited det
    // use exc index data to map into full arrays and create [k,k] tables


    size_t det_offset  = 1;
    size_t data_offset = 1;
    int max_exc_level  = multidiracdet_i.getMaxExcLevel();

    // const std::vector<int>& ndets_per_exc_lvl = *multidiracdet_i.ndets_per_excitation_level_;
    const OffloadVector<int>& excdata = multidiracdet_i.getDetData();

    /// FIXME: do we need signs (parity of perm to normal-order after excitation)?
    ///        everything is a ratio, so signs cancel?
    ///        dmu(OD/D), OD/D, dmu(logD) == (dmu D)/D
    const OffloadVector<RealType>& det_signs = multidiracdet_i.getDetSigns();

    auto update_offsets = [&](size_t ext_level) {
      det_offset += multidiracdet_i.getNdetPerExcLevel(ext_level);
      data_offset += multidiracdet_i.getNdetPerExcLevel(ext_level) * (3 * ext_level + 1);
    };

    // dets ordered as ref, all singles, all doubles, all triples, etc.
    //      shift by ndet_per_exc_level after each exc_level
    // data is ordered in same way, but each det has several contiguous elements
    //      each det of exc_level k has (3*k+1) elements: [k, {pos}, {uno}, {ocp}]
    //      k: exc level
    //      pos : positions of holes in refdet
    //      uno : MO idx of particles
    //      ocp : MO idx of holes

    for (size_t k = 1; k <= max_exc_level; k++)
    {
      size_t ndet_exc = multidiracdet_i.getNdetPerExcLevel(k);
      // calc vals

      // see MultiDiracDeterminant::mw_updateRatios, MultiDiracDeterminant::mw_buildTableMatrix_calculateRatios_impl
      // for examples of excited det idx layout

      // get hole/particle idx
      // take corresponding rows/cols of a, S, mat1, mat2

      // a' = inv({a})
      // dvals_dmu_O[idet] = dval_refdet + tr(-a'.{mat1}.a'.{S} + a'.{mat2})
      // dvals_dmu_O[idet] = dval_refdet + tr(-a'.({mat1}.a'.{S} - {mat2}))

      // a'.{S} also needed for Eq. 29

      /// TODO: do we need to zero these out at each iter?
      ///       should be assigning to all elements for each exc_det, so should be fine
      std::vector<IndexType> hlist(k, 0);
      std::vector<IndexType> plist(k, 0);

      // also ainv (invert in place)
      ValueMatrix a(k, k);
      ValueMatrix m1(k, k);
      // also (m2 - m1.ainv.s) (update in place)
      ValueMatrix m2(k, k);
      ValueMatrix s(k, k);

      ValueMatrix ainv_s(k, k);  // tr -> OD/D
      ValueMatrix ainv_m1(k, k); // tr -> d(logD)


      for (size_t idet = 0; idet < ndet_exc; idet++)
      {
        size_t excdata_offset = data_offset + idet * (3 * k + 1);
        assert(excdata[excdata_offset] == k);
        for (size_t i = 0; i < k; i++)
        {
          plist[i] = excdata[excdata_offset + k + 1 + i];
          hlist[i] = excdata[excdata_offset + 2 * k + 1 + i];
        }

        // construct submatrices
        for (size_t i = 0; i < k; i++)
          for (size_t j = 0; j < k; j++)
          {
            a(i, j)  = a_Ov(hlist[i], plist[j] - virt_offset);
            m1(i, j) = mat1(hlist[i], plist[j] - virt_offset);
            m2(i, j) = mat2(hlist[i], plist[j] - virt_offset);
            s(i, j)  = S_Ov(hlist[i], plist[j] - virt_offset);
          }

        // invert a in place
        invert_matrix(a, false);

        // ainv_s = ainv.s
        BLAS::gemm('n', 'n', k, k, k, 1.0, s.data(), s.cols(), a.data(), a.cols(), 0.0, ainv_s.data(), ainv_s.cols());

        // ainv_m1 = ainv.m1
        BLAS::gemm('n', 'n', k, k, k, 1.0, m1.data(), m1.cols(), a.data(), a.cols(), 0.0, ainv_m1.data(),
                   ainv_m1.cols());

        // m2 = -m1.ainv.s + m2
        BLAS::gemm('n', 'n', k, k, k, -1.0, ainv_s.data(), ainv_s.cols(), m1.data(), m1.cols(), 1.0, m2.data(),
                   m2.cols());

        // dval_dmu_O = dval_dmu_O_refdet + tr(-ainv.m1.ainv.s + ainv.m2)
        //            = dval_dmu_O_refdet + tr(ainv.(m2 - m1.ainv.s))

        // dval_O = dval_O_refdet + tr(ainv.s)

        // dval_dmu = dval_dmu_refdet + tr(ainv.m1)

        // ValueType dval_dmu_O = dval_dmu_O_refdet;
        // ValueType dval_O     = dval_O_refdet;
        // ValueType dval_dmu   = dval_dmu_refdet;
        // debugging: use diff instead of total
        ValueType dval_dmu_O = 0.0;
        ValueType dval_O     = 0.0;
        ValueType dval_dmu   = 0.0;


        for (size_t i = 0; i < k; i++)
        {
          dval_O += ainv_s(i, i);
          dval_dmu += ainv_m1(i, i);
          for (size_t j = 0; j < k; j++)
          {
            dval_dmu_O += m2(i, j) * a(j, i);
          }
        }

        dvals_dmu_O[mdd_id][det_offset + idet] = dval_dmu_O;
        dvals_O[mdd_id][det_offset + idet]     = dval_O;
        dvals_dmu[mdd_id][det_offset + idet]   = dval_dmu;
      }


      // update offsets
      update_offsets(k);
    }
  }

  return;
}

std::tuple<TWFFastDerivWrapper::ValueType, TWFFastDerivWrapper::ValueType, TWFFastDerivWrapper::ValueType>
    TWFFastDerivWrapper::computeMDDerivatives_total(int msd_idx,
                                                    const std::vector<const WaveFunctionComponent*>& mdds,
                                                    const std::vector<ValueVector>& dvals_dmu_O,
                                                    const std::vector<ValueVector>& dvals_O,
                                                    const std::vector<ValueVector>& dvals_dmu) const
{
  /**
   * \f[
   *   \partial_\mu\frac{\hat{O} \Psi}{\Psi} =
   *       \frac{\sum_i c_i D_i \left(\partial_\mu\frac{\hat{O}D_i}{D_i}\right)}{\Psi} 
   *     + \frac{\sum_i c_i D_i \left(\frac{\hat{O}D_i}{D_i}\right)\partial_\mu\log(D_i)}{\Psi} 
   *     - \left(\frac{\sum_i c_i D_i \left(\frac{\hat{O}D_i}{D_i}\right)}{\Psi}\right)
   *     * \left(\frac{\sum_i c_i D_i \partial_\mu \log(D_i)}{\Psi}\right)
   * \f]
   * 
   * with:
   *  w_i = c_i D_i
   *  x_i = d_mu(O D_i/D_i)
   *  y_i = O D_i/D_i
   *  z_i = d_mu(log(D_i))
   * 
   * \f[
   *   \partial_\mu\frac{\hat{O} \Psi}{\Psi} =
   *       \frac{\sum_i w_i * x_i}{\Psi} 
   *     + \frac{\sum_i w_i * y_i * z_i}{\Psi} 
   *     - \left(\frac{\sum_i w_i * y_i}{\Psi}\right)
   *     * \left(\frac{\sum_i w_i * z_i}{\Psi}\right)
   * \f]
   *
   * D is slaterdet composed of spindets
   * D_i = A_{a_i}*B_{b_i} 
   * A and B are alpha/beta spindets
   * a and b map from slaterdet list into unique alpha/beta lists
   * 
   * with D = AB
   *   OD/D = O(AB)/(AB) = ((OA)B + A(OB))/(AB) = OA/A + OB/B
   *   d(OD/D) = d(OA/A + OB/B) = d(OA/A) + d(OB/B)
   *   d(log(D)) = d(AB)/(AB) = ((dA)B + A(dB))/(AB) = dA/A + dB/B = d(log(A)) + d(log(B))
   * 
   */

  // multislaterdet
  const auto& msd = static_cast<const MultiSlaterDetTableMethod&>(getMultiSlaterDet(msd_idx));

  // total slaterdets
  const int num_slaterdets = msd.getNumSlaterDets();

  // Coefs for slaterdets
  const std::vector<ValueType>& C = msd.get_C();

  // mapping from [group_idx][slater_idx] to diracdet_idx
  const std::vector<std::vector<size_t>>& C2node = msd.get_C2node();

  // number of unique diracdets of each group
  std::vector<IndexType> num_diracdets;

  std::vector<const OffloadVector<ValueType>*> diracdet_ratios_to_ref;

  for (size_t i = 0; i < msd.getDetSize(); i++)
  {
    num_diracdets.push_back(msd.getDet(i).getNumDets());
    diracdet_ratios_to_ref.push_back(&msd.getDet(i).getRatiosToRefDet());
  }

  const int num_groups = num_diracdets.size();
  // app_log() << "DEBUG: num_groups = " << num_groups << std::endl;

  ValueType total_psi   = C[0]; // sum_i c_i D_i
  ValueType total_dmu_O = 0.0;  // d_mu(OD/D)
  ValueType total_O     = 0.0;  // OD/D
  ValueType total_dmu   = 0.0;  // d_mu(log(D))
  ValueType total_Odmu  = 0.0;  // (OD/D) * d_mu(log(D))

  /// NOTE: sum doesn't include refdet
  /// total terms are added later; psi C0 term is added above
  for (size_t i_sd = 1; i_sd < num_slaterdets; i_sd++)
  {
    ValueType tmp_psi   = C[i_sd]; // C[i_sd] * prod_i ratio[i][i_sd]
    ValueType tmp_dmu_O = 0.0;     // d_mu(OD/D)
    ValueType tmp_O     = 0.0;     // OD/D
    ValueType tmp_dmu   = 0.0;     // d_mu(log(D))
    ValueType tmp_Odmu  = 0.0;     // (OD/D) * d_mu(log(D))

    // app_log() << "msd::C[" << i_sd << "] = " << tmp_psi << std::endl;


    for (size_t i_group = 0; i_group < num_groups; i_group++)
    {
      size_t i_dd = C2node[i_group][i_sd];
      // std::cout << "DEBUG: local det id " << i_dd << std::endl;


      // app_log() << "tmp_dmu_O  [" << i_group << "][" << i_dd << "]" << dvals_dmu_O[i_group][i_dd] << std::endl;
      // app_log() << "tmp_O  [" << i_group << "][" << i_dd << "]" << dvals_O[i_group][i_dd] << std::endl;
      // app_log() << "tmp_dmu  [" << i_group << "][" << i_dd << "]" << dvals_dmu[i_group][i_dd] << std::endl;
      // app_log() << "tmp_Odmu [" << i_group << "][" << i_dd << "]"
      //           << dvals_dmu[i_group][i_dd] * dvals_O[i_group][i_dd] << std::endl;
      // app_log() << "tmp_psi [" << i_group << "][" << i_dd << "]" << (*diracdet_ratios_to_ref[i_group])[i_dd]
      //           << std::endl;

      tmp_dmu_O += dvals_dmu_O[i_group][i_dd];
      tmp_O += dvals_O[i_group][i_dd];
      tmp_dmu += dvals_dmu[i_group][i_dd];
      tmp_Odmu += dvals_dmu[i_group][i_dd] * dvals_O[i_group][i_dd];
      tmp_psi *= (*diracdet_ratios_to_ref[i_group])[i_dd];
    }

    total_psi += tmp_psi;
    total_dmu_O += tmp_dmu_O * tmp_psi;
    total_O += tmp_O * tmp_psi;
    total_dmu += tmp_dmu * tmp_psi;
    // total_Odmu += tmp_Odmu * tmp_psi;
    total_Odmu += tmp_O * tmp_dmu * tmp_psi;
  }

  ValueType norm_dmu_O = total_dmu_O / total_psi;
  ValueType norm_O     = total_O / total_psi;
  ValueType norm_dmu   = total_dmu / total_psi;
  ValueType norm_Odmu  = total_Odmu / total_psi;

  // d_mu(log(Psi))
  ValueType dmu_psi = norm_dmu;

  // (OPsi/Psi) (this is the same for all spatial derivs, but we get it for free here)
  ValueType Opsi = norm_O;

  // d_mu(OPsi/Psi)
  ValueType dmu_O_psi = norm_Odmu + norm_dmu_O - norm_O * norm_dmu;

  return {dmu_O_psi, dmu_psi, Opsi};
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
void TWFFastDerivWrapper::buildIntermediates(const std::vector<ValueMatrix>& Minv,
                                             const std::vector<ValueMatrix>& B,
                                             const std::vector<ValueMatrix>& M,
                                             std::vector<ValueMatrix>& X,
                                             std::vector<ValueMatrix>& Minv_B,
                                             std::vector<ValueMatrix>& Minv_Mv)
{
  IndexType nspecies = Minv.size();

  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum = Minv[id].rows();
    assert(Minv[id].rows() == Minv[id].cols());
    assert(X[id].rows() == ptclnum);
    assert(X[id].cols() == ptclnum);
    int norb  = M[id].cols();
    int nvirt = norb - ptclnum;

    // Minv_Mv = Minv[e,e].M[e,v]
    BLAS::gemm('n', 'n', nvirt, ptclnum, ptclnum, 1.0, M[id].data() + ptclnum, M[id].cols(), Minv[id].data(),
               Minv[id].cols(), 0.0, Minv_Mv[id].data(), Minv_Mv[id].cols());

    // Minv_B = Minv[e,e].B[e,O]
    BLAS::gemm('n', 'n', norb, ptclnum, ptclnum, 1.0, B[id].data(), B[id].cols(), Minv[id].data(), Minv[id].cols(), 0.0,
               Minv_B[id].data(), Minv_B[id].cols());


    // X = Minv_B[e,e].Minv[e,e]
    BLAS::gemm('n', 'n', ptclnum, ptclnum, ptclnum, 1.0, Minv[id].data(), Minv[id].cols(), Minv_B[id].data(),
               Minv_B[id].cols(), 0.0, X[id].data(), X[id].cols());
  }
}


void TWFFastDerivWrapper::buildIntermediates_dmu(const std::vector<ValueMatrix>& Minv,
                                                 const std::vector<std::vector<ValueMatrix>>& dB,
                                                 const std::vector<std::vector<ValueMatrix>>& dM,
                                                 std::vector<std::vector<ValueMatrix>>& Minv_dB,
                                                 std::vector<std::vector<ValueMatrix>>& Minv_dM)
{
  IndexType nspecies = Minv.size();
  IndexType ndim     = dB.size();

  for (IndexType id = 0; id < nspecies; id++)
  {
    int ptclnum = Minv[id].rows();
    assert(Minv[id].rows() == Minv[id].cols());

    for (IndexType idim = 0; idim < ndim; idim++)
    {
      int norb = dM[idim][id].cols();
      // Minv_dM = Minv[e,e].dM[e,O]
      BLAS::gemm('n', 'n', norb, ptclnum, ptclnum, 1.0, dM[idim][id].data(), dM[idim][id].cols(), Minv[id].data(),
                 Minv[id].cols(), 0.0, Minv_dM[idim][id].data(), Minv_dM[idim][id].cols());
      // Minv_dB = Minv[e,e].dB[e,O]
      BLAS::gemm('n', 'n', norb, ptclnum, ptclnum, 1.0, dB[idim][id].data(), dB[idim][id].cols(), Minv[id].data(),
                 Minv[id].cols(), 0.0, Minv_dB[idim][id].data(), Minv_dB[idim][id].cols());
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

void TWFFastDerivWrapper::wipeVectors(std::vector<ValueVector>& A)
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
