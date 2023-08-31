//////////////////////////////////////////////////////////////////////////////////////
//// This file is distributed under the University of Illinois/NCSA Open Source License.
//// See LICENSE file in top directory for details.
////
//// Copyright (c) QMCPACK developers.
////
//// File developed by: Sergio D. Pineda Flores, sergio_pinedaflores@berkeley.edu, University of California, Berkeley
////                    Eric Neuscamman, eneuscamman@berkeley.edu, University of California, Berkeley
////                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
////
//// File created by: Sergio D. Pineda Flores, sergio_pinedaflores@berkeley.edu, University of California, Berkeley
////////////////////////////////////////////////////////////////////////////////////////
#include "RotatedSPOs.h"
#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
#include "CPU/BLAS.hpp"
#include "io/hdf/hdf_archive.h"


namespace qmcplusplus
{
RotatedSPOs::RotatedSPOs(const std::string& my_name, std::unique_ptr<SPOSet>&& spos)
    : SPOSet(my_name),
      OptimizableObject(my_name),
      Phi(std::move(spos)),
      nel_major_(0),
      params_supplied(false),
      apply_rotation_timer_(createGlobalTimer("RotatedSPOs::apply_rotation", timer_level_fine))
{
  OrbitalSetSize = Phi->getOrbitalSetSize();
}

RotatedSPOs::~RotatedSPOs() {}


void RotatedSPOs::setRotationParameters(const std::vector<RealType>& param_list)
{
  params          = param_list;
  params_supplied = true;
}

void RotatedSPOs::createRotationIndices(int nel, int nmo, RotationIndices& rot_indices)
{
  for (int i = 0; i < nel; i++)
    for (int j = nel; j < nmo; j++)
      rot_indices.emplace_back(i, j);
}

void RotatedSPOs::createRotationIndicesFull(int nel, int nmo, RotationIndices& rot_indices)
{
  rot_indices.reserve(nmo * (nmo - 1) / 2);

  // start with core-active rotations - put them at the beginning of the list
  // so it matches the other list of rotation indices
  for (int i = 0; i < nel; i++)
    for (int j = nel; j < nmo; j++)
      rot_indices.emplace_back(i, j);

  // Add core-core rotations - put them at the end of the list
  for (int i = 0; i < nel; i++)
    for (int j = i + 1; j < nel; j++)
      rot_indices.emplace_back(i, j);

  // Add active-active rotations - put them at the end of the list
  for (int i = nel; i < nmo; i++)
    for (int j = i + 1; j < nmo; j++)
      rot_indices.emplace_back(i, j);
}

void RotatedSPOs::constructAntiSymmetricMatrix(const RotationIndices& rot_indices,
                                               const std::vector<ValueType>& param,
                                               ValueMatrix& rot_mat)
{
  assert(rot_indices.size() == param.size());
  // Assumes rot_mat is of the correct size

  rot_mat = 0.0;

  for (int i = 0; i < rot_indices.size(); i++)
  {
    const int p      = rot_indices[i].first;
    const int q      = rot_indices[i].second;
    const RealType x = param[i];

    rot_mat[q][p] = x;
    rot_mat[p][q] = -x;
  }
}

void RotatedSPOs::extractParamsFromAntiSymmetricMatrix(const RotationIndices& rot_indices,
                                                       const ValueMatrix& rot_mat,
                                                       std::vector<ValueType>& param)
{
  assert(rot_indices.size() == param.size());
  // Assumes rot_mat is of the correct size

  for (int i = 0; i < rot_indices.size(); i++)
  {
    const int p = rot_indices[i].first;
    const int q = rot_indices[i].second;
    param[i]    = rot_mat[q][p];
  }
}

void RotatedSPOs::resetParametersExclusive(const opt_variables_type& active)
{
  std::vector<RealType> delta_param(m_act_rot_inds.size());

  size_t psize = m_act_rot_inds.size();

  if (use_global_rot_)
  {
    psize = m_full_rot_inds.size();
    assert(psize >= m_act_rot_inds.size());
  }

  std::vector<RealType> old_param(psize);
  std::vector<RealType> new_param(psize);

  for (int i = 0; i < m_act_rot_inds.size(); i++)
  {
    int loc        = myVars.where(i);
    delta_param[i] = active[loc] - myVars[i];
    myVars[i]      = active[loc];
  }

  if (use_global_rot_)
  {
    for (int i = 0; i < m_full_rot_inds.size(); i++)
      old_param[i] = myVarsFull[i];

    applyDeltaRotation(delta_param, old_param, new_param);

    // Save the the params
    for (int i = 0; i < m_full_rot_inds.size(); i++)
      myVarsFull[i] = new_param[i];
  }
  else
  {
    apply_rotation(delta_param, false);

    // Save the parameters in the history list
    history_params_.push_back(delta_param);
  }
}

void RotatedSPOs::writeVariationalParameters(hdf_archive& hout)
{
  hout.push("RotatedSPOs");
  if (use_global_rot_)
  {
    hout.push("rotation_global");
    std::string rot_global_name = std::string("rotation_global_") + SPOSet::getName();

    int nparam_full = myVarsFull.size();
    std::vector<RealType> full_params(nparam_full);
    for (int i = 0; i < nparam_full; i++)
      full_params[i] = myVarsFull[i];

    hout.write(full_params, rot_global_name);
    hout.pop();
  }
  else
  {
    hout.push("rotation_history");
    size_t rows = history_params_.size();
    size_t cols = 0;
    if (rows > 0)
      cols = history_params_[0].size();

    Matrix<RealType> tmp(rows, cols);
    for (size_t i = 0; i < rows; i++)
      for (size_t j = 0; j < cols; j++)
        tmp(i, j) = history_params_[i][j];

    std::string rot_hist_name = std::string("rotation_history_") + SPOSet::getName();
    hout.write(tmp, rot_hist_name);
    hout.pop();
  }

  // Save myVars in order to restore object state exactly
  //  The values aren't meaningful, but they need to match those saved in VariableSet
  hout.push("rotation_params");
  std::string rot_params_name = std::string("rotation_params_") + SPOSet::getName();

  int nparam = myVars.size();
  std::vector<RealType> params(nparam);
  for (int i = 0; i < nparam; i++)
    params[i] = myVars[i];

  hout.write(params, rot_params_name);
  hout.pop();

  hout.pop();
}

void RotatedSPOs::readVariationalParameters(hdf_archive& hin)
{
  hin.push("RotatedSPOs", false);

  bool grp_hist_exists   = hin.is_group("rotation_history");
  bool grp_global_exists = hin.is_group("rotation_global");
  if (!grp_hist_exists && !grp_global_exists)
    app_warning() << "Rotation parameters not found in VP file";


  if (grp_global_exists)
  {
    hin.push("rotation_global", false);
    std::string rot_global_name = std::string("rotation_global_") + SPOSet::getName();

    std::vector<int> sizes(1);
    if (!hin.getShape<RealType>(rot_global_name, sizes))
      throw std::runtime_error("Failed to read rotation_global in VP file");

    int nparam_full_actual = sizes[0];
    int nparam_full        = myVarsFull.size();

    if (nparam_full != nparam_full_actual)
    {
      std::ostringstream tmp_err;
      tmp_err << "Expected number of full rotation parameters (" << nparam_full << ") does not match number in file ("
              << nparam_full_actual << ")";
      throw std::runtime_error(tmp_err.str());
    }
    std::vector<RealType> full_params(nparam_full);
    hin.read(full_params, rot_global_name);
    for (int i = 0; i < nparam_full; i++)
      myVarsFull[i] = full_params[i];

    hin.pop();

    applyFullRotation(full_params, true);
  }
  else if (grp_hist_exists)
  {
    hin.push("rotation_history", false);
    std::string rot_hist_name = std::string("rotation_history_") + SPOSet::getName();
    std::vector<int> sizes(2);
    if (!hin.getShape<RealType>(rot_hist_name, sizes))
      throw std::runtime_error("Failed to read rotation history in VP file");

    int rows = sizes[0];
    int cols = sizes[1];
    history_params_.resize(rows);
    Matrix<RealType> tmp(rows, cols);
    hin.read(tmp, rot_hist_name);
    for (size_t i = 0; i < rows; i++)
    {
      history_params_[i].resize(cols);
      for (size_t j = 0; j < cols; j++)
        history_params_[i][j] = tmp(i, j);
    }

    hin.pop();

    applyRotationHistory();
  }

  hin.push("rotation_params", false);
  std::string rot_param_name = std::string("rotation_params_") + SPOSet::getName();

  std::vector<int> sizes(1);
  if (!hin.getShape<RealType>(rot_param_name, sizes))
    throw std::runtime_error("Failed to read rotation_params in VP file");

  int nparam_actual = sizes[0];
  int nparam        = myVars.size();
  if (nparam != nparam_actual)
  {
    std::ostringstream tmp_err;
    tmp_err << "Expected number of rotation parameters (" << nparam << ") does not match number in file ("
            << nparam_actual << ")";
    throw std::runtime_error(tmp_err.str());
  }

  std::vector<RealType> params(nparam);
  hin.read(params, rot_param_name);
  for (int i = 0; i < nparam; i++)
    myVars[i] = params[i];

  hin.pop();

  hin.pop();
}

void RotatedSPOs::buildOptVariables(const size_t nel)
{
#if !defined(QMC_COMPLEX)
  /* Only rebuild optimized variables if more after-rotation orbitals are needed
   * Consider ROHF, there is only one set of SPO for both spin up and down Nup > Ndown.
   * nel_major_ will be set Nup.
   *
   * Use the size of myVars as a flag to avoid building the rotation parameters again
   * when a clone is made (the DiracDeterminant constructor calls buildOptVariables)
   */
  if (nel > nel_major_ && myVars.size() == 0)
  {
    nel_major_ = nel;

    const size_t nmo = Phi->getOrbitalSetSize();

    // create active rotation parameter indices
    RotationIndices created_m_act_rot_inds;

    RotationIndices created_full_rot_inds;
    if (use_global_rot_)
      createRotationIndicesFull(nel, nmo, created_full_rot_inds);

    createRotationIndices(nel, nmo, created_m_act_rot_inds);

    buildOptVariables(created_m_act_rot_inds, created_full_rot_inds);
  }
#endif
}

void RotatedSPOs::buildOptVariables(const RotationIndices& rotations, const RotationIndices& full_rotations)
{
#if !defined(QMC_COMPLEX)
  const size_t nmo = Phi->getOrbitalSetSize();

  // create active rotations
  m_act_rot_inds = rotations;

  if (use_global_rot_)
    m_full_rot_inds = full_rotations;

  if (use_global_rot_)
    app_log() << "Orbital rotation using global rotation" << std::endl;
  else
    app_log() << "Orbital rotation using history" << std::endl;

  // This will add the orbital rotation parameters to myVars
  // and will also read in initial parameter values supplied in input file
  int p, q;
  int nparams_active = m_act_rot_inds.size();

  app_log() << "nparams_active: " << nparams_active << " params2.size(): " << params.size() << std::endl;
  if (params_supplied)
    if (nparams_active != params.size())
      throw std::runtime_error(
          "The number of supplied orbital rotation parameters does not match number prdouced by the slater "
          "expansion. \n");

  myVars.clear();
  for (int i = 0; i < nparams_active; i++)
  {
    p = m_act_rot_inds[i].first;
    q = m_act_rot_inds[i].second;
    std::stringstream sstr;
    sstr << my_name_ << "_orb_rot_" << (p < 10 ? "0" : "") << (p < 100 ? "0" : "") << (p < 1000 ? "0" : "") << p << "_"
         << (q < 10 ? "0" : "") << (q < 100 ? "0" : "") << (q < 1000 ? "0" : "") << q;

    // If the user input parameters, use those. Otherwise, initialize the parameters to zero
    if (params_supplied)
    {
      myVars.insert(sstr.str(), params[i]);
    }
    else
    {
      myVars.insert(sstr.str(), 0.0);
    }
  }

  if (use_global_rot_)
  {
    myVarsFull.clear();
    for (int i = 0; i < m_full_rot_inds.size(); i++)
    {
      p = m_full_rot_inds[i].first;
      q = m_full_rot_inds[i].second;
      std::stringstream sstr;
      sstr << my_name_ << "_orb_rot_" << (p < 10 ? "0" : "") << (p < 100 ? "0" : "") << (p < 1000 ? "0" : "") << p
           << "_" << (q < 10 ? "0" : "") << (q < 100 ? "0" : "") << (q < 1000 ? "0" : "") << q;

      if (params_supplied && i < m_act_rot_inds.size())
        myVarsFull.insert(sstr.str(), params[i]);
      else
        myVarsFull.insert(sstr.str(), 0.0);
    }
  }


  //Printing the parameters
  if (true)
  {
    app_log() << std::string(16, ' ') << "Parameter name" << std::string(15, ' ') << "Value\n";
    myVars.print(app_log());
  }

  if (params_supplied)
  {
    std::vector<RealType> param(m_act_rot_inds.size());
    for (int i = 0; i < m_act_rot_inds.size(); i++)
      param[i] = myVars[i];
    apply_rotation(param, false);
  }
#endif
}

void RotatedSPOs::apply_rotation(const std::vector<RealType>& param, bool use_stored_copy)
{
  assert(param.size() == m_act_rot_inds.size());

  const size_t nmo = Phi->getOrbitalSetSize();
  ValueMatrix rot_mat(nmo, nmo);

  constructAntiSymmetricMatrix(m_act_rot_inds, param, rot_mat);

  /*
    rot_mat is now an anti-hermitian matrix. Now we convert
    it into a unitary matrix via rot_mat = exp(-rot_mat). 
    Finally, apply unitary matrix to orbs.
  */
  exponentiate_antisym_matrix(rot_mat);
  {
    ScopedTimer local(apply_rotation_timer_);
    Phi->applyRotation(rot_mat, use_stored_copy);
  }
}

void RotatedSPOs::applyDeltaRotation(const std::vector<RealType>& delta_param,
                                     const std::vector<RealType>& old_param,
                                     std::vector<RealType>& new_param)
{
  const size_t nmo = Phi->getOrbitalSetSize();
  ValueMatrix new_rot_mat(nmo, nmo);
  constructDeltaRotation(delta_param, old_param, m_act_rot_inds, m_full_rot_inds, new_param, new_rot_mat);

  {
    ScopedTimer local(apply_rotation_timer_);
    Phi->applyRotation(new_rot_mat, true);
  }
}

void RotatedSPOs::constructDeltaRotation(const std::vector<RealType>& delta_param,
                                         const std::vector<RealType>& old_param,
                                         const RotationIndices& act_rot_inds,
                                         const RotationIndices& full_rot_inds,
                                         std::vector<RealType>& new_param,
                                         ValueMatrix& new_rot_mat)
{
  assert(delta_param.size() == act_rot_inds.size());
  assert(old_param.size() == full_rot_inds.size());
  assert(new_param.size() == full_rot_inds.size());

  const size_t nmo = new_rot_mat.rows();
  assert(new_rot_mat.rows() == new_rot_mat.cols());

  ValueMatrix old_rot_mat(nmo, nmo);

  constructAntiSymmetricMatrix(full_rot_inds, old_param, old_rot_mat);
  exponentiate_antisym_matrix(old_rot_mat);

  ValueMatrix delta_rot_mat(nmo, nmo);

  constructAntiSymmetricMatrix(act_rot_inds, delta_param, delta_rot_mat);
  exponentiate_antisym_matrix(delta_rot_mat);

  // Apply delta rotation to old rotation.
  BLAS::gemm('N', 'N', nmo, nmo, nmo, 1.0, delta_rot_mat.data(), nmo, old_rot_mat.data(), nmo, 0.0, new_rot_mat.data(),
             nmo);

  ValueMatrix log_rot_mat(nmo, nmo);
  log_antisym_matrix(new_rot_mat, log_rot_mat);
  extractParamsFromAntiSymmetricMatrix(full_rot_inds, log_rot_mat, new_param);
}

void RotatedSPOs::applyFullRotation(const std::vector<RealType>& full_param, bool use_stored_copy)
{
  assert(full_param.size() == m_full_rot_inds.size());

  const size_t nmo = Phi->getOrbitalSetSize();
  ValueMatrix rot_mat(nmo, nmo);
  rot_mat = ValueType(0);

  constructAntiSymmetricMatrix(m_full_rot_inds, full_param, rot_mat);

  /*
    rot_mat is now an anti-hermitian matrix. Now we convert
    it into a unitary matrix via rot_mat = exp(-rot_mat).
    Finally, apply unitary matrix to orbs.
  */
  exponentiate_antisym_matrix(rot_mat);
  Phi->applyRotation(rot_mat, use_stored_copy);
}

void RotatedSPOs::applyRotationHistory()
{
  for (auto delta_param : history_params_)
  {
    apply_rotation(delta_param, false);
  }
}

// compute exponential of a real, antisymmetric matrix by diagonalizing and exponentiating eigenvalues
void RotatedSPOs::exponentiate_antisym_matrix(ValueMatrix& mat)
{
  const int n = mat.rows();
  std::vector<std::complex<RealType>> mat_h(n * n, 0);
  std::vector<RealType> eval(n, 0);
  std::vector<std::complex<RealType>> work(2 * n, 0);
  std::vector<RealType> rwork(3 * n, 0);
  std::vector<std::complex<RealType>> mat_d(n * n, 0);
  std::vector<std::complex<RealType>> mat_t(n * n, 0);
  // exponentiating e^X = e^iY (Y hermitian)
  // i(-iX) = X, so -iX is hermitian
  // diagonalize -iX = UDU^T, exponentiate e^iD, and return U e^iD U^T
  // construct hermitian analogue of mat by multiplying by -i
  for (int i = 0; i < n; ++i)
  {
    for (int j = i; j < n; ++j)
    {
      mat_h[i + n * j] = std::complex<RealType>(0, -1.0 * mat[j][i]);
      mat_h[j + n * i] = std::complex<RealType>(0, 1.0 * mat[j][i]);
    }
  }
  // diagonalize the matrix
  char JOBZ('V');
  char UPLO('U');
  int N(n);
  int LDA(n);
  int LWORK(2 * n);
  int info = 0;
  LAPACK::heev(JOBZ, UPLO, N, &mat_h.at(0), LDA, &eval.at(0), &work.at(0), LWORK, &rwork.at(0), info);
  if (info != 0)
  {
    std::ostringstream msg;
    msg << "heev failed with info = " << info << " in RotatedSPOs::exponentiate_antisym_matrix";
    throw std::runtime_error(msg.str());
  }
  // iterate through diagonal matrix, exponentiate terms
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      mat_d[i + j * n] = (i == j) ? std::exp(std::complex<RealType>(0.0, eval[i])) : std::complex<RealType>(0.0, 0.0);
    }
  }
  // perform matrix multiplication
  // assume row major
  BLAS::gemm('N', 'C', n, n, n, std::complex<RealType>(1.0, 0), &mat_d.at(0), n, &mat_h.at(0), n,
             std::complex<RealType>(0.0, 0.0), &mat_t.at(0), n);
  BLAS::gemm('N', 'N', n, n, n, std::complex<RealType>(1.0, 0), &mat_h.at(0), n, &mat_t.at(0), n,
             std::complex<RealType>(0.0, 0.0), &mat_d.at(0), n);
  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    {
      if (mat_d[i + n * j].imag() > 1e-12)
      {
        app_log() << "warning: large imaginary value in orbital rotation matrix: (i,j) = (" << i << "," << j
                  << "), im = " << mat_d[i + n * j].imag() << std::endl;
      }
      mat[j][i] = mat_d[i + n * j].real();
    }
}

void RotatedSPOs::log_antisym_matrix(const ValueMatrix& mat, ValueMatrix& output)
{
  const int n = mat.rows();
  std::vector<RealType> mat_h(n * n, 0);
  std::vector<RealType> eval_r(n, 0);
  std::vector<RealType> eval_i(n, 0);
  std::vector<RealType> mat_l(n * n, 0);
  std::vector<RealType> work(4 * n, 0);

  std::vector<std::complex<RealType>> mat_cd(n * n, 0);
  std::vector<std::complex<RealType>> mat_cl(n * n, 0);
  std::vector<std::complex<RealType>> mat_ch(n * n, 0);

  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
      mat_h[i + n * j] = mat[i][j];

  // diagonalize the matrix
  char JOBL('V');
  char JOBR('N');
  int N(n);
  int LDA(n);
  int LWORK(4 * n);
  int info = 0;
  LAPACK::geev(&JOBL, &JOBR, &N, &mat_h.at(0), &LDA, &eval_r.at(0), &eval_i.at(0), &mat_l.at(0), &LDA, nullptr, &LDA,
               &work.at(0), &LWORK, &info);
  if (info != 0)
  {
    std::ostringstream msg;
    msg << "heev failed with info = " << info << " in RotatedSPOs::log_antisym_matrix";
    throw std::runtime_error(msg.str());
  }

  // iterate through diagonal matrix, take log
  for (int i = 0; i < n; ++i)
  {
    for (int j = 0; j < n; ++j)
    {
      auto tmp = (i == j) ? std::log(std::complex<RealType>(eval_r[i], eval_i[i])) : std::complex<RealType>(0.0, 0.0);
      mat_cd[i + j * n] = tmp;

      if (eval_i[j] > 0.0)
      {
        mat_cl[i + j * n]       = std::complex<RealType>(mat_l[i + j * n], mat_l[i + (j + 1) * n]);
        mat_cl[i + (j + 1) * n] = std::complex<RealType>(mat_l[i + j * n], -mat_l[i + (j + 1) * n]);
      }
      else if (!(eval_i[j] < 0.0))
      {
        mat_cl[i + j * n] = std::complex<RealType>(mat_l[i + j * n], 0.0);
      }
    }
  }

  RealType one(1.0);
  RealType zero(0.0);
  BLAS::gemm('N', 'N', n, n, n, one, &mat_cl.at(0), n, &mat_cd.at(0), n, zero, &mat_ch.at(0), n);
  BLAS::gemm('N', 'C', n, n, n, one, &mat_ch.at(0), n, &mat_cl.at(0), n, zero, &mat_cd.at(0), n);


  for (int i = 0; i < n; ++i)
    for (int j = 0; j < n; ++j)
    {
      if (mat_cd[i + n * j].imag() > 1e-12)
      {
        app_log() << "warning: large imaginary value in antisymmetric matrix: (i,j) = (" << i << "," << j
                  << "), im = " << mat_cd[i + n * j].imag() << std::endl;
      }
      output[i][j] = mat_cd[i + n * j].real();
    }
}

void RotatedSPOs::evaluateDerivRatios(const VirtualParticleSet& VP,
                                      const opt_variables_type& optvars,
                                      ValueVector& psi,
                                      const ValueVector& psiinv,
                                      std::vector<ValueType>& ratios,
                                      Matrix<ValueType>& dratios,
                                      int FirstIndex,
                                      int LastIndex)
{
  Phi->evaluateDetRatios(VP, psi, psiinv, ratios);

  const size_t nel = LastIndex - FirstIndex;
  const size_t nmo = Phi->getOrbitalSetSize();

  psiM_inv.resize(nel, nel);
  psiM_all.resize(nel, nmo);
  dpsiM_all.resize(nel, nmo);
  d2psiM_all.resize(nel, nmo);

  psiM_inv   = 0;
  psiM_all   = 0;
  dpsiM_all  = 0;
  d2psiM_all = 0;

  const ParticleSet& P = VP.getRefPS();
  int iel              = VP.refPtcl;

  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_all, dpsiM_all, d2psiM_all);

  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nel; j++)
      psiM_inv(i, j) = psiM_all(i, j);

  Invert(psiM_inv.data(), nel, nel);

  const ValueType* const A(psiM_all.data());
  const ValueType* const Ainv(psiM_inv.data());
  SPOSet::ValueMatrix T_orig;
  T_orig.resize(nel, nmo);

  BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), A, nmo, Ainv, nel, ValueType(0.0), T_orig.data(), nmo);

  SPOSet::ValueMatrix T;
  T.resize(nel, nmo);

  ValueVector tmp_psi;
  tmp_psi.resize(nmo);

  for (int iat = 0; iat < VP.getTotalNum(); iat++)
  {
    Phi->evaluateValue(VP, iat, tmp_psi);

    for (int j = 0; j < nmo; j++)
      psiM_all(iel - FirstIndex, j) = tmp_psi[j];

    for (int i = 0; i < nel; i++)
      for (int j = 0; j < nel; j++)
        psiM_inv(i, j) = psiM_all(i, j);

    Invert(psiM_inv.data(), nel, nel);

    const ValueType* const A(psiM_all.data());
    const ValueType* const Ainv(psiM_inv.data());

    // The matrix A is rectangular.  Ainv is the inverse of the square part of the matrix.
    // The multiply of Ainv and the square part of A is just the identity.
    // This multiply could be reduced to Ainv and the non-square part of A.
    BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), A, nmo, Ainv, nel, ValueType(0.0), T.data(), nmo);

    for (int i = 0; i < m_act_rot_inds.size(); i++)
    {
      int kk = myVars.where(i);
      if (kk >= 0)
      {
        const int p      = m_act_rot_inds.at(i).first;
        const int q      = m_act_rot_inds.at(i).second;
        dratios(iat, kk) = T(p, q) - T_orig(p, q); // dratio size is (nknot, num_vars)
      }
    }
  }
}

void RotatedSPOs::evaluateDerivativesWF(ParticleSet& P,
                                        const opt_variables_type& optvars,
                                        Vector<ValueType>& dlogpsi,
                                        int FirstIndex,
                                        int LastIndex)
{
  const size_t nel = LastIndex - FirstIndex;
  const size_t nmo = Phi->getOrbitalSetSize();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PART1

  psiM_inv.resize(nel, nel);
  psiM_all.resize(nel, nmo);
  dpsiM_all.resize(nel, nmo);
  d2psiM_all.resize(nel, nmo);

  psiM_inv   = 0;
  psiM_all   = 0;
  dpsiM_all  = 0;
  d2psiM_all = 0;

  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_all, dpsiM_all, d2psiM_all);

  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nel; j++)
      psiM_inv(i, j) = psiM_all(i, j);

  Invert(psiM_inv.data(), nel, nel);

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PART2
  const ValueType* const A(psiM_all.data());
  const ValueType* const Ainv(psiM_inv.data());
  SPOSet::ValueMatrix T;
  T.resize(nel, nmo);

  BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), A, nmo, Ainv, nel, ValueType(0.0), T.data(), nmo);

  for (int i = 0; i < m_act_rot_inds.size(); i++)
  {
    int kk = myVars.where(i);
    if (kk >= 0)
    {
      const int p = m_act_rot_inds.at(i).first;
      const int q = m_act_rot_inds.at(i).second;
      dlogpsi[kk] = T(p, q);
    }
  }
}

void RotatedSPOs::evaluateDerivatives(ParticleSet& P,
                                      const opt_variables_type& optvars,
                                      Vector<ValueType>& dlogpsi,
                                      Vector<ValueType>& dhpsioverpsi,
                                      const int& FirstIndex,
                                      const int& LastIndex)
{
  const size_t nel = LastIndex - FirstIndex;
  const size_t nmo = Phi->getOrbitalSetSize();

  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PART1
  myG_temp.resize(nel);
  myG_J.resize(nel);
  myL_temp.resize(nel);
  myL_J.resize(nel);

  myG_temp = 0;
  myG_J    = 0;
  myL_temp = 0;
  myL_J    = 0;

  Bbar.resize(nel, nmo);
  psiM_inv.resize(nel, nel);
  psiM_all.resize(nel, nmo);
  dpsiM_all.resize(nel, nmo);
  d2psiM_all.resize(nel, nmo);

  Bbar       = 0;
  psiM_inv   = 0;
  psiM_all   = 0;
  dpsiM_all  = 0;
  d2psiM_all = 0;


  Phi->evaluate_notranspose(P, FirstIndex, LastIndex, psiM_all, dpsiM_all, d2psiM_all);

  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nel; j++)
      psiM_inv(i, j) = psiM_all(i, j);

  Invert(psiM_inv.data(), nel, nel);

  //current value of Gradient and Laplacian
  // gradient components
  for (int a = 0; a < nel; a++)
    for (int i = 0; i < nel; i++)
      for (int k = 0; k < 3; k++)
        myG_temp[a][k] += psiM_inv(i, a) * dpsiM_all(a, i)[k];
  // laplacian components
  for (int a = 0; a < nel; a++)
  {
    for (int i = 0; i < nel; i++)
      myL_temp[a] += psiM_inv(i, a) * d2psiM_all(a, i);
  }

  // calculation of myG_J which will be used to represent \frac{\nabla\psi_{J}}{\psi_{J}}
  // calculation of myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}}
  // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold
  for (int a = 0, iat = FirstIndex; a < nel; a++, iat++)
  {
    myG_J[a] = (P.G[iat] - myG_temp[a]);
    myL_J[a] = (P.L[iat] + dot(P.G[iat], P.G[iat]) - myL_temp[a]);
  }
  //possibly replace wit BLAS calls
  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nmo; j++)
      Bbar(i, j) = d2psiM_all(i, j) + 2 * dot(myG_J[i], dpsiM_all(i, j)) + myL_J[i] * psiM_all(i, j);


  //~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~PART2
  const ValueType* const A(psiM_all.data());
  const ValueType* const Ainv(psiM_inv.data());
  const ValueType* const B(Bbar.data());
  SPOSet::ValueMatrix T;
  SPOSet::ValueMatrix Y1;
  SPOSet::ValueMatrix Y2;
  SPOSet::ValueMatrix Y3;
  SPOSet::ValueMatrix Y4;
  T.resize(nel, nmo);
  Y1.resize(nel, nel);
  Y2.resize(nel, nmo);
  Y3.resize(nel, nmo);
  Y4.resize(nel, nmo);


  BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), A, nmo, Ainv, nel, ValueType(0.0), T.data(), nmo);
  BLAS::gemm('N', 'N', nel, nel, nel, ValueType(1.0), B, nmo, Ainv, nel, ValueType(0.0), Y1.data(), nel);
  BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), T.data(), nmo, Y1.data(), nel, ValueType(0.0), Y2.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nel, ValueType(1.0), B, nmo, Ainv, nel, ValueType(0.0), Y3.data(), nmo);

  //possibly replace with BLAS call
  Y4 = Y3 - Y2;

  for (int i = 0; i < m_act_rot_inds.size(); i++)
  {
    int kk = myVars.where(i);
    if (kk >= 0)
    {
      const int p = m_act_rot_inds.at(i).first;
      const int q = m_act_rot_inds.at(i).second;
      dlogpsi[kk] += T(p, q);
      dhpsioverpsi[kk] += ValueType(-0.5) * Y4(p, q);
    }
  }
}

void RotatedSPOs::evaluateDerivatives(ParticleSet& P,
                                      const opt_variables_type& optvars,
                                      Vector<ValueType>& dlogpsi,
                                      Vector<ValueType>& dhpsioverpsi,
                                      const ValueType& psiCurrent,
                                      const std::vector<ValueType>& Coeff,
                                      const std::vector<size_t>& C2node_up,
                                      const std::vector<size_t>& C2node_dn,
                                      const ValueVector& detValues_up,
                                      const ValueVector& detValues_dn,
                                      const GradMatrix& grads_up,
                                      const GradMatrix& grads_dn,
                                      const ValueMatrix& lapls_up,
                                      const ValueMatrix& lapls_dn,
                                      const ValueMatrix& M_up,
                                      const ValueMatrix& M_dn,
                                      const ValueMatrix& Minv_up,
                                      const ValueMatrix& Minv_dn,
                                      const GradMatrix& B_grad,
                                      const ValueMatrix& B_lapl,
                                      const std::vector<int>& detData_up,
                                      const size_t N1,
                                      const size_t N2,
                                      const size_t NP1,
                                      const size_t NP2,
                                      const std::vector<std::vector<int>>& lookup_tbl)
{
  bool recalculate(false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (optvars.recompute(kk))
      recalculate = true;
  }
  if (recalculate)
  {
    ParticleSet::ParticleGradient myG_temp, myG_J;
    ParticleSet::ParticleLaplacian myL_temp, myL_J;
    const int NP = P.getTotalNum();
    myG_temp.resize(NP);
    myG_temp = 0.0;
    myL_temp.resize(NP);
    myL_temp = 0.0;
    myG_J.resize(NP);
    myG_J = 0.0;
    myL_J.resize(NP);
    myL_J            = 0.0;
    const size_t nmo = Phi->getOrbitalSetSize();
    const size_t nel = P.last(0) - P.first(0);

    const RealType* restrict C_p = Coeff.data();
    for (int i = 0; i < Coeff.size(); i++)
    {
      const size_t upC     = C2node_up[i];
      const size_t dnC     = C2node_dn[i];
      const ValueType tmp1 = C_p[i] * detValues_dn[dnC];
      const ValueType tmp2 = C_p[i] * detValues_up[upC];
      for (size_t k = 0, j = N1; k < NP1; k++, j++)
      {
        myG_temp[j] += tmp1 * grads_up(upC, k);
        myL_temp[j] += tmp1 * lapls_up(upC, k);
      }
      for (size_t k = 0, j = N2; k < NP2; k++, j++)
      {
        myG_temp[j] += tmp2 * grads_dn(dnC, k);
        myL_temp[j] += tmp2 * lapls_dn(dnC, k);
      }
    }

    myG_temp *= (1 / psiCurrent);
    myL_temp *= (1 / psiCurrent);

    // calculation of myG_J which will be used to represent \frac{\nabla\psi_{J}}{\psi_{J}}
    // calculation of myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}}
    // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold
    for (int iat = 0; iat < (myL_temp.size()); iat++)
    {
      myG_J[iat] = (P.G[iat] - myG_temp[iat]);
      myL_J[iat] = (P.L[iat] + dot(P.G[iat], P.G[iat]) - myL_temp[iat]);
    }


    table_method_eval(dlogpsi, dhpsioverpsi, myL_J, myG_J, nel, nmo, psiCurrent, Coeff, C2node_up, C2node_dn,
                      detValues_up, detValues_dn, grads_up, grads_dn, lapls_up, lapls_dn, M_up, M_dn, Minv_up, Minv_dn,
                      B_grad, B_lapl, detData_up, N1, N2, NP1, NP2, lookup_tbl);
  }
}


void RotatedSPOs::evaluateDerivativesWF(ParticleSet& P,
                                        const opt_variables_type& optvars,
                                        Vector<ValueType>& dlogpsi,
                                        const QTFull::ValueType& psiCurrent,
                                        const std::vector<ValueType>& Coeff,
                                        const std::vector<size_t>& C2node_up,
                                        const std::vector<size_t>& C2node_dn,
                                        const ValueVector& detValues_up,
                                        const ValueVector& detValues_dn,
                                        const ValueMatrix& M_up,
                                        const ValueMatrix& M_dn,
                                        const ValueMatrix& Minv_up,
                                        const ValueMatrix& Minv_dn,
                                        const std::vector<int>& detData_up,
                                        const std::vector<std::vector<int>>& lookup_tbl)
{
  bool recalculate(false);
  for (int k = 0; k < myVars.size(); ++k)
  {
    int kk = myVars.where(k);
    if (kk < 0)
      continue;
    if (optvars.recompute(kk))
      recalculate = true;
  }
  if (recalculate)
  {
    const size_t nmo = Phi->getOrbitalSetSize();
    const size_t nel = P.last(0) - P.first(0);

    table_method_evalWF(dlogpsi, nel, nmo, psiCurrent, Coeff, C2node_up, C2node_dn, detValues_up, detValues_dn, M_up,
                        M_dn, Minv_up, Minv_dn, detData_up, lookup_tbl);
  }
}

void RotatedSPOs::table_method_eval(Vector<ValueType>& dlogpsi,
                                    Vector<ValueType>& dhpsioverpsi,
                                    const ParticleSet::ParticleLaplacian& myL_J,
                                    const ParticleSet::ParticleGradient& myG_J,
                                    const size_t nel,
                                    const size_t nmo,
                                    const ValueType& psiCurrent,
                                    const std::vector<RealType>& Coeff,
                                    const std::vector<size_t>& C2node_up,
                                    const std::vector<size_t>& C2node_dn,
                                    const ValueVector& detValues_up,
                                    const ValueVector& detValues_dn,
                                    const GradMatrix& grads_up,
                                    const GradMatrix& grads_dn,
                                    const ValueMatrix& lapls_up,
                                    const ValueMatrix& lapls_dn,
                                    const ValueMatrix& M_up,
                                    const ValueMatrix& M_dn,
                                    const ValueMatrix& Minv_up,
                                    const ValueMatrix& Minv_dn,
                                    const GradMatrix& B_grad,
                                    const ValueMatrix& B_lapl,
                                    const std::vector<int>& detData_up,
                                    const size_t N1,
                                    const size_t N2,
                                    const size_t NP1,
                                    const size_t NP2,
                                    const std::vector<std::vector<int>>& lookup_tbl)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUIDE TO THE MATICES BEING BUILT
----------------------------------------------
The idea here is that there is a loop over all unique determinants. For each determiant the table method is employed to calculate the contributions to the parameter derivatives (dhpsioverpsi/dlogpsi)

  loop through unquie determinants 
    loop through parameters
      evaluate contributaion to dlogpsi and dhpsioverpsi
\noindent 

  BLAS GUIDE  for matrix multiplication of  [  alpha * A.B + beta * C = C ]
  Matrix A is of dimensions a1,a2 and Matrix B is b1,b2   in which a2=b1
  The BLAS command is as follows...

 BLAS::gemm('N','N', b2, a1, a2 ,alpha, B, b2, A, a2, beta, C, b2);

Below is a human readable format for the matrix multiplications performed below...

This notation is inspired by http://dx.doi.org/10.1063/1.4948778
\newline
\hfill\break
$
    A_{i,j}=\phi_j(r_{i}) \\
    T = A^{-1} \widetilde{A} \\
    B_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i}) \\
    \hat{O_{I}} = \hat{O}D_{I} \\
    D_{I}=det(A_{I}) \newline 
    \psi_{MS} = \sum_{I=0} C_{I} D_{I\uparrow}D_{I\downarrow} \\
    \Psi_{total} = \psi_{J}\psi_{MS} \\
    \alpha_{I} = P^{T}_{I}TQ_{I} \\
    M_{I} = P^{T}_{I} \widetilde{M} Q_{I} = P^{T}_{I} (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )Q_{I} \\
$
\newline
There are three constants I use in the expressions for dhpsioverpsi and dlogpsi
\newline
\hfill\break
$
  const0 = C_{0}*det(A_{0\downarrow})+\sum_{I=1} C_{I}*det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{I=1} C_{I}*\hat{O}det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const2 = \sum_{I=1} C_{I}*det(A_{I\downarrow})* Tr[\alpha_{I}^{-1}M_{I}]*det(\alpha_{I}) \\
$
\newline
Below is a translation of the shorthand I use to represent matrices independent of ``excitation matrix".
\newline
\hfill\break
$
    Y_{1} =  A^{-1}B   \\
    Y_{2} = A^{-1}BA^{-1}\widetilde{A} \\
    Y_{3} = A^{-1}\widetilde{B} \\
    Y_{4} = \widetilde{M} = (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )\\
$
\newline
Below is a translation of the shorthand I use to represent matrices dependent on ``excitation" with respect to the reference Matrix and sums of matrices. Above this line I have represented these excitation matrices with a subscript ``I" but from this point on The subscript will be omitted and it is clear that whenever a matrix depends on $P^{T}_I$ and $Q_{I}$ that this is an excitation matrix. The reference matrix is always $A_{0}$ and is always the Hartree Fock Matrix.
\newline
\hfill\break
$
    Y_{5} = TQ \\
    Y_{6} = (P^{T}TQ)^{-1} = \alpha_{I}^{-1}\\
    Y_{7} = \alpha_{I}^{-1} P^{T} \\
    Y_{11} = \widetilde{M}Q \\
    Y_{23} = P^{T}\widetilde{M}Q \\
    Y_{24} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q \\
    Y_{25} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1} \\
    Y_{26} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1}P^{T}\\
$
\newline
So far you will notice that I have not included up or down arrows to specify what spin the matrices are of. This is because we are calculating the derivative of all up or all down spin orbital rotation parameters at a time. If we are finding the up spin derivatives then any term that is down spin will be constant. The following assumes that we are taking up-spin MO rotation parameter derivatives. Of course the down spin expression can be retrieved by swapping the up and down arrows. I have dubbed any expression with lowercase p prefix as a "precursor" to an expression actually used...
\newline
\hfill\break
$
    \dot{C_{I}} = C_{I}*det(A_{I\downarrow})\\
    \ddot{C_{I}} = C_{I}*\hat{O}det(A_{I\downarrow}) \\
    pK1 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}) \\
    pK2 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK3 = \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK4 = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK5 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T}) \\
$
\newline
Now these p matrices will be used to make various expressions via BLAS commands.
\newline
\hfill\break
$
    K1T = const0^{-1}*pK1.T =const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}T) \\
    TK1T = T.K1T = const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    TK2AiB = T.K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    K2XA =  const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}X\widetilde{A})\\
    TK2XA = T.K2XA = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}X\widetilde{A})\\ \\
    K2T = \frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK2T = T.K2T =\frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\
    MK2T = \frac{const0}{const1} Y_{4}.K2T= const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (\widetilde{M}Q\alpha_{I}^{-1}P^{T}T)\\ \\
    K3T = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK3T = T.K3T  = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK4T = T.K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\ \\
    K5T =  const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
    TK5T = T.K5T  = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (T Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
$
\newline
Now with all these matrices and constants the expressions of dhpsioverpsi and dlogpsi can be created.




In addition I will be using a special generalization of the kinetic operator which I will denote as O. Our Slater matrix with the special O operator applied to each element will be called B_bar

$
``Bbar"_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i})
$
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
  ValueMatrix Table;
  ValueMatrix Bbar;
  ValueMatrix Y1, Y2, Y3, Y4, Y5, Y6, Y7, Y11, Y23, Y24, Y25, Y26;
  ValueMatrix pK1, K1T, TK1T, pK2, K2AiB, TK2AiB, K2XA, TK2XA, K2T, TK2T, MK2T, pK3, K3T, TK3T, pK5, K5T, TK5T;

  Table.resize(nel, nmo);

  Bbar.resize(nel, nmo);

  Y1.resize(nel, nel);
  Y2.resize(nel, nmo);
  Y3.resize(nel, nmo);
  Y4.resize(nel, nmo);

  pK1.resize(nmo, nel);
  K1T.resize(nmo, nmo);
  TK1T.resize(nel, nmo);

  pK2.resize(nmo, nel);
  K2AiB.resize(nmo, nmo);
  TK2AiB.resize(nel, nmo);
  K2XA.resize(nmo, nmo);
  TK2XA.resize(nel, nmo);
  K2T.resize(nmo, nmo);
  TK2T.resize(nel, nmo);
  MK2T.resize(nel, nmo);

  pK3.resize(nmo, nel);
  K3T.resize(nmo, nmo);
  TK3T.resize(nel, nmo);

  pK5.resize(nmo, nel);
  K5T.resize(nmo, nmo);
  TK5T.resize(nel, nmo);

  const int parameters_size(m_act_rot_inds.size());
  const int parameter_start_index(0);

  const size_t num_unique_up_dets(detValues_up.size());
  const size_t num_unique_dn_dets(detValues_dn.size());

  const RealType* restrict cptr = Coeff.data();
  const size_t nc               = Coeff.size();
  const size_t* restrict upC(C2node_up.data());
  const size_t* restrict dnC(C2node_dn.data());
  //B_grad holds the gradient operator
  //B_lapl holds the laplacian operator
  //B_bar will hold our special O operator

  const int offset1(N1);
  const int offset2(N2);
  const int NPother(NP2);

  RealType* T(Table.data());

  //possibly replace wit BLAS calls
  for (int i = 0; i < nel; i++)
    for (int j = 0; j < nmo; j++)
      Bbar(i, j) = B_lapl(i, j) + 2 * dot(myG_J[i + offset1], B_grad(i, j)) + myL_J[i + offset1] * M_up(i, j);

  const RealType* restrict B(Bbar.data());
  const RealType* restrict A(M_up.data());
  const RealType* restrict Ainv(Minv_up.data());
  //IMPORTANT NOTE: THE Dets[0]->psiMinv OBJECT DOES NOT HOLD THE INVERSE IF THE MULTIDIRACDETERMINANTBASE ONLY CONTAINS ONE ELECTRON. NEED A FIX FOR THIS CASE
  // The T matrix should be calculated and stored for use
  // T = A^{-1} \widetilde A
  //REMINDER: that the ValueMatrix "matrix" stores data in a row major order and that BLAS commands assume column major
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), A, nmo, Ainv, nel, RealType(0.0), T, nmo);

  BLAS::gemm('N', 'N', nel, nel, nel, RealType(1.0), B, nmo, Ainv, nel, RealType(0.0), Y1.data(), nel);
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), T, nmo, Y1.data(), nel, RealType(0.0), Y2.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), B, nmo, Ainv, nel, RealType(0.0), Y3.data(), nmo);

  //possibly replace with BLAS call
  Y4 = Y3 - Y2;

  //Need to create the constants: (Oi, const0, const1, const2)to take advantage of minimal BLAS commands;
  //Oi is the special operator applied to the slater matrix "A subscript i" from the total CI expansion
  //\hat{O_{i}} = \hat{O}D_{i} with D_{i}=det(A_{i}) and Multi-Slater component defined as \sum_{i=0} C_{i} D_{i\uparrow}D_{i\downarrow}
  std::vector<RealType> Oi(num_unique_dn_dets);

  for (int index = 0; index < num_unique_dn_dets; index++)
    for (int iat = 0; iat < NPother; iat++)
      Oi[index] += lapls_dn(index, iat) + 2 * dot(grads_dn(index, iat), myG_J[offset2 + iat]) +
          myL_J[offset2 + iat] * detValues_dn[index];

  //const0 = C_{0}*det(A_{0\downarrow})+\sum_{i=1} C_{i}*det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{i=1} C_{i}*\hat{O}det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const2 = \sum_{i=1} C_{i}*det(A_{i\downarrow})* Tr[\alpha_{i}^{-1}M_{i}]*det(\alpha_{i})
  RealType const0(0.0), const1(0.0), const2(0.0);
  for (size_t i = 0; i < nc; ++i)
  {
    const RealType c  = cptr[i];
    const size_t up   = upC[i];
    const size_t down = dnC[i];

    const0 += c * detValues_dn[down] * (detValues_up[up] / detValues_up[0]);
    const1 += c * Oi[down] * (detValues_up[up] / detValues_up[0]);
  }

  std::fill(pK1.begin(), pK1.end(), 0.0);
  std::fill(pK2.begin(), pK2.end(), 0.0);
  std::fill(pK3.begin(), pK3.end(), 0.0);
  std::fill(pK5.begin(), pK5.end(), 0.0);

  //Now we are going to loop through all unique determinants.
  //The few lines above are for the reference matrix contribution.
  //Although I start the loop below from index 0, the loop only performs actions when the index is >= 1
  //the detData object contains all the information about the P^T and Q matrices (projection matrices) needed in the table method
  const int* restrict data_it = detData_up.data();
  for (int index = 0, datum = 0; index < num_unique_up_dets; index++)
  {
    const int k = data_it[datum];

    if (k == 0)
    {
      datum += 3 * k + 1;
    }

    else
    {
      //Number of rows and cols of P^T
      const int prows = k;
      const int pcols = nel;
      //Number of rows and cols of Q
      const int qrows = nmo;
      const int qcols = k;

      Y5.resize(nel, k);
      Y6.resize(k, k);

      //Any matrix multiplication of P^T or Q is simply a projection
      //Explicit matrix multiplication can be avoided; instead column or row copying can be done
      //BlAS::copy(size of col/row being copied,
      //           Matrix pointer + place to begin copying,
      //           storage spacing (number of elements btw next row/col element),
      //           Pointer to resultant matrix + place to begin pasting,
      //           storage spacing of resultant matrix)
      //For example the next 4 lines is the matrix multiplication of T*Q = Y5
      std::fill(Y5.begin(), Y5.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(nel, T + data_it[datum + 1 + k + i], nmo, Y5.data() + i, k);
      }

      std::fill(Y6.begin(), Y6.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y5.data() + (data_it[datum + 1 + i]) * k, 1, (Y6.data() + i * k), 1);
      }


      Vector<ValueType> WS;
      Vector<IndexType> Piv;
      WS.resize(k);
      Piv.resize(k);
      std::complex<RealType> logdet = 0.0;
      InvertWithLog(Y6.data(), k, k, WS.data(), Piv.data(), logdet);

      Y11.resize(nel, k);
      Y23.resize(k, k);
      Y24.resize(k, k);
      Y25.resize(k, k);
      Y26.resize(k, nel);

      std::fill(Y11.begin(), Y11.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(nel, Y4.data() + (data_it[datum + 1 + k + i]), nmo, Y11.data() + i, k);
      }

      std::fill(Y23.begin(), Y23.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y11.data() + (data_it[datum + 1 + i]) * k, 1, (Y23.data() + i * k), 1);
      }

      BLAS::gemm('N', 'N', k, k, k, RealType(1.0), Y23.data(), k, Y6.data(), k, RealType(0.0), Y24.data(), k);
      BLAS::gemm('N', 'N', k, k, k, RealType(1.0), Y6.data(), k, Y24.data(), k, RealType(0.0), Y25.data(), k);


      Y26.resize(k, nel);

      std::fill(Y26.begin(), Y26.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y25.data() + i, k, Y26.data() + (data_it[datum + 1 + i]), nel);
      }


      Y7.resize(k, nel);

      std::fill(Y7.begin(), Y7.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y6.data() + i, k, Y7.data() + (data_it[datum + 1 + i]), nel);
      }

      // c_Tr_AlphaI_MI is a constant contributing to constant const2
      // c_Tr_AlphaI_MI = Tr[\alpha_{I}^{-1}(P^{T}\widetilde{M} Q)]
      RealType c_Tr_AlphaI_MI = 0.0;
      for (int i = 0; i < k; i++)
      {
        c_Tr_AlphaI_MI += Y24(i, i);
      }

      for (int p = 0; p < lookup_tbl[index].size(); p++)
      {
        //el_p is the element position that contains information about the CI coefficient, and det up/dn values associated with the current unique determinant
        const int el_p(lookup_tbl[index][p]);
        const RealType c  = cptr[el_p];
        const size_t up   = upC[el_p];
        const size_t down = dnC[el_p];

        const RealType alpha_1(c * detValues_dn[down] * detValues_up[up] / detValues_up[0] * c_Tr_AlphaI_MI);
        const RealType alpha_2(c * detValues_dn[down] * detValues_up[up] / detValues_up[0]);
        const RealType alpha_3(c * Oi[down] * detValues_up[up] / detValues_up[0]);

        const2 += alpha_1;

        for (int i = 0; i < k; i++)
        {
          BLAS::axpy(nel, alpha_1, Y7.data() + i * nel, 1, pK1.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_2, Y7.data() + i * nel, 1, pK2.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_3, Y7.data() + i * nel, 1, pK3.data() + (data_it[datum + 1 + k + i]) * nel, 1);
          BLAS::axpy(nel, alpha_2, Y26.data() + i * nel, 1, pK5.data() + (data_it[datum + 1 + k + i]) * nel, 1);
        }
      }
      datum += 3 * k + 1;
    }
  }


  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK1.data(), nel, RealType(0.0), K1T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K1T.data(), nmo, T, nmo, RealType(0.0), TK1T.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, Y3.data(), nmo, pK2.data(), nel, RealType(0.0), K2AiB.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2AiB.data(), nmo, T, nmo, RealType(0.0), TK2AiB.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, Y2.data(), nmo, pK2.data(), nel, RealType(0.0), K2XA.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2XA.data(), nmo, T, nmo, RealType(0.0), TK2XA.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, const1 / (const0 * const0), T, nmo, pK2.data(), nel, RealType(0.0), K2T.data(),
             nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K2T.data(), nmo, T, nmo, RealType(0.0), TK2T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, const0 / const1, K2T.data(), nmo, Y4.data(), nmo, RealType(0.0), MK2T.data(),
             nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK3.data(), nel, RealType(0.0), K3T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K3T.data(), nmo, T, nmo, RealType(0.0), TK3T.data(), nmo);

  BLAS::gemm('N', 'N', nmo, nmo, nel, 1.0 / const0, T, nmo, pK5.data(), nel, RealType(0.0), K5T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K5T.data(), nmo, T, nmo, RealType(0.0), TK5T.data(), nmo);


  for (int mu = 0, k = parameter_start_index; k < (parameter_start_index + parameters_size); k++, mu++)
  {
    int kk = myVars.where(k);
    if (kk >= 0)
    {
      const int i(m_act_rot_inds[mu].first), j(m_act_rot_inds[mu].second);
      if (i <= nel - 1 && j > nel - 1)
      {
        dhpsioverpsi[kk] +=
            ValueType(-0.5 * Y4(i, j) -
                      0.5 *
                          (-K5T(i, j) + K5T(j, i) + TK5T(i, j) + K2AiB(i, j) - K2AiB(j, i) - TK2AiB(i, j) - K2XA(i, j) +
                           K2XA(j, i) + TK2XA(i, j) - MK2T(i, j) + K1T(i, j) - K1T(j, i) - TK1T(i, j) -
                           const2 / const1 * K2T(i, j) + const2 / const1 * K2T(j, i) + const2 / const1 * TK2T(i, j) +
                           K3T(i, j) - K3T(j, i) - TK3T(i, j) - K2T(i, j) + K2T(j, i) + TK2T(i, j)));
      }
      else if (i <= nel - 1 && j <= nel - 1)
      {
        dhpsioverpsi[kk] += ValueType(
            -0.5 * (Y4(i, j) - Y4(j, i)) -
            0.5 *
                (-K5T(i, j) + K5T(j, i) + TK5T(i, j) - TK5T(j, i) + K2AiB(i, j) - K2AiB(j, i) - TK2AiB(i, j) +
                 TK2AiB(j, i) - K2XA(i, j) + K2XA(j, i) + TK2XA(i, j) - TK2XA(j, i) - MK2T(i, j) + MK2T(j, i) +
                 K1T(i, j) - K1T(j, i) - TK1T(i, j) + TK1T(j, i) - const2 / const1 * K2T(i, j) +
                 const2 / const1 * K2T(j, i) + const2 / const1 * TK2T(i, j) - const2 / const1 * TK2T(j, i) + K3T(i, j) -
                 K3T(j, i) - TK3T(i, j) + TK3T(j, i) - K2T(i, j) + K2T(j, i) + TK2T(i, j) - TK2T(j, i)));
      }
      else
      {
        dhpsioverpsi[kk] += ValueType(-0.5 *
                                      (-K5T(i, j) + K5T(j, i) + K2AiB(i, j) - K2AiB(j, i) - K2XA(i, j) + K2XA(j, i)

                                       + K1T(i, j) - K1T(j, i) - const2 / const1 * K2T(i, j) +
                                       const2 / const1 * K2T(j, i) + K3T(i, j) - K3T(j, i) - K2T(i, j) + K2T(j, i)));
      }
    }
  }
}

void RotatedSPOs::table_method_evalWF(Vector<ValueType>& dlogpsi,
                                      const size_t nel,
                                      const size_t nmo,
                                      const ValueType& psiCurrent,
                                      const std::vector<RealType>& Coeff,
                                      const std::vector<size_t>& C2node_up,
                                      const std::vector<size_t>& C2node_dn,
                                      const ValueVector& detValues_up,
                                      const ValueVector& detValues_dn,
                                      const ValueMatrix& M_up,
                                      const ValueMatrix& M_dn,
                                      const ValueMatrix& Minv_up,
                                      const ValueMatrix& Minv_dn,
                                      const std::vector<int>& detData_up,
                                      const std::vector<std::vector<int>>& lookup_tbl)
{
  ValueMatrix Table;
  ValueMatrix Y5, Y6, Y7;
  ValueMatrix pK4, K4T, TK4T;

  Table.resize(nel, nmo);

  Bbar.resize(nel, nmo);

  pK4.resize(nmo, nel);
  K4T.resize(nmo, nmo);
  TK4T.resize(nel, nmo);

  const int parameters_size(m_act_rot_inds.size());
  const int parameter_start_index(0);

  const size_t num_unique_up_dets(detValues_up.size());
  const size_t num_unique_dn_dets(detValues_dn.size());

  const RealType* restrict cptr = Coeff.data();
  const size_t nc               = Coeff.size();
  const size_t* restrict upC(C2node_up.data());
  const size_t* restrict dnC(C2node_dn.data());

  RealType* T(Table.data());

  const RealType* restrict A(M_up.data());
  const RealType* restrict Ainv(Minv_up.data());
  //IMPORTANT NOTE: THE Dets[0]->psiMinv OBJECT DOES NOT HOLD THE INVERSE IF THE MULTIDIRACDETERMINANTBASE ONLY CONTAINS ONE ELECTRON. NEED A FIX FOR THIS CASE
  // The T matrix should be calculated and stored for use
  // T = A^{-1} \widetilde A
  //REMINDER: that the ValueMatrix "matrix" stores data in a row major order and that BLAS commands assume column major
  BLAS::gemm('N', 'N', nmo, nel, nel, RealType(1.0), A, nmo, Ainv, nel, RealType(0.0), T, nmo);

  //const0 = C_{0}*det(A_{0\downarrow})+\sum_{i=1} C_{i}*det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  RealType const0(0.0), const1(0.0), const2(0.0);
  for (size_t i = 0; i < nc; ++i)
  {
    const RealType c  = cptr[i];
    const size_t up   = upC[i];
    const size_t down = dnC[i];

    const0 += c * detValues_dn[down] * (detValues_up[up] / detValues_up[0]);
  }

  std::fill(pK4.begin(), pK4.end(), 0.0);

  //Now we are going to loop through all unique determinants.
  //The few lines above are for the reference matrix contribution.
  //Although I start the loop below from index 0, the loop only performs actions when the index is >= 1
  //the detData object contains all the information about the P^T and Q matrices (projection matrices) needed in the table method
  const int* restrict data_it = detData_up.data();
  for (int index = 0, datum = 0; index < num_unique_up_dets; index++)
  {
    const int k = data_it[datum];

    if (k == 0)
    {
      datum += 3 * k + 1;
    }

    else
    {
      //Number of rows and cols of P^T
      const int prows = k;
      const int pcols = nel;
      //Number of rows and cols of Q
      const int qrows = nmo;
      const int qcols = k;

      Y5.resize(nel, k);
      Y6.resize(k, k);

      //Any matrix multiplication of P^T or Q is simply a projection
      //Explicit matrix multiplication can be avoided; instead column or row copying can be done
      //BlAS::copy(size of col/row being copied,
      //           Matrix pointer + place to begin copying,
      //           storage spacing (number of elements btw next row/col element),
      //           Pointer to resultant matrix + place to begin pasting,
      //           storage spacing of resultant matrix)
      //For example the next 4 lines is the matrix multiplication of T*Q = Y5
      std::fill(Y5.begin(), Y5.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(nel, T + data_it[datum + 1 + k + i], nmo, Y5.data() + i, k);
      }

      std::fill(Y6.begin(), Y6.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y5.data() + (data_it[datum + 1 + i]) * k, 1, (Y6.data() + i * k), 1);
      }

      Vector<ValueType> WS;
      Vector<IndexType> Piv;
      WS.resize(k);
      Piv.resize(k);
      std::complex<RealType> logdet = 0.0;
      InvertWithLog(Y6.data(), k, k, WS.data(), Piv.data(), logdet);

      Y7.resize(k, nel);

      std::fill(Y7.begin(), Y7.end(), 0.0);
      for (int i = 0; i < k; i++)
      {
        BLAS::copy(k, Y6.data() + i, k, Y7.data() + (data_it[datum + 1 + i]), nel);
      }

      for (int p = 0; p < lookup_tbl[index].size(); p++)
      {
        //el_p is the element position that contains information about the CI coefficient, and det up/dn values associated with the current unique determinant
        const int el_p(lookup_tbl[index][p]);
        const RealType c  = cptr[el_p];
        const size_t up   = upC[el_p];
        const size_t down = dnC[el_p];

        const RealType alpha_4(c * detValues_dn[down] * detValues_up[up] * (1 / psiCurrent));

        for (int i = 0; i < k; i++)
        {
          BLAS::axpy(nel, alpha_4, Y7.data() + i * nel, 1, pK4.data() + (data_it[datum + 1 + k + i]) * nel, 1);
        }
      }
      datum += 3 * k + 1;
    }
  }

  BLAS::gemm('N', 'N', nmo, nmo, nel, RealType(1.0), T, nmo, pK4.data(), nel, RealType(0.0), K4T.data(), nmo);
  BLAS::gemm('N', 'N', nmo, nel, nmo, RealType(1.0), K4T.data(), nmo, T, nmo, RealType(0.0), TK4T.data(), nmo);

  for (int mu = 0, k = parameter_start_index; k < (parameter_start_index + parameters_size); k++, mu++)
  {
    int kk = myVars.where(k);
    if (kk >= 0)
    {
      const int i(m_act_rot_inds[mu].first), j(m_act_rot_inds[mu].second);
      if (i <= nel - 1 && j > nel - 1)
      {
        dlogpsi[kk] += ValueType(detValues_up[0] * (Table(i, j)) * const0 * (1 / psiCurrent) +
                                 (K4T(i, j) - K4T(j, i) - TK4T(i, j)));
      }
      else if (i <= nel - 1 && j <= nel - 1)
      {
        dlogpsi[kk] += ValueType(detValues_up[0] * (Table(i, j) - Table(j, i)) * const0 * (1 / psiCurrent) +
                                 (K4T(i, j) - TK4T(i, j) - K4T(j, i) + TK4T(j, i)));
      }
      else
      {
        dlogpsi[kk] += ValueType((K4T(i, j) - K4T(j, i)));
      }
    }
  }
}


std::unique_ptr<SPOSet> RotatedSPOs::makeClone() const
{
  auto myclone = std::make_unique<RotatedSPOs>(my_name_, std::unique_ptr<SPOSet>(Phi->makeClone()));

  myclone->params          = this->params;
  myclone->params_supplied = this->params_supplied;
  myclone->m_act_rot_inds  = this->m_act_rot_inds;
  myclone->m_full_rot_inds = this->m_full_rot_inds;
  myclone->myVars          = this->myVars;
  myclone->myVarsFull      = this->myVarsFull;
  myclone->history_params_ = this->history_params_;
  myclone->use_global_rot_ = this->use_global_rot_;
  return myclone;
}

void RotatedSPOs::mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                                       const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                       const RefVector<ValueVector>& psi_list,
                                       const std::vector<const ValueType*>& invRow_ptr_list,
                                       std::vector<std::vector<ValueType>>& ratios_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateDetRatios(phi_list, vp_list, psi_list, invRow_ptr_list, ratios_list);
}

void RotatedSPOs::mw_evaluateValue(const RefVectorWithLeader<SPOSet>& spo_list,
                                   const RefVectorWithLeader<ParticleSet>& P_list,
                                   int iat,
                                   const RefVector<ValueVector>& psi_v_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateValue(phi_list, P_list, iat, psi_v_list);
}

void RotatedSPOs::mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& spo_list,
                                 const RefVectorWithLeader<ParticleSet>& P_list,
                                 int iat,
                                 const RefVector<ValueVector>& psi_v_list,
                                 const RefVector<GradVector>& dpsi_v_list,
                                 const RefVector<ValueVector>& d2psi_v_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateVGL(phi_list, P_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list);
}

void RotatedSPOs::mw_evaluateVGLWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                         const RefVectorWithLeader<ParticleSet>& P_list,
                                         int iat,
                                         const RefVector<ValueVector>& psi_v_list,
                                         const RefVector<GradVector>& dpsi_v_list,
                                         const RefVector<ValueVector>& d2psi_v_list,
                                         OffloadMatrix<ComplexType>& mw_dspin) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateVGLWithSpin(phi_list, P_list, iat, psi_v_list, dpsi_v_list, d2psi_v_list, mw_dspin);
}

void RotatedSPOs::mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                                 const RefVectorWithLeader<ParticleSet>& P_list,
                                                 int iat,
                                                 const std::vector<const ValueType*>& invRow_ptr_list,
                                                 OffloadMWVGLArray& phi_vgl_v,
                                                 std::vector<ValueType>& ratios,
                                                 std::vector<GradType>& grads) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateVGLandDetRatioGrads(phi_list, P_list, iat, invRow_ptr_list, phi_vgl_v, ratios, grads);
}

void RotatedSPOs::mw_evaluateVGLandDetRatioGradsWithSpin(const RefVectorWithLeader<SPOSet>& spo_list,
                                                         const RefVectorWithLeader<ParticleSet>& P_list,
                                                         int iat,
                                                         const std::vector<const ValueType*>& invRow_ptr_list,
                                                         OffloadMWVGLArray& phi_vgl_v,
                                                         std::vector<ValueType>& ratios,
                                                         std::vector<GradType>& grads,
                                                         std::vector<ValueType>& spingrads) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluateVGLandDetRatioGradsWithSpin(phi_list, P_list, iat, invRow_ptr_list, phi_vgl_v, ratios, grads,
                                                spingrads);
}

void RotatedSPOs::mw_evaluate_notranspose(const RefVectorWithLeader<SPOSet>& spo_list,
                                          const RefVectorWithLeader<ParticleSet>& P_list,
                                          int first,
                                          int last,
                                          const RefVector<ValueMatrix>& logdet_list,
                                          const RefVector<GradMatrix>& dlogdet_list,
                                          const RefVector<ValueMatrix>& d2logdet_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.mw_evaluate_notranspose(phi_list, P_list, first, last, logdet_list, dlogdet_list, d2logdet_list);
}

void RotatedSPOs::createResource(ResourceCollection& collection) const { Phi->createResource(collection); }

void RotatedSPOs::acquireResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.acquireResource(collection, phi_list);
}

void RotatedSPOs::releaseResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const
{
  auto phi_list = extractPhiRefList(spo_list);
  auto& leader  = phi_list.getLeader();
  leader.releaseResource(collection, phi_list);
}

RefVectorWithLeader<SPOSet> RotatedSPOs::extractPhiRefList(const RefVectorWithLeader<SPOSet>& spo_list)
{
  auto& spo_leader = spo_list.getCastedLeader<RotatedSPOs>();
  const auto nw    = spo_list.size();
  RefVectorWithLeader<SPOSet> phi_list(*spo_leader.Phi);
  phi_list.reserve(nw);
  for (int iw = 0; iw < nw; iw++)
  {
    RotatedSPOs& rot = spo_list.getCastedElement<RotatedSPOs>(iw);
    phi_list.emplace_back(*rot.Phi);
  }
  return phi_list;
}

} // namespace qmcplusplus
