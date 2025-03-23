//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H
#define QMCPLUSPLUS_SOA_LINEARCOMIBINATIONORBITALSET_TEMP_H

#include <memory>
#include "QMCWaveFunctions/SPOSet.h"
#include "QMCWaveFunctions/BasisSetBase.h"

#include "Numerics/MatrixOperators.h"
#include "Numerics/DeterminantOperators.h"
#include <mkl_spblas.h>

namespace qmcplusplus
{
// ---- CSRMatrix struct ----
template<typename T = double, typename I = MKL_INT>
struct CSRMatrix
{
  std::vector<T> values;   // nonzero values
  std::vector<I> rowIndex; // row start idx
  std::vector<I> columns;  // column indices
  int rows, cols, nnz;

  CSRMatrix(int r, int c) : rows(r), cols(c), nnz(0) {}

  void reserve(int nnz_capacity)
  {
    values.reserve(nnz_capacity);
    columns.reserve(nnz_capacity);
    rowIndex.resize(rows + 1, 0);
  }
};

// ---- Conversion functions ----
template<typename T = double, typename I = MKL_INT>
CSRMatrix<T, I> dense_to_sparse(T* A, size_t nrows, size_t ncols, double tol = 1e-12)
{
  CSRMatrix<T, I> Acsr(nrows, ncols);
  Acsr.rowIndex.resize(nrows + 1, 0);
  Acsr.rowIndex[0] = 0;

  for (size_t i = 0; i < nrows; ++i)
  {
    for (size_t j = 0; j < ncols; ++j)
    {
      T Aij = A[i * ncols + j]; // row-major
      // T Aij = A[i + j * nrows]; // col-major
      if (std::abs(Aij) > tol)
      {
        Acsr.values.push_back(Aij);
        Acsr.columns.push_back(static_cast<I>(j));
        Acsr.nnz++;
      }
    }
    Acsr.rowIndex[i + 1] = static_cast<I>(Acsr.nnz);
  }
  app_log() << "converted dense " << nrows << " x " << ncols << " = " << (nrows * ncols) << std::endl;
  app_log() << "to sparse with nnz =  " << Acsr.nnz << std::endl;

  return Acsr;
}


template<typename T = double, typename I = MKL_INT>
sparse_matrix_t make_mkl_sparse(const CSRMatrix<T, I>& A)
{
  sparse_matrix_t A_mkl;
  /// TODO: could just pass by reference, but modifications to A_mkl would change CSRMatrix
  /// (make input non-const)
  std::vector<I> rowIndex = A.rowIndex;
  std::vector<I> columns  = A.columns;
  std::vector<T> values   = A.values;

  mkl_sparse_d_create_csr(&A_mkl, SPARSE_INDEX_BASE_ZERO, A.rows, A.cols, rowIndex.data(), rowIndex.data() + 1,
                          columns.data(), values.data());

  return A_mkl;
}


template<typename T = double, typename I = MKL_INT>
sparse_matrix_t dense_to_mkl(T* A, size_t nrows, size_t ncols, double tol = 1e-12)
{
  CSRMatrix<T, I> Acsr = dense_to_sparse<T, I>(A, nrows, ncols, tol);
  return make_mkl_sparse<T, I>(Acsr);
}

class MklSparseHandle
{
public:
  using default_value_type = double;
  using default_index_type = MKL_INT;

  MklSparseHandle() : descr_({SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_FULL, SPARSE_DIAG_NON_UNIT}) {}

  MklSparseHandle(sparse_matrix_t handle, matrix_descr descr)
      : handle_(handle,
                [](sparse_matrix_t h) {
                  if (h)
                    mkl_sparse_destroy(h);
                }),
        descr_(descr)
  {}

  template<typename T = default_value_type, typename I = default_index_type>
  static MklSparseHandle from_csr(const CSRMatrix<T, I>& csr,
                                  matrix_descr descr = {SPARSE_MATRIX_TYPE_GENERAL, SPARSE_FILL_MODE_FULL,
                                                        SPARSE_DIAG_NON_UNIT})
  {
    auto owned_csr = std::make_shared<CSRMatrix<T, I>>(csr);

    sparse_matrix_t raw;
    mkl_sparse_d_create_csr(&raw, SPARSE_INDEX_BASE_ZERO, csr.rows, csr.cols,
                            const_cast<I*>(owned_csr->rowIndex.data()), const_cast<I*>(owned_csr->rowIndex.data()) + 1,
                            const_cast<I*>(owned_csr->columns.data()), const_cast<T*>(owned_csr->values.data()));

    MklSparseHandle handle(raw, descr);
    handle.csr_owned_ = owned_csr;
    return handle;
  }

  sparse_matrix_t get() const { return handle_.get(); }
  const matrix_descr& descr() const { return descr_; }
  explicit operator bool() const { return static_cast<bool>(handle_); }

  template<typename T = default_value_type, typename I = default_index_type>
  MklSparseHandle view_first_n_rows(I m) const
  {
    auto owned = std::static_pointer_cast<CSRMatrix<T, I>>(csr_owned_);
    assert(owned && "MklSparseHandle does not own CSR data");
    assert(m <= owned->rows);

    sparse_matrix_t view;
    mkl_sparse_d_create_csr(&view, SPARSE_INDEX_BASE_ZERO, m, owned->cols, const_cast<I*>(owned->rowIndex.data()),
                            const_cast<I*>(owned->rowIndex.data()) + 1, const_cast<I*>(owned->columns.data()),
                            const_cast<T*>(owned->values.data()));

    MklSparseHandle sliced(view, descr_);
    sliced.csr_owned_ = owned; // share ownership
    return sliced;
  }

  // accessor for owned CSR data
  std::shared_ptr<void> get_csr() const { return csr_owned_; }
  bool owns_csr() const { return static_cast<bool>(csr_owned_); }

private:
  std::shared_ptr<sparse_matrix> handle_;
  matrix_descr descr_;
  std::shared_ptr<void> csr_owned_; // type-erased ownership

  friend std::ostream& operator<<(std::ostream& os, const MklSparseHandle& h);
};


template<typename T = double, typename I = MKL_INT>
MklSparseHandle from_dense(T* A, size_t nrows, size_t ncols, double tol = 1e-12)
{
  return MklSparseHandle::from_csr<T, I>(dense_to_sparse<T, I>(A, nrows, ncols, tol));
}

/** class to handle linear combinations of basis orbitals used to evaluate the Dirac determinants.
   *
   * SoA verson of LCOrtbitalSet
   * Localized basis set is always real 
   */
struct LCAOrbitalSet : public SPOSet
{
public:
  using basis_type         = SoaBasisSetBase<ValueType>;
  using vgl_type           = basis_type::vgl_type;
  using vgh_type           = basis_type::vgh_type;
  using vghgh_type         = basis_type::vghgh_type;
  using OffloadValueMatrix = OffloadMatrix<ValueType>;

  ///pointer to the basis set
  std::unique_ptr<basis_type> myBasisSet;
  /// pointer to matrix containing the coefficients
  std::shared_ptr<OffloadValueMatrix> C;
  MklSparseHandle C_csr;
  bool use_sparse_coefs;

  /** constructor
     * @param my_name name of the SPOSet object
     * @param bs pointer to the BasisSet
     * @param norb number of orbitals
     * @param identity if true, the MO coefficients matrix is identity
     */
  LCAOrbitalSet(const std::string& my_name,
                std::unique_ptr<basis_type>&& bs,
                size_t norbs,
                bool identity,
                bool use_offload);

  LCAOrbitalSet(const LCAOrbitalSet& in);

  bool isOMPoffload() const override { return useOMPoffload_; }

  std::string getClassName() const final { return "LCAOrbitalSet"; }

  bool isRotationSupported() const final { return true; }

  bool hasIonDerivs() const final { return true; }

  std::unique_ptr<SPOSet> makeClone() const final;

  void storeParamsBeforeRotation() final { C_copy = std::make_shared<OffloadValueMatrix>(*C); }

  void applyRotation(const ValueMatrix& rot_mat, bool use_stored_copy) final;

  /** set the OrbitalSetSize and Identity=false and initialize internal storages
    */
  void setOrbitalSetSize(int norbs) final;

  /** return the size of the basis set
    */
  int getBasisSetSize() const { return (myBasisSet == nullptr) ? 0 : myBasisSet->getBasisSetSize(); }

  bool isIdentity() const { return Identity; };

  /** check consistency between Identity and C
    *
    */
  void checkObject() const final;

  /** update C on device
   */
  void finalizeConstruction() override;

  void evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) final;

  void evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) final;

  void mw_evaluateValue(const RefVectorWithLeader<SPOSet>& spo_list,
                        const RefVectorWithLeader<ParticleSet>& P_list,
                        int iat,
                        const RefVector<ValueVector>& psi_v_list) const final;

  void mw_evaluateVGL(const RefVectorWithLeader<SPOSet>& spo_list,
                      const RefVectorWithLeader<ParticleSet>& P_list,
                      int iat,
                      const RefVector<ValueVector>& psi_v_list,
                      const RefVector<GradVector>& dpsi_v_list,
                      const RefVector<ValueVector>& d2psi_v_list) const final;

  void mw_evaluateDetRatios(const RefVectorWithLeader<SPOSet>& spo_list,
                            const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                            const RefVector<ValueVector>& psi_list,
                            const std::vector<const ValueType*>& invRow_ptr_list,
                            std::vector<std::vector<ValueType>>& ratios_list) const final;

  void evaluateDetRatios(const VirtualParticleSet& VP,
                         ValueVector& psi,
                         const ValueVector& psiinv,
                         std::vector<ValueType>& ratios) final;

  void mw_evaluateVGLandDetRatioGrads(const RefVectorWithLeader<SPOSet>& spo_list,
                                      const RefVectorWithLeader<ParticleSet>& P_list,
                                      int iat,
                                      const std::vector<const ValueType*>& invRow_ptr_list,
                                      OffloadMWVGLArray& phi_vgl_v,
                                      std::vector<ValueType>& ratios,
                                      std::vector<GradType>& grads) const final;

  void evaluateVGH(const ParticleSet& P, int iat, ValueVector& psi, GradVector& dpsi, HessVector& grad_grad_psi) final;

  void evaluateVGHGH(const ParticleSet& P,
                     int iat,
                     ValueVector& psi,
                     GradVector& dpsi,
                     HessVector& grad_grad_psi,
                     GGGVector& grad_grad_grad_psi) final;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            ValueMatrix& d2logdet) final;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet) final;

  void evaluate_notranspose(const ParticleSet& P,
                            int first,
                            int last,
                            ValueMatrix& logdet,
                            GradMatrix& dlogdet,
                            HessMatrix& grad_grad_logdet,
                            GGGMatrix& grad_grad_grad_logdet) final;

  //NOTE:  The data types get complicated here, so here's an overview of the
  //       data types associated with ionic derivatives, and how to get their data.
  //
  //NOTE:  These data structures hold the data for one particular ion, and so the ID is implicit.
  //       It's up to the user to keep track of which ion these derivatives refer to.
  //
  // 1.) GradMatrix grad_phi:  Holds the ionic derivatives of each SPO for each electron.
  //            Example:  grad_phi[iel][iorb][idim].  iel  -- electron index.
  //                                                iorb -- orbital index.
  //                                                idim  -- cartesian index of ionic derivative.
  //                                                        X=0, Y=1, Z=2.
  //
  // 2.) HessMatrix grad_grad_phi:  Holds the ionic derivatives of the electron gradient components
  //                                   for each SPO and each electron.
  //            Example:  grad_grad_phi[iel][iorb](idim,edim)  iel  -- electron index.
  //                                                           iorb -- orbital index.
  //                                                           idim -- ionic derivative's cartesian index.
  //                                                              X=0, Y=1, Z=2
  //                                                           edim -- electron derivative's cartesian index.
  //                                                              x=0, y=1, z=2.
  //
  // 3.) GradMatrix grad_lapl_phi:  Holds the ionic derivatives of the electron laplacian for each SPO and each electron.
  //            Example:  grad_lapl_phi[iel][iorb][idim].  iel  -- electron index.
  //                                                       iorb -- orbital index.
  //                                                       idim -- cartesian index of ionic derivative.
  //                                                           X=0, Y=1, Z=2.

  /**
 * \brief Calculate ion derivatives of SPO's.
 *  
 *  @param P Electron particle set.
 *  @param first index of first electron 
 *  @@param last index of last electron
 *  @param source Ion particle set.
 *  @param iat_src  Index of ion.
 *  @param gradphi Container storing ion gradients for all particles and all orbitals.
 */
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi) final;

  /**
 * \brief Calculate ion derivatives of SPO's, their gradients, and their laplacians.
 *  
 *  @param P Electron particle set.
 *  @param first index of first electron.
 *  @@param last index of last electron
 *  @param source Ion particle set.
 *  @param iat_src  Index of ion.
 *  @param grad_phi Container storing ion gradients for all particles and all orbitals.
 *  @param grad_grad_phi Container storing ion gradients of electron gradients for all particles and all orbitals.
 *  @param grad_lapl_phi Container storing ion gradients of SPO laplacians for all particles and all orbitals.
 */
  void evaluateGradSource(const ParticleSet& P,
                          int first,
                          int last,
                          const ParticleSet& source,
                          int iat_src,
                          GradMatrix& grad_phi,
                          HessMatrix& grad_grad_phi,
                          GradMatrix& grad_lapl_phi) final;

  void evaluateGradSourceRow(const ParticleSet& P,
                             int iel,
                             const ParticleSet& source,
                             int iat_src,
                             GradVector& grad_phi) final;

  void createResource(ResourceCollection& collection) const final;
  void acquireResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const final;
  void releaseResource(ResourceCollection& collection, const RefVectorWithLeader<SPOSet>& spo_list) const final;

protected:
  ///number of Single-particle orbitals
  const IndexType BasisSetSize;
  /// a copy of the original C before orbital rotation is applied;
  std::shared_ptr<OffloadValueMatrix> C_copy;
  ///true if C is an identity matrix
  const bool Identity;
  /// whether offload is on or off at runtime.
  const bool useOMPoffload_;

  ///Temp(BasisSetSize) : Row index=V,Gx,Gy,Gz,L
  vgl_type Temp;
  ///Tempv(OrbitalSetSize) Tempv=C*Temp
  vgl_type Tempv;

  ///These are temporary VectorSoAContainers to hold value, gradient, and hessian for
  ///all basis or SPO functions evaluated at a given point.
  ///Nbasis x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temph;
  ///Norbitals x [1(value)+3(gradient)+6(hessian)]
  vgh_type Temphv;

  ///These are temporary VectorSoAContainers to hold value, gradient, hessian, and
  /// gradient hessian for all basis or SPO functions evaluated at a given point.
  ///Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempgh;
  ///Nbasis x [1(value)+3(gradient)+6(hessian)+10(grad_hessian)]
  vghgh_type Tempghv;

private:
  ///helper functions to handle Identity
  void evaluate_vgl_impl(const vgl_type& temp, ValueVector& psi, GradVector& dpsi, ValueVector& d2psi) const;

  void evaluate_vgl_impl(const vgl_type& temp,
                         int i,
                         ValueMatrix& logdet,
                         GradMatrix& dlogdet,
                         ValueMatrix& d2logdet) const;
  ///These two functions unpack the data in vgh_type temp object into wavefunction friendly data structures.


  ///This unpacks temp into vectors psi, dpsi, and d2psi.
  void evaluate_vgh_impl(const vgh_type& temp, ValueVector& psi, GradVector& dpsi, HessVector& d2psi) const;

  ///Unpacks temp into the ith row (or electron index) of logdet, dlogdet, dhlogdet.
  void evaluate_vgh_impl(const vgh_type& temp,
                         int i,
                         ValueMatrix& logdet,
                         GradMatrix& dlogdet,
                         HessMatrix& dhlogdet) const;
  ///Unpacks data in vghgh_type temp object into wavefunction friendly data structures for value, gradient, hessian
  ///and gradient hessian.
  void evaluate_vghgh_impl(const vghgh_type& temp,
                           ValueVector& psi,
                           GradVector& dpsi,
                           HessVector& d2psi,
                           GGGVector& dghpsi) const;

  void evaluate_vghgh_impl(const vghgh_type& temp,
                           int i,
                           ValueMatrix& logdet,
                           GradMatrix& dlogdet,
                           HessMatrix& dhlogdet,
                           GGGMatrix& dghlogdet) const;


  ///Unpacks data in vgl object and calculates/places ionic gradient result into dlogdet.
  void evaluate_ionderiv_v_impl(const vgl_type& temp, int i, GradMatrix& dlogdet) const;

  ///Unpacks data in vgl object and calculates/places ionic gradient of value,
  ///  electron gradient, and electron laplacian result into dlogdet, dglogdet, and dllogdet respectively.
  void evaluate_ionderiv_vgl_impl(const vghgh_type& temp,
                                  int i,
                                  GradMatrix& dlogdet,
                                  HessMatrix& dglogdet,
                                  GradMatrix& dllogdet) const;

  ///Unpacks data in vgl object and calculates/places ionic gradient of a single row (phi_j(r)) into dlogdet.
  void evaluate_ionderiv_v_row_impl(const vgl_type& temp, GradVector& dlogdet) const;

  void mw_evaluateVGLImplGEMM(const RefVectorWithLeader<SPOSet>& spo_list,
                              const RefVectorWithLeader<ParticleSet>& P_list,
                              int iat,
                              OffloadMWVGLArray& phi_vgl_v) const;

  /// packed walker GEMM implementation
  void mw_evaluateValueImplGEMM(const RefVectorWithLeader<SPOSet>& spo_list,
                                const RefVectorWithLeader<ParticleSet>& P_list,
                                int iat,
                                OffloadMWVArray& phi_v) const;

  /// packed walker GEMM implementation with multi virtual particle sets
  void mw_evaluateValueVPsImplGEMM(const RefVectorWithLeader<SPOSet>& spo_list,
                                   const RefVectorWithLeader<const VirtualParticleSet>& vp_list,
                                   OffloadMWVArray& phi_v) const;

  /// helper function for extracting a list of basis sets from a list of LCAOrbitalSet
  RefVectorWithLeader<basis_type> extractBasisRefList(const RefVectorWithLeader<SPOSet>& spo_list) const;

  struct LCAOMultiWalkerMem;
  ResourceHandle<LCAOMultiWalkerMem> mw_mem_handle_;
  /// timer for basis set
  NewTimer& basis_timer_;
  /// timer for MO
  NewTimer& mo_timer_;
};

inline std::ostream& operator<<(std::ostream& os, const CSRMatrix<double, MKL_INT>& csr)
{
  os << "\nCSRMatrix:\n";
  size_t idx = 0;
  for (size_t i = 0; i < csr.rows; i++)
  {
    idx = csr.rowIndex[i];
    for (size_t j = 0; j < csr.cols; j++)
    {
      if (csr.columns[idx] == j && idx < csr.rowIndex[i + 1])
      {
        os << csr.values[idx++] << " ";
      }
      else
      {
        os << double(0.0) << " ";
      }
    }
    os << "\n";
  }
  os << std::endl;
  return os;
}

inline std::ostream& operator<<(std::ostream& os, const MklSparseHandle& h)
{
  os << "MklSparseHandle:\n";
  /// FIXME: print more descriptive enum names instead of just the associated ints
  os << "  Description: [type = " << h.descr().type << ", mode = " << h.descr().mode << ", diag = " << h.descr().diag
     << "]\n";

  if (h.owns_csr())
  {
    // assuming that csr_owned_ was created via from_csr and holds a CSRMatrix<double, MKL_INT>
    auto csr = std::static_pointer_cast<CSRMatrix<double, MKL_INT>>(h.csr_owned_);
    os << "Printing Matrix (CSR):\n" << std::endl;
    os << *csr << std::endl; // print full CSR mat (including zeros, one matrix row per row)
    os << "  CSRMatrix details:\n";
    os << "    rows: " << csr->rows << ", cols: " << csr->cols << ", nnz: " << csr->nnz << "\n";
    // print internal CSR data (rowidx, columns, values), one array per row
    os << "    rowIndex: ";
    for (auto v : csr->rowIndex)
      os << v << " ";
    os << "\n    columns:  ";
    for (auto v : csr->columns)
      os << v << " ";
    os << "\n    values:   ";
    for (auto v : csr->values)
      os << v << " ";
    os << "\n";
  }
  else
  {
    os << "  No underlying CSRMatrix retained.\n";
  }
  return os;
}

} // namespace qmcplusplus
#endif
