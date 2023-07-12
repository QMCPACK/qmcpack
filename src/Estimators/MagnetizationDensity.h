//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//
// File created by: Raymond Clay, rclay@sandia.gov, Sandia National Laboratories
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_MAGNETIZATION_DENSITY_H
#define QMCPLUSPLUS_MAGNETIZATION_DENSITY_H

#include <vector>
#include <functional>

#include "Estimators/OperatorEstBase.h"
#include "type_traits/complex_help.hpp"
#include "ParticleBase/RandomSeqGenerator.h"
#include <SpeciesSet.h>
#include <StdRandom.h>
#include "Estimators/MagnetizationDensityInput.h"

namespace qmcplusplus
{
namespace testing
{
class MagnetizationDensityTests;
}
/** Magnetization density estimator for non-collinear spin calculations.
 *
 * As documented in the manual, the following formula is computed:  
 *\mathbf{m}_c = \int d\mathbf{X} \left|{\Psi(\mathbf{X})}\right|^2 \int_{\Omega_c}d\mathbf{r} \sum_i\delta(\mathbf{r}-\hat{\mathbf{r}}_i)\int_0^{2\pi} \frac{ds'_i}{2\pi} \frac{\Psi(\ldots \mathbf{r}_i s'_i \ldots )}{\Psi(\ldots \mathbf{r}_i s_i \ldots)}\langle s_i | \hat{\sigma} | s'_i \rangle
 *
 * The way that the above relates to the underlying data structure data_ is as follows.  We grid space up and assign an index
 * for each of the real space bins (identical to SpinDensityNew).  To account for the fact that magnetization is vectorial, we 
 * triple the length of this array.  If grid_i is the index of the real space gridpoint i, then the data is layed out like:
 * [grid_0_x, grid_0_y, grid_0_z, grid_1_x, ..., grid_N_x, grid_N_y, grid_N_z].  This is also the way it is stored in HDF5.  
 *
 */
class MagnetizationDensity : public OperatorEstBase
{
public:
  using Value              = QMCTraits::ValueType;
  using FullPrecValue      = QMCTraits::FullPrecValueType;
  using Real               = RealAlias<Value>;
  using FullPrecReal       = RealAlias<FullPrecValue>;
  using Grad               = TinyVector<Value, OHMMS_DIM>;
  using Lattice            = PtclOnLatticeTraits::ParticleLayout;
  using Position           = QMCTraits::PosType;
  using Integrator         = MagnetizationDensityInput::Integrator;
  static constexpr int DIM = QMCTraits::DIM;

  MagnetizationDensity(MagnetizationDensityInput&& min, const Lattice& lattice);
  MagnetizationDensity(const MagnetizationDensity& magdens, DataLocality dl);

  void startBlock(int steps) override;

  void accumulate(const RefVector<MCPWalker>& walkers,
                  const RefVector<ParticleSet>& psets,
                  const RefVector<TrialWaveFunction>& wfns,
                  RandomBase<FullPrecReal>& rng) override;

  void collect(const RefVector<OperatorEstBase>& operator_estimators) override;
  /**
  * This returns the total size of data object required for this estimator.
  * Right now, it's 3*number_of_realspace_gridpoints
  *
  * @return Size of data.
  */
  size_t getFullDataSize();
  std::unique_ptr<OperatorEstBase> spawnCrowdClone() const override;
  void registerOperatorEstimator(hdf_archive& file) override;

  /**
  * This is a convenience function that handles \int_0^{2\pi} dx f(x).
  * There are two provided methods for integrating this, so this function picks Simpson's
  * or Monte Carlo to perform this integration, while keeping awareness of the integration
  * interval and number of samples.
  *
  * @tparam Value type.  Real or complex in double or single precision.  Return type is same as fgrid type.
  * @param[in] fgrid f(x), the function to integrate.  Assumed to be on a [0-2pi] interval, and if the grid isn't
  *        random, it's assumed to be uniform.
  * @return Value of integral.
  */
  Value integrateMagnetizationDensity(const std::vector<Value>& fgrid) const;

private:
  MagnetizationDensity(const MagnetizationDensity& magdens) = default;


  /**
  * Generates the spin integrand \Psi(s')/Psi(s)* \langle s | \vec{\sigma} | s'\rangle for a specific
  *  electron iat.  Since this is a vectorial quantity, this function returns sx, sy, and sz in their own
  *  arrays.  
  *
  * @param[in] pset ParticleSet  
  * @param[in] wfn TrialWaveFunction
  * @param[in] iat electron index
  * @param[out] sx x component of spin integrand
  * @param[out] sy y component of spin integrand
  * @param[out] sz z component of spin integrand
  */
  void generateSpinIntegrand(ParticleSet& pset,
                             TrialWaveFunction& wfn,
                             const int iat,
                             std::vector<Value>& sx,
                             std::vector<Value>& sy,
                             std::vector<Value>& sz);

  /**
  *  Implementation of Simpson's 1/3 rule to integrate a function on a uniform grid. 
  *
  * @param[in] fgrid f(x), the function to integrate. 
  * @param[in] gridDx, the grid spacing for the uniform grid.  Assumed to be consistent with size of fgrid.
  * @return Value of integral.
  */
  Value integrateBySimpsonsRule(const std::vector<Value>& fgrid, Real gridDx) const;


  /**
  * Convenience function to generate a grid between 0 and 2pi, consistent with nsamples_ and integration method.  
  * Can be uniform or random, depending on choice of integrator.   
  *
  * @param[out] sgrid A grid with nsamples_ points between 0 and 2pi.  Overwritten.
  */
  void generateGrid(std::vector<Real>& sgrid) const;

  /**
  * Generate a uniform grid between [start,stop] for numerical quadrature.  
  *
  * @param[out] sgrid Random number grid between "start" and "stop".  Number of points taken from size of sgrid.
  * @param[in] start start of grid interval 
  * @param[in] stop end of grid interval
  */
  void generateUniformGrid(std::vector<Real>& sgrid, const Real start, const Real stop) const;

  /**
  * Generate random grid between [start,stop] for MC integration.   
  *
  * @tparam RAN_GEN Random number generator type.
  * @param[out] sgrid Random number grid between "start" and "stop".  Number of points taken from size of sgrid.
  * @param[in] rng Random number generator queried to generate random points.   
  * @param[in] start start of grid interval 
  * @param[in] stop end of grid interval
  */

  void generateRandomGrid(std::vector<Real>& sgrid, RandomBase<FullPrecReal>& rng, Real start, Real stop) const;

  /**
  * For a given spatial position r and spin component s, this returns the bin for accumulating the observable.
  * 
  *
  * @param[in]  Position in real space. 3D
  * @param[in]  component, the x=0,y=1,z=2 component of the spin.
  * @return Index of appropriate bin for this position and spin component 
  */
  size_t computeBin(const Position& r, const unsigned int component) const;
  MagnetizationDensityInput input_;
  //These are the same variables as in SpinDensityNew.
  Integrator integrator_;
  Lattice lattice_;
  Position rcorner_;
  Position center_;
  TinyVector<int, DIM> grid_;
  TinyVector<int, DIM> gdims_;
  size_t npoints_;
  int nsamples_;

  friend class testing::MagnetizationDensityTests;
};
} // namespace qmcplusplus
#endif /* QMCPLUSPLUS_MAGNETIZATION_DENSITY_H */
