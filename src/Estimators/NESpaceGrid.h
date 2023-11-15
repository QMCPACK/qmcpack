//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: QMCHamiltonian/SpaceGrid.h
//////////////////////////////////////////////////////////////////////////////////////

/** \file
 *  This is the port of QMCHamiltonian/SpaceGrid to the new Estimator design.
 *  Name clashes were undesirable as the legacy implementation needed to remain.
 *  \todo rename to more obvious Spacegrid once QMCHamiltonian/SpaceGrid is removed
 */
#ifndef QMCPLUSPLUS_NESPACEGRID_H
#define QMCPLUSPLUS_NESPACEGRID_H

#include <Configuration.h>

#include "SpaceGridInput.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Pools/PooledData.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "Particle/DistanceTable.h"
#include "NEReferencePoints.h"

namespace qmcplusplus
{
/** SpaceGrid refactored for use with batched estimator design
 *  NE should be dropped when QMCHamiltonian/SpaceGrid has been deleted.
 *
 *  This class has more going on than just representing a spacial grid
 *  I'm still working out how much of this just because of the Voronoi code that shouldn't be
 *  part of the same class as the simpler grid and how much is particleset contimination etc.
 */
namespace testing
{
template<typename T>
class NESpaceGridTests;
}

class NESpaceGrid
{
public:
  using Real        = QMCTraits::RealType;
  using Point       = typename NEReferencePoints::Point;
  using Points      = typename NEReferencePoints::Points;
  using BufferType  = PooledData<Real>;
  using POLT        = PtclOnLatticeTraits;
  using ParticlePos = POLT::ParticlePos;
  using AxTensor    = Tensor<Real, OHMMS_DIM>;
  enum class ReferenceEnergy
  {
    vacuum,
    neutral,
    noref
  };

  /** This constructor is used for electron only grids
   * \param[in]  sgi            input object for space grid.
   * \param[in]  reference      reference points from which origin and on other reference points referenced in input object are to be found
   * \param[in]  nvalues        number of fields the owning class wants for each grid point.
   * \param[in]  is_period      properly names is what is says
   */
  NESpaceGrid(SpaceGridInput& sgi, const Points& points, const int nvalues, const bool is_periodic);

  /** This is the general constructor
   * \param[in]  sgi            input object for space grid.
   * \param[in]  reference      reference points from which origin and on other reference points referenced in input object are to be found
   * \param[in]  ndp            number of particles that can move
   * \param[in]  nvalues        number of fields the owning class wants for each grid point.
   * \param[in]  is_period      properly names is what is says
   */
  NESpaceGrid(SpaceGridInput& sgi, const Points& points, const int ndp, const int nvalues, const bool is_periodic);

  /** This is the constructor for when PStatic is used.
   */
  NESpaceGrid(SpaceGridInput& sgi,
              const Points& points,
              ParticlePos& static_particle_positions,
              std::vector<Real>& static_particle_charges,
              const int ndp,
              const int nvalues,
              const bool is_periodic);

  NESpaceGrid(const NESpaceGrid& sg) = default;
  NESpaceGrid& operator=(const NESpaceGrid& sg) = default;
  
  void write_description(std::ostream& os, const std::string& indent);

  /** set up Observable helper(s) for this grid
   *  almost unchanged from legacy
   *  \todo uses Observable helper unpleasantly in implementation, remove
   */
  void registerGrid(hdf_archive& file, int grid_index);

  void write(hdf_archive& file) const;
  /// @}

  void accumulate(const ParticlePos& R,
                  const Matrix<Real>& values,
                  std::vector<bool>& particles_outside);

  /** SpaceGridAccumulate not type erased and with its own particular interface.
   *  the composing class needs to provide the following to spave grid.
   *  \param[in]      R                    particle positions
   *  \param[in]      values               matrix indexed particle, value
   *  \param[in/out]  buf                  buffer to accumulating grid to
   *  \param[out]     particles_outside    mask vector of particles falling outside the grid box
   *  \param[in]      dtab                 particle A to Particle B distance table
   *
   *  right now cartesian grids are all accumulated as if they were "periodic" which honestly does not
   *  seem to be well defined with repsect to these grids.  But for the particle cell itself it doesn't make
   *  sense that it be periodic unless it is exactly comenserate with the particle cell (at least IMHO)
   */
  void accumulate(const ParticlePos& R,
                  const Matrix<Real>& values,
                  std::vector<bool>& particles_outside,
		  const DistanceTableAB& dtab);

  bool check_grid(void);
  int nDomains(void) const { return ndomains_; }

  void sum(const BufferType& buf, Real* vals);

  void static collect(NESpaceGrid& reduction_grid, RefVector<NESpaceGrid> grid_for_each_crowd);

  auto& getDataVector() { return data_; }
private:
  /** copy AxisGrid data to SoA layout for evaluation
   */
  void copyToSoA();


  void zero();
  
  /** return actual origin point based on input
   */
  static Point deriveOrigin(const SpaceGridInput& input, const Points& points);

  /** Initialize NESpaceGrid for rectilinear grid
   *  \param[in]  input     SpaceGridInput object
   *  \param[in]  points    ReferencePoints object for grid
   *  Causes side effects updating
   *    origin_    fixed up origin for grid
   *    axes_      axes with scaling applied to it.
   *    axinv_     the inverse of the axes with scaling applied   
   */
  bool initializeRectilinear(const SpaceGridInput& input, const Points& points);

  // /** Initialize NESpaceGrid for voronoi grid
  //  *  \param[in]  input        SpaceGridInput object
  //  *  \param[in]  points       ReferencePoints object for grid
  //  *  \param[in]  rs_static    Initial static particle positions
  //  *  Causes side effects updating
  //  *    origin_    fixed up origin for grid
  //  *    axes_      axes with scaling applied to it.
  //  *    axinv_     the inverse of the axes with scaling applied   
  //  */
  // bool initializeVoronoi(const SpaceGridInput& input, const Points& points, ParticlePos& r_static);

  /** Another function to cut scopes to sort of manageable size.
   *  does nothing but create many side effects
   */
  void someMoreAxisGridStuff();

  /** create axes and axinv tensors
   *  \param[in]  input   space grid input
   *  \param[out] axes
   *  \param[out] axinv
   */
  static void processAxis(const SpaceGridInput& input, const Points& points, AxTensor& axes, AxTensor& axinv);

  static bool checkAxisGridValues(const SpaceGridInput& input_, const AxTensor& axes);

  /** refrence points for the space grid
   *  this reference it to the EstimatorManagers EDE's spacegrid_inputs_
   */
  SpaceGridInput& input_;
  int ndparticles_;
  bool is_periodic_;
  /** refrence points for the space grid
   *  this reference is to the EstimatorManagers EDE's reference points
   */
  const Points& points_;

  // _NO_
  int buffer_start_;
  int buffer_end_;

  /** in legacy used to be the starting index into the collectibles buffer.
   *  Maintained to use more legacy code without modificaiton in the short term.
   *  In the long term its possible the entire way the grid data is structured in memory should be redesigned.
   */
  const int buffer_offset_{0}; 
  int ndomains_{1};
  int nvalues_per_domain_;
  /** @ingroup Calculated by NESpaceGrid
   *  These are calculated by NESpaceGrid before accumulation and possibly belong in SpaceGridInput
   *  as derived inputs.  i.e. they are immutable and only based on the input.
   *  Alternately they would be appropriate to calculate at construction time.
   *  @{ */
  Matrix<Real> domain_volumes_;
  Matrix<Real> domain_centers_;
  //really only used for cartesian-like grids
  Point origin_;
  AxTensor axes_;
  AxTensor axinv_;

  /// @}

  Real volume_;
  Matrix<Real> domain_uwidths_;
  std::string axlabel_[OHMMS_DIM];
  std::array<std::vector<int>, 3> gmap_;
  Real odu_[OHMMS_DIM];
  Real umin_[OHMMS_DIM];
  Real umax_[OHMMS_DIM];
  int dm_[OHMMS_DIM];
  ReferenceEnergy reference_energy_;
  std::vector<Real> data_;
  std::shared_ptr<ObservableHelper> observable_helper_;
  
  struct IRPair
  {
    Real r;
    int i;
  };
  std::vector<IRPair> nearcell_;

public:
  template<typename T>
  friend class testing::NESpaceGridTests;
};


} // namespace qmcplusplus

#endif
