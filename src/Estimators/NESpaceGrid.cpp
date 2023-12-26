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
// Some code refactored from: QMCHamiltonian/SpaceGrid.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include "NESpaceGrid.h"

#include <hdf5.h>

#include "OhmmsData/AttributeSet.h"
#include "Utilities/string_utils.h"
#include <cmath>
#include "OhmmsPETE/OhmmsArray.h"
#include "NEReferencePoints.h"
#include "Concurrency/OpenMP.h"

namespace qmcplusplus
{
using std::acos;
using std::atan2;
using std::cos;
using std::floor;
using std::sin;
using std::sqrt;

template<typename REAL>
NESpaceGrid<REAL>::NESpaceGrid(SpaceGridInput& sgi,
                               const NEReferencePoints::Points& points,
                               const int nvalues,
                               const bool is_periodic)
    : NESpaceGrid(sgi, points, 0, nvalues, is_periodic)
{}

template<typename REAL>
NESpaceGrid<REAL>::NESpaceGrid(SpaceGridInput& sgi,
                               const NEReferencePoints::Points& points,
                               const int ndp,
                               const int nvalues,
                               const bool is_periodic)
    : input_(sgi), ndparticles_(ndp), is_periodic_(is_periodic), points_(points), nvalues_per_domain_(nvalues)
{
  bool init_success = initializeCoordSystem();
  if (!init_success)
    throw std::runtime_error("NESpaceGrid initialization failed");
}

template<typename REAL>
NESpaceGrid<REAL>::NESpaceGrid(SpaceGridInput& sgi,
                               const NEReferencePoints::Points& points,
                               ParticlePos& static_particle_positions,
                               std::vector<Real>& Z,
                               const int ndp,
                               const int nvalues,
                               const bool is_periodic)
    : input_(sgi), ndparticles_(ndp), is_periodic_(is_periodic), points_(points), nvalues_per_domain_(nvalues)
{
  bool init_success = initializeCoordSystem();
  if (!init_success)
    throw std::runtime_error("NESpaceGrid initialization failed");
}

template<typename REAL>
bool NESpaceGrid<REAL>::initializeCoordSystem()
{
  using CoordForm   = SpaceGridInput::CoordForm;
  bool init_success = false;
  switch (input_.get_coord_form())
  {
  case (CoordForm::CARTESIAN):
    init_success = initializeRectilinear(input_, points_);
    break;
  case (CoordForm::CYLINDRICAL):
    init_success = initializeCylindrical(input_, points_);
    break;
  case (CoordForm::SPHERICAL):
    init_success = initializeSpherical(input_, points_);
    break;
  }
  return init_success;
}

template<typename REAL>
void NESpaceGrid<REAL>::processAxis(const SpaceGridInput& input, const Points& points, AxTensor& axes, AxTensor& axinv)
{
  auto& axis_labels = input.get_axis_labels();
  auto& axis_p1s    = input.get_axis_p1s();
  auto& axis_p2s    = input.get_axis_p2s();
  auto& axis_scales = input.get_axis_scales();
  auto& axis_grids  = input.get_axis_grids();

  // SpaceGridInput check particular validity insures the number of axis == OHMMS_DIM
  for (int iaxis = 0; iaxis < OHMMS_DIM; ++iaxis)
  {
    Real frac = axis_scales[iaxis];
    for (int d = 0; d < OHMMS_DIM; d++)
      axes(d, iaxis) = frac * (points.at(axis_p1s[iaxis])[d] - points.at(axis_p2s[iaxis])[d]);
  }
  axinv = inverse(axes);
}

template<typename REAL>
typename NESpaceGrid<REAL>::Point NESpaceGrid<REAL>::deriveOrigin(const SpaceGridInput& input, const Points& points)
{
  const std::string& origin_p1 = input.get_origin_p1();
  const std::string& origin_p2 = input.get_origin_p2();

  if (origin_p1.size() > 0 && origin_p2.size() > 0)
  {
    return points.at(origin_p1) + input.get_origin_fraction() * (points.at(origin_p2) - points.at(origin_p1));
  }
  else if (origin_p1.size() > 0)
    return points.at(origin_p1);
  else
    return points.at("zero");
}

template<typename REAL>
bool NESpaceGrid<REAL>::initializeRectilinear(const SpaceGridInput& input, const Points& points)
{
  // This code should be refactored to SpaceGridInput such that a simple map of
  // axis is available.

  origin_ = deriveOrigin(input, points);
  processAxis(input, points, axes_, axinv_);

  // this should all be done in SpaceGrid input parsing/construction now.
  // bool succeeded = checkAxisGridValues(input, axes_);

  someMoreAxisGridStuff();

  copyToSoA();

  return true;
}

template<typename REAL>
bool NESpaceGrid<REAL>::initializeCylindrical(const SpaceGridInput& input, const Points& points)
{
  // This code should be refactored to SpaceGridInput such that a simple map of
  // axis is available.

  origin_ = deriveOrigin(input, points);
  processAxis(input, points, axes_, axinv_);
  // bool succeeded = checkAxisGridValues(input, axes_);

  someMoreAxisGridStuff();

  copyToSoA();

  return true;
}

template<typename REAL>
bool NESpaceGrid<REAL>::initializeSpherical(const SpaceGridInput& input, const Points& points)
{
  // This code should be refactored to SpaceGridInput such that a simple map of
  // axis is available.

  origin_ = deriveOrigin(input, points);
  processAxis(input, points, axes_, axinv_);
  //bool succeeded = checkAxisGridValues(input, axes_);

  someMoreAxisGridStuff();

  copyToSoA();

  return true;
}

template<typename REAL>
void NESpaceGrid<REAL>::someMoreAxisGridStuff()
{
  auto& axis_grids = input_.get_axis_grids();
  // This dates back to the legacy implementation and I'm not sure why both code blocks are here.
  //set grid dimensions
  // C/Python style indexing
  dm_[0] = axis_grids[1].dimensions * axis_grids[2].dimensions;
  dm_[1] = axis_grids[2].dimensions;
  dm_[2] = 1;
  // Fortran/Matlab style indexing
  //dm[0] = 1;
  //dm[1] = dimensions[0];
  //dm[2] = dimensions[0]*dimensions[1];

  ndomains_ = axis_grids[0].dimensions * axis_grids[1].dimensions * axis_grids[2].dimensions;

  data_.resize(ndomains_ * nvalues_per_domain_);

  volume_ = std::abs(det(axes_)) * 8.0; //axes span only one octant
  //compute domain volumes, centers, and widths

  domain_volumes_.resize(ndomains_, 1);
  domain_centers_.resize(ndomains_, OHMMS_DIM);
  domain_uwidths_.resize(ndomains_, OHMMS_DIM);
  std::vector<Real> interval_centers[OHMMS_DIM];
  std::vector<Real> interval_widths[OHMMS_DIM];

  auto& agr = axis_grids;

  for (int d = 0; d < OHMMS_DIM; d++)
  {
    int nintervals = agr[d].ndu_per_interval.size();
    app_log() << "nintervals " << d << " " << nintervals << std::endl;
    interval_centers[d].resize(nintervals);
    interval_widths[d].resize(nintervals);
    interval_widths[d][0]  = agr[d].ndu_per_interval[0] / agr[d].odu;
    interval_centers[d][0] = interval_widths[d][0] / 2.0 + agr[d].umin;
    for (int i = 1; i < nintervals; i++)
    {
      interval_widths[d][i]  = agr[d].ndu_per_interval[i] / agr[d].odu;
      interval_centers[d][i] = interval_centers[d][i - 1] + .5 * (interval_widths[d][i] + interval_widths[d][i - 1]);
    }
    //app_log()<<"  interval widths"<< std::endl;
    //for(int i=0;i<nintervals;i++)
    //  app_log()<<"    "<<i<<" "<<interval_widths[d][i]<< std::endl;
    //app_log()<<"  interval centers"<< std::endl;
    //for(int i=0;i<nintervals;i++)
    //  app_log()<<"    "<<i<<" "<<interval_centers[d][i]<< std::endl;
  }
  Point du, uc, ubc, rc;
  Real vol        = 0.0;
  Real vscale     = std::abs(det(axes_));
  using CoordForm = SpaceGridInput::CoordForm;
  for (int i = 0; i < agr[0].dimensions; i++)
  {
    for (int j = 0; j < agr[1].dimensions; j++)
    {
      for (int k = 0; k < agr[2].dimensions; k++)
      {
        int idomain = dm_[0] * i + dm_[1] * j + dm_[2] * k;
        du[0]       = interval_widths[0][i];
        du[1]       = interval_widths[1][j];
        du[2]       = interval_widths[2][k];
        uc[0]       = interval_centers[0][i];
        uc[1]       = interval_centers[1][j];
        uc[2]       = interval_centers[2][k];
        switch (input_.get_coord_form())
        {
        case (CoordForm::CARTESIAN):
          vol = du[0] * du[1] * du[2];
          ubc = uc;
          break;
        case (CoordForm::CYLINDRICAL):
          uc[1]  = 2.0 * M_PI * uc[1] - M_PI;
          du[1]  = 2.0 * M_PI * du[1];
          vol    = uc[0] * du[0] * du[1] * du[2];
          ubc[0] = uc[0] * cos(uc[1]);
          ubc[1] = uc[0] * sin(uc[1]);
          ubc[2] = uc[2];
          break;
        case (CoordForm::SPHERICAL):
          uc[1] = 2.0 * M_PI * uc[1] - M_PI;
          du[1] = 2.0 * M_PI * du[1];
          uc[2] = M_PI * uc[2];
          du[2] = M_PI * du[2];
          vol   = (uc[0] * uc[0] + du[0] * du[0] / 12.0) * du[0] //r
              * du[1]                                            //theta
              * 2.0 * sin(uc[2]) * sin(.5 * du[2]);              //phi
          ubc[0] = uc[0] * sin(uc[2]) * cos(uc[1]);
          ubc[1] = uc[0] * sin(uc[2]) * sin(uc[1]);
          ubc[2] = uc[0] * cos(uc[2]);
          break;
        default:
          break;
        }
        vol *= vscale;
        rc = dot(axes_, ubc) + origin_;
        //app_log()<< std::endl;
        //app_log()<<"umin "<<uc-du/2<< std::endl;
        //app_log()<<"uc "<<uc<< std::endl;
        //app_log()<<"umax "<<uc+du/2<< std::endl;
        //app_log()<<"rc "<<rc<< std::endl;
        domain_volumes_(idomain, 0) = vol;
        for (int d = 0; d < OHMMS_DIM; d++)
        {
          domain_uwidths_(idomain, d) = du[d];
          domain_centers_(idomain, d) = rc[d];
        }
      }
    }
  }

  //find the actual volume of the grid
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    du[d] = agr[d].umax - agr[d].umin;
    uc[d] = .5 * (agr[d].umax + agr[d].umin);
  }
  switch (input_.get_coord_form())
  {
  case CoordForm::CARTESIAN:
    vol = du[0] * du[1] * du[2];
    break;
  case CoordForm::CYLINDRICAL:
    uc[1] = 2.0 * M_PI * uc[1] - M_PI;
    du[1] = 2.0 * M_PI * du[1];
    vol   = uc[0] * du[0] * du[1] * du[2];
    break;
  case CoordForm::SPHERICAL:
    uc[1] = 2.0 * M_PI * uc[1] - M_PI;
    du[1] = 2.0 * M_PI * du[1];
    uc[2] = M_PI * uc[2];
    du[2] = M_PI * du[2];
    vol   = (uc[0] * uc[0] + du[0] * du[0] / 12.0) * du[0] //r
        * du[1]                                            //theta
        * 2.0 * sin(uc[2]) * sin(.5 * du[2]);              //phi
    break;
  default:
    break;
  }
  volume_ = vol * det(axes_);

  return;
}

// template<typename REAL>
// bool NESpaceGrid<REAL>::checkAxisGridValues(const SpaceGridInput& input, const AxTensor& axes)
// {
//   auto& axis_labels = input.get_axis_labels();
//   auto& axis_grids  = input.get_axis_grids();

//   //check that all axis grid values fall in the allowed intervals for the coord label
//   bool succeeded = true;
//   for (int d = 0; d < OHMMS_DIM; d++)
//   {
    
//     if (axis_labels[d] == "phi" || axis_labels[d] == "theta" )
//       if (axis_grids[d].umin < 0.0 || axis_grids[d].umax > 1.0)
//       {
//         app_log() << "  grid values for " << axis_labels[d] << " must fall in [0,1]" << std::endl;
//         app_log() << "  interval provided: [" << axis_grids[d].umin << "," << axis_grids[d].umax << "]" << std::endl;
//         succeeded = false;
//       }
//    else
//         if (axis_grids[d].umin < -1.0 || axis_grids[d].umax > 1.0)
//         {
//           app_log() << "  grid values for " << axis_labels[d] << " must fall in [-1,1]" << std::endl;
//           app_log() << "  interval provided: [" << axis_grids[d].umin << "," << axis_grids[d].umax << "]" << std::endl;
//           succeeded = false;
//         }
//   return succeeded;
// }

template<typename REAL>
void NESpaceGrid<REAL>::write_description(std::ostream& os, const std::string& indent)
{
  os << indent + "SpaceGrid" << std::endl;
  std::string s;
  using CoordForm = SpaceGridInput::CoordForm;
  switch (input_.get_coord_form())
  {
  case CoordForm::CARTESIAN:
    s = "cartesian";
    break;
  case CoordForm::CYLINDRICAL:
    s = "cylindrical";
    break;
  case CoordForm::SPHERICAL:
    s = "spherical";
    break;
  default:
    break;
  }
  auto& agr         = input_.get_axis_grids();
  auto& axis_labels = input_.get_axis_labels();
  os << indent + "  SpaceGridInput::lookup_input_ename_value(input_.get_coord_form()  = " + s << std::endl;
  os << indent + "  ndomains_      = " << ndomains_ << std::endl;
  os << indent + "  axes  = " << axes_ << std::endl;
  os << indent + "  axinv = " << axinv_ << std::endl;
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    os << indent + "  axis " << axis_labels[d] << ":" << std::endl;
    os << indent + "    umin = " << agr[d].umin << std::endl;
    os << indent + "    umax = " << agr[d].umax << std::endl;
    os << indent + "    du   = " << 1.0 / static_cast<double>(agr[d].odu) << std::endl;
    os << indent + "    dm   = " << dm_[d] << std::endl;
    os << indent + "    gmap = ";
    for (int g = 0; g < gmap_[d].size(); g++)
    {
      os << gmap_[d][g] << " ";
    }
    os << std::endl;
  }
  os << indent + "end NESpaceGrid" << std::endl;
}

template<typename REAL>
void NESpaceGrid<REAL>::registerGrid(hdf_archive& file, int grid_index)
{
  using iMatrix = Matrix<int>;
  iMatrix imat;
  std::vector<int> ng(1);
  int cshift = 1;
  std::stringstream ss;
  ss << grid_index + cshift;
  hdf_path hdf_name{"spacegrid" + ss.str()};
  observable_helper_ = std::make_shared<ObservableHelper>(hdf_name);
  auto& oh           = *observable_helper_;
  ng[0]              = nvalues_per_domain_ * ndomains_;
  oh.set_dimensions(ng, buffer_offset_);

  // Create a bunch of temporary SoA data from input to write the grid attributes
  auto& agr = input_.get_axis_grids();
  std::vector<int> dimensions(OHMMS_DIM);
  for (int id = 0; id < OHMMS_DIM; ++id)
  {
    dimensions[id] = agr[id].dimensions;
  }

  int coord = static_cast<int>(input_.get_coord_form());
  oh.addProperty(const_cast<int&>(coord), "coordinate", file);
  oh.addProperty(const_cast<int&>(ndomains_), "ndomains", file);
  oh.addProperty(const_cast<int&>(nvalues_per_domain_), "nvalues_per_domain_", file);
  oh.addProperty(const_cast<Real&>(volume_), "volume", file);
  oh.addProperty(const_cast<Matrix<Real>&>(domain_volumes_), "domain_volumes", file);
  oh.addProperty(const_cast<Matrix<Real>&>(domain_centers_), "domain_centers", file);
  oh.addProperty(const_cast<Point&>(origin_), "origin", file);
  oh.addProperty(const_cast<Tensor<Real, OHMMS_DIM>&>(axes_), "axes", file);
  oh.addProperty(const_cast<Tensor<Real, OHMMS_DIM>&>(axinv_), "axinv", file);
  oh.addProperty(const_cast<Matrix<Real>&>(domain_uwidths_), "domain_uwidths", file);
  //add dimensioned quantities
  std::map<std::string, int> axtmap;
  axtmap["x"]     = 0;
  axtmap["y"]     = 1;
  axtmap["z"]     = 2;
  axtmap["r"]     = 0;
  axtmap["phi"]   = 1;
  axtmap["theta"] = 2;
  int axtypes[OHMMS_DIM];
  auto& axis_labels = input_.get_axis_labels();
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    axtypes[d] = axtmap[axis_labels[d]];
  }
  int n;
  const int ni = 3;
  // Now we arrive to a bunch of pointers to pointers that assume that all the axis grid
  // information is laid out soa versus aos
  int* ivar[ni];
  std::string iname[ni];
  n        = 0;
  ivar[n]  = (int*)axtypes;
  iname[n] = "axtypes";
  n++;
  ivar[n]  = dimensions.data();
  iname[n] = "dimensions";
  n++;
  ivar[n]  = (int*)dm_;
  iname[n] = "dm";
  n++;
  const int nr = 3;
  Real* rvar[nr];
  std::string rname[nr];
  n        = 0;
  rvar[n]  = (Real*)odu_;
  rname[n] = "odu";
  n++;
  rvar[n]  = (Real*)umin_;
  rname[n] = "umin";
  n++;
  rvar[n]  = (Real*)umax_;
  rname[n] = "umax";
  n++;
  imat.resize(OHMMS_DIM, 1);
  for (int i = 0; i < ni; i++)
  {
    for (int d = 0; d < OHMMS_DIM; d++)
      imat(d, 0) = ivar[i][d];
    oh.addProperty(const_cast<iMatrix&>(imat), iname[i], file);
  }
  Matrix<Real> rmat;
  rmat.resize(OHMMS_DIM, 1);
  for (int i = 0; i < nr; i++)
  {
    for (int d = 0; d < OHMMS_DIM; d++)
      rmat(d, 0) = rvar[i][d];
    oh.addProperty(const_cast<Matrix<Real>&>(rmat), rname[i], file);
  }
  for (int d = 0; d < OHMMS_DIM; d++)
  {
    int gsize = gmap_[d].size();
    imat.resize(gsize, 1);
    for (int i = 0; i < gsize; i++)
    {
      int gval   = gmap_[d][i];
      imat(i, 0) = gval;
    }
    int ival           = d + 1;
    std::string gmname = "gmap" + int2string(ival);
    oh.addProperty(const_cast<iMatrix&>(imat), gmname, file);
  }
  return;
}

template<typename REAL>
void NESpaceGrid<REAL>::write(hdf_archive& file) const
{
  if (observable_helper_)
  {
    if constexpr (std::is_same<Real, QMCTraits::FullPrecRealType>::value)
    {
      observable_helper_->write(data_.data(), file);
    }
    else
    {
      std::vector<QMCTraits::FullPrecRealType> expanded_data(data_.size(), 0.0);
      std::copy_n(data_.begin(), data_.size(), expanded_data.begin());
      assert(!data_.empty());
      // auto total = std::accumulate(data_->begin(), data_->end(), 0.0);
      // std::cout << "data size: " << data_->size() << " : " << total << '\n';
      observable_helper_->write(expanded_data.data(), file);
    }
    file.pop();
  }
  return;
}

#define NESpaceGrid_CHECK

template<typename REAL>
void NESpaceGrid<REAL>::copyToSoA()
{
  auto& agr = input_.get_axis_grids();
  for (int id = 0; id < OHMMS_DIM; ++id)
  {
    odu_[id]  = agr[id].odu;
    umin_[id] = agr[id].umin;
    umax_[id] = agr[id].umax;
    gmap_[id].resize(floor((umax_[id] - umin_[id]) * odu_[id] + .5));
    gmap_[id] = agr[id].gmap;
  }
}

template<typename REAL>
void NESpaceGrid<REAL>::accumulate(const ParticlePos& R,
                                   const Matrix<Real>& values,
                                   std::vector<bool>& particles_outside,
                                   const DistanceTableAB& dtab)
{
  // //find cell center nearest to each dynamic particle
  // int nd, nn;
  // RealType dist;
  // for (p = 0; p < ndparticles; p++)
  // {
  //   const auto& dist = dtab.getDistRow(p);
  //   for (nd = 0; nd < ndomains; nd++)
  //     if (dist[nd] < nearcell[p].r)
  //     {
  //       nearcell[p].r = dist[nd];
  //       nearcell[p].i = nd;
  //     }
  // }
  // //accumulate values for each dynamic particle
  // for (p = 0; p < ndparticles; p++)
  // {
  //   buf_index = buffer_offset + nvalues * nearcell[p].i;
  //   for (v = 0; v < nvalues; v++, buf_index++)
  //     buf[buf_index] += values(p, v);
  // }
  // //accumulate values for static particles (static particles == cell centers)
  // buf_index = buffer_offset;
  // for (p = ndparticles; p < nparticles; p++)
  //   for (v = 0; v < nvalues; v++, buf_index++)
  //     buf[buf_index] += values(p, v);
  // //each particle belongs to some voronoi cell
  // for (p = 0; p < nparticles; p++)
  //   particles_outside[p] = false;
  // //reset distances
  // for (p = 0; p < ndparticles; p++)
  //   nearcell[p].r = std::numeric_limits<RealType>::max();
  accumulate(R, values, particles_outside);
}

template<typename REAL>
void NESpaceGrid<REAL>::accumulate(const ParticlePos& R,
                                   const Matrix<Real>& values,
                                   std::vector<bool>& particles_outside)
{
  int p, v;
  int nparticles = values.size1();
  int nvalues    = values.size2();
  int iu[OHMMS_DIM];
  int buf_index;
  const Real o2pi = 1.0 / (2.0 * M_PI);
  using CoordForm = SpaceGridInput::CoordForm;
  auto& agr       = input_.get_axis_grids();
  std::fill(particles_outside.begin(),particles_outside.end(), true);
  
  switch (input_.get_coord_form())
  {
  case CoordForm::CARTESIAN:
    if (input_.isPeriodic())
    {
      for (p = 0; p < nparticles; p++)
      {
        particles_outside[p] = false;
        Point u              = dot(axinv_, (R[p] - origin_));
        try
        {
          for (int d = 0; d < OHMMS_DIM; ++d)
            iu[d] = gmap_[d].at(floor((u[d] - umin_[d]) * odu_[d]));
        }
        catch (const std::exception& exc)
        {
          std::ostringstream error;
          error << "NESpaceGrid: particle: " << p << " position: " << R[p]
                << " falls outside of the cell, for a period system all particle positions must be in the cell!\n";
          std::throw_with_nested(std::runtime_error(error.str()));
        }
        buf_index = buffer_offset_;
        for (int d = 0; d < OHMMS_DIM; ++d)
          buf_index += nvalues * dm_[d] * iu[d];
        for (v = 0; v < nvalues; v++, buf_index++)
          data_[buf_index] += values(p, v);
      }
    }
    else
    {
      for (p = 0; p < nparticles; p++)
      {
        Point u = dot(axinv_, (R[p] - origin_));
        if (u[0] > umin_[0] && u[0] < umax_[0] && u[1] > umin_[1] && u[1] < umax_[1] && u[2] > umin_[2] &&
            u[2] < umax_[2])
        {
          particles_outside[p] = false;
          iu[0]                = gmap_[0][floor((u[0] - umin_[0]) * odu_[0])];
          iu[1]                = gmap_[1][floor((u[1] - umin_[1]) * odu_[1])];
          iu[2]                = gmap_[2][floor((u[2] - umin_[2]) * odu_[2])];
          buf_index            = buffer_offset_ + nvalues * (dm_[0] * iu[0] + dm_[1] * iu[1] + dm_[2] * iu[2]);
          for (v = 0; v < nvalues; v++, buf_index++)
            data_[buf_index] += values(p, v);
        }
      }
    }
    break;
  case CoordForm::CYLINDRICAL:
    for (p = 0; p < nparticles; p++)
    {
      Point ub = dot(axinv_, (R[p] - origin_));
      Point u{sqrt(ub[0] * ub[0] + ub[1] * ub[1]), static_cast<Real>(atan2(ub[1], ub[0]) * o2pi + .5), ub[2]};
      if (u[0] > umin_[0] && u[0] < umax_[0] && u[1] > umin_[1] && u[1] < umax_[1] && u[2] > umin_[2] &&
          u[2] < umax_[2])
      {
        particles_outside[p] = false;
        iu[0]                = gmap_[0][floor((u[0] - umin_[0]) * odu_[0])];
        iu[1]                = gmap_[1][floor((u[1] - umin_[1]) * odu_[1])];
        iu[2]                = gmap_[2][floor((u[2] - umin_[2]) * odu_[2])];
        buf_index            = buffer_offset_ + nvalues * (dm_[0] * iu[0] + dm_[1] * iu[1] + dm_[2] * iu[2]);
        for (v = 0; v < nvalues; v++, buf_index++)
          data_[buf_index] += values(p, v);
      }
    }
    break;
  case CoordForm::SPHERICAL:
    for (p = 0; p < nparticles; p++)
    {
      Point ub = dot(axinv_, (R[p] - origin_));
      Point u;
      u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2]);
      u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
      u[2] = acos(ub[2] / u[0]) * o2pi * 2.0;
      if (u[0] > umin_[0] && u[0] < umax_[0] && u[1] > umin_[1] && u[1] < umax_[1] && u[2] > umin_[2] &&
          u[2] < umax_[2])
      {
        particles_outside[p] = false;
        iu[0]                = gmap_[0][floor((u[0] - umin_[0]) * odu_[0])];
        iu[1]                = gmap_[1][floor((u[1] - umin_[1]) * odu_[1])];
        iu[2]                = gmap_[2][floor((u[2] - umin_[2]) * odu_[2])];
        buf_index            = buffer_offset_ + nvalues * (dm_[0] * iu[0] + dm_[1] * iu[1] + dm_[2] * iu[2]);
        for (v = 0; v < nvalues; v++, buf_index++)
          data_[buf_index] += values(p, v);
      }
    }
    break;
  default:
    app_log() << "  coordinate type must be cartesian, cylindrical, spherical" << std::endl;
    throw std::runtime_error("SpaceGrid<REAL>::evaluate received an invalid coordinate type");
  }
}

template<typename REAL>
void NESpaceGrid<REAL>::sum(const BufferType& buf, Real* vals)
{
  for (int v = 0; v < nvalues_per_domain_; v++)
  {
    vals[v] = 0.0;
  }
  for (int i = 0, n = buffer_offset_; i < ndomains_; i++, n += nvalues_per_domain_)
  {
    for (int v = 0; v < nvalues_per_domain_; v++)
    {
      vals[v] += data_[n + v];
    }
  }
}

template<typename REAL>
void NESpaceGrid<REAL>::collect(NESpaceGrid& reduction_grid, RefVector<NESpaceGrid> grid_for_each_crowd)
{
  for (NESpaceGrid& crowd_grid : grid_for_each_crowd)
  {
    std::transform(reduction_grid.data_.begin(), reduction_grid.data_.end(), crowd_grid.data_.begin(),
                   reduction_grid.data_.begin(), std::plus<>{});
    crowd_grid.zero();
  }
}

template<typename REAL>
void NESpaceGrid<REAL>::zero()
{
  data_.clear();
}

template<typename REAL>
bool NESpaceGrid<REAL>::check_grid(void)
{
  app_log() << "SpaceGrid<REAL>::check_grid" << std::endl;
  const Real o2pi = 1.0 / (2.0 * M_PI);
  int iu[OHMMS_DIM];
  int idomain;
  bool ok = true;
  Point dc;
  auto& agr = input_.get_axis_grids();
  for (int i = 0; i < ndomains_; i++)
  {
    for (int d = 0; d < OHMMS_DIM; d++)
      dc[d] = domain_centers_(i, d);
    Point ub = dot(axinv_, (dc - origin_));
    Point u;
    using CoordForm = SpaceGridInput::CoordForm;
    switch (input_.get_coord_form())
    {
    case CoordForm::CARTESIAN:
      u = ub;
      break;
    case CoordForm::CYLINDRICAL:
      u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1]);
      u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
      u[2] = ub[2];
      break;
    case CoordForm::SPHERICAL:
      u[0] = sqrt(ub[0] * ub[0] + ub[1] * ub[1] + ub[2] * ub[2]);
      u[1] = atan2(ub[1], ub[0]) * o2pi + .5;
      u[2] = acos(ub[2] / u[0]) * o2pi * 2.0;
      break;
    default:
      break;
    }
    iu[0]   = gmap_[0][floor((u[0] - agr[0].umin) * odu_[0])];
    iu[1]   = gmap_[1][floor((u[1] - agr[1].umin) * odu_[1])];
    iu[2]   = gmap_[2][floor((u[2] - agr[2].umin) * odu_[2])];
    idomain = dm_[0] * iu[0] + dm_[1] * iu[1] + dm_[2] * iu[2];
    if (idomain != i)
    {
      app_log() << "  cell mismatch " << i << " " << idomain << std::endl;
      ok = false;
    }
  }
  if (!ok)
  {
    app_log() << "  NESpaceGrid cells do not map onto themselves" << std::endl;
  }
  app_log() << "end NESpaceGrid<REAL>::check_grid" << std::endl;
  return ok;
}

template class NESpaceGrid<float>;
template class NESpaceGrid<double>;
} // namespace qmcplusplus
