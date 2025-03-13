//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter  W. Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "ParseGridInput.hpp"

#include <cmath>

#include <ModernStringUtils.hpp>
#include "Message/UniformCommunicateError.h"

namespace qmcplusplus
{

template<typename REAL>
AxisGrid<REAL> parseGridInput(std::istringstream& grid_input_stream)
{
  using namespace modernstrutil;
  const REAL utol = 1e-5;

  std::string grid_line;
  std::getline(grid_input_stream, grid_line);
  // This refactored from QMCHamiltonian/SpaceGrid::initializeRectilinear
  std::vector<std::string_view> tokens = split(grid_line, " ");

  //  count the number of intervals
  int nintervals = 0;
  for (int i = 0; i < tokens.size(); i++)
    if (tokens[i][0] != '(')
      nintervals++;
  nintervals--;
  //  allocate temporary interval variables
  AxisGrid<REAL> agr;
  agr.ndom_int.resize(nintervals);
  agr.du_int.resize(nintervals);
  agr.ndu_int.resize(nintervals);
  //  determine number of domains in each interval and the width of each domain
  REAL u1 = string2Real<REAL>(tokens[0]);

  agr.umin = u1;
  if (std::abs(u1) > 1.0000001)
  {
    std::ostringstream msg;
    msg << "  interval endparticles cannot be greater than " << 1 << '\n' << "  endpoint provided: " << u1;
    throw UniformCommunicateError(msg.str());
  }
  bool is_int        = false;
  bool has_paren_val = false;
  REAL du_i;
  int ndom_i   = 1;
  int interval = -1;
  for (int i = 1; i < tokens.size(); i++)
  {
    if (tokens[i][0] != '(')
    {
      REAL u2  = string2Real<REAL>(tokens[i]);
      agr.umax = u2;
      if (!has_paren_val)
        du_i = u2 - u1;
      has_paren_val = false;
      interval++;
      if (u2 < u1)
      {
        std::ostringstream msg;

        msg << "  interval (" << u1 << "," << u2 << ") is negative" << std::endl;
        throw UniformCommunicateError(msg.str());
      }
      if (std::abs(u2) > 1.0000001)
      {
        std::ostringstream msg;

        msg << "  interval endparticles cannot be greater than " << 1 << std::endl;
        msg << "  endpoint provided: " << u2 << std::endl;
        throw UniformCommunicateError(msg.str());
      }
      if (is_int)
      {
        agr.du_int[interval]   = (u2 - u1) / ndom_i;
        agr.ndom_int[interval] = ndom_i;
      }
      else
      {
        agr.du_int[interval]   = du_i;
        agr.ndom_int[interval] = std::floor((u2 - u1) / du_i + .5);
        if (std::abs(u2 - u1 - du_i * agr.ndom_int[interval]) > utol)
        {
          std::ostringstream msg;

          msg << "  interval (" << u1 << "," << u2 << ") not divisible by du=" << du_i << std::endl;
          throw UniformCommunicateError(msg.str());
        }
      }
      u1 = u2;
    }
    else
    {
      has_paren_val = true;
      std::string paren_val{tokens[i].substr(1, tokens[i].length() - 2)};
      is_int = tokens[i].find(".") == std::string::npos;
      if (is_int)
      {
        ndom_i = string2Int<int>(paren_val);
        du_i   = -1.0;
      }
      else
      {
        ndom_i = 0;
        du_i   = string2Real<REAL>(paren_val);
      }
    }
  }
  // find the smallest domain width
  REAL du_min = 1.0;
  for (int i = 0; i < agr.du_int.size(); i++)
    du_min = std::min(du_min, agr.du_int[i]);
  agr.odu = 1.0 / du_min;
  // make sure it divides into all other domain widths
  for (int i = 0; i < agr.du_int.size(); i++)
  {
    agr.ndu_int[i] = floor(agr.du_int[i] / du_min + .5);
    if (std::abs(agr.du_int[i] - agr.ndu_int[i] * du_min) > utol)
    {
      std::ostringstream msg;
      msg << "interval " << i + 1 << " of axis is not divisible by smallest subinterval " << du_min;
      throw UniformCommunicateError(msg.str());
    }
  }
  // set up the interval map such that gmap[u/du]==domain index
  agr.gmap.resize(floor((agr.umax - agr.umin) * agr.odu + .5));
  int n  = 0;
  int nd = -1;
  for (int i = 0; i < agr.ndom_int.size(); i++)
    for (int j = 0; j < agr.ndom_int[i]; j++)
    {
      nd++;
      for (int k = 0; k < agr.ndu_int[i]; k++)
      {
        //app_log()<<"        accessing gmap "<<iaxis<<" "<<n<< std::endl;
        agr.gmap[n] = nd;
        n++;
      }
    }
  agr.dimensions = nd + 1;
  //end read in the grid contents
  //save interval width information
  int ndom_tot = 0;
  for (int i = 0; i < agr.ndom_int.size(); i++)
    ndom_tot += agr.ndom_int[i];
  agr.ndu_per_interval.resize(ndom_tot);
  int idom = 0;
  for (int i = 0; i < agr.ndom_int.size(); i++)
    for (int ii = 0; ii < agr.ndom_int[i]; ii++)
    {
      agr.ndu_per_interval[idom] = agr.ndu_int[i];
      idom++;
    }
  return agr;
}

template<typename REAL>
AxisGrid<REAL>::AxisGrid(std::vector<int>&& rhs_ndom_int,
                         std::vector<int>&& rhs_ndu_int,
                         std::vector<REAL>&& rhs_du_int,
                         REAL rhs_umin,
                         REAL rhs_umax,
                         REAL rhs_odu,
                         std::vector<int>&& rhs_gmap,
                         std::vector<int>&& rhs_ndu_per_interval,
                         int rhs_dimensions)
{
  ndom_int         = rhs_ndom_int;
  ndu_int          = rhs_ndu_int;
  du_int           = rhs_du_int;
  umin             = rhs_umin;
  umax             = rhs_umax;
  odu              = rhs_odu;
  gmap             = rhs_gmap;
  ndu_per_interval = rhs_ndu_per_interval;
  dimensions       = rhs_dimensions;
}

template<typename REAL>
AxisGrid<REAL>::AxisGrid(const AxisGrid& rhs)
{
  ndom_int         = rhs.ndom_int;
  ndu_int          = rhs.ndu_int;
  du_int           = rhs.du_int;
  umin             = rhs.umin;
  umax             = rhs.umax;
  odu              = rhs.odu;
  gmap             = rhs.gmap;
  ndu_per_interval = rhs.ndu_per_interval;
  dimensions       = rhs.dimensions;
}

template<typename REAL>
AxisGrid<REAL>& AxisGrid<REAL>::operator=(AxisGrid<REAL> rhs)
{
  std::swap(ndom_int, rhs.ndom_int);
  std::swap(ndu_int, rhs.ndu_int);
  std::swap(du_int, rhs.du_int);
  std::swap(umin, rhs.umin);
  std::swap(umax, rhs.umax);
  std::swap(odu, rhs.odu);
  std::swap(gmap, rhs.gmap);
  std::swap(ndu_per_interval, rhs.ndu_per_interval);
  std::swap(dimensions, rhs.dimensions);
  return *this;
}


// explicit instantiation
template struct AxisGrid<float>;
template struct AxisGrid<double>;
template AxisGrid<double> parseGridInput<double>(std::istringstream& grid_input_stream);
template AxisGrid<float> parseGridInput<float>(std::istringstream& grid_input_stream);
} // namespace qmcplusplus
