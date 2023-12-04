//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file BsplineReader.h
 *
 * base class to read data and manage spline tables
 */
#ifndef QMCPLUSPLUS_BSPLINE_READER_H
#define QMCPLUSPLUS_BSPLINE_READER_H

#include <vector>
#include <einspline/bspline_base.h>
#include <BandInfo.h>
#include "EinsplineSetBuilder.h"

namespace qmcplusplus
{
struct SPOSetInputInfo;

/**
 * Each SplineC2X needs a reader derived from BsplineReader.
 * This base class handles common chores
 * - check_twists : read gvectors, set twists for folded bands if needed, and set the phase for the special K
 * - set_grid : create the basic grid and boundary conditions for einspline
 * Note that template is abused but it works.
 */
struct BsplineReader
{
  ///pointer to the EinsplineSetBuilder
  EinsplineSetBuilder* mybuilder;
  ///communicator
  Communicate* myComm;
  ///mesh size
  ///check the norm of orbitals
  bool checkNorm;
  ///save spline coefficients to storage
  bool saveSplineCoefs;
  ///apply orbital rotations
  bool rotate;
  ///map from spo index to band index
  std::vector<std::vector<int>> spo2band;

  BsplineReader(EinsplineSetBuilder* e);

  virtual ~BsplineReader();

  std::string getSplineDumpFileName(const BandInfoGroup& bandgroup) const
  {
    auto& MeshSize = mybuilder->MeshSize;
    std::ostringstream oo;
    oo << bandgroup.myName << ".g" << MeshSize[0] << "x" << MeshSize[1] << "x" << MeshSize[2] << ".h5";
    return oo.str();
  }

  /** read gvectors and set the mesh, and prepare for einspline
   */
  template<typename GT, typename BCT>
  inline bool set_grid(const TinyVector<int, 3>& halfg, GT* xyz_grid, BCT* xyz_bc) const
  {
    //This sets MeshSize from the input file
    bool havePsig = mybuilder->ReadGvectors_ESHDF();

    for (int j = 0; j < 3; ++j)
    {
      xyz_grid[j].start = 0.0;
      xyz_grid[j].end   = 1.0;
      xyz_grid[j].num   = mybuilder->MeshSize[j];

      if (halfg[j])
      {
        xyz_bc[j].lCode = ANTIPERIODIC;
        xyz_bc[j].rCode = ANTIPERIODIC;
      }
      else
      {
        xyz_bc[j].lCode = PERIODIC;
        xyz_bc[j].rCode = PERIODIC;
      }

      xyz_bc[j].lVal = 0.0;
      xyz_bc[j].rVal = 0.0;
    }
    return havePsig;
  }

  /** initialize twist-related data for N orbitals
   */
  template<typename SPE>
  inline void check_twists(SPE& bspline, const BandInfoGroup& bandgroup) const
  {
    //init(orbitalSet,bspline);
    bspline.PrimLattice = mybuilder->PrimCell;
    bspline.GGt         = dot(transpose(bspline.PrimLattice.G), bspline.PrimLattice.G);

    int N       = bandgroup.getNumDistinctOrbitals();
    int numOrbs = bandgroup.getNumSPOs();

    bspline.setOrbitalSetSize(numOrbs);
    bspline.resizeStorage(N, N);

    bspline.first_spo = bandgroup.getFirstSPO();
    bspline.last_spo  = bandgroup.getLastSPO();

    int num                                = 0;
    const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
    for (int iorb = 0; iorb < N; iorb++)
    {
      int ti                      = cur_bands[iorb].TwistIndex;
      bspline.kPoints[iorb]       = mybuilder->PrimCell.k_cart(-mybuilder->primcell_kpoints[ti]);
      bspline.MakeTwoCopies[iorb] = (num < (numOrbs - 1)) && cur_bands[iorb].MakeTwoCopies;
      num += bspline.MakeTwoCopies[iorb] ? 2 : 1;
    }

    app_log() << "NumDistinctOrbitals " << N << " numOrbs = " << numOrbs << std::endl;

    bspline.HalfG             = 0;
    TinyVector<int, 3> bconds = mybuilder->TargetPtcl.getLattice().BoxBConds;
    if (!bspline.isComplex())
    {
      //no k-point folding, single special k point (G, L ...)
      TinyVector<double, 3> twist0 = mybuilder->primcell_kpoints[bandgroup.TwistIndex];
      for (int i = 0; i < 3; i++)
        if (bconds[i] && ((std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)))
          bspline.HalfG[i] = 1;
        else
          bspline.HalfG[i] = 0;
      app_log() << "  TwistIndex = " << cur_bands[0].TwistIndex << " TwistAngle " << twist0 << std::endl;
      app_log() << "   HalfG = " << bspline.HalfG << std::endl;
    }
    app_log().flush();
  }

  /** return the path name in hdf5
   * @param ti twist index
   * @param spin spin index
   * @param ib band index
   */
  inline std::string psi_g_path(int ti, int spin, int ib) const
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << ib << "/psi_g";
    return path.str();
  }

  /** return the path name in hdf5
   * @param ti twist index
   * @param spin spin index
   * @param ib band index
   */
  inline std::string psi_r_path(int ti, int spin, int ib) const
  {
    std::ostringstream path;
    path << "/electrons/kpoint_" << ti << "/spin_" << spin << "/state_" << ib << "/psi_r";
    return path.str();
  }

  /** create the actual spline sets
   */
  virtual std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                                    int spin,
                                                    const BandInfoGroup& bandgroup) = 0;

  /** setting common parameters
   */
  void setCommon(xmlNodePtr cur);

  /** create the spline after one of the kind is created */
  std::unique_ptr<SPOSet> create_spline_set(int spin, xmlNodePtr cur, SPOSetInputInfo& input_info);

  /** create the spline set */
  std::unique_ptr<SPOSet> create_spline_set(int spin, xmlNodePtr cur);

  /** Set the checkNorm variable */
  inline void setCheckNorm(bool new_checknorm) { checkNorm = new_checknorm; };

  /** Set the orbital rotation flag. Rotations are applied to balance the real/imaginary components. */
  inline void setRotate(bool new_rotate) { rotate = new_rotate; };

  void initialize_spo2band(int spin,
                           const std::vector<BandInfo>& bigspace,
                           SPOSetInfo& sposet,
                           std::vector<int>& band2spo);
};

} // namespace qmcplusplus
#endif
