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
#include "BsplineSet.h"

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
class BsplineReader
{
  /** create the actual spline sets
   */
  virtual std::unique_ptr<SPOSet> create_spline_set(const std::string& my_name,
                                                    int spin,
                                                    const BandInfoGroup& bandgroup) = 0;

  void initialize_spo2band(const std::string& spo_name,
                           int spin,
                           const std::vector<BandInfo>& bigspace,
                           SPOSetInfo& sposet,
                           std::vector<int>& band2spo);

public:
  BsplineReader(EinsplineSetBuilder* e, bool use_duplex_splines);

  virtual ~BsplineReader();

  /** setting common parameters
   */
  void setCommon(xmlNodePtr cur);

  /** create the spline after one of the kind is created */
  std::unique_ptr<SPOSet> create_spline_set(const std::string& spo_name, int spin, SPOSetInputInfo& input_info);

  /** create the spline set */
  std::unique_ptr<SPOSet> create_spline_set(const std::string& spo_name, int spin, const size_t size);

  /** Set the checkNorm variable */
  inline void setCheckNorm(bool new_checknorm) { checkNorm = new_checknorm; };

  /** Set the orbital rotation flag. Rotations are applied to balance the real/imaginary components. */
  inline void setRotate(bool new_rotate) { rotate = new_rotate; };

  ///pointer to the EinsplineSetBuilder
  EinsplineSetBuilder* mybuilder;
  ///communicator
  Communicate* myComm;

protected:
  ///check the norm of orbitals
  bool checkNorm;
  ///save spline coefficients to storage
  bool saveSplineCoefs;
  ///apply orbital rotations
  bool rotate;
  ///map from spo index to band index
  std::vector<std::vector<int>> spo2band;
  /// if true, use two real-valued splines for one complex-valued DFT orbital.
  const bool use_duplex_splines_;
  /// if true, use offload
  bool use_offload;

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
  static void set_grid(const TinyVector<int, 3>& mesh_sizes, const TinyVector<int, 3>& halfg, GT* xyz_grid, BCT* xyz_bc)
  {
    for (int j = 0; j < 3; ++j)
    {
      xyz_grid[j].start = 0.0;
      xyz_grid[j].end   = 1.0;
      xyz_grid[j].num   = mesh_sizes[j];

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
  }

  /** initialize twist-related data for N orbitals
   */
  inline void check_twists(BsplineSet& bspline, const BandInfoGroup& bandgroup) const
  {
    const int N       = bandgroup.getNumDistinctOrbitals();
    const int numOrbs = bandgroup.getNumSPOs();

    const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
    for (int iorb = 0, num = 0; iorb < N; iorb++)
    {
      int ti                      = cur_bands[iorb].TwistIndex;
      bspline.kPoints[iorb]       = mybuilder->PrimCell.k_cart(-mybuilder->primcell_kpoints[ti]);
      bspline.MakeTwoCopies[iorb] = (num < (numOrbs - 1)) && cur_bands[iorb].MakeTwoCopies;
      num += bspline.MakeTwoCopies[iorb] ? 2 : 1;
    }

    bspline.resize_kpoints();
    app_log() << "NumDistinctOrbitals " << N << " numOrbs = " << numOrbs << std::endl;
  }

  /** compute sign bits at the G/2 boundaries
   * no supercell, no k-point folding, single special k point (G, L ...)
   */
  static TinyVector<int, 3> computeHalfG(const TinyVector<int, OHMMS_DIM>& bconds,
                                         const std::vector<TinyVector<double, OHMMS_DIM>>& primcell_kpoints,
                                         size_t twist0_index)
  {
    TinyVector<int, 3> halfG;
    const auto& twist0 = primcell_kpoints[twist0_index];
    app_log() << "  TwistIndex = " << twist0_index << " TwistAngle " << twist0 << std::endl;
    for (int i = 0; i < 3; i++)
      if (bconds[i] && ((std::abs(std::abs(twist0[i]) - 0.5) < 1.0e-8)))
        halfG[i] = 1;
      else
        halfG[i] = 0;
    app_log() << "   HalfG = " << halfG << std::endl;
    return halfG;
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

  /** create data space in the spline object and try open spline dump files.
   * @param bandgroup band info
   * @param bspline the spline object being worked on
   * @return true if dumpfile pass class name and data type size check
   */
  bool lookforSplineDataDumpFile(const BandInfoGroup& bandgroup,
                                 const std::string& keyword,
                                 size_t datatype_size) const;

  /** read planewave coefficients from h5 file
   * @param s data set full path in h5
   * @param h5f hdf5 file handle
   * @param cG vector to store coefficients
   */
  void readOneOrbitalCoefs(const std::string& s, hdf_archive& h5f, Vector<std::complex<double>>& cG) const;
};

} // namespace qmcplusplus
#endif
