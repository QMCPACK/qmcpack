//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "Message/CommOperators.h"
#include <fftw3.h>
#include "QMCWaveFunctions/einspline_helper.hpp"

namespace qmcplusplus
{

bool EinsplineSetBuilder::ReadGvectors_ESHDF()
{
  bool root = myComm->rank() == 0;
  //this is always ugly
  MeshSize    = 0;
  int hasPsig = 1;
  if (root)
  {
    H5File.readEntry(MeshSize, "/electrons/psi_r_mesh");
    H5File.readEntry(MeshSize, "/electrons/mesh");
  }
  myComm->bcast(MeshSize);
  hasPsig = (MeshSize[0] == 0);
  if (hasPsig)
  {
    int nallowed  = 257;
    int allowed[] = {72,    75,    80,    81,    90,    96,    100,   108,   120,   125,   128,   135,   144,   150,
                     160,   162,   180,   192,   200,   216,   225,   240,   243,   250,   256,   270,   288,   300,
                     320,   324,   360,   375,   384,   400,   405,   432,   450,   480,   486,   500,   512,   540,
                     576,   600,   625,   640,   648,   675,   720,   729,   750,   768,   800,   810,   864,   900,
                     960,   972,   1000,  1024,  1080,  1125,  1152,  1200,  1215,  1250,  1280,  1296,  1350,  1440,
                     1458,  1500,  1536,  1600,  1620,  1728,  1800,  1875,  1920,  1944,  2000,  2025,  2048,  2160,
                     2187,  2250,  2304,  2400,  2430,  2500,  2560,  2592,  2700,  2880,  2916,  3000,  3072,  3125,
                     3200,  3240,  3375,  3456,  3600,  3645,  3750,  3840,  3888,  4000,  4050,  4096,  4320,  4374,
                     4500,  4608,  4800,  4860,  5000,  5120,  5184,  5400,  5625,  5760,  5832,  6000,  6075,  6144,
                     6250,  6400,  6480,  6561,  6750,  6912,  7200,  7290,  7500,  7680,  7776,  8000,  8100,  8192,
                     8640,  8748,  9000,  9216,  9375,  9600,  9720,  10000, 10125, 10240, 10368, 10800, 10935, 11250,
                     11520, 11664, 12000, 12150, 12288, 12500, 12800, 12960, 13122, 13500, 13824, 14400, 14580, 15000,
                     15360, 15552, 15625, 16000, 16200, 16384, 16875, 17280, 17496, 18000, 18225, 18432, 18750, 19200,
                     19440, 19683, 20000, 20250, 20480, 20736, 21600, 21870, 22500, 23040, 23328, 24000, 24300, 24576,
                     25000, 25600, 25920, 26244, 27000, 27648, 28125, 28800, 29160, 30000, 30375, 30720, 31104, 31250,
                     32000, 32400, 32768, 32805, 33750, 34560, 34992, 36000, 36450, 36864, 37500, 38400, 38880, 39366,
                     40000, 40500, 40960, 41472, 43200, 43740, 45000, 46080, 46656, 46875, 48000, 48600, 49152, 50000,
                     50625, 51200, 51840, 52488, 54000, 54675, 55296, 56250, 57600, 58320, 59049, 60000, 60750, 61440,
                     62208, 62500, 64000, 64800, 65536};
    MaxNumGvecs   = 0;
    //    std::set<TinyVector<int,3> > Gset;
    // Read k-points for all G-vectors and take the union
    TinyVector<int, 3> maxIndex(0, 0, 0);
    Gvecs.resize(NumTwists);
    {
      int numg = 0;
      if (root)
      {
        std::ostringstream Gpath;
        Gpath << "/electrons/kpoint_0/gvectors";
        H5File.read(Gvecs[0], Gpath.str());
        numg = Gvecs[0].size();
      }
      myComm->bcast(numg);
      if (!root)
        Gvecs[0].resize(numg);
      myComm->bcast(Gvecs[0]);
      MaxNumGvecs = Gvecs[0].size();
      for (int ig = 0; ig < Gvecs[0].size(); ig++)
      {
        maxIndex[0] = std::max(maxIndex[0], std::abs(Gvecs[0][ig][0]));
        maxIndex[1] = std::max(maxIndex[1], std::abs(Gvecs[0][ig][1]));
        maxIndex[2] = std::max(maxIndex[2], std::abs(Gvecs[0][ig][2]));
      }
      // for (int ig=0; ig<Gvecs.size(); ig++)
      // 	if (Gset.find(Gvecs[ig]) == Gset.end())
      // 	  Gset.insert(Gvecs[ig]);
    } //done with kpoint_0
    MeshSize[0] = (int)std::ceil(4.0 * MeshFactor * maxIndex[0]);
    MeshSize[1] = (int)std::ceil(4.0 * MeshFactor * maxIndex[1]);
    MeshSize[2] = (int)std::ceil(4.0 * MeshFactor * maxIndex[2]);
    //only use 2^a 3^b 5^c where a>=2  up to 65536
    int* ix     = std::lower_bound(allowed, allowed + nallowed, MeshSize[0]);
    int* iy     = std::lower_bound(allowed, allowed + nallowed, MeshSize[1]);
    int* iz     = std::lower_bound(allowed, allowed + nallowed, MeshSize[2]);
    MeshSize[0] = (MeshSize[0] > 128) ? *ix : (MeshSize[0] + MeshSize[0] % 2);
    MeshSize[1] = (MeshSize[1] > 128) ? *iy : (MeshSize[1] + MeshSize[1] % 2);
    MeshSize[2] = (MeshSize[2] > 128) ? *iz : (MeshSize[2] + MeshSize[2] % 2);
    if (Version[0] < 2)
    {
      //get the map for each twist, but use the MeshSize from kpoint_0
      app_log() << "  ESHDF::Version " << Version << std::endl;
      app_log() << "  Assumes distinct Gvecs set for different twists. Regenerate orbital files using updated QE."
                << std::endl;
      for (int k = 0; k < DistinctTwists.size(); ++k)
      {
        int ik = DistinctTwists[k];
        if (ik == 0)
          continue; //already done
        int numg = 0;
        if (root)
        {
          std::ostringstream Gpath;
          Gpath << "/electrons/kpoint_" << ik << "/gvectors";
          H5File.read(Gvecs[ik], Gpath.str());
          numg = Gvecs[ik].size();
        }
        myComm->bcast(numg);
        if (numg == 0)
        {
          //copy kpoint_0, default
          Gvecs[ik] = Gvecs[0];
        }
        else
        {
          if (numg != MaxNumGvecs)
          {
            std::ostringstream o;
            o << "Twist " << ik << ": The number of Gvecs is different from kpoint_0."
              << " This is not supported anymore. Rerun pw2qmcpack.x or equivalent";
            APP_ABORT(o.str());
          }
          if (!root)
            Gvecs[ik].resize(numg);
          myComm->bcast(Gvecs[ik]);
        }
      }
    }
  }
  app_log() << "B-spline mesh factor is " << MeshFactor << std::endl;
  app_log() << "B-spline mesh size is (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";
  app_log() << "Maxmimum number of Gvecs " << MaxNumGvecs << std::endl;
  app_log().flush();
  return hasPsig;
}


} // namespace qmcplusplus
