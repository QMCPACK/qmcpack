//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "ParticleIOUtility.h"
#include "Utilities/ProgressReportEngine.h"
namespace qmcplusplus
{

Lattice createSuperLattice(const Lattice& in, const Tensor<int, OHMMS_DIM>& tmat)
{
  Lattice super_lattice;
  super_lattice.set(dot(tmat, in.R));
  super_lattice.BoxBConds = in.BoxBConds;
  return super_lattice;
}

void expandSuperCell(const ParticleSet& in, const Tensor<int, OHMMS_DIM>& tmat, ParticleSet& out)
{
  app_log() << "  Expanding a simulation cell from " << in.getName() << " to " << out.getName() << std::endl;
  {
    char buff[500];
    snprintf(buff, 500, "   tilematrix= %4d %4d %4d %4d %4d %4d %4d %4d %4d\n", tmat[0], tmat[1], tmat[2], tmat[3],
             tmat[4], tmat[5], tmat[6], tmat[7], tmat[8]);
    app_log() << buff << std::endl;
  }

  if (in.R.InUnit != PosUnit::Lattice)
    throw std::runtime_error("The coordinates of input particle set was not in PosUnit::Lattice!");

  const auto& prim_cell(in.getLattice());
  const int natoms    = in.getTotalNum();
  const int numCopies = std::abs(det(tmat));

  const auto& primPos(in.R);
  const auto& primTypes(in.GroupID);
  out.create(natoms * numCopies);
  const int maxCopies = 10;
  int index           = 0;
  //set the unit to the Cartesian!
  out.R.InUnit = PosUnit::Cartesian;
  app_log() << "  Reduced coord    Cartesion coord    species.\n";
  for (int ns = 0; ns < in.getSpeciesSet().getTotalNum(); ++ns)
  {
    for (int i0 = -maxCopies; i0 <= maxCopies; i0++)
      for (int i1 = -maxCopies; i1 <= maxCopies; i1++)
        for (int i2 = -maxCopies; i2 <= maxCopies; i2++)
          for (int iat = 0; iat < primPos.size(); iat++)
          {
            if (primTypes[iat] != ns)
              continue;
            auto uPrim = primPos[iat];
            for (int i = 0; i < 3; i++)
              uPrim[i] -= std::floor(uPrim[i]);
            auto r = prim_cell.toCart(uPrim) + (double)i0 * prim_cell.a(0) + (double)i1 * prim_cell.a(1) +
                (double)i2 * prim_cell.a(2);
            auto uSuper = out.getLattice().toUnit(r);
            if ((uSuper[0] >= -1.0e-6) && (uSuper[0] < 0.9999) && (uSuper[1] >= -1.0e-6) && (uSuper[1] < 0.9999) &&
                (uSuper[2] >= -1.0e-6) && (uSuper[2] < 0.9999))
            {
              char buff[500];
              snprintf(buff, 500, "  %10.4f  %10.4f %10.4f   %12.6f %12.6f %12.6f %d\n", uSuper[0], uSuper[1],
                       uSuper[2], r[0], r[1], r[2], ns);
              app_log() << buff;
              out.R[index]       = r;
              out.GroupID[index] = ns; //primTypes[iat];
              index++;
            }
          }
  }

  app_log() << "  Simulationcell after tiling" << std::endl;
  out.getLattice().print(app_log());
  app_log() << std::endl;

  out.resetGroups();
}
} // namespace qmcplusplus
