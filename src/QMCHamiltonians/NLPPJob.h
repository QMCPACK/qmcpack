//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_NLPPJob_H
#define QMCPLUSPLUS_NLPPJob_H

#include "OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{
/** meta data for NLPP calculation of a pair of ion and electron
   */
template<typename T>
struct NLPPJob
{
  using RealType = T;
  using PosType  = TinyVector<RealType, 3>;
  const int ion_id;
  const int electron_id;
  const PosType elec_pos;
  const RealType ion_elec_dist;
  const PosType ion_elec_displ;

  NLPPJob(const int ion_id_in,
          const int electron_id_in,
          const PosType& elec_pos_in,
          const RealType ion_elec_dist_in,
          const PosType& ion_elec_displ_in)
      : ion_id(ion_id_in),
        electron_id(electron_id_in),
        elec_pos(elec_pos_in),
        ion_elec_dist(ion_elec_dist_in),
        ion_elec_displ(ion_elec_displ_in)
  {}
};
} // namespace qmcplusplus
#endif
