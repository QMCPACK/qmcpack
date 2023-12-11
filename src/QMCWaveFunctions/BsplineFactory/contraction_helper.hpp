//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CONTRACTION_HELPER_HPP
#define QMCPLUSPLUS_CONTRACTION_HELPER_HPP

namespace qmcplusplus
{
enum SoAFields3D
{
  VAL = 0,
  GRAD0,
  GRAD1,
  GRAD2,
  HESS00,
  HESS01,
  HESS02,
  HESS11,
  HESS12,
  HESS22,
  LAPL,
  NUM_FIELDS
};

/** compute Trace(H*G)
   *
   * gg is symmetrized as
   *  gg[0]=GG(0,0)
   *  gg[1]=GG(0,1)+GG(1,0)
   *  gg[2]=GG(0,2)+GG(2,0)
   *  gg[3]=GG(1,1)
   *  gg[4]=GG(1,2)+GG(2,1)
   *  gg[5]=GG(2,2)
   */
template<typename T>
inline T SymTrace(T h00, T h01, T h02, T h11, T h12, T h22, const T gg[6])
{
  return h00 * gg[0] + h01 * gg[1] + h02 * gg[2] + h11 * gg[3] + h12 * gg[4] + h22 * gg[5];
}

/** compute vector[3]^T x matrix[3][3] x vector[3]
   *
   */
template<typename T>
T v_m_v(T h00, T h01, T h02, T h11, T h12, T h22, T g1x, T g1y, T g1z, T g2x, T g2y, T g2z)
{
  return g1x * g2x * h00 + g1x * g2y * h01 + g1x * g2z * h02 + g1y * g2x * h01 + g1y * g2y * h11 + g1y * g2z * h12 +
      g1z * g2x * h02 + g1z * g2y * h12 + g1z * g2z * h22;
}
/** Coordinate transform for a 3rd rank symmetric tensor representing coordinate derivatives
   *  (hence t3_contract, for contraction with vectors).
   *
   * hijk are the symmetry inequivalent tensor elements, i,j,k range from 0 to 2 for x to z.
   * (gix,giy,giz) are vectors, labelled 1, 2 and 3.  g1 is contracted with the first tensor index,
   * g2 with the second, g3 with the third. 
   *
   * This would be easier with a for loop, but I'm sticking with the convention in this section.
   */
template<typename T>
T t3_contract(T h000,
              T h001,
              T h002,
              T h011,
              T h012,
              T h022,
              T h111,
              T h112,
              T h122,
              T h222,
              T g1x,
              T g1y,
              T g1z,
              T g2x,
              T g2y,
              T g2z,
              T g3x,
              T g3y,
              T g3z)
{
  return h000 * (g1x * g2x * g3x) + h001 * (g1x * g2x * g3y + g1x * g2y * g3x + g1y * g2x * g3x) +
      h002 * (g1x * g2x * g3z + g1x * g2z * g3x + g1z * g2x * g3x) +
      h011 * (g1x * g2y * g3y + g1y * g2x * g3y + g1y * g2y * g3x) +
      h012 *
      (g1x * g2y * g3z + g1x * g2z * g3y + g1y * g2x * g3z + g1y * g2z * g3x + g1z * g2x * g3y + g1z * g2y * g3x) +
      h022 * (g1x * g2z * g3z + g1z * g2x * g3z + g1z * g2z * g3x) + h111 * (g1y * g2y * g3y) +
      h112 * (g1y * g2y * g3z + g1y * g2z * g3y + g1z * g2y * g3y) +
      h122 * (g1y * g2z * g3z + g1z * g2y * g3z + g1z * g2z * g3y) + h222 * (g1z * g2z * g3z);
}

} // namespace qmcplusplus

#endif
