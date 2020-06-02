//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_POSUNIT_H
#define QMCPLUSPLUS_POSUNIT_H

namespace qmcplusplus
{
/** enum class to assist copy and unit conversion operations on position vectors
*/
enum class PosUnit
{
  Cartesian = 0, /*!< indicates that the values are in Cartesian units*/
  Lattice        /*!< indicates that the values are in Lattice units*/
};
/** write unit type in human readable format
 *
 *  This could break tools if they rely on parsing log.
 */
inline std::ostream& operator<<(std::ostream& o_stream, PosUnit pos_unit)
{
  switch (pos_unit)
  {
  case PosUnit::Cartesian:
    o_stream << "Cartesian";
    break;
  case PosUnit::Lattice:
    o_stream << "Lattice";
    break;
  }
  return o_stream;
}

/** Read unit type recorded in int
 *
 *  This should really be human readable
 *  TODO: support both until going to string only.
 */
inline std::istream& operator>>(std::istream& i_stream, PosUnit& pos_unit)
{
  int unit_read;
  i_stream >> unit_read;
  pos_unit = static_cast<PosUnit>(unit_read);
  return i_stream;
}

} // namespace qmcplusplus
#endif
