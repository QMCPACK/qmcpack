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


#ifndef OHMMS_RANDOMNUMBERCONTROL_H__
#define OHMMS_RANDOMNUMBERCONTROL_H__
#include "OhmmsData/OhmmsElementBase.h"
#include "Utilities/RandomGenerator.h"
#include "Utilities/PrimeNumberSet.h"
#include <io/hdf_archive.h>

class Communicate;

namespace qmcplusplus
{

/**class RandomNumberControl
 *\brief Encapsulate data to initialize and save the status of the random number generator
 *
 * Default:  myName = "random"
 * 2007-12-01
 *   Use PrimeNumbers to generate random seeds.
 */
class RandomNumberControl : public OhmmsElementBase
{

public:

  typedef RandomGenerator_t::uint_type uint_type;
  static PrimeNumberSet<uint_type> PrimeNumbers;
  //children random number generator
  static std::vector<RandomGenerator_t*>  Children;

  /// constructors and destructors
  RandomNumberControl(const char* aname="random");

  bool get(std::ostream& os) const;
  bool put(std::istream& is);
  bool put(xmlNodePtr cur);
  void reset();
  static void test();

  static void make_seeds();
  static void make_children();

  xmlNodePtr initialize(xmlXPathContextPtr);

  /** read in parallel or serial
   * @param fname file name
   * @param comm communicator
   */
  static void read(const std::string& fname, Communicate* comm);
  /** write in parallel or serial
   * @param fname file name
   * @param comm communicator
   */
  static void write(const std::string& fname, Communicate* comm);
  /** read random state from a hdf file in parallel
   * @param hdf_archive set to parallel
   * @param comm communicator
   */
  static void read_parallel(hdf_archive& hin, Communicate* comm);
  /** write random state to a hdf file in parallel
   * @param hdf_archive set to parallel
   * @param comm communicator
   */
  static void write_parallel(hdf_archive& hout, Communicate* comm);
  /** rank 0 reads random states from a hdf file
   * and distributes them to all the other ranks
   * @param hdf_archive set to serial
   * @param comm communicator
   */
  static void read_rank_0(hdf_archive& hin, Communicate* comm);
  /** rank 0 gathers the random states from all the other ranks
   * and write them to a hdf file
   * @param hin hdf_archive object set to serial
   * @param comm communicator
   */
  static void write_rank_0(hdf_archive& hout, Communicate* comm);
   /** read random state from a xml file
   * @param fname file name
   * @param comm communicator
   */
  static void read_old(const std::string& fname, Communicate* comm);
   /** write random state to a xml file
   * @param fname file name
   * @param comm communicator
   */
  static void write_old(const std::string& fname, Communicate* comm);
private:

  bool NeverBeenInitialized;
  xmlNodePtr myCur;
  static uint_type Offset;
};
}

#endif

