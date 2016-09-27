//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file BandInfo.h
 *
 * Define helper class to sort bands according to the band and twist and functions
 */
#ifndef QMCPLUSPLUS_BANDINFO_H
#define QMCPLUSPLUS_BANDINFO_H
#include <Configuration.h>

namespace qmcplusplus
{

class SPOSetInfo;

struct BandInfo
{
  ///twist index
  int TwistIndex;
  ///band index 
  int BandIndex;
  /// band group index
  int BandGroup;
  ///spin index 
  int Spin;
  ///energy associated with this band
  double Energy;
  ///  This is true if we should make distinct copies represeninting a +k, -k pair
  bool MakeTwoCopies;
  /// True if this state is a core state
  bool IsCoreState;
  ///default constructor
  BandInfo()
    :TwistIndex(0),BandIndex(-1),BandGroup(0),Spin(0),Energy(1e9),MakeTwoCopies(false),IsCoreState(false)
  { }

  /** operator to determine the order of any band
   * 
   * - energy
   * - twist index if the energies are similar
   */
  inline bool operator<(BandInfo other) const
  {
    if  ((Energy < other.Energy+1e-6)&&(Energy > other.Energy-1e-6))
    {
      if (TwistIndex == other.TwistIndex)
        return BandIndex<other.BandIndex;
      else
        return TwistIndex < other.TwistIndex;
    }
    else
      return Energy < other.Energy;
  }
};

/** a group of bands 
 */
struct BandInfoGroup
{
  ///index of the group
  int GroupID;
  ///starting SPO
  int FirstSPO;
  ///number of SPOs handled by this object
  int NumSPOs;
  ///starting band
  int FirstBand;
  ///twist index set by the full band not by the subset
  int TwistIndex;
  ///Bands that belong to this group
  std::vector<BandInfo> myBands;
  ///name of this band
  std::string myName;
  ///constructor
  BandInfoGroup();
  ///return the size of this band
  inline int getNumDistinctOrbitals() const { return myBands.size();}
  ///return the indext of the first SPO set
  inline int getFirstSPO() const { return FirstSPO;}
  ///return the indext of the last SPO set
  inline int getLastSPO() const { return NumSPOs+FirstSPO;}
  ///return the number of SPOs
  inline int getNumSPOs() const { return NumSPOs;}

  /** select the bands within an energy range [emin,emax)
   *
   * @param bigspace a set of sorted bands
   * @param emin minimum energy
   * @param emax maxmimum energy
   */
  void selectBands(const std::vector<BandInfo>& bigspace, double emin, double emax);

  /** get the bands within [first_spo,first_spo+num_spos)
   * @param bigspace a set of sorted bands
   * @param first_orb index of the first uniquie orbitals
   * @param num_spos number of SPOs to be created
   * @param relative if(relative) FirstSPO is set to any valid state index  \f$[0,\infty)\f$
   */
  void selectBands(const std::vector<BandInfo>& bigspace, int first_orb, int num_spos, bool relative);
};

}

#endif
