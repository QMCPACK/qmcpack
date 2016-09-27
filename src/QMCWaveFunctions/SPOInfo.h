//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_SPO_INFO_H
#define QMCPLUSPLUS_SPO_INFO_H

#include <Configuration.h>

namespace qmcplusplus
{

  /** base class to describe a single orbital in an SPOSet
   */
  struct SPOInfo
  {
    typedef QMCTraits::RealType RealType;

    enum{DIM=OHMMS_DIM};

    static const int      no_index     ;
    static const int      no_degeneracy;
    static const RealType no_energy    ;
      
    /// original orbital index in the maximal basis
    int index;

    /// energetic degeneracy of the orbital
    int degeneracy;

    /// energy of the orbital (in Hartree units)
    RealType energy;

    SPOInfo();

    SPOInfo(int orb_index,RealType en=no_energy);
  
    virtual ~SPOInfo() { }


    inline bool has_index()
    {
      return index!=no_index;
    }

    inline bool has_energy()
    {
      return energy!=no_energy;
    }

    inline bool has_degeneracy()
    {
      return degeneracy!=no_degeneracy;
    }

    inline SPOInfo* copy()
    {
      return new SPOInfo(*this);
    }

    /// write orbital info to stdout
    void report(const std::string& pad="");
  };


  namespace spoinfo
  {
    /// enumeration of possible orbital info orderings
    enum orderings{index_ordered=0,energy_ordered,energy_and_index_ordered,unordered,no_order};

    /// comparison function for sorting SPOInfo based on orbital index
    inline bool index_order(const SPOInfo* left, const SPOInfo* right)
    {
      return left->index < right->index;
    }

    /// comparison functor for sorting SPOInfo based on energy
    struct EnergyOrder
    {
      typedef QMCTraits::RealType RealType;
      RealType energy_tol;
      EnergyOrder(RealType tol=1e-6) : energy_tol(tol) {};
      ~EnergyOrder() { };
      inline bool operator()(const SPOInfo* left, const SPOInfo* right)
      {
        if(std::abs(left->energy-right->energy)<energy_tol)
          return index_order(left,right);
        else
          return left->energy < right->energy;
      }
    };
  }


}

#endif
