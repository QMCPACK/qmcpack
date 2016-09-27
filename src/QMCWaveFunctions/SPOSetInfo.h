//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#ifndef QMCPLUSPLUS_SPOSET_INFO_H
#define QMCPLUSPLUS_SPOSET_INFO_H

#include <QMCWaveFunctions/SPOInfo.h>
#include <Utilities/IteratorUtility.h>


namespace qmcplusplus
{

  /** collection of orbital info for SPOSet instance or builder
   */
  class SPOSetInfo
  {
  public:
    typedef QMCTraits::RealType RealType;
    typedef std::vector<SPOInfo*> states_t;
    typedef spoinfo::orderings orderings;

    // construction/destruction
    SPOSetInfo(); 
    ~SPOSetInfo(); 
    
    // initialization
    void add(SPOInfo& state);
    void add(SPOInfo* state);
    void add(std::vector<SPOInfo*>& state_vector);
    void add(SPOSetInfo& other);
    /// renders collection immutable, must be called at end of initialization
    void finish(orderings ord=spoinfo::no_order,RealType tol=1e-6);
    
    // initialization queries
    bool complete() const;
    bool partial() const;
    bool has_indices() const;
    bool has_energies() const;

    // array-like read-only access
    int size() const;
    const SPOInfo* operator[](int s) const;
    const SPOInfo* operator[](int s);

    // data access
    int min_index() const;
    int max_index() const;
    RealType energy_tolerance() const;

    // state properties
    bool contiguous() const;
    bool unordered() const;
    bool index_ordered() const;
    bool energy_ordered() const;

    // printing
    void report(const std::string& pad="");

    // templated versions of finish to work with arbitrary vectors
    template<typename SPOI>
    inline void finish(std::vector<SPOI*>& state_vector,orderings ord=spoinfo::no_order,RealType  tol=1e-6)
    {
      for(int i=0;i<state_vector.size();++i)
        add(state_vector[i]);
      finish(ord,tol);
    }

    template<typename SPOI>
    inline void finish(std::vector<int>& subset,std::vector<SPOI*>& state_vector,orderings ord=spoinfo::no_order,RealType  tol=1e-6)
    {
      for(int i=0;i<subset.size();++i)
        add(state_vector[subset[i]]);
      finish(ord,tol);
    }

  private:
    /// whether initialization is complete and SPOSetInfo is ready for use
    bool is_complete;

    /// whether all states have an index assigned
    bool indices_present;

    /// whether all states have an energy assigned
    bool energies_present;

    /// enum for how states are ordered
    orderings order;

    /// tolerance used to sort energies
    RealType energy_tol;

    /// minimum orbital index in the set (w.r.t the full set)
    int index_min;

    /// maximum orbital index in the set (w.r.t the full set)
    int index_max;

    /// collection of SPOInfo
    std::vector<SPOInfo*> states;

    /// sort states by index
    void index_sort();

    /// sort states by energy
    void energy_sort(RealType tol);

    /// count energetic degeneracy of states
    void count_degeneracies();

    /// determine the ordering of the states, if any
    void determine_order(RealType tol);

    /// render collection mutable
    void modify();

    /// empty collection and render mutable
    void clear();

    friend class BasisSetBuilder;
  };


  template<typename SPOI>
  struct SPOSetInfoSimple
  {
    typedef QMCTraits::RealType RealType;
    std::vector<SPOI*> states; //SPOI should derive from SPOInfo

    SPOSetInfoSimple() { }

    ~SPOSetInfoSimple() 
    { 
      delete_iter(states.begin(),states.end());
    }

    inline void add(SPOI* state)
    {
      states.push_back(state);
      state=0;
    }

    inline void clear()
    {
      delete_iter(states.begin(),states.end());
      states.clear();
    }

    int size() const
    {
      return states.size();
    }

    SPOI* operator[](int s) const
    {
      return states[s];
    }

    SPOI*& operator[](int s)
    {
      return states[s];
    }

    void index_sort()
    {
      sort(states.begin(),states.end(),spoinfo::index_order);
    }

    void energy_sort(RealType tol=1e-6,bool assign_indices=false)
    {
      spoinfo::EnergyOrder energy_order(tol);
      sort(states.begin(),states.end(),energy_order);
      if(assign_indices)
        for(int i=0;i<size();++i)
          states[i]->index = i;
    }
  };

}

#endif
