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
    
    

#include <QMCWaveFunctions/SPOSetInfo.h>
#include <limits>


namespace qmcplusplus
{

  typedef QMCTraits::RealType RealType;
  using namespace spoinfo;

  const RealType no_energy_tol = std::numeric_limits<int>::max();
  const int      no_index      = -1;


  // construction/destruction
  SPOSetInfo::SPOSetInfo() 
  { 
    modify();
  }

  SPOSetInfo::~SPOSetInfo() 
  {
    delete_iter(states.begin(),states.end());
  }


  // initialization
  void SPOSetInfo::add(SPOInfo& state)
  {
    states.push_back(state.copy());
  }
    
  void SPOSetInfo::add(SPOInfo* state)
  {
    states.push_back(state->copy());
  }

  void SPOSetInfo::add(std::vector<SPOInfo*>& state_vector)
  {
    for(int i=0;i<state_vector.size();++i)
      add(state_vector[i]);
  }

  void SPOSetInfo::add(SPOSetInfo& other)
  {
    states_t::iterator so=other.states.begin(),so_end=other.states.end();
    for(;so!=so_end;++so)
      add((*so));
  }

  void SPOSetInfo::finish(orderings ord,RealType tol)
  {
    if(!complete())
    {
      index_min = std::numeric_limits<int>::max();
      index_max = std::numeric_limits<int>::min();
      indices_present  = true;
      energies_present  = true;
      states_t::iterator i=states.begin(),i_end=states.end();
      for(;i!=i_end;++i)
      {
        SPOInfo& s = **i;
        index_min = std::min(s.index,index_min);
        index_max = std::max(s.index,index_max);
        indices_present  &= s.has_index();
        energies_present &= s.has_energy();
      }
      if(!has_indices())
        APP_ABORT("SPOSetInfo::finish\n  indices have not been assigned to states as required\n  this is a developer error");
      if(ord==spoinfo::index_ordered)
        index_sort();
      else if(ord==spoinfo::energy_ordered)
        energy_sort(tol);
      else
        determine_order(tol);
      is_complete = true;
    }
  }


  // initialization queries
  bool SPOSetInfo::complete() const
  {
    return is_complete;
  }

  bool SPOSetInfo::partial() const
  {
    return size()>0 && !complete();
  }

  bool SPOSetInfo::has_indices() const
  {
    return indices_present;
  }

  bool SPOSetInfo::has_energies() const
  {
    return energies_present;
  }


  // array-like read-only access
  int SPOSetInfo::size() const
  {
    return states.size();
  }

  const SPOInfo* SPOSetInfo::operator[](int s) const
  {
    return states[s];
  }

  const SPOInfo* SPOSetInfo::operator[](int s)
  {
    return states[s];
  }


  // data access
  int SPOSetInfo::min_index() const
  {
    return index_min;
  }

  int SPOSetInfo::max_index() const
  {
    return index_max;
  }

  RealType SPOSetInfo::energy_tolerance() const
  {
    return energy_tol;
  }


  // state properties
  bool SPOSetInfo::contiguous() const
  {
    return index_max-index_min==size()-1;
  }

  bool SPOSetInfo::unordered() const
  {
    return order==spoinfo::unordered;
  }

  bool SPOSetInfo::index_ordered() const
  {
    return order==spoinfo::index_ordered || energy_and_index_ordered;
  }

  bool SPOSetInfo::energy_ordered() const
  {
    return order==spoinfo::energy_ordered || energy_and_index_ordered;
  }


  // printing
  void SPOSetInfo::report(const std::string& pad)
  {
    app_log()<<pad<<"complete           = "<< complete()           << std::endl;
    if(complete())
    {
      app_log()<<pad<<"min_index  = "<< min_index()  << std::endl;
      app_log()<<pad<<"max_index  = "<< max_index()  << std::endl;
      app_log()<<pad<<"contiguous         = "<< contiguous()         << std::endl;
      app_log()<<pad<<"has_energies       = "<< has_energies()       << std::endl;
      app_log()<<pad<<"energy_ordered     = "<< energy_ordered()     << std::endl;
      if(energy_ordered())
        app_log()<<pad<<"energy_tolerance   = "<< energy_tolerance()   << std::endl;
      app_log()<<pad<<"# of states        = "<<size()<< std::endl;
      app_log()<<pad<<"state information:"<< std::endl;
      states_t::iterator s=states.begin(),s_end=states.end();
      std::string pad2 = pad+"    ";
      int ns=0;
      for(;s!=s_end;++s,++ns)
      {
        app_log()<<pad<<"  state "<<ns<< std::endl;
        (*s)->report(pad2);
      }
    }
    app_log().flush();
  }


  void SPOSetInfo::index_sort()
  {
    sort(states.begin(),states.end(),index_order);
    order = spoinfo::index_ordered;
  }

  void SPOSetInfo::energy_sort(RealType tol)
  {
    if(!has_energies())
      APP_ABORT("SHOSetInfo::energy_sort  not all states have an energy assigned");
    energy_tol = tol;
    EnergyOrder energy_order(energy_tol);
    sort(states.begin(),states.end(),energy_order);
    order = spoinfo::energy_ordered;
    count_degeneracies();
  }

  void SPOSetInfo::count_degeneracies()
  {
    if(energy_ordered())
    {
      states_t::iterator stmp,s=states.begin(),s_end=states.end();
      while(s!=s_end)
      {
        int g=1;
        stmp=s;
        ++stmp;
        // look ahead to count
        while(stmp!=s_end && std::abs((*stmp)->energy-(*s)->energy)<energy_tol)
        {
          g++;
          ++stmp;
        }
        (*s)->degeneracy = g;
        //(*s)->degeneracy_index = g-1;
        // run over degenerate states to assign
        for(int n=0;n<g;++n)
        {
          (*s)->degeneracy       = g;
          //(*s)->degeneracy_index = g-1-n;
          ++s;
        }
      }
    }
  }

  void SPOSetInfo::determine_order(RealType tol)
  {
    bool index_ord  = false;
    bool energy_ord = false;
    states_t::iterator s_prev,s;
    s=states.begin();
    for(int i=1;i<size();++i)
    {
      s_prev=s;
      ++s;
      index_ord &= index_order(*s_prev,*s);
    }
    if(has_energies())
    {
      energy_tol = tol;
      EnergyOrder energy_order(energy_tol);
      s=states.begin();
      for(int i=1;i<size();++i)
      {
        s_prev=s;
        ++s;
        energy_ord &= energy_order(*s_prev,*s);
      }
    }
    if(index_ord && energy_ord)
      order = energy_and_index_ordered;
    else if(index_ord)
      order = spoinfo::index_ordered;
    else if(energy_ord)
      order = spoinfo::energy_ordered;
    else
      order = spoinfo::unordered;
    count_degeneracies();
  }

  void SPOSetInfo::modify()
  {
    is_complete      = false;
    order            = no_order;
    indices_present  = false;
    energies_present = false;
    energy_tol       = no_energy_tol;
    index_min        = no_index;
    index_max        = no_index;
  }

  void SPOSetInfo::clear()
  {
    delete_iter(states.begin(),states.end());
    states.clear();
    modify();
  }

}
