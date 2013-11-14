//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_SPOSET_INPUT_INFO_H
#define QMCPLUSPLUS_SPOSET_INPUT_INFO_H

#include <QMCWaveFunctions/SPOSetInfo.h>


namespace qmcplusplus
{
  /** class to read state range information from sposet input
   *  
   *    typically just a temporary object to read and get info from xml
   *
   */
  struct SPOSetInputInfo
  {
    typedef QMCTraits::RealType RealType;
    typedef vector<int> indices_t;
    typedef vector<RealType> energies_t;

    int group;
    int size;
    int index_min;
    int index_max;
    string occ;
    indices_t indices;
    RealType ecut;
    RealType energy_min;
    RealType energy_max;
    energies_t energies;
    RealType matching_tol;

    bool has_size;
    bool has_index_range;
    bool has_occ;
    bool has_ecut;
    bool has_energy_range;
    bool has_indices;
    bool has_energies;

    bool has_index_info;
    bool has_energy_info;

    bool legacy_request;

    int lowest_index;
    int highest_index;
    RealType lowest_energy;
    RealType highest_energy;

    bool all_indices_computed;
    vector<bool> occupations;
    indices_t    all_indices;

    SPOSetInputInfo(xmlNodePtr cur) 
    { 
      put(cur);
    }
    
    SPOSetInputInfo()
    {
      reset();
    }

    ~SPOSetInputInfo() { }

    void reset();

    void put(xmlNodePtr cur);

    void report(const string& pad="");

    inline int min_index()
    {
      return lowest_index;
    }

    inline int max_index()
    {
      return highest_index;
    }

    inline RealType min_energy()
    {
      return lowest_energy;
    }

    inline RealType max_energy()
    {
      return highest_energy;
    }

    indices_t& get_indices(const vector<SPOSetInfo*>& states_vec);

    inline indices_t& get_indices(xmlNodePtr cur,const vector<SPOSetInfo*>& states_vec)
    {
      put(cur);
      return get_indices(states_vec);
    }


  private:
    void find_index_extrema();
    void find_energy_extrema();

    void occupy_size();
    void occupy_index_range();
    void occupy_occ();
    void occupy_indices();
    void occupy_ecut(const SPOSetInfo& states);
    void occupy_energy_range(const SPOSetInfo& states);
    void occupy_energies(const SPOSetInfo& states);

    void occupy(const string& loc,const indices_t& ind);

  };

}

#endif
