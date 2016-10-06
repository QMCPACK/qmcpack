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
    
    

#include <QMCWaveFunctions/SPOSetInputInfo.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/unit_conversion.h>
#include <limits>

namespace qmcplusplus
{
  typedef QMCTraits::RealType RealType;

  const int      inone = std::numeric_limits<int>::min();
  const RealType rnone = std::numeric_limits<RealType>::max();
  const std::string&  snone = "none";


  void SPOSetInputInfo::reset()
  {
    group        = 0;
    size         = inone;
    index_min    = inone;
    index_max    = inone;
    occ          = snone;
    ecut         = rnone;
    energy_min   = rnone;
    energy_max   = rnone;
    matching_tol = rnone;

    lowest_index   = inone;
    highest_index  = inone;
    lowest_energy  = rnone;
    highest_energy = rnone;
    
    has_size         = false;
    has_index_range  = false;
    has_occ          = false;
    has_ecut         = false;
    has_energy_range = false;
    has_indices      = false;
    has_energies     = false;
    has_index_info   = false;
    has_energy_info  = false;

    legacy_request   = false;

    all_indices_computed = false;
  }


  void SPOSetInputInfo::put(xmlNodePtr cur)
  {
    using namespace Units;

    reset();

    indices.clear();
    energies.clear();

    std::string units_in = snone;

    OhmmsAttributeSet attrib;
    attrib.add(group,     "group"     );
    attrib.add(size,      "size"      );
    attrib.add(index_min, "index_min" );
    attrib.add(index_max, "index_max" );
    attrib.add(occ,       "occ"       );
    attrib.add(ecut,      "ecut"      );
    attrib.add(energy_min,"energy_min");
    attrib.add(energy_max,"energy_max");
    attrib.add(units_in,  "units"     );
    attrib.put(cur);

    has_size         = size!=inone;
    has_index_range  = index_min!=inone && index_max!=inone;
    has_occ          = occ!=snone;
    has_ecut         = ecut!=rnone;
    has_energy_range = energy_min!=rnone && energy_max!=rnone;
    has_indices      = false;
    has_energies     = false;

    if(group<0)
      APP_ABORT("SPOSetInputInfo::put  group must be positive");

    if(has_ecut || has_energy_range)
    {
      report();
      if(units_in==snone)
        APP_ABORT("SPOSetInputInfo::put  ecut or energy range present, but units have not been provided");
      units eunits = energy_unit(units_in);
      if(has_ecut)
        ecut = convert(ecut,eunits,Ha);
      else if(has_energy_range)
      {
        energy_min = convert(energy_min,eunits,Ha);
        energy_max = convert(energy_max,eunits,Ha);
      }
    }

    xmlNodePtr element = cur->xmlChildrenNode;
    while(element!=NULL)
    {
      std::string ename((const char*)element->name);
      if(ename=="indices")
      {
        has_indices = true;
        putContent(indices,element);
      }
      else if(ename=="energies")
      {
        has_energies = true;
        units_in = snone;
        OhmmsAttributeSet attrib;
        attrib.add(matching_tol,"matching_tol");
        attrib.add(units_in,    "units"       );
        attrib.put(element);
        putContent(energies,element);
        if(units_in==snone)
          APP_ABORT("SPOSetInputInfo::put  energies present, but units have not been provided");
        units eunits = energy_unit(units_in);
        if(matching_tol==rnone)
          matching_tol = 1e-6;
        else
          matching_tol = convert(matching_tol,eunits,Ha);

        //convert_array(energies,eunits,Ha);
        std::vector<double> entmp;
        convert_array(entmp,eunits,Ha);

        sort(energies.begin(),energies.end());
      }
      element = element->next;
    }

    has_index_info  = has_size || has_index_range || has_occ || has_indices;
    has_energy_info = has_ecut || has_energy_range || has_energies; 

    legacy_request   = !(has_index_info && !has_size) && !has_energy_info;

    find_index_extrema();

    find_energy_extrema();
  }


  void SPOSetInputInfo::find_index_extrema()
  {
    if(has_index_info)
    {
      lowest_index  = std::numeric_limits<int>::max();
      highest_index = std::numeric_limits<int>::min();
      if(has_size)
      {
        lowest_index  = std::min(lowest_index,0);
        highest_index = std::max(highest_index,size-1);
      }
      if(has_index_range)
      {
        lowest_index  = std::min(lowest_index, index_min);
        highest_index = std::max(highest_index,index_max);
      }
      if(has_occ)
      {
        int imin=-1;
        int imax=-1;
        for(int i=0;i<occ.size();++i)
          if(occ[i]=='1')
          {
            if(imin==-1)
              imin = i;
            imax = i;
          }
        if(imin!=-1)
          lowest_index  = std::min(lowest_index,imin);
        if(imax!=-1)
          highest_index = std::max(highest_index,imax);
      }
      if(has_indices)
        for(int i=0;i<indices.size();++i)
        {
          int ind = indices[i];
          lowest_index  = std::min(lowest_index, ind);
          highest_index = std::max(highest_index,ind);
        }
    }
  }


  void SPOSetInputInfo::find_energy_extrema()
  {
    if(has_energy_info)
    {
      lowest_energy  =  rnone;
      highest_energy = -rnone;
      if(has_ecut)
      {
        lowest_energy  = std::min(lowest_energy,-rnone);
        highest_energy = std::max(highest_energy,ecut);
      }
      if(has_energy_range)
      {
        lowest_energy  = std::min(lowest_energy, energy_min);
        highest_energy = std::max(highest_energy,energy_max);
      }
      if(has_energies)
        for(int i=0;i<energies.size();++i)
        {
          RealType en = energies[i];
          lowest_energy  = std::min(lowest_energy, en);
          highest_energy = std::max(highest_energy,en);
        }
    }
  }


  void SPOSetInputInfo::report(const std::string& pad)
  {
    app_log()<<pad<<"SPOSetInput report"<< std::endl;
    app_log()<<pad<<"  has_size         = "<< has_size << std::endl;
    app_log()<<pad<<"  has_index_range  = "<< has_index_range << std::endl;
    app_log()<<pad<<"  has_occ          = "<< has_occ << std::endl;
    app_log()<<pad<<"  has_ecut         = "<< has_ecut << std::endl;
    app_log()<<pad<<"  has_energy_range = "<< has_energy_range << std::endl;
    app_log()<<pad<<"  has_indices      = "<< has_indices << std::endl;
    app_log()<<pad<<"  has_energies     = "<< has_energies << std::endl;
    app_log()<<pad<<"  group            = "<< group << std::endl;
    app_log()<<pad<<"  size             = "<< size << std::endl;
    app_log()<<pad<<"  index_min        = "<< index_min << std::endl;
    app_log()<<pad<<"  index_max        = "<< index_max << std::endl;
    app_log()<<pad<<"  occ              = "<< occ << std::endl;
    app_log()<<pad<<"  ecut             = "<< ecut << std::endl;
    app_log()<<pad<<"  energy_min       = "<< energy_min << std::endl;
    app_log()<<pad<<"  energy_max       = "<< energy_max << std::endl;
    app_log()<<pad<<"  # of indices     = "<<indices.size() << std::endl;
    app_log()<<pad<<"  indices          = \n    ";
    for(int i=0;i<indices.size();++i)
      app_log()<<indices[i]<<" ";
    app_log()<< std::endl;
    app_log()<<pad<<"  # of energies    = "<<energies.size() << std::endl;
    app_log()<<pad<<"  energies         = \n    ";
    for(int i=0;i<energies.size();++i)
      app_log()<<energies[i]<<" ";
    app_log()<< std::endl;
    app_log()<<pad<<"  matching_tol     = "<< matching_tol << std::endl;
    app_log()<<pad<<"  lowest_index     = "<<lowest_index<< std::endl;
    app_log()<<pad<<"  highest_index    = "<<highest_index<< std::endl;
    app_log()<<pad<<"  lowest_energy    = "<<lowest_energy<< std::endl;
    app_log()<<pad<<"  highest_energy   = "<<highest_energy<< std::endl;
    app_log()<<pad<<"end SPOSetInput report"<< std::endl;
    app_log().flush();
  }



  SPOSetInputInfo::indices_t& SPOSetInputInfo::get_indices(const std::vector<SPOSetInfo*>& states_vec)
  {
    if(!all_indices_computed)
    {
      all_indices.clear();

      if(group>=states_vec.size())
        APP_ABORT("SPOSetInputInfo::get_indices  orbital group index is out of range");

      const SPOSetInfo& states = *states_vec[group];

      // ensure that state info has been properly intialized
      bool energy_incomplete = states.partial() || !states.has_indices() || !states.has_energies();
      if(has_energy_info && energy_incomplete)
        APP_ABORT("SPOSetInputInfo::get_indices\n  energies requested for sposet\n  but state info is incomplete");

      occupations.clear();

      occupy_size();  
      occupy_index_range();
      occupy_occ();
      occupy_indices();
      occupy_ecut(states);
      occupy_energy_range(states);
      occupy_energies(states);

      sort(all_indices.begin(),all_indices.end());

      if(all_indices[all_indices.size()-1] >= states.size())
        APP_ABORT("SPOSetInputInfo::get_indices  requested state indices outside the range of states");

      if(all_indices.size()==0)
        APP_ABORT("SPOSetInputInfo::get_indices  no states matched request");

      all_indices_computed = true;
    }

    return all_indices;
  }

  
  void SPOSetInputInfo::occupy_size()
  {
    if(has_size)
    {
      indices_t ind;
      for(int i=0;i<size;++i)
        ind.push_back(i);
      occupy("size",ind);
    }
  }

  void SPOSetInputInfo::occupy_index_range()
  {
    if(has_index_range)
    {
      indices_t ind;
      for(int i=index_min;i<index_max;++i)
        ind.push_back(i);
      occupy("index_range",ind);
    }
  }

  void SPOSetInputInfo::occupy_indices()
  {
    if(has_indices)
      occupy("indices",indices);
  }

  void SPOSetInputInfo::occupy_occ()
  {
    if(has_occ)
    {
      indices_t ind;
      for(int i=0;i<occ.size();++i)
        if(occ[i]=='1')
          ind.push_back(i);
      occupy("occ",ind);
    }
  }

  void SPOSetInputInfo::occupy_ecut(const SPOSetInfo& states)
  {
    if(has_ecut)
    {
      indices_t ind;
      for(int i=0;i<states.size();++i)
      {
        const SPOInfo& state = *states[i];
        if(state.energy < ecut)
          ind.push_back(state.index);
      }
      occupy("ecut",ind);
    }
  }

  void SPOSetInputInfo::occupy_energy_range(const SPOSetInfo& states)
  {
    if(has_energy_range)
    {
      indices_t ind;
      for(int i=0;i<states.size();++i)
      {
        const SPOInfo& state = *states[i];
        if(state.energy<energy_max && state.energy>=energy_min)
          ind.push_back(state.index);
      }
      occupy("energy_range",ind);
    }
  }

  void SPOSetInputInfo::occupy_energies(const SPOSetInfo& states)
  {
    if(has_energies)
    {
      if(!states.energy_ordered())
        APP_ABORT("SPOSetInputInfo::load_indices(energies)\n  states are not energy ordered\n  this is a developer error");
      indices_t ind;
      int i=0;
      for(int n=0;n<energies.size();++n)
      {
        RealType e = energies[n];
        bool found = false;
        while(i<states.size())
        {
          const SPOInfo& state = *states[i];
          while(std::abs(e-state.energy)<matching_tol)
          {
            ind.push_back(state.index);
            i++;
            found = true;
          }
          if(found)
            break;
          else
            i++;
        }
        if(!found)
          APP_ABORT("SPOSetInputInfo::load_indices(energies)\n  energy eigenvalue not found");
      }
      occupy("energies",ind);
    }
  }


  void SPOSetInputInfo::occupy(const std::string& loc,const indices_t& ind)
  {
    int imin = std::numeric_limits<int>::max();
    int imax = std::numeric_limits<int>::min();
    for(int i=0;i<ind.size();++i)
    {
      int ival = ind[i];
      imin = std::min(imin,ival);
      imax = std::max(imax,ival);
    }
    if(imin<0)
      APP_ABORT("SPOSetInputInfo::occupy("+loc+")\n  indices are negative");
    for(int i=0;i<ind.size();++i)
    {
      int iocc = ind[i];
      if(iocc >= occupations.size())
      {
        int old_size = occupations.size();
        occupations.resize(iocc+1);
        for(int j=old_size;j<occupations.size();++j)
          occupations[j] = false;
      }
      if(occupations[iocc])
        APP_ABORT("SPOSetInputInfo::occupy("+loc+")\n  sposet request has overlapping index ranges");
      all_indices.push_back(iocc);
      occupations[iocc] = true;
    }
  }

}
