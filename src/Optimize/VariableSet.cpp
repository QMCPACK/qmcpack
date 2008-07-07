//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file VariableSet.cpp
 * @brief Define VariableSet functions.
 */
#include "Optimize/VariableSet.h"
#include <map>

namespace optimize
{
  VariableSet::VariableSet(variable_map_type& input):num_active_vars(0)
  {
    Index.resize(input.size(),-1);
    NameAndValue.resize(input.size());
    std::copy(input.begin(),input.end(),NameAndValue.begin());
    for(int i=0; i<Index.size(); ++i) Index[i]=i;
  }

  void VariableSet::clear()
  {
    num_active_vars=0;
    Index.clear();
    NameAndValue.clear();
  }

  /** insert name-value pairs  of this object to output
   * @param output parameters to be added
   */
  void VariableSet::insertTo(variable_map_type& output) const
  {
    for(int i=0; i<Index.size(); ++i)
    {
      if(Index[i]>-1)
      {
        output[NameAndValue[i].first]=NameAndValue[i].second;
      }
    }
  }

  void VariableSet::insertFrom(const VariableSet& input)
  {
    for(int i=0; i<input.size(); ++i)
    {
      iterator loc=find(input.name(i));
      if(loc==NameAndValue.end())
      {
        Index.push_back(input.Index[i]);
        NameAndValue.push_back(input.NameAndValue[i]);
      }
    }
  }

  void VariableSet::activate(const variable_map_type& selected)
  {
    //activate the variables 
    variable_map_type::const_iterator it(selected.begin()),it_end(selected.end());
    while(it != it_end)
    {
      iterator loc=find((*it++).first);
      if(loc != NameAndValue.end())
      {
        int i=loc-NameAndValue.begin();
        if(Index[i]<0)
          Index[i]=num_active_vars++;
      }
    }
  }

  void VariableSet::disable(const variable_map_type& selected)
  {
    variable_map_type::const_iterator it(selected.begin()),it_end(selected.end());
    while(it != it_end)
    {
      int loc=find((*it++).first)-NameAndValue.begin();
      if(loc<NameAndValue.size()) Index[loc]=-1;
    }
  }


  void VariableSet::removeInactive()
  {
    std::vector<int> valid(Index);
    std::vector<pair_type> acopy(NameAndValue);
    num_active_vars=0;
    Index.clear();
    NameAndValue.clear();
    for(int i=0; i<valid.size(); ++i)
    {
      if(valid[i]>-1)
      {
        Index.push_back(num_active_vars++);
        NameAndValue.push_back(acopy[i]);
      }
    }
  }

  void VariableSet::resetIndex()
  {
    num_active_vars=0;
    for(int i=0; i<Index.size(); ++i)
    {
      Index[i]=(Index[i]<0)?-1:num_active_vars++;
    }
  }

  void VariableSet::getIndex(const VariableSet& selected)
  {
    for(int i=0; i<NameAndValue.size(); ++i)
    {
      Index[i]=selected.getIndex(NameAndValue[i].first);
      if(Index[i]>=0) num_active_vars++; 
    }
  }

  void VariableSet::setDefaults(bool optimize_all)
  {
    for(int i=0; i<Index.size(); ++i) Index[i]=optimize_all?i:-1;
  }

  void VariableSet::print(std::ostream& os)
  {
    for(int i=0; i<NameAndValue.size(); ++i)
    {
      os <<NameAndValue[i].first << " " << NameAndValue[i].second << " ";
      if(Index[i]<0) 
        os << " OFF" << std::endl;
      else 
        os << " ON " << Index[i] << std::endl;
    }
  }

}
/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2550 $   $Date: 2008-03-26 15:17:43 -0500 (Wed, 26 Mar 2008) $
 * $Id: VarList.h 2550 2008-03-26 20:17:43Z jnkim $ 
 ***************************************************************************/
