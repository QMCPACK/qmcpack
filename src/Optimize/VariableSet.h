//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_OPTIMIZE_VARIABLESET_H
#define QMCPLUSPLUS_OPTIMIZE_VARIABLESET_H
#include "config.h"
#include <map>
#include <vector>
#include <iostream>

namespace optimize
{
/** An enum useful for determing the type of parameter is being optimized.
*   knowing this in the opt routine can reduce the computational load.
*/
enum
{
  OTHER_P=0,
  LOGLINEAR_P, //B-spline Jastrows
  LOGLINEAR_K, //K space Jastrows
  LINEAR_P,    //Multi-determinant coefficients
  SPO_P,       //SPO set Parameters
  BACKFLOW_P   //Backflow parameters
};

/** class to handle a set of variables that can be modified during optimizations
 *
 * A serialized container of named variables.
 */
struct VariableSet
{
  typedef OHMMS_PRECISION                         real_type;
  typedef std::pair<std::string,real_type>        pair_type;
  typedef std::pair<std::string,int>              indx_pair_type;
  typedef std::vector<pair_type>::iterator        iterator;
  typedef std::vector<pair_type>::const_iterator  const_iterator;
  typedef std::vector<pair_type>::size_type       size_type;
  typedef std::map<std::string,real_type>         variable_map_type;

  ///number of active variables
  int num_active_vars;
  /** store locator of the named variable
   *
   * if(Index[i]  == -1), the named variable is not active
   */
  std::vector<int>        Index;
  std::vector<pair_type>  NameAndValue;
  std::vector<indx_pair_type>        ParameterType;
  std::vector<indx_pair_type>        Recompute;

  ///default constructor
  inline VariableSet():num_active_vars(0) {}
  ///constructor using map
//       VariableSet(variable_map_type& input);
  ///viturval destructor for safety
  virtual ~VariableSet() {}
  /** if any of Index value is not zero, return true
   */
  inline bool is_optimizable() const
  {
    return num_active_vars>0;
  }
  ///return the number of active variables
  inline int size_of_active() const
  {
    return num_active_vars;
  }
  ///return the first const_iterator
  inline const_iterator begin() const
  {
    return NameAndValue.begin();
  }
  ///return the last const_iterator
  inline const_iterator end() const
  {
    return NameAndValue.end();
  }
  ///return the first iterator
  inline iterator begin()
  {
    return NameAndValue.begin();
  }
  ///return the last iterator
  inline iterator end()
  {
    return NameAndValue.end();
  }
  ///return the size
  inline size_type size() const
  {
    return NameAndValue.size();
  }
  ///return the locator of the i-th Index
  inline int where(int i) const
  {
    return Index[i];
  }
  /** return the iterator of a named parameter
   * @param vname name of a parameter
   * @return the locator of vname
   *
   * If vname is not found among the Names, return NameAndValue.end()
   * so that ::end() member function can be used to validate the iterator.
   */
  inline iterator find(const std::string& vname)
  {
    iterator it(NameAndValue.begin());
    while(it != NameAndValue.end())
    {
      if((*it).first == vname)
        return it;
      ++it;
    }
    return NameAndValue.end();
  }

  /** return the Index vaule for the named parameter
   * @param vname name of the variable
   *
   * If vname is not found in this variables, return -1;
   */
  inline int getIndex(const std::string& vname) const
  {
    int loc=0;
    while(loc != NameAndValue.size())
    {
      if(NameAndValue[loc].first == vname)
        return Index[loc];
      ++loc;
    }
    return -1;
  }

  inline void insert(const std::string& vname, real_type v, bool enable=true, int type=OTHER_P)
  {
    iterator loc=find(vname);
    int ind_loc=loc-NameAndValue.begin();
    if(loc==NameAndValue.end()) //  && enable==true)
    {
      Index.push_back(ind_loc);
      NameAndValue.push_back(pair_type(vname,v));
      ParameterType.push_back(indx_pair_type(vname,type));
      Recompute.push_back(indx_pair_type(vname,1));
    }
    //disable it if enable == false
    if(!enable)
      Index[ind_loc]=-1;
  }

  inline void setParameterType(int type)
  {
    std::vector<indx_pair_type>::iterator PTit(ParameterType.begin()), PTend(ParameterType.end());
    while (PTit!=PTend)
    {
      (*PTit).second=type;
      PTit++;
    }
  }

  inline void getParameterTypeList(std::vector<int>& types)
  {
    std::vector<indx_pair_type>::iterator PTit(ParameterType.begin()), PTend(ParameterType.end());
    types.resize(PTend-PTit);
    std::vector<int>::iterator tit(types.begin());
    while (PTit!=PTend)
      (*tit++) = (*PTit++).second;
  }


  /** equivalent to std::map<std::string,T>[string] operator
   */
  inline real_type& operator[](const std::string& vname)
  {
    iterator loc=find(vname);
    if(loc==NameAndValue.end())
    {
      Index.push_back(-1);
      NameAndValue.push_back(pair_type(vname,0));
      ParameterType.push_back(indx_pair_type(vname,0));
      Recompute.push_back(indx_pair_type(vname,1));
      return NameAndValue.back().second;
    }
    return (*loc).second;
  }


  /** return the name of i-th variable
   * @param i index
   */
  inline std::string name(int i) const
  {
    return NameAndValue[i].first;
  }

  /** return the i-th value
   * @param i index
   */
  inline real_type operator[](int i) const
  {
    return NameAndValue[i].second;
  }

  /** assign the i-th value
   * @param i index
   */
  inline real_type& operator[](int i)
  {
    return NameAndValue[i].second;
  }

  /** get the i-th parameter's type
  * @param i index
  */
  inline int getType(int i)
  {
    return ParameterType[i].second;
  }

  inline bool recompute(int i) const
  {
    return (Recompute[i].second==1);
  }

  inline int& recompute(int i)
  {
    return Recompute[i].second;
  }

  inline void setComputed()
  {
    for(int i=0; i<Recompute.size(); i++)
    {
      if      (ParameterType[i].second==LOGLINEAR_P)
        Recompute[i].second=0;
      else
        if (ParameterType[i].second==LOGLINEAR_K)
          Recompute[i].second=0;
        else
          Recompute[i].second=1;
    }
  }

  inline void setRecompute()
  {
    for(int i=0; i<Recompute.size(); i++)
      Recompute[i].second=1;
  }

  /** clear the variable set
   *
   * Remove all the data.
   */
  void clear();

  /** insert local variables to output
   */
//       void insertTo(variable_map_type& output) const;

  /** insert a VariableSet to the list
   * @param input varaibles
   */
  void insertFrom(const VariableSet& input);

  /** sum together the values of the optimizable parameter values in
   *  two VariableSet objects, and set this object's values to equal them.
   *  @param first set of input variables
   *  @param second set of input variables
   */
  void insertFromSum(const VariableSet& input_1, const VariableSet& input_2);

  /** take the difference (input_1-input_2) of values of the optimizable
   *  parameter values in two VariableSet objects, and set this object's
   *  values to equal them.
   *  @param first set of input variables
   *  @param second set of input variables
   */
  void insertFromDiff(const VariableSet& input_1, const VariableSet& input_2);

  /** activate variables for optimization
   * @param first iterator of the first name
   * @param last iterator of the last name
   * @param reindex if true, Index is updated
   *
   * The status of a variable that is not included in the [first,last)
   * remains the same.
   */
  template<typename ForwardIterator>
  void activate(ForwardIterator first, ForwardIterator last, bool reindex)
  {
    while(first != last)
    {
      iterator loc=find(*first++);
      if(loc != NameAndValue.end())
      {
        int i=loc-NameAndValue.begin();
        if(Index[i]<0)
          Index[i]=num_active_vars++;
      }
    }
    if(reindex)
    {
      removeInactive();
      resetIndex();
    }
  }

  /** make the selected variables active
   * @param selected input variables that are set to be varied
   */
  void activate(const variable_map_type& selected);


  /** deactivate variables for optimization
   * @param first iterator of the first name
   * @param last iterator of the last name
   * @param reindex if true, the variales are removed and Index is updated
   */
  template<typename ForwardIterator>
  void disable(ForwardIterator first, ForwardIterator last, bool reindex)
  {
    while(first != last)
    {
      int loc=find(*first++)-NameAndValue.begin();
      if(loc<NameAndValue.size())
        Index[loc]=-1;
    }
    if(reindex)
    {
      removeInactive();
      resetIndex();
    }
  }

  ///** make the selected variables active
  // * @param selected input variables that are set to be varied
  // */
  //void activate(const std::vector<std::string>& selected, bool reindex);

  /** exlcude variables
   * @param selected name-value pairs that should be dropped from the set
   */
  void disable(const variable_map_type& selected);

  /** reset Index
   */
  void resetIndex();
  /** remove inactive variables and trim the internal data
   */
  void removeInactive();

  /** set the index table of this VariableSet
   * @param selected input variables
   *
   * This VariableSet is a subset of selected.
   */
  void getIndex(const VariableSet& selected);

  /** set default Indices
   * @param optimize_all if true, all the variables are active
   */
  void setDefaults(bool optimize_all);

  void print(std::ostream& os);
};
}

#endif
