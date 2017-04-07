//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef OHMMS_DISTRIBUTEDINDEX_H
#define OHMMS_DISTRIBUTEDINDEX_H
#include "Utilities/UtilityFunctions.h"
#include <vector>
#include <list>
#include <map>
/**@file DistributedIndex.h
 *@brief Declaration of class DistributedIndex
 *
 */
/**@class DistributedIndex
 *@brief Container and manager of distributed data
 *
 *DistributedIndex are implemented by using std::map and std::list for
 *the container of distributed indexes. The choice is made at compiler
 *time by -DDISTRIBUTEDINDEX_MAP. The performance of the two
 *implementations can vary with respect to the compilers, STL
 *implementations, and the problem (size, homogeneity). List will work
 *in general.  However, for large-scale systems that are very
 *inhomogeneous, using map implementation can perform better.
 */
#ifdef DISTRIBUTEDINDEX_MAP

struct DistributedIndex
{

  typedef std::list<int> List_t;
  typedef std::list<int>::iterator iterator;
  typedef std::list<int>::const_iterator const_iterator;

  std::map<int,List_t*> I;
  std::vector<int> M;
  int NumData;

  inline DistributedIndex():NumData(0) { }
  inline ~DistributedIndex()
  {
    std::map<int,List_t*>::iterator it=I.begin();
    while(it != I.end())
    {
      delete (*it).second;
      ++it;
    }
  }

  void create(int npart)
  {
    ///first clean up
    std::map<int,List_t*>::iterator it=I.begin();
    while(it != I.end())
    {
      delete (*it).second;
      ++it;
    }
    for(int i=0; i<npart; i++)
    {
      I[i] = new std::list<int>;
    }
  }

  template <class _InputIterator>
  inline
  void distribute(_InputIterator first, _InputIterator last)
  {
    int n=0,i=0;
    _InputIterator first_s = first;
    while(first_s != last)
    {
      n++;
      first_s++;
    }
    M.resize(n+1);
    M[0] = 0;
    while(first != last)
    {
      M[i+1]=M[i]+(*first);
      first++;
    }
    NumData = M[n];
  }

  inline void distribute(int ntot, int npart)
  {
    FairDivide(ntot,npart,M);
    NumData = M[npart];
  }

  inline int size() const
  {
    return NumData;
  }
  inline int capacity() const
  {
    if(I.empty())
    {
      return M[1];
    }
    else
    {
      std::map<int,List_t*>::const_iterator it=I.begin();
      std::list<int>::size_type n=0;
      while(it != I.end())
      {
        n = std::max(n,(*it).second->size());
        ++it;
      }
      return int(n);
    }
  }

  inline int size(int i) const
  {
    return M[i+1]-M[i];
  }
  inline int first_data(int i) const
  {
    return M[i];
  }
  inline int last_data(int i) const
  {
    return M[i];
  }
  inline iterator begin(int i)
  {
    return I[i]->begin();
  }
  inline iterator end(int i)
  {
    return I[i]->end();
  }
  inline const_iterator begin(int i)
  {
    return I[i]->begin();
  }
  inline const_iterator end(int i)
  {
    return I[i]->end();
  }

  inline void clear()
  {
    std::map<int,List_t*>::iterator it=I.begin();
    while(it != I.end())
    {
      (*it).second->clear();
      ++it;
    }
    NumData = 0;
  }

  inline void add(int i, int d)
  {
    I[i]->push_back(d);
    NumData++;
  }

};
#else

struct DistributedIndex
{

  typedef std::list<int> List_t;
  typedef std::list<int>::iterator iterator;
  typedef std::list<int>::const_iterator const_iterator;

  /**container of data lists
   *
   *A dataset is distributed over a number of partitions or groups.
   *The objects that constitute a dataset can be anything as far as an
   *object in the dataset can be accessed by an index (integer type).
   *
   *The memebers of each group can be accessed by usual STL::list
   *iterators and member functions.  For instance, I[0]->begin() will
   *return the iterator for the first element of the 0th group.
   */
  std::vector<List_t*> I;

  ///the offsets to retrieve the size of a group effectively : M[i+1]-M[i] = I[i].size()
  std::vector<int> M;

  ///the total number of data: NumData = M.back()
  int NumData;

  ///default constructor
  inline DistributedIndex():NumData(0) { }

  inline DistributedIndex(const DistributedIndex& old):
    M(old.M),NumData(old.NumData)
  {
    for(int i=0; i<old.I.size(); ++i)
      I.push_back(new List_t(*I[i]));
  }

  ///destructor
  inline ~DistributedIndex()
  {
    for(int i=0; i<I.size(); i++)
      delete I[i];
  }

  /** create lists to store data
   *@param npart the number of lists
   *
   *A new list is added to the existing lists.
   */
  void create(int npart)
  {
    int dn = npart-I.size();
    while(dn)
    {
      I.push_back(new std::list<int>);
      --dn;
    }
  }


  /** distribute the dataset
   *@param ntot the size of the dataset to be partitioned
   *@param npart the number of groups or partitions
   *
   *Partition ntot data over npart groups so that the size of each
   *group can only vary by 1 using FairDivide(int,int,std::vector<int>&)
   *(see UtilityFunctions.h).
   */
  inline void distribute(int ntot, int npart)
  {
    FairDivide(ntot,npart,M);
    NumData = M[npart];
  }

  /** distribute the dataset
   *@param first the iterator of the first group
   *@param last the iterator of the last group
   *
   *The number of groups are determined by last-first.
   *The containers to store data lists and aux counters are resized.
   */
  template <class _InputIterator>
  inline
  void distribute(_InputIterator first, _InputIterator last)
  {
    int n=0,i=0;
    _InputIterator first_s = first;
    while(first_s != last)
    {
      n++;
      first_s++;
    }
    M.resize(n+1);
    M[0] = 0;
    while(first != last)
    {
      M[i+1]=M[i]+(*first);
      first++;
      i++;
    }
    NumData = M[n];
  }

  ///return the capacity of containers
  inline int capacity() const
  {
    if(I.empty())
    {
      return M[1];
    }
    else
    {
      std::list<int>::size_type n=0;
      for(int i=0; i<I.size(); i++)
      {
        n = std::max(n,I[i]->size());
      }
      return int(n);
    }
  }

  ///return the size of groups
  inline int getNumDataSets() const
  {
    return I.size();
  }

  ///return the size of the data
  inline int getNumData() const
  {
    return NumData;
  }

  /**return the size of the group i
   *@param i the group index
   */
  inline int getNumData(int i) const
  {
    return M[i+1]-M[i];
  }

  /**return the first data of the group i
   *@param i the group index
   */
  inline int first_data(int i) const
  {
    return M[i];
  }

  /**return the last data of the group i
   *@param i the group index
   */
  inline int last_data(int i) const
  {
    return M[i+1];
  }

  ///return the size of the group i
  inline int size(int i) const
  {
    return I[i]->size();
  }
  ///return the starting iterator of the group i
  inline iterator begin(int i)
  {
    return I[i]->begin();
  }
  ///return the ending iterator of the group i
  inline iterator end(int i)
  {
    return I[i]->end();
  }
  ///return the starting const_iterator of the group i
  inline const_iterator begin(int i) const
  {
    return I[i]->begin();
  }
  ///return the ending const_iterator of the group i
  inline const_iterator end(int i) const
  {
    return I[i]->end();
  }

  ///print groupings for test/debug
  inline void print(std::ostream& os) const
  {
    os  << "The number of data sets " << I.size() << std::endl;
    for(int ig=0; ig<I.size(); ig++)
    {
      os << "Data set " << ig << std::endl;
      copy(I[ig]->begin(), I[ig]->end(),std::ostream_iterator<int>(os, " "));
      os << std::endl;
    }
  }

  /**clear the list for each group.
   *
   *The groups themselves are not affected but the data that belong
   *to each group are removed so that a new list can be built.
   */
  inline void clear()
  {
    for(int i=0; i<I.size(); i++)
      I[i]->clear();
    NumData = 0;
  }

  /**add a dat to a list
   *@param i the group index
   *@param d the data index to be added to the group i
   */
  inline void add(int i, int d)
  {
    I[i]->push_back(d);
    NumData++;
  }

};
#endif
#endif
