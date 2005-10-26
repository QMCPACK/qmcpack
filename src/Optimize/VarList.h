//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_QMC_OPTIMIZE_VARREGISTRY_H
#define OHMMS_QMC_OPTIMIZE_VARREGISTRY_H

/**class VarRegistry
 *\brief register variables for optimization. 
 *
 For a set of variables, register the name, pointer and size.  
 *\note a set can consist of one variable or an array of variables
*/

template<class T>
struct VarRegistry {

  ///vector containing the names of all the sets of optimizable variables
  vector<string> Names;
  ///vector containing the pointers to all the sets of optimizable variables 
  vector<T*> Pointers;
  ///vector containing the sizes of the sets of optimizable variables 
  vector<int> Sizes;

  /** return the current size of variables */
  inline int size() const { 
    return Names.size();
  }

  /** Add a new optimizable variable set to the list of sets.
   * \param aname name of the optimizable set
   * \param ptr pointer to the set
   * \param size size of the set 
   */
  void add(const string& aname, T* ptr, int size=1) {
    vector<string>::iterator it = Names.begin();
    while(it != Names.end()) {
      if((*it) == aname) {
	return;
      }
      it++;
    }
    Names.push_back(aname);
    Pointers.push_back(ptr);
    Sizes.push_back(size);
  }

  ///return the index of aname
  int find(const string& aname) {
    int id = 0;
    vector<string>::iterator it = Names.begin();
    while(it != Names.end()) {
      if((*it) == aname) {
	return id;
      }
      it++;id++;
    }
    return -1;
  }

  void print(ostream& os) {
    os << "Optimizable variable list " << std::endl;
    for(int i=0; i<Names.size(); i++)
      os << Names[i] << " " << *(Pointers[i]) << std::endl;
  }
};
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
