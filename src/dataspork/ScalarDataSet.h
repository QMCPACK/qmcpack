//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_SCALARDATAFILE_H
#define QMCPLUSPLUS_SCALARDATAFILE_H
#include <map>
///forward declation of ScalarObservable
class ScalarObservable;

/** container for an input file
*/
struct ScalarDataSet
{

  typedef std::map<std::string, ScalarObservable*> RawDataSetType;
  bool Grouped;
  int FirstRow;
  int LastRow;
  int NumRows;
  RawDataSetType Data;

  ScalarDataSet():Grouped(false) {}
  ScalarDataSet(const std::string& fname, int first, int last);
  ~ScalarDataSet();

  inline int rows() const
  {
    return NumRows;
  }

  /** Clone itself
   * @param nc number of columns
   * @return ScalarDataSet*
   */
  ScalarDataSet* clone(int nc=-1);

  void read(const std::string& fname, int first, int last);

  void accumulate(const RawDataSetType& raw);
  void group(const std::map<std::string,ScalarDataSet*>& dlist);

  void writeCollect(std::ostream& os, const std::string& cname);
  void writeObservable(std::ostream& os, const std::string& cname);
};

#endif
