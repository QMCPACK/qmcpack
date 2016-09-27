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
    
    


#ifndef QMCPLUSPLUS_DATASETTYPE_H
#define QMCPLUSPLUS_DATASETTYPE_H

struct ScalarObservable
{

  typedef double                  value_type;
  typedef std::vector<value_type> Base_t;
  typedef Base_t::iterator        iterator;
  typedef Base_t::const_iterator  const_iterator;
  typedef Base_t::size_type       size_type;

  /** typedef of the return type of get_covariance
   *
   * first : correlation time
   * second : covariance
   */
  typedef std::pair<value_type,value_type> CoVarType;

  bool Modified;
  int NumCuts;
  Base_t Self;
  value_type Average;
  value_type Variance;
  value_type CorrTime;
  value_type NumEffData;

  ScalarObservable(int ncut=1):
    Modified(true), NumCuts(ncut), Average(0.0),Variance(-1),CorrTime(-1) {}

  template<class IT>
  inline void copy(IT first, IT last)
  {
    Self.clear();
    Self.insert(Self.end(), first, last);
    Modified=true;
  }

  template<class IT>
  inline void append(IT first, IT last)
  {
    Self.insert(Self.end(), first, last);
    Modified=true;
  }


  inline void push_back(value_type x)
  {
    Self.push_back(x);
  }

  inline void resize(int n, value_type val=0.0)
  {
    Self.resize(n,val);
  }
  inline iterator begin()
  {
    return Self.begin();
  }
  inline iterator end()
  {
    return Self.end();
  }
  inline const_iterator begin() const
  {
    return Self.begin();
  }
  inline const_iterator end() const
  {
    return Self.end();
  }
  inline size_type size() const
  {
    return Self.size();
  }


  /** get statistic data of this observable
   *
   * Equivalent to get_corrtime() and reserved for future use
   */
  void get_stat();

  /** get the average and variance
   * @return average
   */
  value_type get_average();

  /** error corrected by a blocking
   * @param bl block size
   */
  value_type get_blockerr(int bl);

  /** error corrected by a blocking
   * @param other ScalarObservable* whose covariance is evaluated
   * @param bl block size
   * @param return the error
   */
  value_type get_blockerr(ScalarObservable* other, int bl);

  /** correlation time of Self
   */
  value_type get_corrtime();

  /** correlation time of Self using the original algorithm of dataspork
   * @param imin first index
   * @param imax last index
   *
   * Correlation evaluated in [imin,imax]. Note that both ends are closed.
   */
  CoVarType get_corrtime_dataspork(int imin, int imax);

  /** correlation time of (Self,*other)
   * @param other ScalarObservable*
   */
  value_type get_corrtime(ScalarObservable* other);

  /** covariance of (Self,*other)
   * @param other ScalarObservable*
   * @return std::pair<doulbe,double>
   *
   * CoVarType::first correlation time
   * CoVarType::second covariance
   */
  CoVarType get_covariance(ScalarObservable* other);

  /** correlation of (Self,*other+ic)
   * @param other ScalarObservable*
   * @param ic offset for the correlation function
   */
  value_type get_correlation(ScalarObservable* other, int ic);

  void write(std::ostream& os);
};

#endif
