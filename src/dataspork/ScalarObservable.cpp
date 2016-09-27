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
    
    


#include <vector>
#include <map>
#include <numeric>
#include <cmath>
#include <iostream>
#include <iomanip>
#include "ScalarObservable.h"

void
ScalarObservable::get_stat()
{
  get_average();
  get_corrtime();
}

ScalarObservable::value_type
ScalarObservable::get_average()
{
  if(Modified)
  {
    Average = std::accumulate(Self.begin(),Self.end(),0.0)
             /static_cast<value_type>(Self.size());
    //iterator first(Self.begin()),last(Self.end());
    //value_type s(0),s2(0);
    //while(first != last) {
    //  s2+= (*first)*(*first);
    //  s += (*first++);
    //}
    //Average = s/static_cast<value_type>(Self.size());
    //Variance = s2/static_cast<value_type>(Self.size())-Average*Average;
    Modified=false;
  }
  return Average;
}

ScalarObservable::value_type
ScalarObservable::get_corrtime()
{
  value_type avg=get_average(); //always call the average
  if(CorrTime<0)
  {
    //CoVarType res=get_covariance(this);
    CoVarType res=get_corrtime_dataspork(0,Self.size()-1);
    CorrTime=res.first;
    Variance=res.second;
  }
  return CorrTime;
}

ScalarObservable::value_type
ScalarObservable::get_corrtime(ScalarObservable* other)
{
  get_average();
  other->get_average();
  return get_covariance(other).first;
}

ScalarObservable::CoVarType
ScalarObservable::get_covariance(ScalarObservable* other)
{
  int NumData=size();
  int ic=0;
  value_type covar=get_correlation(other,ic);
  value_type idealcovar(covar),corrt(1.0);
  do
  {
    ic++;
    //for(ic=1; NumCuts*corrt >= ic; ic++){
    covar += get_correlation(other,ic);
    covar += get_correlation(other,-ic);
    corrt = (covar/idealcovar)*(value_type(NumData-1)/value_type(NumData-2*ic-1));
  }
  while(NumCuts*corrt>ic);
  NumEffData=static_cast<value_type>(NumData-2*(ic-1)-1);
  return CoVarType(corrt,covar/NumEffData);
  //return CoVarType(corrt,covar/static_cast<value_type>(NumData-2*(--ic)-1));
}

ScalarObservable::value_type
ScalarObservable::get_correlation(ScalarObservable* other, int ic)
{
  int NumData=size();
  iterator data_x(Self.begin());
  iterator data_y(other->begin()+ic);
  value_type corr(0.0);
  int istart(0),iend(NumData-ic);
  if(ic<0)
  {
    istart=-ic;
    iend=NumData;
    data_x=Self.begin()-ic;
    data_y=other->begin();
  }
  value_type avgI(Average);
  value_type avgJ(other->Average);
  for(int k=istart; k< iend; k++)
  {
    corr+=((*data_x-avgI)*(*data_y-avgJ));
    ++data_x;
    ++data_y;
  }
  return corr;
  //return corr/NumEffData;
}

ScalarObservable::value_type
ScalarObservable::get_blockerr(int bl)
{
  return get_blockerr(this,bl);
}

ScalarObservable::value_type
ScalarObservable::get_blockerr(ScalarObservable* other, int bl)
{
  int ndata=Self.size();
  iterator data_x(Self.begin());
  iterator data_y(other->begin());
  value_type avgI(Average);
  value_type avgJ(other->Average);
  value_type invbl=1.0/static_cast<value_type>(bl);
  value_type blockcov=0.0;
  int nbl(0);
  for (int k=(bl-1); k<ndata; k+=bl)
  {
    nbl++;
    double avebl_x(0.0),avebl_y(0.0);
    for(int kb=0; kb<bl; kb++)
    {
      avebl_x+=(*data_x++);
      avebl_y+=(*data_y++);
    }
    blockcov+=((avebl_x*invbl-avgI)*(avebl_y*invbl-avgJ));
  }
  return sqrt(blockcov/static_cast<value_type>(nbl*(nbl-1)));
}

void ScalarObservable::write(std::ostream& os)
{
  value_type sig=sqrt(Variance);
  value_type err=sig*sqrt (CorrTime/static_cast<value_type>(Self.size()));
  /*
  value_type err=get_blockerr(static_cast<int>(CorrTime+1.0));
    */
  os << "      <average>"<<Average << "</average>\n";
  streamsize p=os.precision();
  os.precision(1);
  os  << "      <error>"<<err<<"</error>\n";
  os.precision(p);
  os << "      <sigma>"<<sig<<"</sigma>\n"
     << "      <corrtime>"<<CorrTime<<"</corrtime>\n";
  std::cout <<setw(16) <<  Average
       << std::setw(16) << err
       << std::setw(16) << sig
       << std::setw(16) << CorrTime << std::endl;
}

/** Original dataspork method to evaluate the autocorrelation */
ScalarObservable::CoVarType
ScalarObservable::get_corrtime_dataspork(int imin, int imax)
{
  int n=imax-imin+1;
  int cutoff=n;
  value_type auto_min=1;
  value_type auto_max=1;
  value_type mean=Average;
  value_type var=0.0;
  iterator first(Self.begin()+imin),last(Self.begin()+imax+1);
  while(first != last)
  {
    var += (*first)*(*first);
    ++first;
  }
  var = var/static_cast<value_type>(n)-Average*Average;
  std::vector<value_type> autoCorr(cutoff,0.0);
  value_type corrtime;
  for(int j=0; j<cutoff; j++)
  {
    for(int i=imin; i+j<=imax; i++)
    {
      autoCorr[j]+=(Self[i]-mean)*(Self[i+j]-mean);
    }
    value_type dnj=static_cast<value_type>(n-j);
    autoCorr[j]/= (var*dnj);
    if(autoCorr[j]<sqrt(2.0/dnj))
    {
      if(cutoff>5*j)
        cutoff=5*j;
    }
    auto_max=std::max(autoCorr[j],auto_max);
    auto_min=std::min(autoCorr[j],auto_min);
  }
  int autoCorrlen=autoCorr.size();
  if(autoCorrlen > cutoff)
  {
    autoCorr.erase(autoCorr.begin()+cutoff,autoCorr.end());
  }
  autoCorrlen=autoCorr.size();
  if(autoCorrlen)
  {
    corrtime=autoCorr[0];
    for (int j=1; j*5< 2*cutoff; j++)
    {
      corrtime += 2.0*autoCorr[j];
    }
  }
  else
  {
    corrtime=1.0;
  }
  return CoVarType(corrtime,var);
}
