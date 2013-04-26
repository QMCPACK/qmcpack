#ifndef QMCPLUSPLUS_REPTILESTATISTICS_H
#define QMCPLUSPLUS_REPTILESTATISTICS_H

#include "ReptileEstimator.h"

namespace qmcplusplus
{
struct ReptileStatistics: public ReptileEstimator
{

  ReptileStatistics(MultiChain* polymer)
  {
    nvals=2;
    Values.resize(nvals);
  }

  ~ReptileStatistics() {}

  void put(xmlNodePtr cur) {}

  void put(xmlNodePtr cur, MCWalkerConfiguration& refWalker)
  {
    put(cur);
  }

  void evaluate(MultiChain::iterator first, MultiChain::iterator last, int ipsi)
  {
    int maxMade=max((*first)->stepmade,(*(last-1))->stepmade);
    int minMade=min((*first)->stepmade,(*(last-1))->stepmade);
    int maxtouch(0);
    while(first!=last)
    {
      if ((*first)->stepmade < minMade)
        minMade = (*first)->stepmade;
      if ((*first)->timesTouched > maxtouch)
        maxtouch = (*first)->timesTouched;
      first++;
    }
    Values[0]=maxtouch;
    Values[1]=maxMade-minMade;
  };

  void addNames(vector<string>& names)
  {
    myIndex=names.size();
    names.push_back("maxTouch");
    names.push_back("maxAge");
  };

  void setValues(vector<accumulator_type>& data, RealType weight)
  {
    data[myIndex](Values[0],1);
    data[myIndex+1](Values[1],1);
  };


private:
  int nvals;

};
}
#endif