
#ifndef QMCPLUSPLUS_COUNTING_FUNCTOR_H
#define QMCPLUSPLUS_COUNTING_FUNCTOR_H

#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

template <class T> class SigmoidFunctor
{
  public:

  SigmoidFunctor() {}
  bool put(xmlNodePtr cur);



};

template <class T> class GaussianFunctor
{
  public: 

  GaussianFunctor() {}
  bool put(xmlNodePtr cur);



};

}
#endif
