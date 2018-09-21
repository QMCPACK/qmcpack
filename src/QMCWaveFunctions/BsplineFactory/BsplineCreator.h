
#ifndef QMCPLUSPLUS_BSPLINE_CREATOR_H
#define QMCPLUSPLUS_BSPLINE_CREATOR_H

#include "QMCWaveFunctions/BsplineFactory/HybridRealAdoptor.h"
#include "QMCWaveFunctions/BsplineFactory/SplineR2RAdoptor.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineReaderBase.h"
#include "QMCWaveFunctions/BsplineFactory/SplineAdoptorReaderP.h"
#include "QMCWaveFunctions/BsplineFactory/SplineHybridAdoptorReaderP.h"

namespace qmcplusplus
{
template<Batching batching = Batching::SINGLE>
class BsplineCreator
{
public:
  static BsplineReaderBase* createBsplineRealDouble(EinsplineSetBuilder* e, bool hybrid_rep)
  {
    BsplineReaderBase* aReader=nullptr;
    // if(hybrid_rep)
    //   aReader= new SplineHybridAdoptorReader<HybridRealSoA<SplineR2RAdoptor<double,OHMMS_PRECISION>>,
    // 					     batching>(e);
    // else
    aReader= new SplineAdoptorReader<SplineR2RAdoptor<double,OHMMS_PRECISION>, batching>(e);
    return aReader;
  }
};

}

#endif
