
#ifndef QMCPLUSPLUS_BSPLINE_CREATOR_H
#define QMCPLUSPLUS_BSPLINE_CREATOR_H
#include "qmc_common.h"
#include <Utilities/ProgressReportEngine.h>
#include "QMCWaveFunctions/EinsplineSetBuilder.h"
#include "QMCWaveFunctions/BsplineFactory/BsplineSet.h"
#include <fftw3.h>

#include <QMCWaveFunctions/einspline_helper.hpp>
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
    if(hybrid_rep)
      aReader= new SplineHybridAdoptorReader<HybridRealSoA<SplineR2RAdoptor<double,OHMMS_PRECISION>>,
    					     batching>(e);
    else 
      aReader= new SplineAdoptorReader<SplineR2RAdoptor<double,OHMMS_PRECISION>, batching>(e);
    return aReader;
  }
#ifdef QMC_COMPLEX
  static   BsplineReaderBase* createBsplineComplexDouble(EinsplineSetBuilder* e, bool hybrid_rep)
  {
    typedef OHMMS_PRECISION RealType;
    BsplineReaderBase* aReader=nullptr;

#if defined(QMC_COMPLEX)
    if(hybrid_rep)
      aReader= new SplineHybridAdoptorReader<HybridCplxSoA<SplineC2CAdoptor<double,RealType> >, batching >(e);
    else
      aReader= new SplineAdoptorReader<SplineC2CAdoptor<double,RealType>, batching >(e);
#else //QMC_COMPLEX
    if(hybrid_rep)
      aReader= new SplineHybridAdoptorReader<HybridCplxSoA<SplineC2RAdoptor<double,RealType> >, batching >(e);
    else
      aReader= new SplineAdoptorReader<SplineC2RAdoptor<double,RealType>, batching >(e);
#endif

    return aReader;
  }
#endif
};

}

#endif
