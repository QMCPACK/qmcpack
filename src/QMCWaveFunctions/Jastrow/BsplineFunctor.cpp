
#include "Configuration.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"

namespace qmcplusplus
{

template<>
const double  BsplineFunctor<double>::A[16] =
{
  -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
  3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
  -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
  1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0
};

template<>
const double BsplineFunctor<double>::dA[16] =
{
  0.0, -0.5,  1.0, -0.5,
  0.0,  1.5, -2.0,  0.0,
  0.0, -1.5,  1.0,  0.5,
  0.0,  0.5,  0.0,  0.0
};

template<>
const double BsplineFunctor<double>::d2A[16] =
{
  0.0, 0.0, -1.0,  1.0,
  0.0, 0.0,  3.0, -2.0,
  0.0, 0.0, -3.0,  1.0,
  0.0, 0.0,  1.0,  0.0
};



template<>
const double BsplineFunctor<double>::d3A[16] =
{
  0.0, 0.0, 0.0, -1.0,
  0.0, 0.0, 0.0,  3.0,
  0.0, 0.0, 0.0, -3.0,
  0.0, 0.0, 0.0,  1.0
};

//   template<>
//   const float  BsplineFunctor<float>::A[16] =
//     { -1.0/6.0,  3.0/6.0, -3.0/6.0, 1.0/6.0,
//        3.0/6.0, -6.0/6.0,  0.0/6.0, 4.0/6.0,
//       -3.0/6.0,  3.0/6.0,  3.0/6.0, 1.0/6.0,
//        1.0/6.0,  0.0/6.0,  0.0/6.0, 0.0/6.0 };

//   template<>
//   const float BsplineFunctor<float>::dA[16] =
//     {  0.0, -0.5,  1.0, -0.5,
//        0.0,  1.5, -2.0,  0.0,
//        0.0, -1.5,  1.0,  0.5,
//        0.0,  0.5,  0.0,  0.0 };

//   template<>
//   const float BsplineFunctor<float>::d2A[16] =
//     {  0.0, 0.0, -1.0,  1.0,
//        0.0, 0.0,  3.0, -2.0,
//        0.0, 0.0, -3.0,  1.0,
//        0.0, 0.0,  1.0,  0.0 };

}
