#ifndef QMC_FINITE_SIZE_H
#define QMC_FINITE_SIZE_H

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Particle/ParticleSetPool.h"
#include "LongRange/LRCoulombSingleton.h"
#include "einspline/bspline_structs.h"
#include "einspline/nubspline_structs.h"

namespace qmcplusplus
{
class QMCFiniteSize : public QMCAppBase, QMCTraits
{
public:
  typedef LRCoulombSingleton::LRHandlerType LRHandlerType;
  typedef LRCoulombSingleton::GridType GridType;
  typedef LRCoulombSingleton::RadFunctorType RadFunctorType;
  typedef LRHandlerType::mRealType mRealType;
  typedef SkParserBase::Grid_t Grid_t;
  typedef QMCTraits::RealType RealType;
  typedef QMCTraits::FullPrecRealType FullPrecRealType;
  typedef QMCTraits::PosType PosType;
  QMCFiniteSize();
  QMCFiniteSize(SkParserBase* skparser_i);
  ~QMCFiniteSize(){};


  inline void setSkParser(SkParserBase* skparser_i) { skparser = skparser_i; };
  bool validateXML();
  bool execute();

  void build_spherical_grid(IndexType mtheta, IndexType mphi);
  void getSkInfo(UBspline_3d_d* spline, vector<RealType>& symmatelem);
  UBspline_3d_d* getSkSpline(vector<RealType> sk, RealType limit = 1.0);
  RealType sphericalAvgSk(UBspline_3d_d* spline, RealType k);

  RealType integrate_spline(NUBspline_1d_d* spline, RealType a, RealType b, IndexType N);
  NUBspline_1d_d* spline_clamped(vector<RealType>& grid, vector<RealType>& vals, RealType lVal, RealType rVal);

  void initialize();
  void calcPotentialCorrection();
  void calcLeadingOrderCorrections();
  void summary();
  RealType calcPotentialDiscrete(vector<RealType> sk);
  RealType calcPotentialInt(vector<RealType> sk);

private:
  SkParserBase* skparser;
  ParticleSetPool ptclPool;
  RealType myRcut;
  RealType myConst;
  ParticleSet* P;
  RealType h; //this is for finite differencing.
  vector<PosType> sphericalgrid;
  GridType* myGrid;
  LRHandlerType* AA;
  RadFunctorType* rVs;
  bool processPWH(xmlNodePtr cur);
  void wfnPut(xmlNodePtr cur);
  void initBreakup();
  Grid_t gridx;
  Grid_t gridy;
  Grid_t gridz;
  void printSkRawSphAvg(const vector<RealType>& sk);
  void printSkSplineSphAvg(UBspline_3d_d* spline);
  KContainer Klist;
  vector<TinyVector<int, OHMMS_DIM>> kpts;
  vector<RealType> SK_raw;
  vector<RealType> SKerr_raw;
  vector<RealType> SK;
  vector<RealType> SKerr;
  IndexType mtheta;
  IndexType mphi;
  IndexType NumSamples;
  RealType Ne, Vol, rs, rho;
  RealType tlo, tloerr, vlo, vloerr, Vfs, Vfserr;
};
} // namespace qmcplusplus

#endif
