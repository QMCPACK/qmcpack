#ifndef QMC_FINITE_SIZE_H
#define QMC_FINITE_SIZE_H

#include "QMCApp/QMCAppBase.h"
#include "QMCTools/QMCFiniteSize/SkParserBase.h"
#include "Particle/ParticleSetPool.h"
#include "LongRange/LRCoulombSingleton.h"
#include "einspline/bspline_structs.h"

namespace qmcplusplus
{
/** Class to handle FS corrections
 *
 * Implements finite size corrections from Holzmann et al., PRB (2016)
 * Currently implements Eqn. (30), using a long-rang break up of the 
 * Coulomb interaction and a spline representation of S(k). 
 * S(k) is obtained from SkParserBase
 */
class QMCFiniteSize : public QMCAppBase, QMCTraits
{
public:
  using LRHandlerType    = LRCoulombSingleton::LRHandlerType;
  using GridType         = LRCoulombSingleton::GridType;
  using RadFunctorType   = LRCoulombSingleton::RadFunctorType;
  using mRealType        = LRHandlerType::mRealType;
  using Grid_t           = SkParserBase::Grid_t;
  using RealType         = QMCTraits::RealType;
  using FullPrecRealType = QMCTraits::FullPrecRealType;
  using PosType          = QMCTraits::PosType;
  QMCFiniteSize();
  QMCFiniteSize(SkParserBase* skparser_i);
  ~QMCFiniteSize(){};


  inline void setSkParser(SkParserBase* skparser_i) { skparser = skparser_i; };
  bool validateXML() override;
  bool execute() override;

  void build_spherical_grid(IndexType mtheta, IndexType mphi);
  void getSkInfo(UBspline_3d_d* spline, std::vector<RealType>& symmatelem);
  UBspline_3d_d* getSkSpline(std::vector<RealType> sk, RealType limit = 1.0);
  RealType sphericalAvgSk(UBspline_3d_d* spline, RealType k);

  RealType integrate_spline(UBspline_1d_d* spline, RealType a, RealType b, IndexType N);
  UBspline_1d_d* spline_clamped(std::vector<RealType>& grid, std::vector<RealType>& vals, RealType lVal, RealType rVal);

  void initialize();
  void calcPotentialCorrection();
  void calcLeadingOrderCorrections();
  void summary();
  RealType calcPotentialDiscrete(std::vector<RealType> sk);
  RealType calcPotentialInt(std::vector<RealType> sk);

private:
  SkParserBase* skparser;
  ParticleSetPool ptclPool;
  RealType myRcut;
  RealType myConst;
  ParticleSet* P;
  RealType h; //this is for finite differencing.
  std::vector<PosType> sphericalgrid;
  GridType* myGrid;
  std::unique_ptr<LRHandlerType> AA;
  std::unique_ptr<RadFunctorType> rVs;
  bool processPWH(xmlNodePtr cur);
  void wfnPut(xmlNodePtr cur);
  void initBreakup();
  Grid_t gridx;
  Grid_t gridy;
  Grid_t gridz;
  void printSkRawSphAvg(const std::vector<RealType>& sk);
  void printSkSplineSphAvg(UBspline_3d_d* spline);
  KContainer Klist;
  std::vector<TinyVector<int, OHMMS_DIM>> kpts;
  std::vector<RealType> SK_raw;
  std::vector<RealType> SKerr_raw;
  std::vector<RealType> SK;
  std::vector<RealType> SKerr;
  IndexType mtheta;
  IndexType mphi;
  IndexType NumSamples;
  RealType Ne, Vol, rs, rho;
  RealType tlo, tloerr, vlo, vloerr, Vfs, Vfserr;
};
} // namespace qmcplusplus

#endif
