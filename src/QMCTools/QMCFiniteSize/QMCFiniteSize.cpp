#include "QMCFiniteSize.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/WaveFunctionComponentBuilder.h"
#include <iostream>
#include <cmath>
#include "Configuration.h"
#include "einspline/bspline_eval_d.h"
#include "QMCTools/QMCFiniteSize/FSUtilities.h"
#include "Utilities/RandomGenerator.h"

namespace qmcplusplus
{
QMCFiniteSize::QMCFiniteSize()
    : skparser(NULL), ptclPool(NULL), myRcut(0.0), myConst(0.0), P(NULL), h(0.0), sphericalgrid(0), myGrid(NULL)
{
  IndexType mtheta = 80;
  IndexType mphi   = 80;
  app_log() << "Building spherical grid. n_theta x n_phi = " << mtheta << " x " << mphi << std::endl;
  build_spherical_grid(mtheta, mphi);
  h = 0.1;
}

QMCFiniteSize::QMCFiniteSize(SkParserBase* skparser_i)
    : skparser(skparser_i), ptclPool(NULL), myRcut(0.0), myConst(0.0), P(NULL), h(0.0), sphericalgrid(0), myGrid(NULL)
{
  mtheta     = 80;
  mphi       = 80;
  h          = 0.1;
  NumSamples = 1000;
  build_spherical_grid(mtheta, mphi);
}

void QMCFiniteSize::build_spherical_grid(IndexType mtheta, IndexType mphi)
{
  //Spherical grid from https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
  RealType alpha = 4.0 * M_PI / (mtheta * mphi);
  RealType d     = std::sqrt(alpha);
  RealType Mt    = int(std::round(M_PI / d));
  RealType Dt    = M_PI / Mt;
  RealType Dp    = alpha / Dt;
  int count      = 0;
  for (int m = 0; m < Mt; m++)
  {
    RealType theta = M_PI * (m + 0.5) / Mt;
    RealType Mp    = int(std::round(2 * M_PI * std::sin(theta) / Dp));
    for (int n = 0; n < Mp; n++)
    {
      IndexType gindex = m * mtheta + n;
      RealType phi     = 2 * M_PI * n / Mp;
      PosType tmp;
      tmp[0] = std::sin(theta) * std::cos(phi);
      tmp[1] = std::sin(theta) * std::sin(phi);
      tmp[2] = std::cos(theta);
      sphericalgrid.push_back(tmp);
    }
  }
}

bool QMCFiniteSize::validateXML()
{
  xmlXPathContextPtr m_context = xml_doc_stack_.top()->getXPathContext();
  xmlNodePtr cur               = xml_doc_stack_.top()->getRoot()->children;

  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "particleset")
    {
      ptclPool.put(cur);
    }
    else if (cname == "wavefunction")
    {
      wfnPut(cur);
    }
    else if (cname == "include")
    {
      //file is provided
      const xmlChar* a = xmlGetProp(cur, (const xmlChar*)"href");
      if (a)
      {
        pushDocument((const char*)a);
        processPWH(xml_doc_stack_.top()->getRoot());
        popDocument();
      }
    }
    else if (cname == "qmcsystem")
    {
      processPWH(cur);
    }
    else {}
    cur = cur->next;
  }

  app_log() << "=========================================================\n";
  app_log() << " Summary of QMC systems \n";
  app_log() << "=========================================================\n";
  ptclPool.get(app_log());
  return true;
}


void QMCFiniteSize::wfnPut(xmlNodePtr cur)
{
  std::string id("psi0"), target("e"), role("extra");
  OhmmsAttributeSet pAttrib;
  pAttrib.add(id, "id");
  pAttrib.add(id, "name");
  pAttrib.add(target, "target");
  pAttrib.add(target, "ref");
  pAttrib.add(role, "role");
  pAttrib.put(cur);
  ParticleSet* qp = ptclPool.getParticleSet(target);

  if (qp == nullptr)
    throw std::runtime_error("target particle set named '" + target + "' not found");
}

bool QMCFiniteSize::processPWH(xmlNodePtr cur)
{
  //return true and will be ignored
  if (cur == NULL)
    return true;
  bool inputnode = true;
  //save the root to grep @tilematrix
  xmlNodePtr cur_root = cur;
  cur                 = cur->children;
  while (cur != NULL)
  {
    std::string cname((const char*)cur->name);
    if (cname == "simulationcell")
      ptclPool.readSimulationCellXML(cur);
    else if (cname == "particleset")
      ptclPool.put(cur);
    else if (cname == "wavefunction")
      wfnPut(cur);

    cur = cur->next;
  }
  return inputnode;
}

void QMCFiniteSize::initBreakup()
{
  app_log() << "=========================================================\n";
  app_log() << " Initializing Long Range Breakup (Esler) \n";
  app_log() << "=========================================================\n";
  P      = ptclPool.getParticleSet("e");
  AA     = LRCoulombSingleton::getHandler(*P);
  myRcut = AA->get_rc();
  if (rVs == nullptr)
  {
    rVs = LRCoulombSingleton::createSpline4RbyVs(AA.get(), myRcut, myGrid);
  }
}

UBspline_3d_d* QMCFiniteSize::getSkSpline(std::vector<RealType> sk, RealType limit)
{
  //get the einspline grids.
  Ugrid esgridx = gridx.einspline_grid();
  Ugrid esgridy = gridy.einspline_grid();
  Ugrid esgridz = gridz.einspline_grid();

  //setup the einspline boundary conditions.
  BCtype_d bcx;
  BCtype_d bcy;
  BCtype_d bcz;

  //This piece iterates through S(k) and sets
  //pieces beyond the k-cutoff equal to 1.
  //A violent approximation if S(k) is not converged, but
  //better than S(k)=0.
  double kc     = AA->get_kc();
  double kcutsq = kc * kc;

  for (int i = int(gridx.lower_bound), skindex = 0; i <= int(gridx.upper_bound); i++)
    for (int j = int(gridy.lower_bound); j <= int(gridy.upper_bound); j++)
      for (int k = int(gridz.lower_bound); k <= int(gridz.upper_bound); k++)
      {
        PosType v;
        v[0]         = i;
        v[1]         = j;
        v[2]         = k;
        RealType ksq = P->getLattice().ksq(v);

        if (ksq > kcutsq)
          sk[skindex] = limit;
        skindex++;
      }
  //No particular BC's on the edge of S(k).

  bcx.lCode = NATURAL;
  bcx.rCode = NATURAL;
  bcx.lVal  = 1.0;
  bcx.rVal  = 1.0;

  bcy.lCode = NATURAL;
  bcy.rCode = NATURAL;
  bcy.lVal  = 1.0;
  bcy.rVal  = 1.0;

  bcz.lCode = NATURAL;
  bcz.rCode = NATURAL;
  bcz.lVal  = 1.0;
  bcz.rVal  = 1.0;

  //hack for QMC_MIXED_PRECISION to interface to UBspline_3d_d
  std::vector<FullPrecRealType> sk_fp(sk.begin(), sk.end());
  UBspline_3d_d* spline = create_UBspline_3d_d(esgridx, esgridy, esgridz, bcx, bcy, bcz, sk_fp.data());

  return spline;
}

void QMCFiniteSize::getSkInfo(UBspline_3d_d* spline, std::vector<RealType>& symmatelem)
{
  symmatelem.resize(6);
  FullPrecRealType sx(0), sy(0), sz(0), sxy(0), sxz(0), syz(0);
  RealType h2 = h * h;

  PosType disp;
  PosType disp_lat;

  disp[0]  = h;
  disp[1]  = 0;
  disp[2]  = 0;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &sx);

  disp[0]  = 0;
  disp[1]  = h;
  disp[2]  = 0;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &sy);

  disp[0]  = 0;
  disp[1]  = 0;
  disp[2]  = h;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &sz);

  disp[0]  = h;
  disp[1]  = h;
  disp[2]  = 0;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &sxy);

  disp[0]  = h;
  disp[1]  = 0;
  disp[2]  = h;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &sxz);

  disp[0]  = 0;
  disp[1]  = h;
  disp[2]  = h;
  disp_lat = P->getLattice().k_unit(disp);
  eval_UBspline_3d_d(spline, disp_lat[0], disp_lat[1], disp_lat[2], &syz);

  symmatelem[0] = RealType(sx) / h2;
  symmatelem[1] = RealType(sy) / h2;
  symmatelem[2] = RealType(sz) / h2;
  symmatelem[3] = 0.5 * RealType(sxy - sx - sy) / h2;
  symmatelem[4] = 0.5 * RealType(sxz - sx - sz) / h2;
  symmatelem[5] = 0.5 * RealType(syz - sy - sz) / h2;
}

QMCFiniteSize::RealType QMCFiniteSize::sphericalAvgSk(UBspline_3d_d* spline, RealType k)
{
  RealType sum         = 0.0;
  FullPrecRealType val = 0.0;
  PosType kvec(0);
  IndexType ngrid = sphericalgrid.size();
  for (IndexType i = 0; i < ngrid; i++)
  {
    kvec     = P->getLattice().k_unit(k * sphericalgrid[i]); // to reduced coordinates
    bool inx = true;
    bool iny = true;
    bool inz = true;
    if (kvec[0] <= gridx.lower_bound || kvec[0] >= gridx.upper_bound)
      inx = false;
    if (kvec[1] <= gridy.lower_bound || kvec[1] >= gridy.upper_bound)
      iny = false;
    if (kvec[2] <= gridz.lower_bound || kvec[2] >= gridz.upper_bound)
      inz = false;
    if (!(inx & iny & inz))
      sum += 1;
    else
    {
      eval_UBspline_3d_d(spline, kvec[0], kvec[1], kvec[2], &val);
      sum += RealType(val);
    }
  }

  return sum / RealType(ngrid);
}

UBspline_1d_d* QMCFiniteSize::spline_clamped(std::vector<RealType>& grid,
                                             std::vector<RealType>& vals,
                                             RealType lVal,
                                             RealType rVal)
{
  //hack to interface to NUgrid stuff in double prec for MIXED build
  std::vector<FullPrecRealType> grid_fp(grid.begin(), grid.end());

  Grid_t lingrid;
  lingrid.set(grid_fp[0], grid_fp.back(), grid_fp.size());
  Ugrid esgrid = lingrid.einspline_grid();

  BCtype_d xBC;
  xBC.lVal  = lVal;
  xBC.rVal  = rVal;
  xBC.lCode = DERIV1;
  xBC.rCode = DERIV1;
  //hack to interface to NUgrid stuff in double prec for MIXED build
  std::vector<FullPrecRealType> vals_fp(vals.begin(), vals.end());
  return create_UBspline_1d_d(esgrid, xBC, vals_fp.data());
}

//Integrate the spline using Simpson's 5/8 rule.  For Bsplines, this should be exact
//provided your delta is smaller than the smallest bspline mesh spacing.
// JPT 13/03/2018 - Fixed an intermittant segfault that occurred b/c
//                  eval_NUB_spline_1d_d sometimes went out of bounds.
// #3677 changed NUBspline to UBspline.
QMCFiniteSize::RealType QMCFiniteSize::integrate_spline(UBspline_1d_d* spline, RealType a, RealType b, IndexType N)
{
  if (N % 2 != 0) // if N odd, warn that destruction is imminent
  {
    std::cerr << "Warning in integrate_spline: N must be even!\n";
    N = N - 1; // No risk of overflow
  }

  RealType eps         = (b - a) / RealType(N);
  RealType sum         = 0.0;
  FullPrecRealType tmp = 0.0; //hack to interface to UBspline_1d_d
  RealType xi          = 0.0;
  for (int i = 1; i < N / 2; i++)
  {
    xi = a + (2 * i - 2) * eps;
    eval_UBspline_1d_d(spline, xi, &tmp);
    sum += RealType(tmp);

    xi = a + (2 * i - 1) * eps;
    eval_UBspline_1d_d(spline, xi, &tmp);
    sum += 4 * tmp;

    xi = a + (2 * i) * eps;
    eval_UBspline_1d_d(spline, xi, &tmp);
    sum += tmp;
  }

  return (eps / 3.0) * sum;
}

void QMCFiniteSize::initialize()
{
  //Initialize the long range breakup. Chosen in input xml
  initBreakup();
  Ne    = P->getTotalNum();
  Vol   = P->getLattice().Volume;
  rs    = std::pow(3.0 / (4 * M_PI) * Vol / RealType(Ne), 1.0 / 3.0);
  rho   = RealType(Ne) / Vol;
  Klist = P->getSimulationCell().getKLists();
  kpts  = Klist.kpts; //These are in reduced coordinates.
                      //Easier to spline, but will have to convert
                      //for real space integration.

  if (!skparser->has_grid())
    skparser->set_grid(kpts);
  std::cout << "Grid computed.\n";

  skparser->get_grid(gridx, gridy, gridz);
}

void QMCFiniteSize::printSkRawSphAvg(const std::vector<RealType>& sk)
{
  std::vector<RealType> vsk_1d(Klist.kshell.size());

  // Average within each shell
  for (int ks = 0; ks < Klist.kshell.size() - 1; ks++)
  {
    RealType u = 0;
    RealType n = 0;
    for (int ki = Klist.kshell[ks]; ki < Klist.kshell[ks + 1]; ki++)
    {
      u += sk[ki];
      n++;
    }
    if (n != 0)
    {
      vsk_1d[ks] = u / n;
    }
    else
    {
      vsk_1d[ks] = 0;
    }
  }

  app_log() << std::fixed;
  app_log() << "\nSpherically averaged raw S(k):\n";
  app_log() << std::setw(12) << "k" << std::setw(12) << "S(k)" << std::setw(12) << "vk"
            << "\n";
  for (int ks = 0; ks < Klist.kshell.size() - 1; ks++)
  {
    app_log() << std::setw(12) << std::setprecision(8) << std::sqrt(Klist.ksq[Klist.kshell[ks]]) << std::setw(12)
              << std::setprecision(8) << vsk_1d[ks] << std::setw(12) << std::setprecision(8) << AA->Fk_symm[ks] << '\n';
  }

  if (vsk_1d[Klist.kshell.size() - 2] < 0.99)
  {
    app_log() << "####################################################################\n";
    app_log() << "WARNING: The S(k) in the largest kshell is less than 0.99\n";
    app_log() << "         This code assumes the S(k) is converged to 1.0 at large k\n";
    app_log() << "         You may need to rerun with a larger LR_dim_cutoff\n";
    app_log() << "####################################################################\n";
  }
}

void QMCFiniteSize::printSkSplineSphAvg(UBspline_3d_d* spline)
{
  std::vector<RealType> Amat;
  getSkInfo(spline, Amat);

  app_log() << "\n=========================================================\n";
  app_log() << " S(k) Info \n";
  app_log() << "=========================================================\n";
  app_log() << "S(k) anisotropy near k=0\n";
  app_log() << "------------------------\n";
  app_log() << "  a_xx = " << Amat[0] << std::endl;
  app_log() << "  a_yy = " << Amat[1] << std::endl;
  app_log() << "  a_zz = " << Amat[2] << std::endl;
  app_log() << "  a_xy = " << Amat[3] << std::endl;
  app_log() << "  a_xz = " << Amat[4] << std::endl;
  app_log() << "  a_yz = " << Amat[5] << std::endl;
  app_log() << "------------------------\n";

  RealType b = (Amat[0] + Amat[1] + Amat[2]) / 3.0;

  app_log() << "Spherically averaged S(k) near k=0\n";
  app_log() << "S(k)=b*k^2   b = " << b << std::endl;
  app_log() << "------------------------\n";
  app_log() << std::endl;

  RealType kmax = AA->get_kc();
  RealType nk   = 100;
  RealType kdel = kmax / (nk - 1.0);

  app_log() << "\nSpherically averaged splined S(k):\n";
  app_log() << std::setw(12) << "k" << std::setw(12) << "S(k)"
            << "\n";
  for (int k = 0; k < nk; k++)
  {
    RealType kval = kdel * k;
    app_log() << std::setw(12) << std::setprecision(8) << kval << std::setw(12) << std::setprecision(8)
              << sphericalAvgSk(spline, kval) << "\n";
  }
}

QMCFiniteSize::RealType QMCFiniteSize::calcPotentialDiscrete(std::vector<RealType> sk)
{
  //This is the \frac{1}{Omega} \sum_{\mathbf{k}} \frac{v_k}{2} S(\mathbf{k}) term.
  return 0.5 * AA->evaluate_w_sk(Klist.kshell, sk.data());
}

QMCFiniteSize::RealType QMCFiniteSize::calcPotentialInt(std::vector<RealType> sk)
{
  auto spline = std::unique_ptr<UBspline_3d_d, void (*)(void*)>{getSkSpline(sk), destroy_Bspline};

  RealType kmax   = AA->get_kc();
  IndexType ngrid = 2 * Klist.kshell.size() - 1; //make a lager kmesh

  std::vector<RealType> unigrid1d, k2vksk;
  RealType dk = kmax / ngrid;

  unigrid1d.push_back(0.0);
  k2vksk.push_back(0.0);
  for (int i = 1; i < ngrid; i++)
  {
    RealType kval = i * dk;
    unigrid1d.push_back(kval);
    RealType skavg = sphericalAvgSk(spline.get(), kval);
    RealType k2vk  = kval * kval * AA->evaluate_vlr_k(kval); //evaluation for arbitrary kshell for any LRHandler
    k2vksk.push_back(0.5 * k2vk * skavg);
  }

  k2vksk.push_back(0.0);
  unigrid1d.push_back(kmax);

  auto integrand =
      std::unique_ptr<UBspline_1d_d, void (*)(void*)>{spline_clamped(unigrid1d, k2vksk, 0.0, 0.0), destroy_Bspline};

  //Integrate the spline and compute the thermodynamic limit.
  RealType integratedval = integrate_spline(integrand.get(), 0.0, kmax, 200);
  RealType intnorm       = Vol / 2.0 / M_PI / M_PI; //The volume factor here is because 1/Vol is
                                                    //included in QMCPACK's v_k.  See CoulombFunctor.

  return intnorm * integratedval;
}

void QMCFiniteSize::calcPotentialCorrection()
{
  //resample vsums and vints
  std::vector<RealType> vsums, vints;
  vsums.resize(NumSamples);
  vints.resize(NumSamples);

  RandomGenerator rng;
#pragma omp parallel for
  for (int i = 0; i < NumSamples; i++)
  {
    std::vector<RealType> newSK_raw(SK_raw.size());
    for (int j = 0; j < SK_raw.size(); j++)
    {
      FullPrecRealType chi;
      chi          = rng();
      newSK_raw[j] = SK_raw[j] + SKerr_raw[j] * chi;
    }
    vsums[i] = calcPotentialDiscrete(newSK_raw);

    std::vector<RealType> newSK(SK.size());
    for (int j = 0; j < SK.size(); j++)
    {
      FullPrecRealType chi;
      chi      = rng();
      newSK[j] = SK[j] + SKerr[j] * chi;
    }
    vints[i] = calcPotentialInt(newSK);
  }

  RealType vint, vinterr;
  getStats(vints, vint, vinterr);

  RealType vsum, vsumerr;
  getStats(vsums, vsum, vsumerr);

  Vfs    = vint - vsum;
  Vfserr = std::sqrt(vinterr * vinterr + vsumerr * vsumerr);
}

void QMCFiniteSize::calcLeadingOrderCorrections()
{
  RandomGenerator rng;

  std::vector<RealType> bs(NumSamples);
#pragma omp parallel for
  for (int i = 0; i < NumSamples; i++)
  {
    std::vector<RealType> newSK(SK.size());
    for (int j = 0; j < SK.size(); j++)
    {
      FullPrecRealType chi;
      chi      = rng();
      newSK[j] = SK[j] + SKerr[j] * chi;
    }
    UBspline_3d_d* spline = getSkSpline(newSK);
    std::vector<RealType> Amat;
    getSkInfo(spline, Amat);
    bs[i] = (Amat[0] + Amat[1] + Amat[2]) / 3.0;
  }

  RealType b, berr;
  getStats(bs, b, berr);

  vlo    = 2 * M_PI * rho * b / RealType(Ne);
  vloerr = (2 * M_PI * rho / RealType(Ne)) * berr;
  tlo    = 1.0 / RealType(Ne * b * 8);
  tloerr = berr / (8 * RealType(Ne) * b * b);
}

void QMCFiniteSize::summary()
{
  // Here are the fsc corrections to potential
  app_log() << "\n=========================================================\n";
  app_log() << " Finite Size Corrections:\n";
  app_log() << "=========================================================\n";
  app_log() << " System summary:\n";
  app_log() << std::fixed;
  app_log() << "  Nelec = " << std::setw(12) << Ne << "\n";
  app_log() << "  Vol   = " << std::setw(12) << std::setprecision(8) << Vol << " [a0^3]\n";
  app_log() << "  Ne/V  = " << std::setw(12) << std::setprecision(8) << rho << " [1/a0^3]\n";
  app_log() << "  rs/a0 = " << std::setw(12) << std::setprecision(8) << rs << "\n";
  app_log() << "\n";
  app_log() << " Leading Order Corrections:\n";
  app_log() << "  V_LO / electron = " << std::setw(12) << std::setprecision(8) << vlo << " +/- " << vloerr
            << " [Ha/electron]\n";
  app_log() << "  V_LO            = " << std::setw(12) << std::setprecision(8) << vlo * Ne << " +/- " << vloerr * Ne
            << " [Ha]\n";
  app_log() << "  T_LO / electron = " << std::setw(12) << std::setprecision(8) << tlo << " +/- " << tloerr
            << " [Ha/electron]\n";
  app_log() << "  T_LO            = " << std::setw(12) << std::setprecision(8) << tlo * Ne << " +/- " << tloerr * Ne
            << " [Ha]\n";
  app_log() << "  NB: This is a crude estimate of the kinetic energy correction!\n";
  app_log() << "\n";
  app_log() << " Beyond Leading Order (Integrated corrections):\n";
  app_log() << "  V_Int / electron = " << std::setw(12) << std::setprecision(8) << Vfs << " +/- " << Vfserr
            << " [Ha/electron]\n";
  app_log() << "  V_Int            = " << std::setw(12) << std::setprecision(8) << Vfs * Ne << " +/- " << Vfserr * Ne
            << " [Ha]\n";
}

bool QMCFiniteSize::execute()
{
  initialize();
  //Print Spherical Avg from data
  SK_raw    = skparser->get_sk_raw();
  SKerr_raw = skparser->get_skerr_raw();
  if (skparser->is_normalized() == false)
  {
    for (int i = 0; i < SK_raw.size(); i++)
    {
      SK_raw[i] /= RealType(Ne);
      SKerr_raw[i] /= RealType(Ne);
    }
  }
  printSkRawSphAvg(SK_raw);

  //Print Spherical Avg from spline
  skparser->get_sk(SK, SKerr); //now have SK on full grid
  if (skparser->is_normalized() == false)
  {
    for (IndexType i = 0; i < SK.size(); i++)
    {
      SK[i] /= RealType(Ne);
      SKerr[i] /= RealType(Ne);
    }
  }
  UBspline_3d_d* sk3d_spline = getSkSpline(SK);
  printSkSplineSphAvg(sk3d_spline);

  calcLeadingOrderCorrections();
  calcPotentialCorrection();

  summary();

  return true;
}

} // namespace qmcplusplus
