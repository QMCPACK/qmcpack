//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCHamiltonians/MPC.h"
#include "Lattice/ParticleBConds.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsData/AttributeSet.h"
#include "Particle/DistanceTable.h"
#include "Particle/DistanceTableData.h"
#include "Particle/MCWalkerConfiguration.h"
#include "Utilities/IteratorUtility.h"

#if defined(HAVE_LIBFFTW)
#include <fftw3.h>
#endif

namespace qmcplusplus
{

void
MPC::resetTargetParticleSet(ParticleSet& ptcl)
{ }

MPC::MPC(ParticleSet& ptcl, double cutoff) :
  PtclRef(&ptcl), Ecut(cutoff), FirstTime(true),
  VlongSpline(0), DensitySpline(0)
{
  int it=ptcl.addTable(ptcl,DT_AOS);
  initBreakup();
}

MPC::~MPC()
{
  if (VlongSpline)
    destroy_Bspline(VlongSpline);
  if (DensitySpline)
    destroy_Bspline(DensitySpline);
}

void
MPC::init_gvecs()
{
  TinyVector<int,OHMMS_DIM> maxIndex;
  PosType b[OHMMS_DIM];
  for (int j=0; j < OHMMS_DIM; j++)
  {
    maxIndex[j] = 0;
    b[j] = 2.0*M_PI*PtclRef->Lattice.b(j);
  }
  int numG1 = PtclRef->Density_G.size();
  int numG2 = PtclRef->DensityReducedGvecs.size();
  assert (PtclRef->Density_G.size() == PtclRef->DensityReducedGvecs.size());
  // Loop through all the G-vectors, and find the largest
  // indices in each direction with energy less than the cutoff
  for (int iG=0; iG < PtclRef->DensityReducedGvecs.size(); iG++)
  {
    TinyVector<int,OHMMS_DIM> gint = PtclRef->DensityReducedGvecs[iG];
    PosType G = (double)gint[0] * b[0];
    for (int j=1; j<OHMMS_DIM; j++)
      G += (double)gint[j]*b[j];
    if (0.5*dot(G,G) < Ecut)
    {
      for (int j=0; j<OHMMS_DIM; j++)
        maxIndex[j] = std::max(maxIndex[j], std::abs(gint[j]));
      Gvecs.push_back(G);
      Gints.push_back(gint);
      Rho_G.push_back(PtclRef->Density_G[iG]);
    }
  }
  SplineDim = 4 * maxIndex;
  MaxDim = std::max(maxIndex[0], std::max(maxIndex[1], maxIndex[2]));
  app_log() << "  Using " << Gvecs.size()
            << " G-vectors for MPC interaction.\n";
  app_log() << "   Using real-space box of size [" << SplineDim[0]
            << "," << SplineDim[1] << "," << SplineDim[2]
            << "] for MPC spline.\n";
}



void
MPC::compute_g_G(double &g_0, std::vector<double> &g_G, int N)
{
  double L = PtclRef->Lattice.WignerSeitzRadius;
  double Linv = 1.0/L;
  double Linv3 = Linv*Linv*Linv;
  // create an FFTW plan
  Array<std::complex<double>,3> rBox(N,N,N);
  Array<std::complex<double>,3> GBox(N,N,N);
  // app_log() << "Doing " << N << " x " << N << " x " << N << " FFT.\n";
  //create BC handler
  DTD_BConds<RealType,3,SUPERCELL_BULK> mybc(PtclRef->Lattice);
  // Fill the real-space array with f(r)
  double Ninv = 1.0/(double)N;
  TinyVector<RealType,3> u, r;
  for (int ix=0; ix<N; ix++)
  {
    u[0] = Ninv*ix;
    for (int iy=0; iy<N; iy++)
    {
      u[1] = Ninv*iy;
      for (int iz=0; iz<N; iz++)
      {
        u[2] = Ninv*iz;
        r = PtclRef->Lattice.toCart (u);
        //DTD_BConds<double,3,SUPERCELL_BULK>::apply (PtclRef->Lattice, r);
        //double rmag = std::sqrt(dot(r,r));
        double rmag = std::sqrt(mybc.apply_bc(r));
        if (rmag < L)
          rBox(ix,iy,iz) = -0.5*rmag*rmag*Linv3 + 1.5*Linv;
        else
          rBox(ix,iy,iz) = 1.0/rmag;
      }
    }
  }
  fftw_plan fft = fftw_plan_dft_3d
                  (N, N, N, (fftw_complex*)rBox.data(),
                   (fftw_complex*) GBox.data(), 1, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);
  // Now, copy data into output, and add on analytic part
  double norm = Ninv*Ninv*Ninv;
  int numG = Gints.size();
  for (int iG=0; iG < numG; iG++)
  {
    TinyVector<int,OHMMS_DIM> gint = Gints[iG];
    for (int j=0; j<OHMMS_DIM; j++)
      gint[j] = (gint[j] + N)%N;
    g_G[iG] = norm * real(GBox(gint[0], gint[1], gint[2]));
  }
  g_0 = norm * real(GBox(0,0,0));
}



inline double
extrap (int N, TinyVector<double,3> g_124)
{
  Tensor<double,3> A, Ainv;
  double N1inv2 = 1.0/(double)(N*N);
  double N2inv2 = 1.0/(double)((2*N)*(2*N));
  double N4inv2 = 1.0/(double)((4*N)*(4*N));
  A(0,0) = 1.0;
  A(0,1) = N1inv2;
  A(0,2) = N1inv2*N1inv2;
  A(1,0) = 1.0;
  A(1,1) = N2inv2;
  A(1,2) = N2inv2*N2inv2;
  A(2,0) = 1.0;
  A(2,1) = N4inv2;
  A(2,2) = N4inv2*N4inv2;
  Ainv = inverse(A);
  TinyVector<double,3> abc = dot(Ainv,g_124);
  return abc[0];
}


inline double
extrap (int N, TinyVector<double,2> g_12)
{
  Tensor<double,2> A, Ainv;
  double N1inv2 = 1.0/(double)(N*N);
  double N2inv2 = 1.0/(double)((2*N)*(2*N));
  A(0,0) = 1.0;
  A(0,1) = N1inv2;
  A(1,0) = 1.0;
  A(1,1) = N2inv2;
  Ainv = inverse(A);
  TinyVector<double,2> abc = dot(Ainv,g_12);
  return abc[0];
}



void
MPC::init_f_G()
{
  int numG = Gints.size();
  f_G.resize(numG);
  int N = std::max(64, 2*MaxDim+1);
  std::vector<double> g_G_N(numG), g_G_2N(numG), g_G_4N(numG);
  double g_0_N, g_0_2N, g_0_4N;
  compute_g_G(g_0_N , g_G_N,  1*N);
  compute_g_G(g_0_2N, g_G_2N, 2*N);
  compute_g_G(g_0_4N, g_G_4N, 4*N);
  // fprintf (stderr, "g_G_1N[0]      = %18.14e\n", g_G_N[0]);
  // fprintf (stderr, "g_G_2N[0]      = %18.14e\n", g_G_2N[0]);
  // fprintf (stderr, "g_G_4N[0]      = %18.14e\n", g_G_4N[0]);
  double volInv = 1.0/PtclRef->Lattice.Volume;
  double L = PtclRef->Lattice.WignerSeitzRadius;
  TinyVector<double,2> g0_12  (g_0_2N, g_0_4N);
  TinyVector<double,3> g0_124 (g_0_N, g_0_2N, g_0_4N);
  f_0 = extrap (N, g0_124);
  app_log().precision(12);
  app_log() << "    Linear extrap    = " << std::scientific
            << extrap (2*N, g0_12) << std::endl;
  app_log() << "    Quadratic extrap = " << std::scientific
            << f_0 << std::endl;
  f_0 += 0.4*M_PI*L*L*volInv;
  // std::cerr << "f_0 = " << f_0/volInv << std::endl;
  double worst = 0.0, worstLin, worstQuad;
  int iworst = 0;
  for (int iG=0; iG<numG; iG++)
  {
    TinyVector<double,2> g_12  (g_G_2N[iG], g_G_4N[iG]);
    TinyVector<double,3> g_124 (g_G_N [iG], g_G_2N[iG], g_G_4N[iG]);
    double linearExtrap = extrap (2*N, g_12);
    double quadExtrap   = extrap (N, g_124);
    double diff = std::abs(linearExtrap - quadExtrap);
    if (diff > worst)
    {
      worst = diff;
      iworst = iG;
      worstLin  = linearExtrap;
      worstQuad = quadExtrap;
    }
    f_G[iG] = quadExtrap;
    double G2 = dot(Gvecs[iG], Gvecs[iG]);
    double G  = std::sqrt(G2);
    f_G[iG] += volInv*M_PI*(4.0/G2 + 12.0/(L*L*G2*G2)*
                            (std::cos(G*L) - std::sin(G*L)/(G*L)));
    // std::cerr << "f_G = " << f_G[iG]/volInv << std::endl;
    // std::cerr << "f_G - 4*pi/G2= " << f_G[iG]/volInv - 4.0*M_PI/G2 << std::endl;
  }
  char buff[1000];
  snprintf (buff, 1000 ,
            "    Worst MPC discrepancy:\n"
            "      Linear Extrap   : %18.14e\n"
            "      Quadratic Extrap: %18.14e\n", worstLin, worstQuad);
  app_log() << buff;
}



void
MPC::init_spline()
{
  Array<std::complex<double>,3> rBox(SplineDim[0], SplineDim[1], SplineDim[2]),
        GBox(SplineDim[0], SplineDim[1], SplineDim[2]);
  Array<double,3> splineData(SplineDim[0], SplineDim[1], SplineDim[2]);
  GBox = std::complex<double>();
  Vconst = 0.0;
  // Now fill in elements of GBox
  const RealType vol = PtclRef->Lattice.Volume;
  const RealType volInv = 1.0/PtclRef->Lattice.Volume;
  const RealType halfvol=vol/2.0;
  for (int iG=0; iG < Gvecs.size(); iG++)
  {
    TinyVector<int,OHMMS_DIM> gint = Gints[iG];
    PosType G = Gvecs[iG];
    double G2 = dot(G,G);
    TinyVector<int,OHMMS_DIM> index;
    for (int j=0; j<OHMMS_DIM; j++)
      index[j] = (gint[j] + SplineDim[j]) % SplineDim[j];

    const RealType xxx(vol*(4.0*M_PI*volInv/G2-f_G[iG]));
    if (!(index[0]==0 && index[1]==0 && index[2]==0))
    {
      GBox(index[0], index[1], index[2]) = xxx*Rho_G[iG];
      Vconst -= halfvol*xxx*norm(Rho_G[iG]);
      //GBox(index[0], index[1], index[2]) = vol *
      //                                     Rho_G[iG] * (4.0*M_PI*volInv/G2 - f_G[iG]);
      //Vconst -= 0.5 * vol * vol * norm(Rho_G[iG])
      //          * (4.0*M_PI*volInv/G2 - f_G[iG]);
    }
  }
  // G=0 component calculated seperately
  //GBox(0,0,0) = -vol * f_0 * Rho_G[0];
  //Vconst += 0.5 * vol * vol * f_0 * norm(Rho_G[0]);
  const RealType volf=vol*f_0;
  GBox(0,0,0) = -volf * Rho_G[0];
  Vconst += halfvol * volf * norm(Rho_G[0]);

  app_log() << "  Constant potential = " << Vconst << std::endl;
  fftw_plan fft = fftw_plan_dft_3d
                  (SplineDim[0], SplineDim[1], SplineDim[2], (fftw_complex*)GBox.data(),
                   (fftw_complex*) rBox.data(), -1, FFTW_ESTIMATE);
  fftw_execute (fft);
  fftw_destroy_plan (fft);
  for (int i0=0; i0<SplineDim[0]; i0++)
    for (int i1=0; i1<SplineDim[1]; i1++)
      for (int i2=0; i2<SplineDim[2]; i2++)
        splineData(i0, i1, i2) = real(rBox(i0,i1,i2));
  BCtype_d bc0, bc1, bc2;
  Ugrid grid0, grid1, grid2;
  grid0.start=0.0;
  grid0.end=1.0;
  grid0.num = SplineDim[0];
  grid1.start=0.0;
  grid1.end=1.0;
  grid1.num = SplineDim[1];
  grid2.start=0.0;
  grid2.end=1.0;
  grid2.num = SplineDim[2];
  bc0.lCode = bc0.rCode = PERIODIC;
  bc1.lCode = bc1.rCode = PERIODIC;
  bc2.lCode = bc2.rCode = PERIODIC;
  VlongSpline = create_UBspline_3d_d (grid0, grid1, grid2, bc0, bc1, bc2,
                                      splineData.data());
//     grid0.num = PtclRef->Density_r.size(0);
//     grid1.num = PtclRef->Density_r.size(1);
//     grid2.num = PtclRef->Density_r.size(2);
//     DensitySpline = create_UBspline_3d_d (grid0, grid1, grid2, bc0, bc1, bc2,
// 					  PtclRef->Density_r.data());
}

void
MPC::initBreakup()
{
  NParticles = PtclRef->getTotalNum();
  app_log() << "\n  === Initializing MPC interaction === " << std::endl;
  if (PtclRef->Density_G.size() == 0)
  {
    app_error() << "************************\n"
                << "** Error in MPC setup **\n"
                << "************************\n"
                << "    The electron density was not setup by the "
                << "wave function builder.\n";
    abort();
  }
  init_gvecs();
  init_f_G();
  init_spline();
  // FILE *fout = fopen ("MPC.dat", "w");
  // double vol = PtclRef->Lattice.Volume;
  // PosType r0 (0.0, 0.0, 0.0);
  // PosType r1 (10.26499236, 10.26499236, 10.26499236);
  // int nPoints=1001;
  // for (int i=0; i<nPoints; i++) {
  //   double s = (double)i/(double)(nPoints-1);
  //   PosType r = (1.0-s)*r0 + s*r1;
  //   PosType u = PtclRef->Lattice.toUnit(r);
  //   double V, rho(0.0);
  //   eval_UBspline_3d_d (VlongSpline, u[0], u[1], u[2], &V);
  //   // eval_UBspline_3d_d (DensitySpline, u[0], u[1], u[2], &rho);
  //   fprintf (fout, "%6.4f %14.10e %14.10e\n", s, V, rho);
  // }
  // fclose(fout);
  app_log() << "  === MPC interaction initialized === \n\n";
}

QMCHamiltonianBase*
MPC::makeClone(ParticleSet& qp, TrialWaveFunction& psi)
{
  // return new MPC(qp, Ecut);
  MPC* newMPC = new MPC(*this);
  newMPC->resetTargetParticleSet(qp);
  return newMPC;
}

MPC::Return_t
MPC::evalSR(ParticleSet& P) const
{
  const DistanceTableData& d_aa=(*P.DistTables[0]);
  RealType SR=0.0;
  if(d_aa.DTType == DT_SOA)
  {
    const RealType cone(1);
    for(size_t ipart=0; ipart<NParticles; ipart++)
    {
      RealType esum(0);
      const RealType* restrict dist=d_aa.Distances[ipart];
      for(size_t j=0; j<ipart; ++j)
        esum += cone/dist[j];
      SR += esum;
    }
  }
  else
  {
    for(int ipart=0; ipart<NParticles; ipart++)
    {
      RealType esum = 0.0;
      for(int nn=d_aa.M[ipart],jpart=ipart+1;
          nn<d_aa.M[ipart+1]; nn++,jpart++)
        esum += d_aa.rinv(nn);
      //Accumulate pair sums...species charge for atom i.
      SR += esum;
    }
  }
  return SR;
}

MPC::Return_t
MPC::evalLR(ParticleSet& P) const
{
  RealType LR=0.0;
  double val;
  for(int i=0; i<NParticles; i++)
  {
    //PosType r = P.R[i];
    //PosType u = P.Lattice.toUnit(r);
    PosType u = P.Lattice.toUnit(P.R[i]);
    for (int j=0; j<OHMMS_DIM; j++)
      u[j] -= std::floor(u[j]);
    eval_UBspline_3d_d (VlongSpline, u[0], u[1], u[2], &val);
    LR += val;
  }
  return LR;
}

MPC::Return_t
MPC::evaluate (ParticleSet& P)
{
  //if (FirstTime || P.tag() == PtclRef->tag())
  Value = evalSR(P) + evalLR(P) + Vconst;
  return Value;
}

void MPC::addEnergy(MCWalkerConfiguration &W, std::vector<RealType> &LocalEnergy)
{
//only used for debugging
//  const int nw=W.getActiveWalkers();
//  const int indx=NUMPROPERTIES+myIndex;
//#ifdef QMC_CUDA
//#pragma omp parallel for
//  for(int iw=0; iw<nw; ++iw)
//  {
//    ParticleSet& p(*myPtcl[omp_get_thread_num()]);
//    p.R=W[iw]->R;
//    p.update();
//    RealType e=evalSR(p)+evalLR(p)+Vconst;
//    W[iw]->getPropertyBase()[indx] = e;
//    LocalEnergy[iw]+=e;
//  }
//#endif
}

bool MPC::put(xmlNodePtr cur)
{
  Ecut = -1.0;
  OhmmsAttributeSet attribs;
  attribs.add (Ecut, "cutoff");
  attribs.put (cur);
  if (Ecut < 0.0)
  {
    Ecut = 30.0;
    app_log() << "    MPC cutoff not found.  Set using \"cutoff\" attribute.\n"
              << "    Setting to default value of " << Ecut << std::endl;
  }
  return true;
}

}
