#include "QMCHamiltonians/MPC.h"
#include "Lattice/ParticleBConds.h"
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsData/AttributeSet.h"

#include <fftw3.h>

namespace qmcplusplus {

  MPC::MPC(ParticleSet& ptcl, double cutoff) :
    PtclRef(ptcl), Ecut(cutoff)
  {
    initBreakup();
 
  }

  MPC::~MPC() {


  }
  
  

  
  void
  MPC::compute_g_G(double &g_0, vector<double> &g_G, int N)
  {
    double L = PtclRef.Lattice.WignerSeitzRadius;
    double Linv = 1.0/L;
    double Linv3 = Linv*Linv*Linv;
    // create an FFTW plan
    Array<complex<double>,3> rBox(N,N,N);
    Array<complex<double>,3> GBox(N,N,N);
    // app_log() << "Doing " << N << " x " << N << " x " << N << " FFT.\n";
    
    // Fill the real-space array with f(r)
    double Ninv = 1.0/(double)N;
    TinyVector<double,3> u, r;
    for (int ix=0; ix<N; ix++) {
      u[0] = Ninv*ix;
      for (int iy=0; iy<N; iy++) {
	u[1] = Ninv*iy;
	for (int iz=0; iz<N; iz++) {
	  u[2] = Ninv*iz;
	  r = PtclRef.Lattice.toCart (u);
	  DTD_BConds<double,3,SUPERCELL_BULK>::apply (PtclRef.Lattice, r);
	  double rmag = std::sqrt(dot(r,r));
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
    KContainer &Glist = PtclRef.SK->KLists;
    int numG = Glist.kpts.size();
    for (int iG=0; iG < numG; iG++) {
      TinyVector<int,OHMMS_DIM> gint = Glist.kpts[iG];
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
    A(0,0) = 1.0; A(0,1) = N1inv2; A(0,2) = N1inv2*N1inv2;
    A(1,0) = 1.0; A(1,1) = N2inv2; A(1,2) = N2inv2*N2inv2;
    A(2,0) = 1.0; A(2,1) = N4inv2; A(2,2) = N4inv2*N4inv2;
    Ainv = inverse(A);
    // fprintf (stderr, "Ainv = [ %8.5f %8.5f %8.5f\n"
    // 	     "         %8.5f %8.5f %8.5f\n"
    // 	     "         %8.5f %8.5f %8.5f ] \n");

    TinyVector<double,3> abc = dot(Ainv,g_124);
    return abc[0];
  }
      

  inline double 
  extrap (int N, TinyVector<double,2> g_12)
  {
    Tensor<double,2> A, Ainv;
    double N1inv2 = 1.0/(double)(N*N);
    double N2inv2 = 1.0/(double)((2*N)*(2*N));
    A(0,0) = 1.0; A(0,1) = N1inv2; 
    A(1,0) = 1.0; A(1,1) = N2inv2; 
    Ainv = inverse(A);
    TinyVector<double,2> abc = dot(Ainv,g_12);
    return abc[0];
  }



  void
  MPC::init_f_G()
  {
    KContainer &Glist = PtclRef.SK->KLists;
    int numG = Glist.kpts.size();
    f_G.resize(numG);

    int N = max(64, 4*Glist.mmax[OHMMS_DIM]);
    vector<double> g_G_N(numG), g_G_2N(numG), g_G_4N(numG);
    double g_0_N, g_0_2N, g_0_4N;
    compute_g_G(g_0_N , g_G_N,  1*N);
    compute_g_G(g_0_2N, g_G_2N, 2*N);
    compute_g_G(g_0_4N, g_G_4N, 4*N);
    // fprintf (stderr, "g_G_1N[0]      = %18.14e\n", g_G_N[0]);
    // fprintf (stderr, "g_G_2N[0]      = %18.14e\n", g_G_2N[0]);
    // fprintf (stderr, "g_G_4N[0]      = %18.14e\n", g_G_4N[0]);


    double volInv = 1.0/PtclRef.Lattice.Volume;
    double L = PtclRef.Lattice.WignerSeitzRadius;
    
    TinyVector<double,2> g0_12  (g_0_2N, g_0_4N);
    TinyVector<double,3> g0_124 (g_0_N, g_0_2N, g_0_4N);
    f_0 = extrap (N, g0_124);
    
    app_log().precision(12);
    app_log() << "    Linear extrap    = " << std::scientific 
	      << extrap (2*N, g0_12) << endl;
    app_log() << "    Quadratic extrap = " << std::scientific 
	      << f_0 << endl;
    f_0 += 0.4*M_PI*L*L*volInv;
    // cerr << "f_0 = " << f_0/volInv << endl;

    double worst = 0.0, worstLin, worstQuad;
    int iworst = 0;
    for (int iG=0; iG<numG; iG++) {
      TinyVector<double,2> g_12  (g_G_2N[iG], g_G_4N[iG]);
      TinyVector<double,3> g_124 (g_G_N [iG], g_G_2N[iG], g_G_4N[iG]);
      double linearExtrap = extrap (2*N, g_12);
      double quadExtrap   = extrap (N, g_124);
      double diff = std::fabs(linearExtrap - quadExtrap);
      if (diff > worst) {
	worst = diff;
	iworst = iG;
	worstLin  = linearExtrap;
	worstQuad = quadExtrap;
      }
      f_G[iG] = quadExtrap;
      double G2 = Glist.ksq[iG];
      double G  = std::sqrt(G2);
      f_G[iG] += volInv*M_PI*(4.0/G2 + 12.0/(L*L*G2*G2)*
			      (std::cos(G*L) - std::sin(G*L)/(G*L)));
      // cerr << "f_G = " << f_G[iG]/volInv << endl;
      // cerr << "f_G - 4*pi/G2= " << f_G[iG]/volInv - 4.0*M_PI/G2 << endl;
    }
    fprintf (stderr, "    Worst MPC discrepancy:\n"
	     "      Linear Extrap   : %18.14e\n"
	     "      Quadratic Extrap: %18.14e\n", worstLin, worstQuad);
    
  }

  void
  MPC::init_spline() 
  {
    TinyVector<int,OHMMS_DIM> maxIndex;
    PosType b[OHMMS_DIM];
    for (int j=0; j < OHMMS_DIM; j++) {
      maxIndex[j] = 0;
      b[j] = PtclRef.Lattice.b(j);
    }
    // Loop through all the G-vectors, and find the largest
    // indices in each direction with energy less than the cutoff
    
    for (int iG=0; iG < PtclRef.DensityReducedGvecs.size(); iG++) {
      TinyVector<int,OHMMS_DIM> gint = PtclRef.DensityReducedGvecs[iG];
      PosType G = (double)gint[0] * b[0];
      for (int j=1; j<OHMMS_DIM; j++) 
	G += (double)gint[j]*b[j];
      if (0.5*dot(G,G) < Ecut) 
	for (int j=0; j<OHMMS_DIM; j++)
	  maxIndex[j] = max(maxIndex[j], abs(gint[j]));
    }
    TinyVector<int,OHMMS_DIM> boxSize = 4 * maxIndex;
    app_log() << "   Using real-space box of size [" << boxSize[0] 
	      << "," << boxSize[1] << "," << boxSize[2] << "] for MPC spline.\n";

    Array<complex<double>,3> rBox(boxSize[0], boxSize[1], boxSize[2]),
      GBox(boxSize[0], boxSize[1], boxSize[2]);

    GBox = complex<double>();

    // Now fill in elements of GBox
    double volInv = 1.0/PtclRef.Lattice.Volume;
    for (int iG=0; iG < PtclRef.DensityReducedGvecs.size(); iG++) {
      TinyVector<int,OHMMS_DIM> gint = PtclRef.DensityReducedGvecs[iG];
      PosType G = (double)gint[0] * b[0];
      complex<double> rho = PtclRef.Density_G[iG];
      for (int j=1; j<OHMMS_DIM; j++) 
	G += (double)gint[j]*b[j];
      if (0.5*dot(G,G) < Ecut) {
	TinyVector<int,OHMMS_DIM> index;
	for (int j=0; j<OHMMS_DIM; j++)
	  index[j] = (gint[j] + boxSize[j]) % boxSize[j];
	// GBox(index[0], index[1], index[2]) = 
      }
    }
  
    fftw_plan fft = fftw_plan_dft_3d 
      (boxSize[0], boxSize[1], boxSize[2], (fftw_complex*)rBox.data(), 
       (fftw_complex*) GBox.data(), 1, FFTW_ESTIMATE);
    fftw_execute (fft);
    fftw_destroy_plan (fft);

  }

  void
  MPC::initBreakup()
  {
    app_log() << "\n  === Initializing MPC interaction === " << endl;
    if (!PtclRef.Density_G.size()) {
      app_error() << "************************\n"
		  << "** Error in MPC setup **\n"
		  << "************************\n"
		  << "    The electron density was not setup by the wave function builder.\n";
      abort();
    }
    init_f_G();
    init_spline();
    app_log() << "  === MPC interaction initialized === \n\n";
  }

  bool 
  MPC::put(xmlNodePtr cur)
  {
    Ecut = -1.0;
    OhmmsAttributeSet attribs;
    attribs.add (Ecut, "cutoff");
    attribs.put (cur);
    if (Ecut < 0.0) {
      Ecut = 30.0;
      app_log() << "    MPC cutoff not found.  Set using \"cutoff\" attribute.\n"
		<< "    Setting to default value of " << Ecut << endl;
    }
    return true;
  }


  MPC::Return_t
  MPC::evaluate(ParticleSet &P)
  {

    return 0.0;
  }


  QMCHamiltonianBase* 
  MPC::clone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    return new MPC(qp, Ecut);
  }

}
