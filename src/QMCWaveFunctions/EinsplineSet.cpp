//////////////////////////////////////////////////////////////////
// (c) Copyright 2006-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &          //
//   Materials Computation Center                               //
//   University of Illinois, Urbana-Champaign                   //
//   Urbana, IL 61801                                           //
//   e-mail: jnkim@ncsa.uiuc.edu                                //
//                                                              //
// Supported by                                                 //
//   National Center for Supercomputing Applications, UIUC      //
//   Materials Computation Center, UIUC                         //
//////////////////////////////////////////////////////////////////

#include "QMCWaveFunctions/EinsplineSet.h"
#include <einspline/multi_bspline.h>
#include "Configuration.h"
#ifdef HAVE_MKL
  #include "/opt/intel/mkl/10.0.3.020/include/mkl_vml.h"
#endif

namespace qmcplusplus {

  EinsplineSet::UnitCellType
  EinsplineSet::GetLattice()
  {
    return SuperLattice;
  }
  
  void
  EinsplineSet::resetTargetParticleSet(ParticleSet& e)
  {
  }

  void
  EinsplineSet::resetSourceParticleSet(ParticleSet& ions)
  {
  }
  
  void
  EinsplineSet::setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
  }
  
  void 
  EinsplineSet::evaluate (const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    app_error() << "Should never instantiate EinsplineSet.\n";
    abort();
  }

  void 
  EinsplineSet::evaluate (const ParticleSet& P, int iat, 
			  ValueVector_t& psi, GradVector_t& dpsi, 
			  ValueVector_t& d2psi)
  {
    app_error() << "Should never instantiate EinsplineSet.\n";
    abort();
  }

  
  void 
  EinsplineSet::evaluate (const ParticleSet& P, int first, int last,
			  ValueMatrix_t& vals, GradMatrix_t& grads, 
			  ValueMatrix_t& lapls)
  {
    app_error() << "Should never instantiate EinsplineSet.\n";
    abort();
  }


  void 
  EinsplineSetLocal::evaluate (const ParticleSet& P, int iat, 
			       ValueVector_t& psi)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    for(int j=0; j<OrbitalSetSize; j++) {
      complex<double> val;
      Orbitals[j]->evaluate(ru, val);

      double phase = -dot(r, Orbitals[j]->kVec);
      double s,c;
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      val *= e_mikr;
#ifdef QMC_COMPLEX
      psi[j] = val;
#else
      psi[j] = real(val);
#endif
    }
  }
  
  void 
  EinsplineSetLocal::evaluate (const ParticleSet& P, int iat, 
			       ValueVector_t& psi, GradVector_t& dpsi, 
			       ValueVector_t& d2psi)
  {
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    ru[0] -= std::floor (ru[0]);
    ru[1] -= std::floor (ru[1]);
    ru[2] -= std::floor (ru[2]);
    complex<double> val;
    TinyVector<complex<double>,3> gu;
    Tensor<complex<double>,3> hess;
    complex<double> eye (0.0, 1.0);
    for(int j=0; j<OrbitalSetSize; j++) {
      complex<double> u;
      TinyVector<complex<double>,3> gradu;
      complex<double> laplu;

      Orbitals[j]->evaluate(ru, val, gu, hess);
      u  = val;
      // Compute gradient in cartesian coordinates
      gradu = dot(PrimLattice.G, gu);
      laplu = trace(hess, GGt);      
      
      PosType k = Orbitals[j]->kVec;
      TinyVector<complex<double>,3> ck;
      ck[0]=k[0];  ck[1]=k[1];  ck[2]=k[2];
      double s,c;
      double phase = -dot(P.R[iat], k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
#ifdef QMC_COMPLEX
      psi[j]   = e_mikr * u;
      dpsi[j]  = e_mikr*(-eye * ck * u + gradu);
      d2psi[j] = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
#else
      psi[j]   = real(e_mikr * u);
      dpsi[j]  = real(e_mikr*(-eye * ck * u + gradu));
      d2psi[j] = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
#endif

    }
  }
  
  void 
  EinsplineSetLocal::evaluate (const ParticleSet& P, int first, int last,
			       ValueMatrix_t& vals, GradMatrix_t& grads, 
			       ValueMatrix_t& lapls)
  {
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType r (P.R[iat]);
      PosType ru(PrimLattice.toUnit(r));
      ru[0] -= std::floor (ru[0]);
      ru[1] -= std::floor (ru[1]);
      ru[2] -= std::floor (ru[2]);
      complex<double> val;
      TinyVector<complex<double>,3> gu;
      Tensor<complex<double>,3> hess;
      complex<double> eye (0.0, 1.0);
      for(int j=0; j<OrbitalSetSize; j++) {
	complex<double> u;
	TinyVector<complex<double>,3> gradu;
	complex<double> laplu;

	Orbitals[j]->evaluate(ru, val, gu, hess);
	u  = val;
	gradu = dot(PrimLattice.G, gu);
	laplu = trace(hess, GGt);
	
	PosType k = Orbitals[j]->kVec;
	TinyVector<complex<double>,3> ck;
	ck[0]=k[0];  ck[1]=k[1];  ck[2]=k[2];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
#ifdef QMC_COMPLEX
	vals(j,i)  = e_mikr * u;
	grads(i,j) = e_mikr*(-eye*u*ck + gradu);
	lapls(i,j) = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
#else
	vals(j,i)  = real(e_mikr * u);
	grads(i,j) = real(e_mikr*(-eye*u*ck + gradu));
	lapls(i,j) = real(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu));
#endif

      }
    }
  }

  string 
  EinsplineSet::Type()
  {
    return "EinsplineSet";
  }




  ///////////////////////////////////////////
  // EinsplineSetExtended Member functions //
  ///////////////////////////////////////////

  inline void
  convert (complex<double> a, complex<double> &b)
  { b = a;  }

  inline void
  convert (complex<double> a, double &b)
  { b = a.real();  }

  inline void
  convert (double a, complex<double>&b)
  { b = complex<double>(a,0.0); }

  inline void
  convert (double a, double &b)
  { b = a; }

  template<typename T1, typename T2> void
  convertVec (TinyVector<T1,3> a, TinyVector<T2,3> &b)
  {
    for (int i=0; i<3; i++)
      convert (a[i], b[i]);
  }

  // Real evaluation functions
  inline void 
  EinsplineMultiEval (multi_UBspline_3d_d *restrict spline,
		      TinyVector<double,3> r, 
		      Vector<double> &psi)
  {
    eval_multi_UBspline_3d_d (spline, r[0], r[1], r[2], psi.data());
  }

  inline void
  EinsplineMultiEval (multi_UBspline_3d_d *restrict spline,
		      TinyVector<double,3> r,
		      Vector<double> &psi,
		      Vector<TinyVector<double,3> > &grad,
		      Vector<Tensor<double,3> > &hess)
  {
    eval_multi_UBspline_3d_d_vgh (spline, r[0], r[1], r[2],
				  psi.data(), 
				  (double*)grad.data(), 
				  (double*)hess.data());
  }

  //////////////////////////////////
  // Complex evaluation functions //
  //////////////////////////////////
  inline void 
  EinsplineMultiEval (multi_UBspline_3d_z *restrict spline,
		      TinyVector<double,3> r, 
		      Vector<complex<double> > &psi)
  {
    eval_multi_UBspline_3d_z (spline, r[0], r[1], r[2], psi.data());
  }


  inline void
  EinsplineMultiEval (multi_UBspline_3d_z *restrict spline,
		      TinyVector<double,3> r,
		      Vector<complex<double> > &psi,
		      Vector<TinyVector<complex<double>,3> > &grad,
		      Vector<Tensor<complex<double>,3> > &hess)
  {
    eval_multi_UBspline_3d_z_vgh (spline, r[0], r[1], r[2],
				  psi.data(), 
				  (complex<double>*)grad.data(), 
				  (complex<double>*)hess.data());
  }
			   

  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::resetParameters(const opt_variables_type& active)
  {

  }

  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::resetTargetParticleSet(ParticleSet& e)
  {
  }

  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::setOrbitalSetSize(int norbs)
  {
    OrbitalSetSize = norbs;
  }
  
  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int iat, RealValueVector_t& psi)
  {
    ValueTimer.start();
    PosType r (P.R[iat]);

    // Do core states first
    int icore = NumValenceOrbs;
    for (int tin=0; tin<MuffinTins.size(); tin++) {
      MuffinTins[tin].evaluateCore(r, StorageValueVector, icore);
      icore += MuffinTins[tin].get_num_core();
    }
    // Add phase to core orbitals
    for (int j=NumValenceOrbs; j<StorageValueVector.size(); j++) {
      PosType k = kPoints[j];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      StorageValueVector[j] *= e_mikr;
    }

    // Check if we are inside a muffin tin.  If so, compute valence
    // states in the muffin tin.
    bool inTin = false;
    bool need2blend = false;
    double b(0.0);
    for (int tin=0; tin<MuffinTins.size() && !inTin; tin++) {
      MuffinTins[tin].inside(r, inTin, need2blend);
      if (inTin) {
	MuffinTins[tin].evaluate (r, StorageValueVector);
	if (need2blend) {
	  PosType disp = MuffinTins[tin].disp(r);
	  double dr = std::sqrt(dot(disp, disp));
	  MuffinTins[tin].blend_func(dr, b);
	}
	break;
      }
    }

    StorageValueVector_t &valVec = 
      need2blend ? BlendValueVector : StorageValueVector;

    if (!inTin || need2blend) {
      PosType ru(PrimLattice.toUnit(P.R[iat]));
      for (int i=0; i<OHMMS_DIM; i++)
	ru[i] -= std::floor (ru[i]);
      EinsplineTimer.start();
      EinsplineMultiEval (MultiSpline, ru, valVec);
      EinsplineTimer.stop();
      // Add e^ikr phase to B-spline orbitals
      for (int j=0; j<NumValenceOrbs; j++) {
	PosType k = kPoints[j];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
	valVec[j] *= e_mikr;
      }
    }
    int N = StorageValueVector.size();

    // If we are in a muffin tin, don't add the e^ikr term
    // We should add it to the core states, however

    if (need2blend) {
      int psiIndex = 0;
      for (int j=0; j<N; j++) {
	complex<double> psi1 = StorageValueVector[j];
	complex<double> psi2 =   BlendValueVector[j];
	
	complex<double> psi_val = b * psi1 + (1.0-b) * psi2;

	psi[psiIndex] = real(psi_val);
	psiIndex++;
	if (MakeTwoCopies[j]) {
	  psi[psiIndex] = imag(psi_val);
	  psiIndex++;
	}
      }

    }
    else {
      int psiIndex = 0;
      for (int j=0; j<N; j++) {
	complex<double> psi_val = StorageValueVector[j];
	psi[psiIndex] = real(psi_val);
	psiIndex++;
	if (MakeTwoCopies[j]) {
	  psi[psiIndex] = imag(psi_val);
	  psiIndex++;
	}
      }
    }
    ValueTimer.stop();
  }


  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int iat, ComplexValueVector_t& psi)
  {
    ValueTimer.start();
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= std::floor (ru[i]);
    EinsplineTimer.start();
    EinsplineMultiEval (MultiSpline, ru, StorageValueVector);
    EinsplineTimer.stop();
    //computePhaseFactors(r);
    for (int i=0; i<psi.size(); i++) {
      PosType k = kPoints[i];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      convert (e_mikr*StorageValueVector[i], psi[i]);
    }
    ValueTimer.stop();
  }

  // This is an explicit specialization of the above for real orbitals
  // with a real return value, i.e. simulations at the gamma or L 
  // point.
  template<> void
  EinsplineSetExtended<double>::evaluate
  (const ParticleSet &P, int iat, RealValueVector_t& psi)
  {
    ValueTimer.start();
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    int sign=0;
    for (int i=0; i<OHMMS_DIM; i++) {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      sign += HalfG[i] * (int)img;
    }
    EinsplineTimer.start();
    EinsplineMultiEval (MultiSpline, ru, psi);
    EinsplineTimer.stop();
    if (sign & 1) 
      for (int j=0; j<psi.size(); j++)
	psi[j] *= -1.0;
    ValueTimer.stop();
  }

  // Value, gradient, and laplacian
  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int iat, RealValueVector_t& psi, 
   RealGradVector_t& dpsi, RealValueVector_t& d2psi)
  {
    VGLTimer.start();
    PosType r (P.R[iat]);
    complex<double> eye (0.0, 1.0);
    
    // Do core states first
    int icore = NumValenceOrbs;
    for (int tin=0; tin<MuffinTins.size(); tin++) {
      MuffinTins[tin].evaluateCore(r, StorageValueVector, StorageGradVector, 
				   StorageLaplVector, icore);
      icore += MuffinTins[tin].get_num_core();
    }
    
    // Add phase to core orbitals
    for (int j=NumValenceOrbs; j<StorageValueVector.size(); j++) {
      complex<double> u = StorageValueVector[j];
      TinyVector<complex<double>,OHMMS_DIM> gradu = StorageGradVector[j];
      complex<double> laplu = StorageLaplVector[j];
      PosType k = kPoints[j];
      TinyVector<complex<double>,OHMMS_DIM> ck;
      for (int n=0; n<OHMMS_DIM; n++)	  ck[n] = k[n];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      StorageValueVector[j] = e_mikr*u;
      StorageGradVector[j]  = e_mikr*(-eye*u*ck + gradu);
      StorageLaplVector[j]  = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
    }

    // Check muffin tins;  if inside evaluate the orbitals
    bool inTin = false;
    bool need2blend = false;
    PosType disp;
    double b, db, d2b;
    for (int tin=0; tin<MuffinTins.size(); tin++) {
      MuffinTins[tin].inside(r, inTin, need2blend);
      if (inTin) {
	MuffinTins[tin].evaluate (r, StorageValueVector, StorageGradVector, StorageLaplVector);
	if (need2blend) {
	  disp = MuffinTins[tin].disp(r);
	  double dr = std::sqrt(dot(disp, disp));
	  MuffinTins[tin].blend_func(dr, b, db, d2b);
	}
	break;
      }
    }

    StorageValueVector_t &valVec =  
      need2blend ? BlendValueVector : StorageValueVector;
    StorageGradVector_t &gradVec =  
      need2blend ? BlendGradVector : StorageGradVector;
    StorageValueVector_t &laplVec =  
      need2blend ? BlendLaplVector : StorageLaplVector;

    // Otherwise, evaluate the B-splines
    if (!inTin || need2blend) {
      PosType ru(PrimLattice.toUnit(P.R[iat]));
      for (int i=0; i<OHMMS_DIM; i++)
	ru[i] -= std::floor (ru[i]);
      EinsplineTimer.start();
      EinsplineMultiEval (MultiSpline, ru, valVec, gradVec, StorageHessVector);
      EinsplineTimer.stop();
      for (int j=0; j<NumValenceOrbs; j++) {
	gradVec[j] = dot (PrimLattice.G, gradVec[j]);
	laplVec[j] = trace (StorageHessVector[j], GGt);
      }
      // Add e^-ikr phase to B-spline orbitals
      for (int j=0; j<NumValenceOrbs; j++) {
	complex<double> u = valVec[j];
	TinyVector<complex<double>,OHMMS_DIM> gradu = gradVec[j];
	complex<double> laplu = laplVec[j];
	PosType k = kPoints[j];
	TinyVector<complex<double>,OHMMS_DIM> ck;
	for (int n=0; n<OHMMS_DIM; n++)	  ck[n] = k[n];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
	valVec[j]   = e_mikr*u;
	gradVec[j]  = e_mikr*(-eye*u*ck + gradu);
	laplVec[j]  = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
      }
    }

    // Finally, copy into output vectors
    int psiIndex = 0;
    int N = StorageValueVector.size();
    if (need2blend) {
      for (int j=0; j<NumValenceOrbs; j++) {
	complex<double> psi_val, psi_lapl;
	TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	PosType rhat = 1.0/std::sqrt(dot(disp,disp)) * disp;
	complex<double> psi1 = StorageValueVector[j];
	complex<double> psi2 =   BlendValueVector[j];
	TinyVector<complex<double>,OHMMS_DIM> dpsi1 = StorageGradVector[j];
	TinyVector<complex<double>,OHMMS_DIM> dpsi2 = BlendGradVector[j];
	complex<double> d2psi1 = StorageLaplVector[j];
	complex<double> d2psi2 =   BlendLaplVector[j];
	
	TinyVector<complex<double>,OHMMS_DIM> zrhat;
	for (int i=0; i<OHMMS_DIM; i++)
	  zrhat[i] = rhat[i];

	psi_val  = b * psi1 + (1.0-b)*psi2;
	psi_grad = b * dpsi1 + (1.0-b)*dpsi2 + db * (psi1 - psi2)* zrhat;
	psi_lapl = b * d2psi1 + (1.0-b)*d2psi2 +
	  2.0*db * (dot(zrhat,dpsi1) - dot(zrhat, dpsi2)) +
	  d2b * (psi1 - psi2);
	
	psi[psiIndex] = real(psi_val);
	for (int n=0; n<OHMMS_DIM; n++)
	  dpsi[psiIndex][n] = real(psi_grad[n]);
	d2psi[psiIndex] = real(psi_lapl);
	psiIndex++;
	if (MakeTwoCopies[j]) {
	  psi[psiIndex] = imag(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi[psiIndex][n] = imag(psi_grad[n]);
	  d2psi[psiIndex] = imag(psi_lapl);
	  psiIndex++;
	}
      } 
      for (int j=NumValenceOrbs; j<N; j++) {
	complex<double> psi_val, psi_lapl;
	TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	psi_val  = StorageValueVector[j];
	psi_grad = StorageGradVector[j];
	psi_lapl = StorageLaplVector[j];
	
	psi[psiIndex] = real(psi_val);
	for (int n=0; n<OHMMS_DIM; n++)
	  dpsi[psiIndex][n] = real(psi_grad[n]);
	d2psi[psiIndex] = real(psi_lapl);
	psiIndex++;
	if (MakeTwoCopies[j]) {
	  psi[psiIndex] = imag(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi[psiIndex][n] = imag(psi_grad[n]);
	  d2psi[psiIndex] = imag(psi_lapl);
	  psiIndex++;
	}
      }
    }
    else {
      for (int j=0; j<N; j++) {
	complex<double> psi_val, psi_lapl;
	TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	psi_val  = StorageValueVector[j];
	psi_grad = StorageGradVector[j];
	psi_lapl = StorageLaplVector[j];
	
	psi[psiIndex] = real(psi_val);
	for (int n=0; n<OHMMS_DIM; n++)
	  dpsi[psiIndex][n] = real(psi_grad[n]);
	d2psi[psiIndex] = real(psi_lapl);
	psiIndex++;
	if (MakeTwoCopies[j]) {
	  psi[psiIndex] = imag(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi[psiIndex][n] = imag(psi_grad[n]);
	  d2psi[psiIndex] = imag(psi_lapl);
	  psiIndex++;
	}
      }
    }
    VGLTimer.stop();
  }
  
  // Value, gradient, and laplacian
  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int iat, ComplexValueVector_t& psi, 
   ComplexGradVector_t& dpsi, ComplexValueVector_t& d2psi)
  {
    VGLTimer.start();
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
    for (int i=0; i<OHMMS_DIM; i++)
      ru[i] -= std::floor (ru[i]);

    EinsplineTimer.start();
    EinsplineMultiEval (MultiSpline, ru, StorageValueVector,
			StorageGradVector, StorageHessVector);
    EinsplineTimer.stop();
    //computePhaseFactors(r);
    complex<double> eye (0.0, 1.0);
    for (int j=0; j<psi.size(); j++) {
      complex<double> u, laplu;
      TinyVector<complex<double>, OHMMS_DIM> gradu;
      u = StorageValueVector[j];
      gradu = dot(PrimLattice.G, StorageGradVector[j]);
      laplu = trace(StorageHessVector[j], GGt);
      
      PosType k = kPoints[j];
      TinyVector<complex<double>,OHMMS_DIM> ck;
      for (int n=0; n<OHMMS_DIM; n++)	
	ck[n] = k[n];
      double s,c;
      double phase = -dot(r, k);
      sincos (phase, &s, &c);
      complex<double> e_mikr (c,s);
      convert(e_mikr * u, psi[j]);
      convertVec(e_mikr*(-eye*u*ck + gradu), dpsi[j]);
      convert(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu), d2psi[j]);
    }
    VGLTimer.stop();
  }
  
  
  template<> void
  EinsplineSetExtended<double>::evaluate
  (const ParticleSet& P, int iat, RealValueVector_t& psi, 
   RealGradVector_t& dpsi, RealValueVector_t& d2psi)
  {
    VGLTimer.start();
    PosType r (P.R[iat]);
    PosType ru(PrimLattice.toUnit(P.R[iat]));
 
    int sign=0;
    for (int i=0; i<OHMMS_DIM; i++) {
      RealType img = std::floor(ru[i]);
      ru[i] -= img;
      sign += HalfG[i] * (int)img;
    }

    EinsplineTimer.start();
    EinsplineMultiEval (MultiSpline, ru, psi, StorageGradVector, 
			StorageHessVector);
    EinsplineTimer.stop();
    if (sign & 1) 
      for (int j=0; j<psi.size(); j++) {
	psi[j] *= -1.0;
	StorageGradVector[j] *= -1.0;
	StorageHessVector[j] *= -1.0;
      }
    for (int i=0; i<psi.size(); i++) {
      dpsi[i]  = dot(PrimLattice.G, StorageGradVector[i]);
      d2psi[i] = trace(StorageHessVector[i], GGt);
    }
    VGLTimer.stop();
  }
  
  
  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int first, int last, RealValueMatrix_t& psi, 
   RealGradMatrix_t& dpsi, RealValueMatrix_t& d2psi)
  {
    complex<double> eye(0.0,1.0);
    VGLMatTimer.start();
    for (int iat=first,i=0; iat<last; iat++,i++) {
      PosType r (P.R[iat]);
      
      // Do core states first
      int icore = NumValenceOrbs;
      for (int tin=0; tin<MuffinTins.size(); tin++) {
	MuffinTins[tin].evaluateCore(r, StorageValueVector, StorageGradVector,
				     StorageLaplVector, icore);
	icore += MuffinTins[tin].get_num_core();
      }
      
      // Add phase to core orbitals
      for (int j=NumValenceOrbs; j<StorageValueVector.size(); j++) {
	complex<double> u = StorageValueVector[j];
	TinyVector<complex<double>,OHMMS_DIM> gradu = StorageGradVector[j];
	complex<double> laplu = StorageLaplVector[j];
	PosType k = kPoints[j];
	TinyVector<complex<double>,OHMMS_DIM> ck;
	for (int n=0; n<OHMMS_DIM; n++)	  ck[n] = k[n];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
	StorageValueVector[j] = e_mikr*u;
	StorageGradVector[j]  = e_mikr*(-eye*u*ck + gradu);
	StorageLaplVector[j]  = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
      }
      
      // Check if we are in the muffin tin;  if so, evaluate
      bool inTin = false, need2blend = false;
      PosType disp;
      double b, db, d2b;
      for (int tin=0; tin<MuffinTins.size(); tin++) {
	MuffinTins[tin].inside(r, inTin, need2blend);
	if (inTin) {
	  MuffinTins[tin].evaluate (r, StorageValueVector, 
				    StorageGradVector, StorageLaplVector);
	  if (need2blend) {
	    disp = MuffinTins[tin].disp(r);
	    double dr = std::sqrt(dot(disp, disp));
	    MuffinTins[tin].blend_func(dr, b, db, d2b);
	  }
	  break;
	}
      }
      
      StorageValueVector_t &valVec =  
	need2blend ? BlendValueVector : StorageValueVector;
      StorageGradVector_t &gradVec =  
	need2blend ? BlendGradVector : StorageGradVector;
      StorageValueVector_t &laplVec =  
	need2blend ? BlendLaplVector : StorageLaplVector;
      
      // Otherwise, evaluate the B-splines
      if (!inTin || need2blend) {
	PosType ru(PrimLattice.toUnit(P.R[iat]));
	for (int i=0; i<OHMMS_DIM; i++)
	  ru[i] -= std::floor (ru[i]);
	EinsplineTimer.start();
	EinsplineMultiEval (MultiSpline, ru, valVec, gradVec, StorageHessVector);
	EinsplineTimer.stop();
	for (int j=0; j<NumValenceOrbs; j++) {
	  gradVec[j] = dot (PrimLattice.G, gradVec[j]);
	  laplVec[j] = trace (StorageHessVector[j], GGt);
	}
	// Add e^-ikr phase to B-spline orbitals
	for (int j=0; j<NumValenceOrbs; j++) {
	  complex<double> u = valVec[j];
	  TinyVector<complex<double>,OHMMS_DIM> gradu = gradVec[j];
	  complex<double> laplu = laplVec[j];
	  PosType k = kPoints[j];
	  TinyVector<complex<double>,OHMMS_DIM> ck;
	  for (int n=0; n<OHMMS_DIM; n++)	  ck[n] = k[n];
	  double s,c;
	  double phase = -dot(r, k);
	  sincos (phase, &s, &c);
	  complex<double> e_mikr (c,s);
	  valVec[j]   = e_mikr*u;
	  gradVec[j]  = e_mikr*(-eye*u*ck + gradu);
	  laplVec[j]  = e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu);
	}
      }
      
      // Finally, copy into output vectors
      int psiIndex = 0;
      int N = StorageValueVector.size();
      if (need2blend) {
	for (int j=0; j<NumValenceOrbs; j++) {
	  complex<double> psi_val, psi_lapl;
	  TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	  PosType rhat = 1.0/std::sqrt(dot(disp,disp)) * disp;
	  complex<double> psi1 = StorageValueVector[j];
	  complex<double> psi2 =   BlendValueVector[j];
	  TinyVector<complex<double>,OHMMS_DIM> dpsi1 = StorageGradVector[j];
	  TinyVector<complex<double>,OHMMS_DIM> dpsi2 = BlendGradVector[j];
	  complex<double> d2psi1 = StorageLaplVector[j];
	  complex<double> d2psi2 =   BlendLaplVector[j];
	  
	  TinyVector<complex<double>,OHMMS_DIM> zrhat;
	  for (int n=0; n<OHMMS_DIM; n++)
	    zrhat[n] = rhat[n];
	  
	  psi_val  = b * psi1 + (1.0-b)*psi2;
	  psi_grad = b * dpsi1 + (1.0-b)*dpsi2 + db * (psi1 - psi2)* zrhat;
	  psi_lapl = b * d2psi1 + (1.0-b)*d2psi2 +
	    2.0*db * (dot(zrhat,dpsi1) - dot(zrhat, dpsi2)) +
	    d2b * (psi1 - psi2);
	  
	  psi(psiIndex,i) = real(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi(i,psiIndex)[n] = real(psi_grad[n]);
	  d2psi(i,psiIndex) = real(psi_lapl);
	  psiIndex++;
	  if (MakeTwoCopies[j]) {
	    psi(psiIndex,i) = imag(psi_val);
	    for (int n=0; n<OHMMS_DIM; n++)
	      dpsi(i,psiIndex)[n] = imag(psi_grad[n]);
	    d2psi(i,psiIndex) = imag(psi_lapl);
	    psiIndex++;
	  }
	} 
	// Copy core states
	for (int j=NumValenceOrbs; j<N; j++) {
	  complex<double> psi_val, psi_lapl;
	  TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	  psi_val  = StorageValueVector[j];
	  psi_grad = StorageGradVector[j];
	  psi_lapl = StorageLaplVector[j];
	  
	  psi(psiIndex,i) = real(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi(i,psiIndex)[n] = real(psi_grad[n]);
	  d2psi(i,psiIndex) = real(psi_lapl);
	  psiIndex++;
	  if (MakeTwoCopies[j]) {
	    psi(psiIndex,i) = imag(psi_val);
	    for (int n=0; n<OHMMS_DIM; n++)
	      dpsi(i,psiIndex)[n] = imag(psi_grad[n]);
	    d2psi(i,psiIndex) = imag(psi_lapl);
	    psiIndex++;
	  }
	}
      }
      else { // No blending needed
	for (int j=0; j<N; j++) {
	  complex<double> psi_val, psi_lapl;
	  TinyVector<complex<double>, OHMMS_DIM> psi_grad;
	  psi_val  = StorageValueVector[j];
	  psi_grad = StorageGradVector[j];
	  psi_lapl = StorageLaplVector[j];
	  
	  psi(psiIndex,i) = real(psi_val);
	  for (int n=0; n<OHMMS_DIM; n++)
	    dpsi(i,psiIndex)[n] = real(psi_grad[n]);
	  d2psi(i,psiIndex) = real(psi_lapl);
	  psiIndex++;
	  // if (psiIndex >= dpsi.cols()) {
	  //   cerr << "Error:  out of bounds writing in EinsplineSet::evalate.\n"
	  // 	 << "psiIndex = " << psiIndex << "  dpsi.cols() = " << dpsi.cols() << endl;
	  // }
	  if (MakeTwoCopies[j]) {
	    psi(psiIndex,i) = imag(psi_val);
	    for (int n=0; n<OHMMS_DIM; n++)
	      dpsi(i,psiIndex)[n] = imag(psi_grad[n]);
	    d2psi(i,psiIndex) = imag(psi_lapl);
	    psiIndex++;
	  }
	}
      }
    }
    VGLMatTimer.stop();
  }
  
  
  
  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::evaluate
  (const ParticleSet& P, int first, int last, ComplexValueMatrix_t& psi, 
   ComplexGradMatrix_t& dpsi, ComplexValueMatrix_t& d2psi)
  {
    VGLMatTimer.start();
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType r (P.R[iat]);
      PosType ru(PrimLattice.toUnit(P.R[iat]));
      for (int n=0; n<OHMMS_DIM; n++)
	ru[n] -= std::floor (ru[n]);
      EinsplineTimer.start();
      EinsplineMultiEval (MultiSpline, ru, StorageValueVector,
			  StorageGradVector, StorageHessVector);
      EinsplineTimer.stop();
      //computePhaseFactors(r);
      complex<double> eye (0.0, 1.0);
      for (int j=0; j<OrbitalSetSize; j++) {
	complex<double> u, laplu;
	TinyVector<complex<double>, OHMMS_DIM> gradu;
	u = StorageValueVector[j];
	gradu = dot(PrimLattice.G, StorageGradVector[j]);
	laplu = trace(StorageHessVector[j], GGt);
	
	PosType k = kPoints[j];
	TinyVector<complex<double>,OHMMS_DIM> ck;
	for (int n=0; n<OHMMS_DIM; n++)	
	  ck[n] = k[n];
	double s,c;
	double phase = -dot(r, k);
	sincos (phase, &s, &c);
	complex<double> e_mikr (c,s);
	convert(e_mikr * u, psi(j,i));
	convertVec(e_mikr*(-eye*u*ck + gradu), dpsi(i,j));
	convert(e_mikr*(-dot(k,k)*u - 2.0*eye*dot(ck,gradu) + laplu), d2psi(i,j));
      } 
    }
    VGLMatTimer.stop();
  }
  
  
  
  template<> void
  EinsplineSetExtended<double>::evaluate
  (const ParticleSet& P, int first, int last, RealValueMatrix_t& psi, 
   RealGradMatrix_t& dpsi, RealValueMatrix_t& d2psi)
  {
    VGLMatTimer.start();
    for(int iat=first,i=0; iat<last; iat++,i++) {
      PosType r (P.R[iat]);
      PosType ru(PrimLattice.toUnit(P.R[iat]));
      int sign=0;
      for (int n=0; n<OHMMS_DIM; n++) {
	RealType img = std::floor(ru[n]);
	ru[n] -= img;
	sign += HalfG[n] * (int)img;
      }
      for (int n=0; n<OHMMS_DIM; n++)
	ru[n] -= std::floor (ru[n]);
      EinsplineTimer.start();
      EinsplineMultiEval (MultiSpline, ru, StorageValueVector,
			  StorageGradVector, StorageHessVector);
      EinsplineTimer.stop();
      if (sign & 1) 
	for (int j=0; j<OrbitalSetSize; j++) {
	  StorageValueVector[j] *= -1.0;
	  StorageGradVector[j]  *= -1.0;
	  StorageHessVector[j]  *= -1.0;
	}
      complex<double> eye (0.0, 1.0);
      for (int j=0; j<OrbitalSetSize; j++) {
        psi(j,i)   = StorageValueVector[j];
	dpsi(i,j)  = dot(PrimLattice.G, StorageGradVector[j]);
	d2psi(i,j) = trace(StorageHessVector[j], GGt);
      }
    }
    VGLMatTimer.stop();
  }
  
  template<typename StorageType> string
  EinsplineSetExtended<StorageType>::Type()
  {
    return "EinsplineSetExtended";
  }


  template<typename StorageType> void
  EinsplineSetExtended<StorageType>::registerTimers()
  {
    ValueTimer.reset();
    VGLTimer.reset();
    VGLMatTimer.reset();
    EinsplineTimer.reset();
    TimerManager.addTimer (&ValueTimer);
    TimerManager.addTimer (&VGLTimer);
    TimerManager.addTimer (&VGLMatTimer);
    TimerManager.addTimer (&EinsplineTimer);
  }






  SPOSetBase*
  EinsplineSetLocal::makeClone() const 
  {
    return new EinsplineSetLocal(*this);
  }

  void
  EinsplineSetLocal::resetParameters(const opt_variables_type& active)
  {
  }

  template<typename StorageType> SPOSetBase*
  EinsplineSetExtended<StorageType>::makeClone() const
  {
    EinsplineSetExtended<StorageType> *clone = 
      new EinsplineSetExtended<StorageType> (*this);
    clone->registerTimers();
    return clone;
  }

  template class EinsplineSetExtended<complex<double> >;
  template class EinsplineSetExtended<        double  >;
}
