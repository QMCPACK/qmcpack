//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file SplineC2XAdoptor.h
 *
 * Adoptor classes to handle complex-to-(real,complex) with arbitrary precision
 */
#ifndef QMCPLUSPLUS_EINSPLINE_C2X_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_C2X_ADOPTOR_H

#include <spline/einspline_engine.hpp>

namespace qmcplusplus {

  /** adoptor class to match complex<ST> spline stored in a packed array with complex<TT> SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineC2CAdoptorPacked
    {
      static const bool is_complex=true;

      typedef ST                                                      real_type;
      typedef complex<ST>                                             value_type;
      typedef typename einspline_traits<real_type,D>::SplineType      SplineType;
      typedef typename einspline_traits<real_type,D>::BCType          BCType;
      typedef typename OrbitalSetTraits<real_type>::ValueVector_t     StorageValueVector_t;
      typedef typename OrbitalSetTraits<real_type>::GradVector_t      StorageGradVector_t;
      typedef typename OrbitalSetTraits<real_type>::HessVector_t      StorageHessVector_t;
      typedef typename OrbitalSetTraits<real_type>::GradHessVector_t  StorageGradHessVector_t;

      typedef CrystalLattice<ST,D> UnitCellType;
      typedef TinyVector<ST,D>     PointType;

      SplineType          *MultiSpline;
      UnitCellType        SuperLattice;
      UnitCellType        PrimLattice;
      TinyVector<int,D>   HalfG;
      vector<bool>        MakeTwoCopies;
      Tensor<real_type,D> GGt;
      vector<PointType>   kPoints;

      vector<real_type> phase;
      vector<value_type> eikr;

      StorageValueVector_t     myV;
      StorageValueVector_t     myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      inline void resizeStorage(int n, int nvals)
      {
        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);

        kPoints.resize(n);
        MakeTwoCopies.resize(n);

        myV.resize(2*n);
        myL.resize(2*n);
        myG.resize(2*n);
        myH.resize(2*n);
      }

    template<typename GT, typename BCT>
      void create_spline(GT& xyz_g, BCT& xyz_bc)
      {
        MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());
      }

      void set_spline(int ival,  ST* restrict psi_r, ST* restrict psi_i)
      {
        einspline::set(MultiSpline, 2*ival, psi_r);
        einspline::set(MultiSpline, 2*ival+1, psi_i);
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate(MultiSpline,ru,myV);

          //computePhases(r);
          //simd::multiadd(psi.size(),eikr.data(),psi.data());
          register ST s,c;
          TT* restrict t_ptr=reinterpret_cast<TT*>(psi.data());
          for(int j=0,jr=0; j<psi.size(); ++j, jr+=2)
          {
            sincos(-dot(r,kPoints[j]),&s,&c);
            t_ptr[jr]  =c*myV[jr]-s*myV[jr+1];
            t_ptr[jr+1]=s*myV[jr]+c*myV[jr+1];
          }
        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);

          const int N=kPoints.size();
          for (int j=0; j<2*N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
          for (int j=0; j<2*N; j++) myL[j] = trace(myH[j],GGt);

          const ST two=2.0;
          ST s,c;
          TinyVector<ST,D> g_r, g_i;
          //can easily make three independent loops
          for (int j=0,jr=0,ji=1; j<N; j++,jr+=2,ji+=2) 
          {
            g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
            g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
            ST kk=-dot(kPoints[j],kPoints[j]);
            myL[jr]+=kk*myV[jr]+two*dot(kPoints[j],myG[ji]);
            myL[ji]+=kk*myV[ji]-two*dot(kPoints[j],myG[jr]);
            sincos(-dot(r,kPoints[j]),&s,&c); //e-ikr (beware of -1)
            psi[j]=complex<TT>(c*myV[jr]-s*myV[ji],c*myV[ji]+s*myV[jr]);
            d2psi[j]=complex<TT>(c*myL[jr]-s*myL[ji],c*myL[ji]+s*myL[jr]);
            for(int idim=0; idim<D; ++idim)
              dpsi[j][idim]=complex<TT>(c*g_r[idim]-s*g_i[idim], c*g_i[idim]+s*g_r[idim]);
            //complex<ST> e_mikr(c,s);
            //convert(e_mikr * myV[j], psi[j]);
            //convert(e_mikr*(-myV[j]*ck + myG[j]), dpsi[j]);
            //convert(e_mikr*(-myV[j]*kk - two*dot(ck,myG[j]) + myL[j]), d2psi[j]);
          }
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {
        }
    };

  /** adoptor class to match complex<ST> spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   *
   * Requires temporage storage and multiplication of phase vectors
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineC2RAdoptorPacked
    {
      static const bool is_complex=true;

      typedef ST                                                    real_type;
      typedef complex<ST>                                           value_type;

      typedef typename einspline_traits<real_type,D>::SplineType           SplineType;
      typedef typename einspline_traits<real_type,D>::BCType               BCType;
      typedef typename OrbitalSetTraits<real_type>::ValueVector_t          StorageValueVector_t;
      typedef typename OrbitalSetTraits<real_type>::GradVector_t           StorageGradVector_t;
      typedef typename OrbitalSetTraits<real_type>::HessVector_t           StorageHessVector_t;
      typedef typename OrbitalSetTraits<real_type>::GradHessVector_t       StorageGradHessVector_t;

      typedef CrystalLattice<ST,D> UnitCellType;
      typedef TinyVector<ST,D> PointType;

      SplineType          *MultiSpline;
      UnitCellType        SuperLattice;
      UnitCellType        PrimLattice;
      TinyVector<int,D>   HalfG;
      Tensor<real_type,D> GGt;
      vector<PointType>   kPoints;
      vector<bool>        MakeTwoCopies;
      vector<real_type>   CosV;
      vector<real_type>   SinV;
      vector<real_type>   mKK;

      // Temporary storage for Eispline calls
      StorageValueVector_t     myV;
      StorageValueVector_t     myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      SplineC2RAdoptorPacked():MultiSpline(0) { }

      virtual ~SplineC2RAdoptorPacked(){}

      inline void resizeStorage(int n, int nvals)
      {

        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);

        kPoints.resize(n);
        MakeTwoCopies.resize(n);
      
        myV.resize(2*n);
        myL.resize(2*n);
        myG.resize(2*n);
        myH.resize(2*n);
        CosV.resize(n);SinV.resize(n);
      }

    template<typename GT, typename BCT>
      void create_spline(GT& xyz_g, BCT& xyz_bc)
      {
        mKK.resize(kPoints.size());
        for(int i=0; i<kPoints.size(); ++i) mKK[i]=-dot(kPoints[i],kPoints[i]);
        MultiSpline=einspline::create(MultiSpline,xyz_g,xyz_bc,myV.size());
      }


      void set_spline(int ival,  ST* restrict psi_r, ST* restrict psi_i)
      {
        einspline::set(MultiSpline, 2*ival, psi_r);
        einspline::set(MultiSpline, 2*ival+1, psi_i);
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate(MultiSpline,ru,myV);

          int N=kPoints.size();
          
          ST phase[N];
          for(int j=0;j<N; ++j) phase[j]=-dot(r,kPoints[j]);
          eval_e2iphi(N,phase,CosV.data(),SinV.data());
          //for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);

          int psiIndex = 0;
          for (int j=0,jr=0; j<N; j++,jr+=2) 
          {
            psi[psiIndex] = static_cast<TT>(myV[jr]*CosV[j]-myV[jr+1]*SinV[j]);
            psiIndex++;
            if (MakeTwoCopies[j]) 
            {
              psi[psiIndex] = static_cast<TT>(myV[jr+1]*CosV[j]+myV[jr]*SinV[j]);
              psiIndex++;
            }
          }
        }

      template<typename VV, typename GV>
      inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);

          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);

          const int N=kPoints.size();
          for (int j=0; j<2*N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
          for (int j=0; j<2*N; j++) myL[j] = trace(myH[j],GGt);

          const ST zero=0.0;
          const ST two=2.0;
          for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);
          int psiIndex=0;
          TinyVector<ST,D> g_r, g_i;
          for (int j=0,jr=0,ji=1; j<N; j++,jr+=2,ji+=2) 
          {
            g_r=myG[jr]+myV[ji]*kPoints[j]; // \f$\nabla \psi_r + {\bf k}\psi_i\f$
            g_i=myG[ji]-myV[jr]*kPoints[j]; // \f$\nabla \psi_i - {\bf k}\psi_r\f$
            myL[jr]+=mKK[j]*myV[jr]+two*dot(kPoints[j],myG[ji]);
            myL[ji]+=mKK[j]*myV[ji]-two*dot(kPoints[j],myG[jr]);

            psi[psiIndex]=CosV[j]*myV[jr]-SinV[j]*myV[ji];
            dpsi[psiIndex]=CosV[j]*g_r-SinV[j]*g_i; // multiply phase 
            d2psi[psiIndex]=CosV[j]*myL[jr]-SinV[j]*myL[ji];

            ++psiIndex;
            if(MakeTwoCopies[j])
            {
              psi[psiIndex]=CosV[j]*myV[ji]+SinV[j]*myV[jr];
              dpsi[psiIndex]=CosV[j]*g_i+SinV[j]*g_r;
              d2psi[psiIndex]=CosV[j]*myL[ji]+SinV[j]*myL[jr];
              ++psiIndex;
            }
          }
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {
        }
    };

  //inline void computePhases(const PointType& r)
  //{
  //  for (int i=0; i<kPoints.size(); i++) phase[i] = -dot(r, kPoints[i]);
  //  eval_e2iphi(kPoints.size(),phase.data(),eikr.data());
  //}

  /** adoptor class to match complex<ST> spline with complex<TT> SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   *
   * This is a reference implementation but NOT USED. Use Packed version
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineC2CAdoptor
    {
      typedef ST                                                       real_type;
      typedef complex<ST>                                              value_type;
      typedef typename einspline_traits<value_type,D>::SplineType      SplineType;
      typedef typename einspline_traits<value_type,D>::BCType          BCType;
      typedef typename OrbitalSetTraits<value_type>::ValueVector_t     StorageValueVector_t;
      typedef typename OrbitalSetTraits<value_type>::GradVector_t      StorageGradVector_t;
      typedef typename OrbitalSetTraits<value_type>::HessVector_t      StorageHessVector_t;
      typedef typename OrbitalSetTraits<value_type>::GradHessVector_t  StorageGradHessVector_t;

      typedef CrystalLattice<ST,D> UnitCellType;
      typedef TinyVector<ST,D>     PointType;

      SplineType          *MultiSpline;
      UnitCellType        SuperLattice;
      UnitCellType        PrimLattice;
      vector<PointType>   kPoints;
      TinyVector<int,D>   HalfG;
      vector<bool>        MakeTwoCopies;
      Tensor<real_type,D> GGt;

      vector<real_type> phase;
      vector<value_type> eikr;

      StorageValueVector_t     myV;
      StorageValueVector_t     myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      inline void resizeStorage(int n, int nvals)
      {
        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);

        kPoints.resize(n);
        MakeTwoCopies.resize(n);
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
        myH.resize(n);
        myH.resize(n);
      }

      inline bool isready()
      {
        return true;
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate(MultiSpline,ru,myV);

          const int N=myV.size();
          //computePhases(r);
          //simd::multiadd(psi.size(),eikr.data(),psi.data());
          register ST s,c;
          for(int j=0; j<psi.size(); ++j)
          {
            sincos(-dot(r,kPoints[j]),&s,&c);
            psi[j]=complex<TT>(static_cast<TT>(c*myV[j].real()-s*myV[j].imag()),
                static_cast<TT>(c*myV[j].imag()+s*myV[j].real()));
          }

        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);

          const int N=kPoints.size();
          for (int j=0; j<N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
          for (int j=0; j<N; j++) myL[j] = trace(myH[j],GGt);

          const ST two=2.0;
          TinyVector<complex<ST>,D> ck;
          register ST s,c;
          for (int j=0; j<psi.size(); j++) 
          {
            ST kk=dot(kPoints[j],kPoints[j]);
            for(int i=0; i<D; ++i) ck[i]=complex<ST>(0.0,kPoints[j][i]);
            sincos(-dot(r,kPoints[j]),&s,&c);
            complex<ST> e_mikr(c,s);
            convert(e_mikr * myV[j], psi[j]);
            convert(e_mikr*(-myV[j]*ck + myG[j]), dpsi[j]);
            convert(e_mikr*(-myV[j]*kk - two*dot(ck,myG[j]) + myL[j]), d2psi[j]);
          }
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {
        }
    };

  /** adoptor class to match complex<ST> spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   *
   * Requires temporage storage and multiplication of phase vectors
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineC2RAdoptor
    {
      typedef ST                                                    real_type;
      typedef complex<ST>                                           value_type;

      typedef typename einspline_traits<value_type,D>::SplineType  SplineType;
      typedef typename einspline_traits<value_type,D>::BCType      BCType;
      typedef typename OrbitalSetTraits<value_type>::ValueVector_t          StorageValueVector_t;
      typedef typename OrbitalSetTraits<value_type>::GradVector_t           StorageGradVector_t;
      typedef typename OrbitalSetTraits<value_type>::HessVector_t           StorageHessVector_t;
      typedef typename OrbitalSetTraits<value_type>::GradHessVector_t       StorageGradHessVector_t;

      typedef CrystalLattice<ST,D> UnitCellType;
      typedef TinyVector<ST,D> PointType;

      SplineType*         MultiSpline;
      UnitCellType        SuperLattice;
      UnitCellType        PrimLattice;
      TinyVector<int,D>   HalfG;
      Tensor<real_type,D> GGt;
      vector<PointType>   kPoints;
      vector<bool>        MakeTwoCopies;
      vector<real_type>   CosV;
      vector<real_type>   SinV;
      vector<value_type>  mKK;

      // Temporary storage for Eispline calls
      StorageValueVector_t     myV;
      StorageValueVector_t     myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      SplineC2RAdoptor():MultiSpline(0) { }

      virtual ~SplineC2RAdoptor(){}

      inline void resizeStorage(int n, int nvals)
      {

        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);

        kPoints.resize(n);
        MakeTwoCopies.resize(n);
      
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
        myH.resize(n);
        CosV.resize(n);
        SinV.resize(n);
      }

      inline bool isready()
      {
        mKK.resize(kPoints.size());
        for(int i=0; i<kPoints.size(); ++i) mKK[i]=-dot(kPoints[i],kPoints[i]);
        return true;
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);
          einspline::evaluate(MultiSpline,ru,myV);
          const int N=kPoints.size();
          for(int j=0; j<N; ++j) sincos(-dot(r,kPoints[j]),&SinV[j],&CosV[j]);

          int psiIndex = 0;
          for (int j=0; j<N; j++) 
          {
            psi[psiIndex] = static_cast<TT>(myV[j].real()*CosV[j]-myV[j].imag()*SinV[j]);
            //psi[psiIndex] = static_cast<TT>(myV[j].real());
            psiIndex++;
            if (MakeTwoCopies[j]) 
            {
              //psi[psiIndex] = static_cast<TT>(myV[j].imag());
              psi[psiIndex] = static_cast<TT>(myV[j].imag()*CosV[j]+myV[j].real()*SinV[j]);
              psiIndex++;
            }
          }
        }

      template<typename VV, typename GV>
      inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          for (int i=0; i<D; i++) ru[i] -= std::floor (ru[i]);

          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
          const int N=kPoints.size();
          for (int j=0; j<N; j++) myG[j] = dot(PrimLattice.G, myG[j]);
          for (int j=0; j<N; j++) myL[j] = trace(myH[j],GGt);

          const ST zero=0.0;
          const ST two=2.0;
          TinyVector<complex<ST>,D> ck;
          ST s,c;
          for(int j=0; j<N; ++j)
          {
            for (int n=0; n<D; n++)ck[n] = complex<ST>(zero,-kPoints[j][n]);
            sincos (-dot(r,kPoints[j]), &s, &c);
            complex<ST> e_mikr(c,s);
            myV[j]  = e_mikr*myV[j];
            myL[j]  = -dot(kPoints[j],kPoints[j])*myV[j]+e_mikr*(two*dot(ck,myG[j]) + myL[j]);
            myG[j]  = myV[j]*ck+e_mikr*myG[j];
          }

          //const complex<ST> eye (0.0, 1.0);
          //ST s,c;
          //for(int j=0; j<N; ++j)
          //{
          //  complex<ST> u = myV[j];
          //  TinyVector<complex<ST>,D> gradu = myG[j];
          //  complex<ST> laplu = myL[j];
          //  PointType k = kPoints[j];
          //  TinyVector<complex<ST>,D> ck;
          //  for (int n=0; n<D; n++)	  ck[n] = k[n];
          //  sincos (-dot(r,k), &s, &c);
          //  complex<ST> e_mikr (c,s);
          //  myV[j]  = e_mikr*u;
          //  myG[j]  = e_mikr*(-eye*u*ck + gradu);
          //  myL[j]  = e_mikr*(-dot(k,k)*u - two*eye*dot(ck,gradu) + laplu);
          //}

          int psiIndex=0;
          for(int j=0; j<N; ++j)
          {
            psi[psiIndex]=static_cast<TT>(myV[j].real());
            for(int n=0; n<D; ++n) dpsi[psiIndex][n]=static_cast<TT>(myG[j][n].real());
            d2psi[psiIndex]=static_cast<TT>(myL[j].real());
            ++psiIndex;
            if(MakeTwoCopies[j])
            {
              psi[psiIndex]=static_cast<TT>(myV[j].imag());
              for(int n=0; n<D; ++n) dpsi[psiIndex][n]=static_cast<TT>(myG[j][n].imag());
              d2psi[psiIndex]=static_cast<TT>(myL[j].imag());
              ++psiIndex;
            }
          }
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {
        }
    };

  template<typename SA>
    struct SplineC2XAdoptorReader: BsplineReaderBase
  {
    typedef SA adoptor_type;
    typedef typename adoptor_type::real_type  spline_real_type;
    typedef typename adoptor_type::value_type spline_value_type;
    typedef typename adoptor_type::SplineType SplineType;

    SplineC2XAdoptorReader(EinsplineSetBuilder* e): BsplineReaderBase(e)
    {}

    SPOSetBase* create_spline_set(int spin, EinsplineSet* orbitalSet)
    {

      ReportEngine PRE("SplineC2XAdoptorREader","create_spline_set(spin,SPE*)");

      Timer c_prep, c_unpack,c_fft, c_phase, c_spline, c_newphase, c_h5, c_init;
      double t_prep=0.0, t_unpack=0.0, t_fft=0.0, t_phase=0.0, t_spline=0.0, t_newphase=0.0, t_h5=0.0, t_init=0.0;

      BsplineSet<adoptor_type>* bspline=new BsplineSet<adoptor_type>;
      init(orbitalSet,bspline);

      int N = mybuilder->NumDistinctOrbitals;
      int NumValenceOrbs = mybuilder->NumValenceOrbs;

      bspline->resizeStorage(N,NumValenceOrbs);
      int numOrbs=bspline->getOrbitalSetSize();

      // Read in k-points
      //int numOrbs = bspline->getOrbitalSetSize();
      int num = 0;
      const std::vector<BandInfo>& SortBands(mybuilder->SortBands);
      for (int iorb=0; iorb<N; iorb++) 
      {
        int ti = SortBands[iorb].TwistIndex;
        bspline->kPoints[iorb] = orbitalSet->PrimLattice.k_cart(mybuilder->TwistAngles[ti]); //twist);
        bspline->MakeTwoCopies[iorb] = 
          (num < (numOrbs-1)) && SortBands[iorb].MakeTwoCopies;
        num += bspline->MakeTwoCopies[iorb] ? 2 : 1;
      }

      // First, check to see if we have already read this in
      string H5FileName(mybuilder->H5FileName);
      H5OrbSet set(H5FileName, spin, N);

      bool havePsig=mybuilder->ReadGvectors_ESHDF();

      TinyVector<int,3> MeshSize=mybuilder->MeshSize;
      app_log() << "MeshSize = (" << MeshSize[0] << ", " << MeshSize[1] << ", " << MeshSize[2] << ")\n";

      int nx, ny, nz, bi, ti;
      nx=MeshSize[0]; ny=MeshSize[1]; nz=MeshSize[2];

      Ugrid xyz_grid[3];
      typename adoptor_type::BCType xyz_bc[3];
      xyz_grid[0].start = 0.0;  xyz_grid[0].end = 1.0;  xyz_grid[0].num = nx;
      xyz_grid[1].start = 0.0;  xyz_grid[1].end = 1.0;  xyz_grid[1].num = ny;
      xyz_grid[2].start = 0.0;  xyz_grid[2].end = 1.0;  xyz_grid[2].num = nz;

      xyz_bc[0].lCode=PERIODIC; xyz_bc[0].rCode=PERIODIC;
      xyz_bc[1].lCode=PERIODIC; xyz_bc[1].rCode=PERIODIC;
      xyz_bc[2].lCode=PERIODIC; xyz_bc[2].rCode=PERIODIC;

      bspline->create_spline(xyz_grid,xyz_bc);

      int TwistNum = mybuilder->TwistNum;

      cout << "WHAT IS GOING ON " << myComm->rank() << " " << TwistNum << endl;

      //app_log() << "  TEST TWIST " << TargetPtcl.getTwist() << endl;
      //string splinefile=make_spline_filename(H5FileName,TwistNum,MeshSize);
      string splinefile=make_spline_filename(H5FileName,mybuilder->TileMatrix,spin,TwistNum,MeshSize);

      bool root=(myComm->rank() == 0);

      int foundspline=0;
      Timer now;
      if(root)
      {
        hdf_archive h5f;
        foundspline=h5f.open(splinefile,H5F_ACC_RDONLY);
        if(foundspline)
        {
          einspline_engine<SplineType> bigtable(bspline->MultiSpline);
          foundspline=h5f.read(bigtable,"spline_0");
        }
      }
      myComm->bcast(foundspline);
      t_h5 = now.elapsed();

      if(foundspline)
      {
        app_log() << "Use existing bspline tables in " << splinefile << endl;
        chunked_bcast(myComm, bspline->MultiSpline); 
        t_init+=now.elapsed();
      }
      else
      {

        /** For valence orbitals,
         * - extended orbitals either in G or in R
         * - localized orbitals
         */
        Array<spline_real_type,3> splineData_r(nx,ny,nz),splineData_i(ny,ny,nz);

        if(havePsig)//perform FFT using FFTW
        {
          c_init.restart();
          Array<complex<double>,3> FFTbox;
          FFTbox.resize(MeshSize[0], MeshSize[1], MeshSize[2]);
          fftw_plan FFTplan = fftw_plan_dft_3d 
            (MeshSize[0], MeshSize[1], MeshSize[2],
             reinterpret_cast<fftw_complex*>(FFTbox.data()),
             reinterpret_cast<fftw_complex*>(FFTbox.data()),
             +1, FFTW_ESTIMATE);

          Vector<complex<double> > cG(mybuilder->MaxNumGvecs);

          //this will be parallelized with OpenMP
          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
          {
            int ti=SortBands[iorb].TwistIndex;
            get_psi_g(ti,spin,SortBands[iorb].BandIndex,cG);

            c_unpack.restart();
            unpack4fftw(cG,mybuilder->Gvecs[ti],MeshSize,FFTbox);
            t_unpack+= c_unpack.elapsed();

            c_fft.restart();
            fftw_execute (FFTplan);
            t_fft+= c_fft.elapsed();

            c_phase.restart();
            fix_phase_rotate_c2c(FFTbox,splineData_r, splineData_i,mybuilder->TwistAngles[ti]);
            t_phase+= c_phase.elapsed();

            c_spline.restart();
            bspline->set_spline(ival,splineData_r.data(),splineData_i.data());
            t_spline+= c_spline.elapsed();
          }

          fftw_destroy_plan(FFTplan);
          t_init+=c_init.elapsed();
        }
        else
        {
          Array<complex<double>,3> rawData(nx,ny,nz);
          //this will be parallelized with OpenMP
          for(int iorb=0,ival=0; iorb<N; ++iorb, ++ival)
          {
            //check dimension
            if(root)
            {
              string path=psi_r_path(SortBands[iorb].TwistIndex,spin,SortBands[iorb].BandIndex);
              HDFAttribIO<Array<complex<double>,3> >  h_splineData(rawData);
              h_splineData.read(mybuilder->H5FileID, path.c_str());
              simd::copy(splineData_r.data(),splineData_i.data(),rawData.data(),rawData.size());
            }
            myComm->bcast(splineData_r);
            myComm->bcast(splineData_i);
            bspline->set_spline(ival,splineData_r.data(),splineData_i.data());
          }
        }

        if(root)
        {
          hdf_archive h5f;
          h5f.create(splinefile);
          einspline_engine<SplineType> bigtable(bspline->MultiSpline);
          string aname("EinsplineC2XAdoptor");
          h5f.write(aname,"bspline_type");
          h5f.write(bigtable,"spline_0");
        }
      }

      app_log() << "    READBANDS::PREP   = " << t_prep << endl;
      app_log() << "    READBANDS::H5     = " << t_h5 << endl;
      app_log() << "    READBANDS::UNPACK = " << t_unpack << endl;
      app_log() << "    READBANDS::FFT    = " << t_fft << endl;
      app_log() << "    READBANDS::PHASE  = " << t_phase << endl;
      app_log() << "    READBANDS::SPLINE = " << t_spline << endl;
      app_log() << "    READBANDS::SUM    = " << t_init << endl;

      return bspline;
    }
  };

}
#endif
