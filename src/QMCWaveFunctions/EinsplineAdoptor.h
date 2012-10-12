//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file EsinplineAdoptor.h
 *
 * Adoptor classes and BsplineSet<SplineAdoptor> 
 * BsplineSet<SplineAdoptor> is a SPOSetBase class to work with determinant classes
 * SplineAdoptor provides these functions 
 * - evaluate_v    value only
 * - evaluate_vgl  vgl
 * - evaluate_vgh  vgh
 * Implemented 
 * - SplineC2RAdoptor<ST,TT,D> : real wavefunction using complex einspline, tiling
 * - SplineC2CAdoptor<ST,TT,D> : complex wavefunction using complex einspline, tiling
 * - SplineR2RAdoptor<ST,TT,D> : real wavefunction using real einspline, a single twist
 * where ST (TT) is the precision of the einspline (SPOSetBase).
 * 
 * typedefs and data members are duplicated for each adoptor class.
 * @todo Specalization and optimization for orthorhombic cells to use vgl not vgh
 */
#ifndef QMCPLUSPLUS_EINSPLINE_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_ADOPTOR_H

#include <spline/einspline_engine.hpp>
namespace qmcplusplus {

  /** einspline trait class equivalent to  MultiOrbitalTraits
   * @tparam ST spline datatype
   * @tparam D dimension
   * @tparam TT target datatype
   */
  template<typename ST, unsigned D>
    struct einspline_traits { };

  /** specialization for 3D double */
  template<>
    struct einspline_traits<double,3>
    {
      typedef multi_UBspline_3d_d SplineType;  
      typedef BCtype_d            BCType;
      typedef double              DataType;
    };

  /** specialization for 3D complex<double> */
  template<>
    struct einspline_traits<complex<double>,3>
    {
      typedef multi_UBspline_3d_z SplineType;  
      typedef BCtype_z            BCType;
      typedef complex<double>     DataType;
    };

  /** specialization for 3D float */
  template<>
    struct einspline_traits<float,3>
    {
      typedef multi_UBspline_3d_s SplineType;  
      typedef BCtype_s            BCType;
      typedef float               DataType;
    };

  /** specialization for 3D complex<float> */
  template<>
    struct einspline_traits<complex<float>,3>
    {
      typedef multi_UBspline_3d_c SplineType;  
      typedef BCtype_c            BCType;
      typedef complex<float>      DataType;
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

  /** adoptor class to match ST real spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineR2RAdoptor
    {
      typedef typename einspline_traits<ST,D>::SplineType SplineType;
      typedef typename einspline_traits<ST,D>::BCType     BCType;
      typedef typename einspline_traits<ST,D>::DataType   DataType;
      typedef CrystalLattice<ST,D>                        UnitCellType;
      typedef TinyVector<ST,D>                            PointType;

      SplineType *MultiSpline;
      Tensor<ST,D> GGt;
      UnitCellType PrimLattice;
      UnitCellType SuperLattice;
      TinyVector<int,D> HalfG;
      vector<PointType> kPoints;
      vector<bool>      MakeTwoCopies;

      typedef typename OrbitalSetTraits<ST>::ValueVector_t      StorageValueVector_t;
      typedef typename OrbitalSetTraits<ST>::GradVector_t       StorageGradVector_t;
      typedef typename OrbitalSetTraits<ST>::HessVector_t       StorageHessVector_t;
      typedef typename OrbitalSetTraits<ST>::GradHessVector_t   StorageGradHessVector_t;

      // Temporary storage for Eispline calls
      StorageValueVector_t myV, myL;
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
        myGH.resize(n);
      }

      inline bool isready()
      {
        return true;
      }
      template<typename VV>
        void evaluate_v(const PointType& r, VV& psi)
        {
          PointType ru(PrimLattice.toUnit(r));
          int sign=0;
          for (int i=0; i<D; i++) {
            ST img = std::floor(ru[i]);
            ru[i] -= img;
            sign += HalfG[i] * (int)img;
          }

          einspline::evaluate(MultiSpline,ru,myV);

          if (sign & 1) 
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(-myV[j]);
          else
            for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(myV[j]);
        }

      template<typename VV, typename GV>
        void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          PointType ru(PrimLattice.toUnit(r));

          int sign=0;
          for (int i=0; i<D; i++) {
            ST img = std::floor(ru[i]);
            ru[i] -= img;
            sign += HalfG[i] * (int)img;
          }

          einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);
          const ST minus_one=-1.0;
          if (sign & 1) 
            for (int j=0; j<psi.size(); j++) 
            {
              psi[j] = -myV[j];
              dpsi[j]  = minus_one*dot(PrimLattice.G, myG[j]);
              d2psi[j] = -trace(myH[j], GGt);
            }
          else
            for (int j=0; j<psi.size(); j++) 
            {
              psi[j] = myV[j];
              dpsi[j]  = dot(PrimLattice.G, myG[j]);
              d2psi[j] = trace(myH[j], GGt);
            }
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {}
    };

  /** a class to map a memory sequence to a vector
   * @tparam T datatype
   */
  template<typename T>
    struct VectorViewer
    {
      T* data_ptr;
      int data_size;
      inline VectorViewer(T* a, int n):data_ptr(a),data_size(n){}
      inline T* data() { return data_ptr;}
      inline int size() const { return data_size;}
      inline T& operator[](int i) { return data_ptr[i]; }
      inline T operator[](int i) const { return data_ptr[i]; }
    };

  /** BsplineSet 
   * @tparam SplineAdoptor implements evaluation functions that matched the storage requirements.
   *
   * Equivalent to EinsplineSetExtended<Storage>
   * Storage is now handled by SplineAdoptor class.
   */
  template<typename SplineAdoptor>
    struct BsplineSet: public SPOSetBase, public SplineAdoptor
  {
    typedef typename SplineAdoptor::SplineType SplineType;
    typedef typename SplineAdoptor::PointType  PointType;

    using SplineAdoptor::MultiSpline;

    /** default constructor */
    BsplineSet() { }

    /** allocate einspline 
     */
    template<typename GT, typename BCT>
      void allocate(GT& xyz_g, BCT& xyz_bc, int nv)
      {
        SplineType* dummy=0;
        MultiSpline=einspline::create(dummy,xyz_g,xyz_bc,nv);
      }

    SPOSetBase* makeClone() const
    {
      return new BsplineSet<SplineAdoptor>(*this);
    }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
    {
      SplineAdoptor::evaluate_v(P.R[iat],psi);
    }

    inline void evaluate(const ParticleSet& P, int iat, 
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      SplineAdoptor::evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
    }
    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
    {
      SplineAdoptor::evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
    }

    void resetParameters(const opt_variables_type& active)
    { }

    void resetTargetParticleSet(ParticleSet& e)
    { }

    void setOrbitalSetSize(int norbs)
    {
      OrbitalSetSize = norbs;
      BasisSetSize=norbs;
    }

    void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      typedef ValueMatrix_t::value_type value_type;
      typedef GradMatrix_t::value_type grad_type;
      //const int N=last-first;

      for(int iat=first, i=0; iat<last; ++iat,++i)
      {
        VectorViewer<value_type> v(logdet[i],OrbitalSetSize);
        VectorViewer<grad_type> g(dlogdet[i],OrbitalSetSize);
        VectorViewer<value_type> l(d2logdet[i],OrbitalSetSize);
        SplineAdoptor::evaluate_vgl(P.R[iat],v,g,l);
      }
    }
  };


}
#endif
