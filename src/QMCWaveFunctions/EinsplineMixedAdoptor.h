//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
/** @file EsinplineMixedAdoptor.h
 *
 * SplineAdoptor to save memory for the vacuum for 3D systems
 *
 * A template parameter for BsplineSet<SplineAdoptor> 
 * - EinsplineMixedAdoptor: for slabs or wires
 * - EinsplineOpenAdoptor : for molecular systems with an orthorhombic cell
 */
#ifndef QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_MIXED_ADOPTOR_H

#include <einspline/multi_bspline_copy.h>
namespace qmcplusplus {

  template<typename GT>
    void print_grid(GT& grid)
    {
      std::cout << grid.start << " " << grid.end << " " 
        << grid.num << " " 
        << grid.delta <<  " " << grid.delta_inv << std::endl; 
    }

  template<typename ENGT>
    void print_spliner(ENGT* spline)
    {
      cout << "xgrid       "; print_grid(spline->x_grid);
      cout << "ygrid       "; print_grid(spline->y_grid);
      cout << "zgrid       "; print_grid(spline->z_grid);
    }

  //these will go to spline directory
  namespace einspline
  {
    /** interfaces to use UBspline_3d_X 
     *
     * - create
     * - set
     * - evaluate
     */
    template<typename VT, typename IT>
    UBspline_3d_d*  create(UBspline_3d_d* s
        , VT& start , VT& end, IT& ng , bc_code bc, int n=1)
    { 
      Ugrid x_grid, y_grid, z_grid;
      BCtype_d xBC,yBC,zBC;
      x_grid.start=start[0]; x_grid.end=end[0]; x_grid.num=ng[0];
      y_grid.start=start[1]; y_grid.end=end[1]; y_grid.num=ng[1];
      z_grid.start=start[2]; z_grid.end=end[2]; z_grid.num=ng[2];
      xBC.lCode=xBC.rCode=bc;
      yBC.lCode=yBC.rCode=bc;
      zBC.lCode=zBC.rCode=bc;
      return create_UBspline_3d_d(x_grid,y_grid,z_grid, xBC, yBC, zBC,NULL);
    }

    void  set(UBspline_3d_d* s, double* restrict data)
    { 
      recompute_UBspline_3d_d(s,data);
    }

    inline void  set(multi_UBspline_3d_d* spline, int i, UBspline_3d_d* spline_in
        , const int* offset, const int *N)
    { 
      copy_UBspline_3d_d(spline, i, spline_in,offset,N);
    }

    inline void  set(multi_UBspline_3d_s* spline, int i, UBspline_3d_d* spline_in
        , const int* offset, const int *N)
    { 
      copy_UBspline_3d_d_s(spline, i, spline_in,offset,N);
    }

    template<typename PT>
    inline double  evaluate(UBspline_3d_d *restrict spline, const PT& r)
    {
      double res;
      eval_UBspline_3d_d(spline,r[0],r[1],r[2],&res);
      return res;
    }

  }

  template <typename T>
  struct GridConvert
  {
    ///number of points in each direction including BC
    int N[3];
    ///offset
    int Offset[3];
    ///number of points of the original grid
    int BaseN[3];
    ///offset of the original grid, always 0
    int BaseOffset[3];
    ///spacing to be added or removed from the buffer
    TinyVector<T,3> Delta;

    template<typename ENGT1, typename ENGT2, typename PT>
      void create(ENGT1*& out, ENGT2* in, PT& lower, PT& upper, int num)
      {

        typedef typename bspline_engine_traits<ENGT1>::real_type real_type;

        Ugrid agrid[3];
        agrid[0]=in->x_grid;
        agrid[1]=in->y_grid;
        agrid[2]=in->z_grid;

        typename bspline_engine_traits<ENGT1>::BCType xyz_bc[3];
        xyz_bc[0].lCode=in->xBC.lCode;xyz_bc[0].rCode=in->xBC.rCode;
        xyz_bc[1].lCode=in->yBC.lCode;xyz_bc[1].rCode=in->yBC.rCode;
        xyz_bc[2].lCode=in->zBC.lCode;xyz_bc[2].rCode=in->zBC.rCode;
        
        for(int i=0; i<3; ++i)
        {
          int ngi=(int)(lower[i]*agrid[i].delta_inv);
          int ngf=(int)(upper[i]*agrid[i].delta_inv)+1;
          agrid[i].num  =std::min(agrid[i].num,ngf-ngi);
          ngf=agrid[i].num+ngi;
          agrid[i].start=static_cast<real_type>(ngi)*agrid[i].delta;
          agrid[i].end  =static_cast<real_type>(ngf)*agrid[i].delta;

          if (xyz_bc[i].lCode == PERIODIC || xyz_bc[i].lCode == ANTIPERIODIC)
            N[i] = agrid[i].num+3;
          else
            N[i] = agrid[i].num+2;

          Delta[i]=agrid[i].delta;
          Offset[i]=ngi;
        }
        out=einspline::create(out,agrid,xyz_bc,num);
      }
  };

  /** adoptor class to match ST real spline with TT real SPOs
   * @tparam ST precision of spline
   * @tparam TT precision of SPOs
   * @tparam D dimension
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineMixedAdoptor
    {
      typedef typename einspline_traits<ST,D>::SplineType SplineType;
      typedef typename einspline_traits<ST,D>::BCType     BCType;
      typedef typename einspline_traits<ST,D>::DataType   DataType;
      typedef CrystalLattice<ST,D>                        UnitCellType;
      typedef TinyVector<ST,D>                            PointType;
      typedef typename OrbitalSetTraits<ST>::ValueVector_t      StorageValueVector_t;
      typedef typename OrbitalSetTraits<ST>::GradVector_t       StorageGradVector_t;
      typedef typename OrbitalSetTraits<ST>::HessVector_t       StorageHessVector_t;
      typedef typename OrbitalSetTraits<ST>::GradHessVector_t   StorageGradHessVector_t;

      SplineType *MultiSpline;
      SplineType *smallBox;

      TinyVector<ST,D> Lower;
      TinyVector<ST,D> Upper;
      Tensor<ST,D> GGt;
      UnitCellType PrimLattice;
      GridConvert<ST> gTransform;
      //these are not needed for R2R adoptor and will be removed
      UnitCellType SuperLattice;

      ///used for testing only, to be removed
      bool UseFullGrid;

      // Temporary storage for Eispline calls
      StorageValueVector_t myV, myL;
      StorageGradVector_t      myG;
      StorageHessVector_t      myH;
      StorageGradHessVector_t  myGH;

      SplineMixedAdoptor(): MultiSpline(0), smallBox(0)
      {
      }

      void resizeStorage(int n, int nv)
      {
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
        myH.resize(n);
        myGH.resize(n);
      }

      void create_spline(TinyVector<int,D>& mesh, int n, bool samegrid)
      {
        UseFullGrid=samegrid;
        GGt=dot(transpose(PrimLattice.G),PrimLattice.G);
        Ugrid xyz_grid[D];
        BCType xyz_bc[D];
        for(int i=0; i<D; ++i)
        {
          xyz_grid[i].start = 0.0;  xyz_grid[i].end = 1.0;  
          xyz_grid[i].num = (samegrid)?mesh[i]:mesh[i]/2;
          xyz_bc[i].lCode=xyz_bc[i].rCode=PERIODIC;
          gTransform.BaseOffset[i]=0;
          gTransform.BaseN[i]=xyz_grid[i].num+3;
        }
        SplineType* dummy=0;
        MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
      }

      template <typename UBspline, typename PT>
      void add_box(UBspline* dense, PT& lower, PT& upper)
      {
        if(smallBox==0)
          gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
        Lower=lower+2.0*gTransform.Delta;
        Upper=upper-2.0*gTransform.Delta;
      }

      template <typename UBspline>
        void init_spline(UBspline* dense, UBspline* coarse, int ival)
        {
          if(UseFullGrid)
            einspline::set(MultiSpline,ival,dense, gTransform.BaseOffset,  gTransform.BaseN);
          else
            einspline::set(MultiSpline,ival,coarse, gTransform.BaseOffset,  gTransform.BaseN);

          einspline::set(smallBox,   ival,dense,  gTransform.Offset,gTransform.N);
        }

      inline bool isready()
      {
        return true;
      }

      ///** return sign */
      inline void convertPos(const PointType& r, PointType& ru, int& sign)
      {
        ru=PrimLattice.toUnit(r);
        for (int i=0; i<D; i++) ru[i]-=std::floor(ru[i]);
        sign=0; //need to add ANTIPERIODIC
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          int phase;
          TinyVector<ST,D> ru;
          convertPos(r,ru,phase);

          if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
            einspline::evaluate(smallBox,ru,myV);
          else
            einspline::evaluate(MultiSpline,ru,myV);

          for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(myV[j]);
        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          int phase;
          TinyVector<ST,D> ru;
          convertPos(r,ru,phase);

          if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
            einspline::evaluate_vgh(smallBox,ru,myV,myG,myH);
          else
            einspline::evaluate_vgh(MultiSpline,ru,myV,myG,myH);

          const int N=psi.size();
          const Tensor<ST,D> gConv(PrimLattice.G);
          for(int j=0; j<N; ++j) psi[j]=myV[j];
          for(int j=0; j<N; ++j) dpsi[j]=dot(gConv,myG[j]);
          for(int j=0; j<N; ++j) d2psi[j]=trace(myH[j],GGt);
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {}
    };

  /** adoptor class for the non-periodic systems
   *
   * Support two-level grid system for big systems
   * - MultiSpline : full cell with the coarse grid 
   * - smallBox    : small cell with the original grid
   * No coordinate transform is needed for the open BC with a orthorhombic cell
   */
  template<typename ST, typename TT, unsigned D>
    struct SplineOpenAdoptor
    {
      typedef typename einspline_traits<ST,D>::SplineType SplineType;
      typedef typename einspline_traits<ST,D>::BCType     BCType;
      typedef typename einspline_traits<ST,D>::DataType   DataType;
      typedef CrystalLattice<ST,D>                        UnitCellType;
      typedef TinyVector<ST,D>                            PointType;
      typedef typename OrbitalSetTraits<ST>::ValueVector_t      StorageValueVector_t;
      typedef typename OrbitalSetTraits<ST>::GradVector_t       StorageGradVector_t;
      typedef typename OrbitalSetTraits<ST>::HessVector_t       StorageHessVector_t;
      typedef typename OrbitalSetTraits<ST>::GradHessVector_t   StorageGradHessVector_t;

      SplineType *MultiSpline;
      SplineType *smallBox;

      TinyVector<ST,D> Lower;
      TinyVector<ST,D> Upper;
      TinyVector<ST,D> L;
      TinyVector<ST,D> InvL;

      Tensor<ST,D> GGt;
      UnitCellType PrimLattice;
      UnitCellType SuperLattice;

      ///used for testing only 
      bool UseFullGrid;
      ///grid transform
      GridConvert<ST> gTransform;

      // Temporary storage for Eispline calls
      StorageValueVector_t myV, myL;
      StorageGradVector_t      myG;
      StorageGradHessVector_t  myGH;

      SplineOpenAdoptor(): MultiSpline(0), smallBox(0)
      {
      }

      void resizeStorage(int n, int nv)
      {
        for(int i=0; i<D; ++i) L[i]=SuperLattice.R(i,i);
        for(int i=0; i<D; ++i) InvL[i]=1.0/L[i];
        myV.resize(n);
        myL.resize(n);
        myG.resize(n);
      }

      /** create MultiSpline for the full cell with a coarse grid
       * @param mesh original mesh (e.g. FFT)
       * @param n number of states
       * @param samegrid if true, use full grid
       */
      void create_spline(TinyVector<int,D>& mesh, int n, bool samegrid)
      {
        UseFullGrid=samegrid;
        Ugrid xyz_grid[D];
        BCType xyz_bc[D];
        for(int i=0; i<D; ++i)
        {
          L[i]=PrimLattice.R(i,i);
          InvL[i]=1.0/L[i];

          xyz_grid[i].start = 0.0;  xyz_grid[i].end = L[i]; //1.0;  
          xyz_grid[i].num = (samegrid)?mesh[i]:mesh[i]/2; 
          xyz_bc[i].lCode=xyz_bc[i].rCode=PERIODIC;
          gTransform.BaseOffset[i]=0;
          gTransform.BaseN[i]=xyz_grid[i].num+3;
        }
        SplineType* dummy=0;
        MultiSpline=einspline::create(dummy,xyz_grid,xyz_bc,n);
      }

      /** create smallBox
       * @param dense single-spline for the full grid
       * @param lower Cartesian lower bound
       * @param upper Cartesian upper bound
       */
      template <typename UBspline, typename PT>
      void add_box(UBspline* dense, PT& lower, PT& upper)
      {
        if(smallBox==0)
          gTransform.create(smallBox,dense,lower,upper,MultiSpline->num_splines);
        Lower=lower;//+2.0*gTransform.Delta;
        Upper=upper;//-2.0*gTransform.Delta;
      }

      /** set ival-th state 
       * @param dense a single-spline on the full grid
       * @param coarse a single-spline on the half grid
       * @param ival the state to be initialized
       */
      template <typename UBspline>
        void init_spline(UBspline* dense, UBspline* coarse, int ival)
        {
          if(UseFullGrid)
            einspline::set(MultiSpline,ival,dense,  gTransform.BaseOffset,gTransform.BaseN);
          else
            einspline::set(MultiSpline,ival,coarse, gTransform.BaseOffset,gTransform.BaseN);

          einspline::set(smallBox,ival,dense,gTransform.Offset,gTransform.N);
        }

      inline bool isready()
      {
        return true;
      }

      ///not needed for molecular system
      inline void convertPos(const PointType& r, PointType& ru)
      {
        for(int i=0; i<D; ++i)
          ru[i] = r[i]-L[i]*std::floor(r[i]*InvL[i]);
      }

      template<typename VV>
        inline void evaluate_v(const PointType& r, VV& psi)
        {
          TinyVector<ST,D> ru;
          convertPos(r,ru);
          if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
            einspline::evaluate(smallBox,ru,myV);
          else
            einspline::evaluate(MultiSpline,ru,myV);

          for (int j=0; j<psi.size(); j++) psi[j]=static_cast<TT>(myV[j]);
        }

      template<typename VV, typename GV>
        inline void evaluate_vgl(const PointType& r, VV& psi, GV& dpsi, VV& d2psi)
        {
          TinyVector<ST,D> ru;
          convertPos(r,ru);

          if(ru[0]>Lower[0] && ru[0]<Upper[0] && ru[1]>Lower[1] && ru[1]<Upper[1] && ru[2]>Lower[2] && ru[2]<Upper[2])
            einspline::evaluate_vgl(smallBox,ru,myV,myG,myL);
          else
            einspline::evaluate_vgl(MultiSpline,ru,myV,myG,myL);

          const int N=psi.size();
          for(int j=0; j<N; ++j) psi[j]=myV[j];
          for(int j=0; j<N; ++j) dpsi[j]=myG[j];
          for(int j=0; j<N; ++j) d2psi[j]=myL[j];
        }

      template<typename VV, typename GV, typename GGV>
        void evaluate_vgh(const PointType& r, VV& psi, GV& dpsi, GGV& grad_grad_psi)
        {}
    };

}
#endif
