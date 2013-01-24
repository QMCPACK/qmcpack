/** @file einspline_validation.cpp
 *
 * Test code to check the correctness of the einspline
 */
#include <SandBox/SplineTestBase.h>
#include <einspline/multi_bspline_copy.h>
#include <Message/Communicate.h>
#include "Numerics/HDFNumericAttrib.h"
#include <mpi/collectives.h>
#include <getopt.h>

//testing subgrid
qmcplusplus::TinyVector<double,3> lower(0.25,0.25,0.25);
qmcplusplus::TinyVector<double,3> upper(0.75,0.75,0.75);

namespace qmcplusplus
{
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
      //cout << "num_splines " << spline->num_splines << endl;
      //cout << "coefs_size  " << spline->coefs_size << endl;
      //cout << "xstride     " << spline->x_stride << endl;
      //cout << "ystride     " << spline->y_stride << endl;
    }

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

  struct GridConvert
  {
    int N[3];
    int Offset[3];

    template<typename ENGT1, typename ENGT2, typename PT>
      void create(ENGT1*& out, ENGT2* in, PT& lower, PT& upper, int num)
      {

        typedef typename bspline_engine_traits<ENGT1>::real_type real_type;

        Ugrid agrid[3];
        agrid[0]=in->x_grid;
        agrid[1]=in->y_grid;
        agrid[2]=in->z_grid;

        for(int i=0; i<3; ++i)
        {
          int ngi=(int)(lower[i]*agrid[i].delta_inv);
          int ngf=(int)(upper[i]*agrid[i].delta_inv)+1;
          agrid[i].start=(real_type)(ngi)*agrid[i].delta;
          agrid[i].end  =std::min((real_type)(ngf)*agrid[i].delta,agrid[i].end);
          agrid[i].num=std::min(agrid[i].num,ngf-ngi);
          Offset[i]=ngi;
        }

        typename bspline_engine_traits<ENGT1>::BCType xyz_bc[3];
        xyz_bc[0].lCode=in->xBC.lCode; xyz_bc[0].rCode=in->xBC.rCode;
        xyz_bc[1].lCode=in->yBC.lCode; xyz_bc[1].rCode=in->yBC.rCode;
        xyz_bc[2].lCode=in->zBC.lCode; xyz_bc[2].rCode=in->zBC.rCode;

        if (in->xBC.lCode == PERIODIC || in->xBC.lCode == ANTIPERIODIC)
          N[0] = agrid[0].num+3;
        else
          N[0] = agrid[0].num+2;

        if (in->yBC.lCode == PERIODIC || in->yBC.lCode == ANTIPERIODIC)
          N[1] = agrid[1].num+3;
        else
          N[1] = agrid[1].num+2;

        if (in->zBC.lCode == PERIODIC || in->zBC.lCode == ANTIPERIODIC)
          N[2] = agrid[2].num+3;
        else
          N[2] = agrid[2].num+2;

        out=einspline::create(out,agrid,xyz_bc,num);
      }
  };

  template<typename PT>
  inline bool inbox(const PT& pos)
  {
    return (pos[0]>lower[0] && pos[0]<upper[0]
        &&  pos[1]>lower[1] && pos[1]<upper[0] 
        &&  pos[2]>lower[2] && pos[2]<upper[0] 
        );
  }

  template<typename ENGT>
    void SplineTest<ENGT>
    ::test(vector<TinyVector<typename SplineTest<ENGT>::real_type,3> >& coord)
    {
      //read data from a file
      int ngx=160;
      int ngy=160;
      int ngz=160;

      Array<double,3> big(ngx,ngy,ngz),small(ngx/2,ngy/2,ngz/2);
      {
        hid_t h5f = H5Fopen("sample.h5",H5F_ACC_RDONLY,H5P_DEFAULT);
        HDFAttribIO<Array<double,3> >  h_big(big);
        h_big.read(h5f, "/big_40");
        HDFAttribIO<Array<double,3> >  h_small(small);
        h_small.read(h5f, "/small_40");
        H5Fclose(h5f);
      }

      typedef TinyVector<real_type,3> pos_type;
      pos_type start(0.0);
      pos_type end(1.0);
      TinyVector<int,3> ng(ngx,ngy,ngz), ngh(ngx/2,ngy/2,ngz/2);


      //create/set UBspline_3d_X with the dense grid
      typename bspline_engine_traits<ENGT>::SingleSplineType* dense=0;
      dense=einspline::create(dense,start,end,ng,PERIODIC);
      einspline::set(dense,big.data());

      //create/set UBspline_3d_X with the coarse grid
      typename bspline_engine_traits<ENGT>::SingleSplineType* coarse=0;
      coarse=einspline::create(coarse,start,end,ngh,PERIODIC);
      einspline::set(coarse,small.data());

      //create/set multi_UBspline_3d as usual on the dense grid
      einspline_engine<ENGT> einspliner;
      einspliner.create(start,end,ng,PERIODIC,num_splines);
      for(int i=0; i<num_splines; ++i) einspliner.set(i,big);

      //create multi_UBspline_3d as usual on a truncated grid based on the dense grid
      einspline_engine<multi_UBspline_3d_d> newspline;
      GridConvert gconv;
      gconv.create(newspline.spliner,dense,lower,upper,num_splines);
      for(int i=0; i<num_splines; ++i) 
        einspline::set(newspline.spliner,i,dense,gconv.Offset,gconv.N);

      cout << "Original dense grid " << endl;
      print_spliner(einspliner.spliner);
      cout << "Truncated dense grid " << endl;
      print_spliner(newspline.spliner);


      Vector<value_type> psi(num_splines);
      Vector<real_type> psi_c(num_splines);
      accumulator_set<real_type> dv,ds,din,dout;


      int n_in=0;
      int n_out=0;
      for(int i=0; i<nsamples; ++i)
      {
        value_type v_dense, v_mixed;

        einspliner.evaluate(coord[i],psi);
        v_dense = einspline::evaluate(dense,coord[i]);
        ds(diff(v_dense,psi[0]));

        if(inbox(coord[i]))
        {
          n_in++;
          newspline.evaluate(coord[i],psi_c);
          v_mixed=psi_c[0];
          din(diff(v_mixed,v_dense));
        }
        else
        {
          n_out++;
          v_mixed = einspline::evaluate(coarse,coord[i]);
          dout(diff(v_mixed,v_dense));
        }

        dv(diff(v_dense,v_mixed));

        if(verbose)
        {
          cout << "At " << coord[i] << endl;
          cout << "Dense Grid " << v_dense
            <<  "\nMixed Grid " << v_mixed
            << endl << endl;
        }
      }

      app_log() << "Number of coords in the box       " << n_in << endl;
      app_log() << "Number of coords outside the box  " << n_out << endl;
      app_log() << "Average difference org-mixed grid " << dv.mean() 
        << "\n\tInside  " << din.mean()
        << "\n\tOutside " << dout.mean()
        << "\n\tmulti-spline:single-spline ="  << ds.mean() 
        << endl;
    }
}
int main(int argc, char** argv) 
{

  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,-1);

  SplineTestBase param(argc,argv);

  typedef TinyVector<double,3> pos_type;
  vector<pos_type> coord(param.nsamples);
  for(int i=0; i<coord.size(); ++i)
    coord[i]=pos_type(Random(),Random(),Random());
//  for(int i=0; i<coord.size(); ++i)
//    coord[i]=pos_type((upper[0]-lower[0])*Random()+lower[0]
//        , (upper[1]-lower[1])*Random()+lower[1]
//        , (upper[2]-lower[2])*Random()+lower[2]
//        );
//
  vector<TinyVector<float,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i) convert(coord[i],coord_s[i]);

  cout << "\nTesting einspline transformation" << endl;

  SplineTest<multi_UBspline_3d_d> test_d(param);
  test_d.test(coord);
  return 0;
}
