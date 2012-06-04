/** @file einspline_validation.cpp
 *
 * Test code to check the correctness of the einspline
 */
#include <Configuration.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Numerics/OneDimGridBase.h>
#include <SandBox/TestFunc.h>
#include <SandBox/einspline_benchmark.h>
#include <Message/Communicate.h>
#include <Estimators/accumulators.h>
#include <mpi/collectives.h>
#include <getopt.h>

namespace qmcplusplus
{
  /** parameters for the test */
  struct SplineTestBase
  {
    /** grid size */
    int nx,ny,nz;
    /** number of splines */
    int num_splines;
    /** number of samples for testing */
    int nsamples;

    /** default constructor */
    SplineTestBase():nx(10),ny(10),nz(10),num_splines(4),nsamples(8)
    { }

    /** constructor using command line */
    SplineTestBase(int argc, char** argv)
      :nx(10),ny(10),nz(10),num_splines(4),nsamples(8)
    {
      initialize(argc,argv);
    }

    /** initialize the parameters */
    inline int initialize(int argc, char** argv)
    {
      int opt;
      while((opt = getopt(argc, argv, "hg:x:y:z:s:p:")) != -1) 
      {
        switch(opt) {
          case 'h':
            printf("[-g grid| -x grid_x -y grid_y -z grid_z | -s splines | -p particles\n");
            return 1;
          case 'g':
            nx=ny=nz=atoi(optarg);
            break;
          case 'x':
            nx=atoi(optarg);
            break;
          case 'y':
            ny=atoi(optarg);
            break;
          case 'z':
            nz=atoi(optarg);
            break;
          case 's':
            num_splines=atoi(optarg);
            break;
          case 'p':
            nsamples=atoi(optarg);
            break;
        }
      }

      app_log().setf(std::ios::scientific, std::ios::floatfield);
      app_log().precision(6);
      app_log() << "#einspline validation grid = " << nx << " " << ny << " " << nz << " num_splines = " << num_splines << " num_samples = " << nsamples <<endl;
      app_log() << endl;
      app_log().flush();
      return 0;
    }
  };

  /** initialize PW<T> function in TestFunc.h*/
  template<typename T>
    void initialize_pw(PW<T>& infunc)
    {
      typedef typename PW<T>::basis_type basis_type;
      infunc.push_back( 0.50,new basis_type(1,2,1));
      infunc.push_back( 0.45,new basis_type(0,1,0));
      infunc.push_back(-0.25,new basis_type(1,1,1));
      infunc.push_back( 0.10,new basis_type(2,4,2));
      //infunc.push_back( 0.01,new basis_type(3,3,6));
    }

  /**  test class for ENGT 
   * @tparam ENGT type of einspline engine
   */
  template<typename ENGT>
    struct SplineTest:public SplineTestBase
  {
    typedef typename einspline3d_benchmark<ENGT>::real_type real_type;
    typedef typename einspline3d_benchmark<ENGT>::value_type value_type;

    /** analytic function */
    PW<value_type> infunc;

    /** storage for the data to be splined*/
    Array<value_type,3> inData;

    /** constructor 
     * @param p parameter set 
     */
    SplineTest(const SplineTestBase& p):SplineTestBase(p)
    {
      inData.resize(nx,ny,nz);
      initialize_pw(infunc);
      find_coeff();
    }

    /** find bspline coefficient using infunc
     */
    inline void find_coeff()
    {
      real_type ri=0.0;
      real_type rf=1.0;
      typedef qmcplusplus::LinearGrid<real_type> GridType;
      GridType gridX, gridY, gridZ;
      gridX.set(ri,rf,nx+1);
      gridY.set(ri,rf,ny+1);
      gridZ.set(ri,rf,nz+1);

      //Write to an array
      for(int ix=0; ix<nx; ix++) {
        real_type x(gridX(ix));
        for(int iy=0; iy<ny; iy++) {
          real_type y(gridY(iy));
          for(int iz=0; iz<nz; iz++) {
            TinyVector<real_type,3> pos(x,y,gridZ(iz));
            infunc.v(pos,inData(ix,iy,iz));
          }
        }
      }
    }

    inline real_type pow2(real_type x)
    {
      return x*x;
    }

    inline real_type pow2(complex<real_type> x)
    {
      return x.real()*x.real()+x.imag()*x.imag();
    }

    /** return \f$ \sqrt{\sum_{i,j} |a_{ij}-b_{ij}|^2}\f$
     */
    inline real_type diff(const Tensor<value_type,3>& a, const Tensor<value_type,3>& b)
    {
      real_type res=0.0;
      for(int i=0; i<9; ++i) res+=pow2(a(i)-b(i));
      return std::sqrt(res);
    }

    /** return \f$ \sqrt{\sum_{i} |a_{i}-b_{i}|^2}\f$
     */
    inline real_type diff(const TinyVector<value_type,3>& a, const TinyVector<value_type,3>& b)
    {
      return std::sqrt(pow2(a[0]-b[0])+pow2(a[1]-b[1])+pow2(a[2]-b[2]));
    }

    /** return \f$ \sqrt{(a-b)^2}\f$
     */
    inline real_type diff(value_type& a, value_type b)
    {
      return std::sqrt(pow2(a-b));
    }


    /** test the functions of einspline */
    void test(vector<TinyVector<real_type,3> >& coord)
    {

      typedef TinyVector<real_type,3> pos_type;
      typedef TinyVector<value_type,3> grad_type;
      typedef Tensor<value_type,3> hess_type;

      Vector<value_type> psi(num_splines);
      Vector<value_type> lap(num_splines);
      Vector<TinyVector<value_type,3> > grad(num_splines);
      Vector<Tensor<value_type,3> > hess(num_splines);

      for(int k=0; k<hess.size(); ++k) hess[k]=hess_type(0.0);

      accumulator_set<real_type> dv,dvgl_v,dvgl_g,dvgl_l,dvgh_v,dvgh_g,dvgh_h;

      einspline3d_benchmark<ENGT> d_bench;
      d_bench.set(nx,ny,nz,num_splines,false);
      for(int i=0; i<num_splines; ++i) d_bench.assign(i,inData);

      for(int i=0; i<nsamples; ++i)
      {
        value_type v,l;
        grad_type g;
        hess_type h;

        //evaluate analytic solution
        infunc.v(coord[i],v);
        d_bench.einspliner.evaluate(coord[i],psi);
        dv(diff(v,psi[0]));

        infunc.vgl(coord[i],v,g,l);
        d_bench.einspliner.evaluate_vgl(coord[i],psi,grad,lap);
        dvgl_v(diff(v,psi[0]));
        dvgl_g(diff(g,grad[0]));
        dvgl_l(diff(l,lap[0]));

        //for(int k=0; k<hess.size(); ++k)
        //  hess[k]=hess_type(0.0);

        infunc.vgh(coord[i],v,g,h);
        d_bench.einspliner.evaluate_vgh(coord[i],psi,grad,hess);

        dvgh_v(diff(v,psi[0]));
        dvgh_g(diff(g,grad[0]));
        dvgh_h(diff(h,hess[0]));

        for(int k=1; k<hess.size(); ++k)
          if(diff(hess[0],hess[k]) > numeric_limits<real_type>::epsilon())
            app_log() << "Check multi evaluation = " << hess[0] << "\n"  << hess[k] << endl;
      }

      app_log() << "Average difference = " << dv.mean() << endl;
      app_log() << "Average difference VGL = " << dvgl_v.mean() << " " << dvgl_g.mean() << " " << dvgl_l.mean()  <<  endl;
      app_log() << "Average difference VGH = " << dvgh_v.mean() << " " << dvgh_g.mean() << " " << dvgh_h.mean()  <<  endl;
    }
  };
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

  vector<TinyVector<float,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i) convert(coord[i],coord_s[i]);

  cout << "\n Testing double " << endl;
  SplineTest<multi_UBspline_3d_d> test_d(param);
  test_d.test(coord);

  cout << "\n Testing complex<double> " << endl;
  SplineTest<multi_UBspline_3d_z> test_z(param);
  test_z.test(coord);

  cout << "\n Testing float " << endl;
  SplineTest<multi_UBspline_3d_s> test_s(param);
  test_s.test(coord_s);

  cout << "\n Testing complex<float> " << endl;
  SplineTest<multi_UBspline_3d_c> test_c(param);
  test_c.test(coord_s);
  return 0;
}
