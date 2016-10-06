//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/** @file SplineTestBase.h
 *
 * Test code to check the correctness of the einspline
 */
#ifndef QMCPLUSPLUS_SPLINETESTBASE_H
#define QMCPLUSPLUS_SPLINETESTBASE_H
#include <Configuration.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Numerics/OneDimGridBase.h>
#include <SandBox/TestFunc.h>
#include <SandBox/einspline_benchmark.h>
#include <Estimators/accumulators.h>

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
  /** print data */
  bool verbose;

  /** default constructor */
  SplineTestBase():nx(10),ny(10),nz(10),num_splines(4),nsamples(8), verbose(false)
  { }

  /** constructor using command line */
  SplineTestBase(int argc, char** argv)
    :nx(10),ny(10),nz(10),num_splines(4),nsamples(8),verbose(false)
  {
    initialize(argc,argv);
  }

  /** initialize the parameters */
  inline int initialize(int argc, char** argv)
  {
    int opt;
    while((opt = getopt(argc, argv, "hvg:x:y:z:s:p:")) != -1)
    {
      switch(opt)
      {
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
      case 'v':
        verbose=true;
        break;
      }
    }
    app_log().setf(std::ios::scientific, std::ios::floatfield);
    app_log().precision(6);
    app_log() << "#einspline validation grid = " << nx << " " << ny << " " << nz << " num_splines = " << num_splines << " num_samples = " << nsamples << std::endl;
    app_log() << std::endl;
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
  infunc.push_back( 0.01,new basis_type(3,3,6));
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
    for(int ix=0; ix<nx; ix++)
    {
      real_type x(gridX(ix));
      for(int iy=0; iy<ny; iy++)
      {
        real_type y(gridY(iy));
        for(int iz=0; iz<nz; iz++)
        {
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

  inline real_type pow2(std::complex<real_type> x)
  {
    return x.real()*x.real()+x.imag()*x.imag();
  }

  /** return \f$ \sqrt{\sum_{i,j} |a_{ij}-b_{ij}|^2}\f$
   */
  inline real_type diff(const Tensor<value_type,3>& a, const Tensor<value_type,3>& b)
  {
    real_type res=0.0;
    for(int i=0; i<9; ++i)
      res+=pow2(a(i)-b(i));
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


  /** test the functions of einspline
   * Each test code will implement the function, e.g.
   * einspline_validate.cpp
   * einspline_filter.cpp
   */
  void test(std::vector<TinyVector<real_type,3> >& coord);
};
}


#endif
