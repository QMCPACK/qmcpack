#include <Configuration.h>
#include <OhmmsPETE/OhmmsArray.h>
#include <Numerics/OneDimGridBase.h>
#include <SandBox/einspline_benchmark.h>
#include <SandBox/TestFunc.h>
#include <Message/Communicate.h>
#include <mpi/collectives.h>
#include <getopt.h>
using namespace qmcplusplus;


int main(int argc, char** argv) {

  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,-1);
  //Random.init(0,1,11);

  int nx=10,ny=10,nz=10;
  int num_splines=4;
  int nsamples=8;
  int niters=10;
  int opt;

  while((opt = getopt(argc, argv, "hg:x:y:z:i:s:p:")) != -1) {
    switch(opt) {
      case 'h':
        printf("[-g grid| -x grid_x -y grid_y -z grid_z] -s states -p particles -i iterations -t [d|s|z|c] \n");
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
      case 'i':
        niters=atoi(optarg);
        break;
    }
  }

  app_log().setf(std::ios::scientific, std::ios::floatfield);
  app_log().precision(6);
  app_log() << "#einspline benchmark grid = " << nx << " " << ny << " " << nz
    << " num_splines = " << num_splines << " num_samples = " << nsamples 
    <<endl;
  app_log().flush();

  //grid is [0,1)
  double ri=0.0;
  double rf=1.0;
  typedef TinyVector<double,3> pos_type;
  typedef TinyVector<double,3> grad_type;
  typedef Tensor<double,3> hess_type;


  //Create an analytic function for assignment
  ComboFunc<double> infunc;
  infunc.push_back(0.5,new TestFunc<double>(1,0,0));
  infunc.push_back(0.5,new TestFunc<double>(0,1,0));
  infunc.push_back(0.5,new TestFunc<double>(0,0,1));
  infunc.push_back(0.1,new TestFunc<double>(1,1,1));
  infunc.push_back(0.03,new TestFunc<double>(1,1,2));
  //infunc.push_back(0.01,new TestFunc(5,3,2));
  //infunc.push_back(0.01,new TestFunc(5,7,1));

  typedef qmcplusplus::LinearGrid<double> GridType;
  GridType gridX, gridY, gridZ;
  gridX.set(ri,rf,nx+1);
  gridY.set(ri,rf,ny+1);
  gridZ.set(ri,rf,nz+1);

  //Write to an array
  Array<double,3> inData(nx,ny,nz);
  for(int ix=0; ix<nx; ix++) {
    double x(gridX(ix));
    for(int iy=0; iy<ny; iy++) {
      double y(gridY(iy));
      for(int iz=0; iz<nz; iz++) {
        pos_type pos(x,y,gridZ(iz));
        inData(ix,iy,iz)=infunc.f(pos);
      }
    }
  }

  vector<pos_type> coord(nsamples);
  for(int i=0; i<nsamples; ++i)
    coord[i]=pos_type(Random(),Random(),Random());

  Vector<double> psi(num_splines);
  Vector<double> lap(num_splines);
  Vector<TinyVector<double,3> > grad(num_splines);
  Vector<Tensor<double,3> > hess(num_splines);

  for(int k=0; k<hess.size(); ++k) hess[k]=hess_type(0.0);


  app_log() << endl;
  {
    einspline3d_benchmark<multi_UBspline_3d_d> d_bench;
    d_bench.set(nx,ny,nz,num_splines,false);
    for(int i=0; i<num_splines; ++i) d_bench.assign(i,inData);
    double dv=0.0;
    double dvgl_v=0.0, dvgl_g=0.0, dvgl_l=0.0;
    double dvgh_v=0.0, dvgh_g=0.0, dvgh_h=0.0;
    for(int i=0; i<nsamples; ++i)
    {
      //evaluate analytic solution
      double v=infunc.f(coord[i]);
      Tensor<double,3> h;
      grad_type g=infunc.df(coord[i]);
      double l=infunc.d2f(coord[i]);
      infunc.d2f(coord[i],h);

      d_bench.einspliner.evaluate(coord[i],psi);
      dv+=(v-psi[0])*(v-psi[0]);


      d_bench.einspliner.evaluate_vgl(coord[i],psi,grad,lap);
      dvgl_v+=(v-psi[0])*(v-psi[0]);
      dvgl_g += dot(g,grad[0]);
      dvgl_l += (l-lap[0])*(l-lap[0]);

      for(int k=0; k<hess.size(); ++k)
        hess[k]=hess_type(0.0);

      d_bench.einspliner.evaluate_vgh(coord[i],psi,grad,hess);
      dvgh_v +=(v-psi[0])*(v-psi[0]);
      dvgh_g += dot(g,grad[0]);
      dvgh_h += diff(h,hess[0]);
    }

    app_log() << "Average difference = " << std::sqrt(dv)/nsamples << endl;
    app_log() << "Average difference VGL = " << std::sqrt(dvgl_v)/nsamples << " " << std::sqrt(dvgl_g)/nsamples << " " << std::sqrt(dvgl_l)/nsamples  <<  endl;
    app_log() << "Average difference VGH = " << std::sqrt(dvgh_v)/nsamples << " " << std::sqrt(dvgh_g)/nsamples << " " << std::sqrt(dvgh_h)/nsamples  <<  endl;
  }

  app_log() << endl;
  {
    // First, create splines the normal way
    multi_UBspline_3d_d *multi_spline;

    Ugrid x_grid, y_grid, z_grid;
    x_grid.start = ri; x_grid.end = rf; x_grid.num = nx;
    y_grid.start = ri; y_grid.end = rf; y_grid.num = ny;
    z_grid.start = ri; z_grid.end = rf; z_grid.num = nz;

    BCtype_d xBC, yBC, zBC;
    xBC.lCode = xBC.rCode = PERIODIC;
    yBC.lCode = yBC.rCode = PERIODIC;
    zBC.lCode = zBC.rCode = PERIODIC;

    // First, create multispline
    multi_spline = create_multi_UBspline_3d_d (x_grid, y_grid, z_grid, xBC, yBC, zBC, num_splines);

    for(int i=0; i<num_splines; ++i) 
      set_multi_UBspline_3d_d (multi_spline, i, inData.data());

    double dv=0.0;
    double dvgl_v=0.0, dvgl_g=0.0, dvgl_l=0.0;
    double dvgh_v=0.0, dvgh_g=0.0, dvgh_h=0.0;
    for(int i=0; i<nsamples; ++i)
    {
      double v=infunc.f(coord[i]);
      eval_multi_UBspline_3d_d(multi_spline,coord[i][0],coord[i][1],coord[i][2],psi.data());
      dv+=(v-psi[0])*(v-psi[0]);

      Tensor<double,3> h;
      grad_type g=infunc.df(coord[i]);
      double l=infunc.d2f(coord[i]);
      infunc.d2f(coord[i],h);

      eval_multi_UBspline_3d_d_vgl(multi_spline,coord[i][0],coord[i][1],coord[i][2],psi.data(),&(grad[0][0]),lap.data());

      dvgl_v+=(v-psi[0])*(v-psi[0]);
      dvgl_g += dot(g,grad[0]);
      dvgl_l += (l-lap[0])*(l-lap[0]);

      for(int k=0; k<hess.size(); ++k)
        hess[k]=hess_type(0.0);
      eval_multi_UBspline_3d_d_vgh(multi_spline,coord[i][0],coord[i][1],coord[i][2],psi.data(),&(grad[0][0]),&(hess[0][0]));
      dvgh_v +=(v-psi[0])*(v-psi[0]);
      dvgh_g += dot(g,grad[0]);
      dvgh_h += diff(h,hess[0]);

    }

    app_log() << "Average difference = " << std::sqrt(dv)/nsamples << endl;
    app_log() << "Average difference VGL = " << std::sqrt(dvgl_v)/nsamples << " " << std::sqrt(dvgl_g)/nsamples << " " << std::sqrt(dvgl_l)/nsamples  <<  endl;
    app_log() << "Average difference VGH = " << std::sqrt(dvgh_v)/nsamples << " " << std::sqrt(dvgh_g)/nsamples << " " << std::sqrt(dvgh_h)/nsamples  <<  endl;
  }

//  qmcplusplus::UBspline<double,3> einspline;



  return 0;
}
