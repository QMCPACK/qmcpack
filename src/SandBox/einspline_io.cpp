/** @file einspline_validation.cpp
 *
 * Test code to check the correctness of the einspline
 */
#include <SandBox/SplineTestBase.h>
#include <einspline/multi_bspline_copy.h>
#include <Message/Communicate.h>
#include "Numerics/HDFNumericAttrib.h"
#include <io/hdf_archive.h>
#include <mpi/collectives.h>
#include <getopt.h>

//testing subgrid
qmcplusplus::TinyVector<double,3> lower(0.25,0.25,0.25);
qmcplusplus::TinyVector<double,3> upper(0.75,0.75,0.75);

namespace qmcplusplus
{
  template<typename ENGT>
    struct h5data_proxy<einspline_engine<ENGT> >
    : public h5_space_type<typename einspline_engine<ENGT>::value_type,4>
  {
    typedef typename einspline_engine<ENGT>::value_type value_type;
    using h5_space_type<value_type,4>::dims;
    using h5_space_type<value_type,4>::get_address;
    typedef einspline_engine<ENGT> data_type;

    data_type& ref_;

    inline h5data_proxy(data_type& a): ref_(a)
    {
      dims[0]=a.spliner->x_grid.num+3;
      dims[1]=a.spliner->y_grid.num+3;
      dims[2]=a.spliner->z_grid.num+3;
      dims[3]=a.spliner->z_stride;
    }

    inline bool read(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      //if(!h5d_getspace(grp,aname,this->size(),dims)) ref_.resize(dims[0]);
      if(ref_.spliner) 
        return h5d_read(grp,aname,get_address(ref_.spliner->coefs),xfer_plist);
      else
        return false;
    }

    inline bool write(hid_t grp, const std::string& aname, hid_t xfer_plist=H5P_DEFAULT)
    {
      return h5d_write(grp,aname.c_str(),this->size(),dims,get_address(ref_.spliner->coefs),xfer_plist);
    }
  };

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
        hdf_archive h5f;
        h5f.open("sample.h5");
        h5f.read(big,"big_40");
        h5f.read(small,"small_40");
      }

      typedef TinyVector<real_type,3> pos_type;
      pos_type start(0.0);
      pos_type end(1.0);
      TinyVector<int,3> ng(ngx,ngy,ngz), ngh(ngx/2,ngy/2,ngz/2);


      //create/set multi_UBspline_3d as usual on the dense grid
      einspline_engine<ENGT> einspliner;
      einspliner.create(start,end,ng,PERIODIC,num_splines);

      Timer clock;
      clock.restart();
      for(int i=0; i<num_splines; ++i) einspliner.set(i,big);
      double t_spline=clock.elapsed();

//#pragma omp parallel for 
//      for(int i=0; i<num_splines; ++i) 
//      {
//        einspliner.set(i,big);
//      }

      einspline_engine<ENGT> einspliner2;
      einspliner2.create(start,end,ng,PERIODIC,num_splines);

      //dump the spline to a hdf5
      clock.restart();
      {
        hdf_archive h5f;
        h5f.create("spline.h5");
        h5f.write(einspliner,"bspline");
      }
      double t_write=clock.elapsed();

      clock.restart();
      //initialize einspline2 from the file
      {
        hdf_archive h5f;
        h5f.open("spline.h5");
        h5f.read(einspliner2,"bspline");
      }
      double t_read=clock.elapsed();

      app_log() << "Timing spline=" << t_spline << " write = " << t_write << " read = " << t_read 
        << "  saving = " << t_spline/t_read<< endl;
    }
}
int main(int argc, char** argv) 
{

  using namespace qmcplusplus;
  OHMMS::Controller->initialize(argc,argv);
  Communicate* mycomm=OHMMS::Controller;
  OhmmsInfo Welcome("einspline",mycomm->rank());
  Random.init(0,1,11);

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
  vector<TinyVector<double,3> > coord_s(coord.size());
  for(int i=0; i<coord.size(); ++i) convert(coord[i],coord_s[i]);

  cout << "\nTesting einspline transformation" << endl;

  SplineTest<multi_UBspline_3d_d> test_d(param);
  test_d.test(coord);
  return 0;
}
