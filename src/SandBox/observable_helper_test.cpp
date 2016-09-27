//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



#include <QMCHamiltonians/observable_helper.h>
#include <Utilities/RandomGenerator.h>
#include <Message/Communicate.h>
using namespace qmcplusplus;

int main(int argc, char **argv)
{
  OHMMS::Controller->initialize(argc,argv);
  std::vector<observable_helper> observables;
  std::vector<int> den_dims(3,4);
  observables.push_back(observable_helper("density"));
  observables[0].set_dimensions(den_dims,0);
  std::vector<int> dims(1,8);
  observables.push_back(observable_helper("gofr_u_u"));
  observables[1].set_dimensions(dims,0);
  observables.push_back(observable_helper("gofr_u_d"));
  observables[2].set_dimensions(dims,0);
  std::vector<int> pos_dims(2);
  pos_dims[0]=32;
  pos_dims[1]=3;
  observables.push_back(observable_helper("force"));
  observables[3].set_dimensions(pos_dims,0);
  std::vector<int> onedim(1,1);
  observables.push_back(observable_helper("energy"));
  observables[4].set_dimensions(onedim,0);
  std::string fname("observables.h5");
  hid_t fid= H5Fcreate(fname.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
  hid_t gid=H5Gcreate(fid,"observables",0);
  for(int o=0; o<observables.size(); ++o)
    observables[o].open(gid);
  {
    std::vector<double> bound(6);
    for(int i=0; i<3; ++i)
    {
      bound[2*i]=0.0;
      bound[2*i+1]=1.0;
    }
    observables[0].addProperty(bound,"extends");
  }
  {
    std::vector<double> bound(2);
    bound[0]=0.0;
    bound[1]=10.0;
    observables[1].addProperty(bound,"extends");
    observables[2].addProperty(bound,"extends");
  }
  {
    std::string units("hartree/bohr");
    observables[3].addProperty(units,"units");
    std::string blob("bare estimator");
    observables[3].addProperty(blob,"comments");
  }
  {
    std::string blob("having fun");
    observables[4].addProperty(blob,"comments");
  }
  std::vector<double> data(1000);
  std::vector<double> data2(data.size());
  for(int i=0; i<4; ++i)
  {
    for(int i=0; i<data.size(); ++i)
      data[i]=i+Random();
    for(int i=0; i<data.size(); ++i)
      data2[i]=data[i]*data[i];
    for(int o=0; o<observables.size(); ++o)
      observables[o].write(&data[0],&data2[0]);
  }
  H5Gclose(gid);
  H5Fclose(fid);
  OHMMS::Controller->finalize();
  return 0;
}
