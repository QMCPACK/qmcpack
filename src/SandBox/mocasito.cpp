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
    
    



#include <alps/hdf5.hpp>
int main(int argc, char** argv)
{
  int length=8;
  {
    alps::hdf5::oarchive h5ar("bela.h5");
    int d = 44;
    h5ar << alps::make_pvp("/int", d);
  }
  {
    alps::hdf5::oarchive h5ar("bela.h5");
    std::vector<int> d(length, 42);
    h5ar << alps::make_pvp("/foo/bar2", d);
  }
  {
    alps::hdf5::iarchive h5ar("bela.h5");
    std::vector<int> d;
    h5ar >> alps::make_pvp("/foo/bar2", d);
    copy (d.begin(), d.end(), std::ostream_iterator<int, char, std::char_traits<char> >(std::cout, " "));
    std::cout << std::endl;
  }
  {
    alps::hdf5::oarchive h5ar("bela.h5");
    std::complex<double> *d = new std::complex<double>[length];
    h5ar << alps::make_pvp("test/data", d, length);
    //std::complex<double>* d
    //std::vector<std::complex<double> > d(length, std::complex<double>(0.3,0.1));
    //h5ar << alps::make_pvp("test/data", d, length);
  }
}
