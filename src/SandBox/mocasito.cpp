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
    std::copy (d.begin(), d.end(), std::ostream_iterator<int, char, std::char_traits<char> >(std::cout, " "));
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
