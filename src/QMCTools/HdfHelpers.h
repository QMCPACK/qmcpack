#ifndef HDF_HELPERS_H
#define HDF_HELPERS_H
#include <string>
#include <vector>
#include "hdf5.h"

namespace hdfHelper {
  template<typename T>
  hid_t getHdfDataType(const T& val);

  template<>
    hid_t getHdfDataType<int>(const int& val) {
    return H5T_NATIVE_INT;
  }
  template<>
    hid_t getHdfDataType<float>(const float& val) {
    return H5T_NATIVE_FLOAT;
  }
  template<>
    hid_t getHdfDataType<double>(const double& val) {
    return H5T_NATIVE_DOUBLE;
  }
    
  // note currently only support int,float and double
  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, const std::vector<T>& data, hid_t loc,
			int rank, hsize_t* dimensionality) {
    hid_t type = H5Tcopy(getHdfDataType(data[0]));
    hid_t dspace = H5Screate_simple(rank, dimensionality, NULL);
    hid_t dset = H5Dcreate2(loc, fieldName.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(data.front()));
    H5Sclose(dspace);
    H5Tclose(type);
    H5Dclose(dset);
    return ret;
  }

  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, const std::vector<T>& data, hid_t loc) {
    hsize_t len = data.size();
    return writeNumsToHDF(fieldName, data, loc, 1, &len);
  }

  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, T data, hid_t loc) {
    std::vector<T> tmp{data};
    return writeNumsToHDF(fieldName, tmp, loc);
  }

  herr_t writeStringToHDF(const std::string& fieldName, const std::string& str, hid_t loc) {
    hsize_t ns = 1;
    hid_t strtype = H5Tcopy(H5T_C_S1);
    herr_t ret = H5Tset_size(strtype, str.size());
    hid_t dspace = H5Screate_simple(1, &ns, NULL);
    hid_t dset = H5Dcreate2(loc, fieldName.c_str(), strtype, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    ret = H5Dwrite(dset, strtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, str.c_str());
    H5Sclose(dspace);
    H5Tclose(strtype);
    H5Dclose(dset);
    return ret;
  }
  
  hid_t makeHDFGroup(const std::string& gName, hid_t location) {
    return H5Gcreate2(location, gName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

};

#endif
