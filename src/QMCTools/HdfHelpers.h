#ifndef HDF_HELPERS_H
#define HDF_HELPERS_H
#include <string>
#include <vector>
#include <iostream>
#include "hdf5.h"

namespace hdfhelper {
  template<typename T>
  hid_t getHdfDataType(const T& val);

  template<>
  hid_t getHdfDataType<int>(const int& val) 
  {
    return H5T_NATIVE_INT;
  }
  
  template<>
  hid_t getHdfDataType<float>(const float& val) 
  {
    return H5T_NATIVE_FLOAT;
  }
  
  template<>
  hid_t getHdfDataType<double>(const double& val) 
  {
    return H5T_NATIVE_DOUBLE;
  }
    
  // note currently only support int,float and double  
  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName,
			T* const data,
			hid_t loc,
			int rank,
			hsize_t* dimensionality)
  {
    hid_t type = H5Tcopy(getHdfDataType(*data));
    hid_t dspace = H5Screate_simple(rank, dimensionality, NULL);
    hid_t dset = H5Dcreate2(loc, fieldName.c_str(), type, dspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    H5Sclose(dspace);
    H5Tclose(type);
    H5Dclose(dset);
    return ret;
  }

  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, 
			const std::vector<T>& data, 
			hid_t loc,
			int rank, 
			hsize_t* dimensionality) 
  {
    return writeNumsToHDF(fieldName, &(data.front()), loc, rank, dimensionality);
  }
    
  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, const std::vector<T>& data, hid_t loc) 
  {
    hsize_t len = data.size();
    return writeNumsToHDF(fieldName, data, loc, 1, &len);
  }

  // this assumes a 1d arrangement
  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, T *const data, int len, hid_t loc) 
  {
    hsize_t lena = len;
    return writeNumsToHDF(fieldName, data, loc, 1, &lena);
  }

  
  template<typename T>
  herr_t writeNumsToHDF(const std::string& fieldName, T data, hid_t loc) 
  {
    std::vector<T> tmp{data};
    return writeNumsToHDF(fieldName, tmp, loc);
  }

  herr_t writeStringToHDF(const std::string& fieldName, const std::string& str, hid_t loc) 
  {
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
  
  hid_t makeHDFGroup(const std::string& gName, hid_t location) 
  {
    return H5Gcreate2(location, gName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  // we are assuming that enough memory to hold all the data has already been
  // allocated an is associated with the pointer that is handed in
  template<typename T>
  hid_t readNumsFromHDF(const std::string& fieldName, T* const data, hid_t file)
  {
    hid_t type = H5Tcopy(getHdfDataType(*data));
    hid_t dset = H5Dopen1(file, fieldName.c_str());
    herr_t status = H5Dread(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    status = H5Dclose(dset);
    return status;
  }

  // also do a version to read a single number
  template<typename T>
  hid_t readNumsFromHDF(const std::string& fieldName, T& data, hid_t file)
  {
    return readNumsFromHDF(fieldName, &data, file);
  }

  void getDatasetDimensionality(const std::string& fieldName,
				 std::vector<int>& dims, hid_t file)
  {
    hid_t dset = H5Dopen1(file, fieldName.c_str());
    hid_t dspace = H5Dget_space(dset);
    const int dimensionality = H5Sget_simple_extent_ndims(dspace);
    dims.resize(dimensionality);
    hsize_t maxdims[dimensionality];
    hsize_t intdims[dimensionality];
    H5Sget_simple_extent_dims(dspace, intdims, maxdims);
    for (int i = 0; i < dimensionality; i++)
    {
      dims[i] = intdims[i];
    }

  }

  // hand in name of dataset, a vector whose values specify dimensions to read
  // in from dataset, a pointer to memory that will hold the data that is read out
  // and the handle for the file
  //
  // note the vector must have dimensionality corresponding to the dataset,
  // values for a dimension must be [0,num_entries-1] for that dimension to specify
  // which value to hold and a -1 to say grab all elements from that dimension
  // for example, if the dataset was [5,2,6] and the vector contained (2,1,-1),
  // this would grab 6 elements corresponding to [2,1,:]
  template<typename T>
  void readPortionFromHDF(const std::string& fieldName,
			  const std::vector<int>& readSpec,
			  T* outdata, hid_t file)
  {
    hid_t dset = H5Dopen1(file, fieldName.c_str());
    hid_t dspace = H5Dget_space(dset);
    const int dimensionality = H5Sget_simple_extent_ndims(dspace);
    hsize_t maxdims[dimensionality];
    hsize_t intdims[dimensionality];
    H5Sget_simple_extent_dims(dspace, intdims, maxdims);

    // check that what we have requested is within bounds
    for (int i = 0; i < dimensionality; i++) {
      if (readSpec[i] < -1 || readSpec[i] >= static_cast<int>(intdims[i])) {
	std::cerr << "Requested an invalid slice of the dataset named " << fieldName << std::endl;
	std::cerr << "dataset has dimensionality: ";
	for (int j = 0; j < dimensionality; j++) std::cerr << maxdims[j] << "  ";
	std::cerr << std::endl;
	std::cerr << "requested slice defined by: ";
	for (int j = 0; j < dimensionality; j++) std::cerr << readSpec[j] << "  ";
	std::cerr << std::endl;
	std::cerr << "failed on dim " << i << std::endl;
	if (readSpec[i] < -1) {
	  std::cerr << "  failed because < -1" << std::endl;
	} else {
	  std::cerr << "  failed because " << readSpec[i] << " >= " << intdims[i] << std::endl;
	}
	exit(1);
      }
    }
    hid_t type = H5Tcopy(getHdfDataType(*outdata));

    // set up hyperslab parameters to match what we are asking for
    hsize_t offset[dimensionality];
    hsize_t block[dimensionality];
    hsize_t stride[dimensionality];
    hsize_t count[dimensionality];
    for (int i = 0; i < dimensionality; i++) {
      offset[i] = readSpec[i];
      count[i] = 1;
      if (readSpec[i] == -1) {
	offset[i] = 0;
	count[i] = intdims[i];
      }
      block[i] = 1;
      stride[i] = 1;
    }
    H5Sselect_hyperslab(dspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    
    // make memory dataspace
    hsize_t dimsm[] = {1};
    for (int i = 0; i < dimensionality; i++) {
      dimsm[0] *= count[i];
    }
    //std::cerr << "memory dataspace contains " << dimsm[0] << " elements" << std::endl;

    hid_t memspace = H5Screate_simple(1,dimsm,NULL);  
    hsize_t offset_out[] = {0};
    hsize_t count_out[] = {dimsm[0]};
    H5Sselect_hyperslab(memspace, H5S_SELECT_SET, offset_out, NULL, count_out, NULL);
    
    H5Dread(dset, type, memspace, dspace, H5P_DEFAULT, outdata);
  }

};

#endif
