//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF_STL_INTERFACE_H
#define QMCPLUSPLUS_HDF_STL_INTERFACE_H

#include <vector>
#include <sstream>
#include <bitset>

namespace qmcplusplus
{
/** specialization for std::vector<T>
 *
 * Used with any T with a proper h5_space_type, e.g., intrinsic, TinyVector<T,D>, Tensor<T,D>
 */
template<typename T>
struct h5data_proxy<std::vector<T>> : public h5_space_type<T, 1>
{
  using FileSpace = h5_space_type<T, 1>;
  using FileSpace::dims;
  using FileSpace::get_address;
  using data_type = std::vector<T>;

  inline h5data_proxy(const data_type& a) { dims[0] = a.size(); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    if (!checkShapeConsistency<T>(grp, aname, FileSpace::rank, dims))
      ref.resize(dims[0]);
    return h5d_read(grp, aname, get_address(&ref[0]), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    return h5d_write(grp, aname.c_str(), FileSpace::rank, dims, get_address(&ref[0]), xfer_plist);
  }

  inline bool write(const data_type& ref,
                    hid_t grp,
                    const std::string& aname,
                    const std::vector<hsize_t>& dvec,
                    hid_t xfer_plist) const
  {
    return h5d_write(grp, aname.c_str(), dvec.size(), dvec.data(), get_address(&ref[0]), xfer_plist);
  }
};

/** specialization for std::bitset<N>
 */
template<std::size_t N>
struct h5data_proxy<std::bitset<N>>
{
  using data_type = std::bitset<N>;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    unsigned long c = ref.to_ulong();
    h5data_proxy<unsigned long> hc(c);
    if (hc.read(ref, grp, aname, xfer_plist))
    {
      ref = c;
      return true;
    }
    else
      return false;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    unsigned long c = ref.to_ulong();
    h5data_proxy<unsigned long> hc(c);
    return hc.write(ref, grp, aname, xfer_plist);
  }
};


/** Specialization for std::string */
template<>
struct h5data_proxy<std::string>
{
  using data_type = std::string;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    hid_t dataset = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
    if (dataset > -1)
    {
      hid_t datatype = H5Dget_type(dataset);
      hsize_t dim_out;
      if (datatype == H5T_NATIVE_CHAR)
      {
        hid_t dataspace = H5Dget_space(dataset);
        hid_t status    = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);
        H5Sclose(dataspace);
      }
      else
      {
        dim_out = H5Tget_size(datatype);
      }
      ref.resize(dim_out);
      herr_t ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, xfer_plist, &(ref[0]));
      H5Tclose(datatype);
      H5Dclose(dataset);
      return ret != -1;
    }
    return false;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    hid_t str80 = H5Tcopy(H5T_C_S1);
    H5Tset_size(str80, ref.size());
    hsize_t dim = 1;

    herr_t ret = -1;
    hid_t h1   = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
    if (h1 < 0) // missing create one
    {
      hid_t dataspace = H5Screate_simple(1, &dim, NULL);
      hid_t dataset   = H5Dcreate(grp, aname.c_str(), str80, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret             = H5Dwrite(dataset, str80, H5S_ALL, H5S_ALL, xfer_plist, ref.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
    else
    {
      ret = H5Dwrite(h1, str80, H5S_ALL, H5S_ALL, xfer_plist, ref.data());
    }
    H5Dclose(h1);
    return ret != -1;
  }
};

/// Specialization for vector of strings
template<>
struct h5data_proxy<std::vector<std::string>>
{
  using data_type = std::vector<std::string>;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);
    hid_t dataset = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
    std::vector<char*> char_list;
    herr_t ret = -1;
    if (dataset > -1)
    {
      hsize_t dim_out;
      hid_t dataspace = H5Dget_space(dataset);
      hid_t status    = H5Sget_simple_extent_dims(dataspace, &dim_out, NULL);

      char_list.resize(dim_out);
      ret = H5Dread(dataset, datatype, H5S_ALL, H5S_ALL, xfer_plist, char_list.data());

      for (int i = 0; i < dim_out; i++)
        ref.push_back(char_list[i]);

      H5Dvlen_reclaim(datatype, dataspace, xfer_plist, char_list.data());

      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
    H5Tclose(datatype);

    return ret >= 0;
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    // See the section in the HDF user's manual on datatypes,
    // particularly the subsection on strings.
    // (e.g. http://davis.lbl.gov/Manuals/HDF5-1.8.7/UG/11_Datatypes.html)
    // and stackoverflow
    // https://stackoverflow.com/questions/6184817/hdf5-inserting-a-set-of-strings-in-a-dataset
    hid_t datatype = H5Tcopy(H5T_C_S1);
    H5Tset_size(datatype, H5T_VARIABLE);
    hsize_t dim = ref.size();

    // Create vector of pointers to the actual string data
    std::vector<const char*> char_list;
    for (int i = 0; i < ref.size(); i++)
      char_list.push_back(ref[i].data());

    hid_t h1   = H5Dopen(grp, aname.c_str(), H5P_DEFAULT);
    herr_t ret = -1;
    if (h1 < 0) // missing create one
    {
      hid_t dataspace = H5Screate_simple(1, &dim, NULL);
      hid_t dataset   = H5Dcreate(grp, aname.c_str(), datatype, dataspace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      ret             = H5Dwrite(dataset, datatype, H5S_ALL, H5S_ALL, xfer_plist, char_list.data());
      H5Sclose(dataspace);
      H5Dclose(dataset);
    }
    else
      ret = H5Dwrite(h1, datatype, H5S_ALL, H5S_ALL, xfer_plist, char_list.data());

    H5Dclose(h1);
    return ret >= 0;
  }
};

template<>
struct h5data_proxy<std::ostringstream>
{
  using data_type = std::ostringstream;

  h5data_proxy(const data_type& a) {}

  inline bool read(data_type& ref, hid_t grp, char* name, hid_t xfer_plist = H5P_DEFAULT) { return false; }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  {
    std::string clone(ref.str());
    h5data_proxy<std::string> proxy(clone);
    return proxy.write(clone, grp, aname);
  }
};

} // namespace qmcplusplus
#endif
