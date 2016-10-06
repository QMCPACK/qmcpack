//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by:  Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_HDF5_FW_H
#define QMCPLUSPLUS_HDF5_FW_H
#include <cstring>
#include <HDFVersion.h>
#include <io/hdf.h>
namespace qmcplusplus
{
class HDF5_FW_float
{
public:
  HDF5_FW_float():RANK(1) {}
  ~HDF5_FW_float() {}

  void openFile( std::string filename)
  {
    c_file = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  }

  void closeFile()
  {
    if (H5Fclose(c_file)>-1)
      c_file=-1;
  }


  void setStep(int step)
  {
    std::stringstream gname("");
    gname<<"Block_"<< step;
    d_file = H5Gopen(c_file,gname.str().c_str());
    dataset = H5Dopen(d_file, "Positions");
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
//         cparms = H5Dget_create_plist(dataset);
//         if (H5D_CHUNKED == H5Pget_layout(cparms))  {
//           rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
//         }
  }

  int getFloat(int first, int last, std::vector<float>& data_out)
  {
    if (dims[0]< first)
      return 0;
    if (dims[0]< last)
      data_out.resize(dims[0]-first);
    offset[1] = 0;
    offset[0] = first;
    count[1]  = 1;
    count[0]  = data_out.size();
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    dims_m[0]=data_out.size();
    memspace = H5Screate_simple (1, dims_m, NULL);
    offset_m[0] = 0;
    count_m[0]  = data_out.size();
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_m, NULL, count_m, NULL);
    status = H5Dread (dataset, H5T_NATIVE_FLOAT, memspace, dataspace, H5P_DEFAULT, &(data_out[0]));
    H5Sclose(memspace);
    return data_out.size();
  }

  void endStep()
  {
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Gclose(d_file);
  }

private:

  int RANK;

  hid_t       dataset;
  hid_t       dataspace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[2],dims_m[1];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[2],count_m[1];
  hsize_t     offset[2],offset_m[1];

  herr_t      status, status_n;

  int         rank, rank_chunk;
  hsize_t hi, hj;


  hid_t c_file;
  hid_t d_file;
};

class HDF5_FW_long
{
public:
  HDF5_FW_long():RANK(1) {}
  ~HDF5_FW_long() {}

  void openFile( std::string filename)
  {
    c_file = H5Fopen(filename.c_str(),H5F_ACC_RDONLY,H5P_DEFAULT);
  }

  void closeFile()
  {
    if (H5Fclose(c_file)>-1)
      c_file=-1;
  }

  void setID( std::string ID)
  {
    IDstring=ID;
  }


  void setStep(int step)
  {
    std::stringstream gname("");
    gname<<"Block_"<< step;
    d_file = H5Gopen(c_file,gname.str().c_str());
    dataset = H5Dopen(d_file, "Positions");
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))
    {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    }
  }

  int getFloat(int first, int last, std::vector<float>& data_out)
  {
    if (dims[0]< last)
      data_out.resize(dims[0]-first);
    offset[0] = 0;
    offset[1] = first;
    count[0]  = 1;
    count[1]  = data_out.size();
    status = H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, offset, NULL, count, NULL);
    dims_m[0]=data_out.size();
    memspace = H5Screate_simple (RANK, dims_m, NULL);
    offset_m[0] = 0;
    count_m[0]  = data_out.size();
    status = H5Sselect_hyperslab (memspace, H5S_SELECT_SET, offset_m, NULL, count_m, NULL);
    status = H5Dread (dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT, &(data_out[0]));
    return data_out.size();
  }

  void endStep()
  {
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Gclose(d_file);
  }

  void readAll(int step, std::vector<long>& data_out)
  {
    std::stringstream gname("");
    gname<<"Block_"<< step;
    hid_t d_file = H5Gopen(c_file,gname.str().c_str());
    dataset = H5Dopen(d_file, IDstring.c_str());
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    data_out.resize(dims[0]);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))
    {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    }
    memspace = H5Screate_simple(RANK,dims,NULL);
    status = H5Dread(dataset, H5T_NATIVE_LONG, memspace, dataspace, H5P_DEFAULT, &(data_out[0]));
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
    H5Gclose(d_file);
  }

private:

  int RANK;

  hid_t       dataset;
  hid_t       dataspace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[2],dims_m[1];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[2],count_m[1];
  hsize_t     offset[2],offset_m[1];

  herr_t      status, status_n;

  int         rank, rank_chunk;
  hsize_t hi, hj;


  hid_t c_file;
  hid_t d_file;
  std::string IDstring;
};

class HDF5_FW_observables
{
public:
  HDF5_FW_observables() {}
  ~HDF5_FW_observables() {}

  void setFileName( std::string fn)
  {
    std::stringstream sstr("");
    sstr<<fn<<".storedOBS.h5";
    filename=sstr.str();
  }

  std::string getFileName()
  {
    return filename;
  }

  void makeFile()
  {
    hid_t d1= H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    HDFVersion cur_version;
    cur_version.write(d1,hdf::version);
    if (H5Fclose(d1) > -1)
      d1 = -1;
  }

  void openFile()
  {
    c_file = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }

  void closeFile()
  {
    if (H5Fclose(c_file)>-1)
      c_file=-1;
  }

  void addStep(int step, std::vector<double>& Observables)
  {
    std::stringstream sstr("");
    sstr<<"Block_"<<step;
//         d_file = H5Gcreate(c_file,sstr.str().c_str(),0);
    dims[0] = Observables.size();
    maxdims[0] = H5S_UNLIMITED;
    rank=1;
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    hid_t p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
//
//         sstr.str("OBS");
//         std::string groupName = sstr.str();
    dataset =  H5Dcreate(c_file, sstr.str().c_str(), H5T_NATIVE_DOUBLE, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,&(Observables[0]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
//         if (H5Gclose(d_file) > -1) d_file = -1;
  }

  void readStep(int step, std::vector<double>& data_out)
  {
    std::stringstream gname("");
    gname<<"Block_"<< step;
//         d_file = H5Gopen(c_file,gname.str().c_str());
//         gname.str("OBS");
    dataset = H5Dopen(c_file, gname.str().c_str());
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    data_out.resize(dims[0]);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))
    {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    }
    memspace = H5Screate_simple(rank,dims,NULL);
    status = H5Dread(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT, &(data_out[0]));
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
//         H5Gclose(d_file);
  }

  int numObsStep(int step)
  {
    std::stringstream gname("");
    gname<<"Block_"<< step;
    dataset = H5Dopen(c_file, gname.str().c_str());
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return dims[0];
  }

private:

  hid_t       dataset;
  hid_t       dataspace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[2],maxdims[2];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[3],count_m[1];
  hsize_t     offset[3],offset_m[1];

  herr_t      status, status_n;

  int         rank, rank_chunk;
  hsize_t hi, hj;


  hid_t c_file;
  hid_t d_file;
  std::string filename;
};

class HDF5_FW_weights
{
public:
  HDF5_FW_weights() {}
  ~HDF5_FW_weights() {}

  void setFileName( std::string fn)
  {
    std::stringstream sstr("");
    sstr<<fn<<".storedWeights.h5";
    filename=sstr.str();
  }

  std::string getFileName()
  {
    return filename;
  }

  void makeFile()
  {
    hid_t d1= H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    HDFVersion cur_version;
    cur_version.write(d1,hdf::version);
    if (H5Fclose(d1) > -1)
      d1 = -1;
  }

  void openFile()
  {
    c_file = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }

  void closeFile()
  {
    if (H5Fclose(c_file)>-1)
      c_file=-1;
  }

  void addFW(int gap)
  {
    std::stringstream sstr("");
    sstr<<"Age_"<<gap;
    d_file = H5Gcreate(c_file,sstr.str().c_str(),0);
  }

  void closeFW()
  {
    if (H5Gclose(d_file) > -1)
      d_file = -1;
  }

  void addStep(int step, std::vector<int>& Observables)
  {
    dims[0] = Observables.size();
    maxdims[0] = H5S_UNLIMITED;
    rank=1;
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    hid_t p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream gname("");
    gname<<"Block_"<<step;
    dataset =  H5Dcreate(d_file, gname.str().c_str(), H5T_NATIVE_INT, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(Observables[0]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }

  void readStep(int age, int step, std::vector<int>& data_out)
  {
    std::stringstream gname("");
    gname<<"Age_"<<age<<"/Block_"<< step;
//         hid_t d_file = H5Dopen(c_file,gname.str().c_str());
//         gname.str("OBS");
    dataset = H5Dopen(c_file, gname.str().c_str());
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    data_out.resize(dims[0]);
    cparms = H5Dget_create_plist(dataset);
    if (H5D_CHUNKED == H5Pget_layout(cparms))
    {
      rank_chunk = H5Pget_chunk(cparms, 2, chunk_dims);
    }
    memspace = H5Screate_simple(rank,dims,NULL);
    status = H5Dread(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT, &(data_out[0]));
    H5Pclose(cparms);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    H5Sclose(memspace);
//         H5Gclose(d_file);
  }

  int numWgtStep(int step)
  {
    std::stringstream gname("");
    gname<<"Age_0/Block_"<< step;
    dataset = H5Dopen(c_file, gname.str().c_str());
    dataspace = H5Dget_space(dataset);
    rank = H5Sget_simple_extent_ndims(dataspace);
    status_n  = H5Sget_simple_extent_dims(dataspace, dims, NULL);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    return dims[0];
  }


private:

  hid_t       dataset;
  hid_t       dataspace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[1],maxdims[1];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[3],count_m[1];
  hsize_t     offset[3],offset_m[1];

  herr_t      status, status_n;

  int         rank, rank_chunk;
  hsize_t hi, hj;


  hid_t c_file;
  hid_t d_file;
  std::string filename;
};

class HDF5_FW_density
{
public:
  HDF5_FW_density() {}
  ~HDF5_FW_density() {}

  void setFileName( std::string fn)
  {
    std::stringstream sstr("");
    sstr<<fn<<".storedDensity.h5";
    filename=sstr.str();
  }

  std::string getFileName()
  {
    return filename;
  }

  void makeFile()
  {
    hid_t d1= H5Fcreate(filename.c_str(),H5F_ACC_TRUNC,H5P_DEFAULT,H5P_DEFAULT);
    HDFVersion cur_version;
    cur_version.write(d1,hdf::version);
    if (H5Fclose(d1) > -1)
      d1 = -1;
  }

  void openFile()
  {
    c_file = H5Fopen(filename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
  }

  void closeFile()
  {
    if (H5Fclose(c_file)>-1)
      c_file=-1;
  }

  void writeDensity(std::vector<int> info, std::vector<int>& Observables)
  {
    rank=1;
    dims[0]=Observables.size();
    maxdims[0] = Observables.size();
    int age=info[1];
//         for (int x=0;x<DD;x++)
//           for (int y=0;y<DD;y++)
//             for (int z=0;z<DD;z++)
//             {
//               dataspace  = H5Screate_simple(rank, dims, maxdims);
//               hid_t p = H5Pcreate (H5P_DATASET_CREATE);
//               H5Pset_chunk(p,rank,dims);
//
//               std::stringstream gname("");
//               gname<<"Block_"<<age<<"/X_"<<x<<"/Y_"<<y<<"/Z_"<<z;
//               dataset =  H5Dcreate(c_file, gname.str().c_str(), H5T_NATIVE_INT, dataspace, p);
//               memspace = H5Screate_simple( rank, dims, NULL);
//               status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(Observables[DD2*x+DD*y+z]));
//
//               H5Sclose(memspace);
//               H5Dclose(dataset);
//               H5Sclose(dataspace);
//             }
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    hid_t p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream gname("");
    gname<<"Block_"<<age;
    dataset =  H5Dcreate(c_file, gname.str().c_str(), H5T_NATIVE_INT, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(Observables[0]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    dims[0]=1;
    maxdims[0] = 1;
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    gname<<"_totalWeight";
    dataset =  H5Dcreate(c_file, gname.str().c_str(), H5T_NATIVE_INT, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(info[2]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }

  void writeIons(std::vector<std::string> NMS, std::vector<double>& Pos, double& latticeConstant, int& nElectrons, int& ngrids)
  {
    rank=1;
    int sze(-1);
    for(int i=0; i<NMS.size(); i++)
      sze+=NMS[i].size()+1;
    char names[sze];
    strcpy (names,NMS[0].c_str());
    for(int i=1; i<NMS.size(); i++)
    {
      strncat(names,",",1);
      strncat(names,NMS[i].c_str(),NMS[i].size());
    }
    sze = strlen(names);
    dims[0]=sze;
    maxdims[0] = sze;
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    hid_t p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream gname("");
    gname<<"Ion_Str";
    dataset =  H5Dcreate(c_file, gname.str().c_str(), H5T_C_S1, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_C_S1, memspace, dataspace, H5P_DEFAULT,&(names[0]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    dims[0]=Pos.size();
    maxdims[0] = Pos.size();
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream pname("");
    pname<<"Ion_Pos";
    dataset =  H5Dcreate(c_file, pname.str().c_str(), H5T_NATIVE_DOUBLE, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,&(Pos[0]));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    dims[0]=1;
    maxdims[0] = 1;
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream lname("");
    lname<<"LatticeConstant";
    dataset =  H5Dcreate(c_file, lname.str().c_str(), H5T_NATIVE_DOUBLE, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_DOUBLE, memspace, dataspace, H5P_DEFAULT,&(latticeConstant));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream ename("");
    ename<<"nElectrons";
    dataset =  H5Dcreate(c_file, ename.str().c_str(), H5T_NATIVE_INT, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(nElectrons));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
    dataspace  = H5Screate_simple(rank, dims, maxdims);
    p = H5Pcreate (H5P_DATASET_CREATE);
    H5Pset_chunk(p,rank,dims);
    std::stringstream dname("");
    dname<<"nbins";
    dataset =  H5Dcreate(c_file, dname.str().c_str(), H5T_NATIVE_INT, dataspace, p);
    memspace = H5Screate_simple( rank, dims, NULL);
    status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, dataspace, H5P_DEFAULT,&(ngrids));
    H5Sclose(memspace);
    H5Dclose(dataset);
    H5Sclose(dataspace);
  }

private:

  hid_t       dataset;
  hid_t       dataspace;
  hid_t       memspace;
  hid_t       cparms;
  hsize_t     dims[1],maxdims[1];                     /* dataset and chunk dimensions*/
  hsize_t     chunk_dims[1];
  hsize_t     col_dims[1];
  hsize_t     count[3],count_m[1];
  hsize_t     offset[3],offset_m[1];

  herr_t      status, status_n;

  int         rank, rank_chunk;
  hsize_t hi, hj;


  hid_t c_file;
  hid_t d_file;
  std::string filename;
};

}
#endif
