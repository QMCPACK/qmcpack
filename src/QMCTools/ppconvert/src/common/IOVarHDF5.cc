//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "IOVarHDF5.h"

namespace IO {

  template<> bool
  IOVarHDF5<string,0>::VarRead(string &val)
  {
    hid_t type = H5Dget_type(DatasetID);
    size_t length = H5Tget_size(type);
    blitz::Array<char,1> charArray(length);
    herr_t status = H5Dread(DatasetID, type, MemSpaceID,
			    DiskSpaceID, H5P_DEFAULT, charArray.data());
    val = charArray.data();
    H5Tclose(type);
    return (status == 0);
  }

  template<> bool
  IOVarHDF5<string,1>::VarRead(Array<string,1> &val)
  {
    hid_t type = H5Dget_type(DatasetID);
    TinyVector<hsize_t,1> hdims;
    TinyVector<int,1> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    val.resize(dims);
    hsize_t length = H5Tget_size(type);
    blitz::Array<char,2> charArray(dims[0],length);
    herr_t status = H5Dread(DatasetID, type, MemSpaceID,
			    DiskSpaceID, H5P_DEFAULT, charArray.data());
    for (int i=0; i<dims[0]; i++)
      val(i) = &(charArray(i,0));
    H5Tclose(type);
    return (status == 0);
  }
  
  template<> bool
  IOVarHDF5<string,2>::VarRead(Array<string,2> &val)
  {
    hid_t type = H5Dget_type(DatasetID);
    TinyVector<hsize_t,2> hdims;
    TinyVector<int,2> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    val.resize(dims);
    hsize_t length = H5Tget_size(type);
    blitz::Array<char,3> charArray(dims[0],dims[1],length);
    herr_t status = H5Dread(DatasetID, type, MemSpaceID,
			    DiskSpaceID, H5P_DEFAULT, charArray.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	val(i,j) = &(charArray(i,j,0));
    H5Tclose(type);
    return (status == 0);
  }
  
  template<> bool
  IOVarHDF5<string,3>::VarRead(Array<string,3> &val)
  {
    hid_t type = H5Dget_type(DatasetID);
    TinyVector<hsize_t,3> hdims;
    TinyVector<int,3> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    val.resize(dims);
    hsize_t length = H5Tget_size(type);
    blitz::Array<char,4> charArray(dims[0],dims[1],dims[2],length);
    herr_t status = H5Dread(DatasetID, type, MemSpaceID,
			    DiskSpaceID, H5P_DEFAULT, charArray.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	for (int k=0; k<dims[2]; k++)
	  val(i,j,k) = &(charArray(i,j,k,0));
    H5Tclose(type);
    return (status == 0);
  }
  
  template<> bool
  IOVarHDF5<string,4>::VarRead(Array<string,4> &val)
  {
    hid_t type = H5Dget_type(DatasetID);
    TinyVector<hsize_t,4> hdims;
    TinyVector<int,4> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    val.resize(dims);
    hsize_t length = H5Tget_size(type);
    blitz::Array<char,5> charArray(dims[0],dims[1],dims[2],dims[3],length);
    herr_t status = H5Dread(DatasetID, type, MemSpaceID,
			    DiskSpaceID, H5P_DEFAULT, charArray.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	for (int k=0; k<dims[2]; k++)
	  for (int l=0; l<dims[3]; l++)
	    val(i,j,k,l) = &(charArray(i,j,k,l,0));
    H5Tclose(type);
    return (status == 0);    
  }


  template<> bool
  IOVarHDF5<bool,0>::VarRead(bool &val)
  {
    unsigned char cval;
    herr_t status = H5Dread(DatasetID, BoolTypeID, MemSpaceID, 
			    DiskSpaceID, H5P_DEFAULT, &cval);
    val = (cval == (unsigned char)1);
    return (status == 0);
  }  
  
  template<> bool
  IOVarHDF5<bool,1>::VarRead(Array<bool,1> &val)
  {
    TinyVector<hsize_t,1> hdims;
    TinyVector<int,1> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    Array<unsigned char,1> cval(dims);
    val.resize(dims);
    herr_t status = H5Dread(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID,
			     H5P_DEFAULT, cval.data());
    for (int i=0; i<dims[0]; i++)
      val(i) = (cval(i) == 1);
    return (status == 0);
  }
  
  template<> bool
  IOVarHDF5<bool,2>::VarRead(Array<bool,2> &val)
  {
    TinyVector<hsize_t,2> hdims;
    TinyVector<int,2> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    Array<unsigned char,2> cval(dims);
    val.resize(dims);
    herr_t status = H5Dread(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID,
			     H5P_DEFAULT, cval.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	val(i,j) = (cval(i,j) == 1);
    return (status == 0); 
  }
  
  template<> bool
  IOVarHDF5<bool,3>::VarRead(Array<bool,3> &val)
  {
    TinyVector<hsize_t,3> hdims;
    TinyVector<int,3> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    Array<unsigned char,3> cval(dims);
    val.resize(dims);
    herr_t status = H5Dread(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID,
			     H5P_DEFAULT, cval.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	for (int k=0; k<dims[2]; k++)
	  val(i,j,k) = (cval(i,j,k) == 1);
    return (status == 0); 
  }
  
  template<> bool
  IOVarHDF5<bool,4>::VarRead(Array<bool,4> &val)
  {
    TinyVector<hsize_t,4> hdims;
    TinyVector<int,4> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
    dims = hdims;
    Array<unsigned char,4> cval(dims);
    val.resize(dims);
    herr_t status = H5Dread(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID,
			     H5P_DEFAULT, cval.data());
    for (int i=0; i<dims[0]; i++)
      for (int j=0; j<dims[1]; j++)
	for (int k=0; k<dims[2]; k++)
	  for (int l=0; l<dims[3]; l++)
	    val(i,j,k,k) = (cval(i,j,k,l) == 1);
    return (status == 0);     
  }


  ///////////////////////////
  // Complex type VarReads //
  ///////////////////////////
//   template<> bool
//   IOVarHDF5<complex<double>,0>::VarRead(complex<double> &val)
//   {
//     herr_t status = H5Dread(DatasetID, ComplexTypeID, MemSpaceID, 
// 			    DiskSpaceID, H5P_DEFAULT, &val);
//     return (status == 0);
//   }  
  
//   template<> bool
//   IOVarHDF5<complex<double>,1>::VarRead(Array<complex<double>,1> &val)
//   {
//     TinyVector<hsize_t,1> hdims;
//     TinyVector<int,1> dims;
//     H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
//     dims = hdims;
//     Array<unsigned char,1> cval(dims);
//     val.resize(dims);
//     herr_t status = H5Dread(DatasetID, ComplexTypeID, MemSpaceID, DiskSpaceID,
// 			     H5P_DEFAULT, val.data());
//     return (status == 0);
//   }
  
//   template<> bool
//   IOVarHDF5<complex<double>,2>::VarRead(Array<complex<double>,2> &val)
//   {
//     TinyVector<hsize_t,2> hdims;
//     TinyVector<int,2> dims;
//     H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
//     dims = hdims;
//     Array<unsigned char,2> cval(dims);
//     val.resize(dims);
//     herr_t status = H5Dread(DatasetID, ComplexTypeID, MemSpaceID, DiskSpaceID,
// 			     H5P_DEFAULT, val.data());
//     return (status == 0); 
//   }
  
//   template<> bool
//   IOVarHDF5<complex<double>,3>::VarRead(Array<complex<double>,3> &val)
//   {
//     TinyVector<hsize_t,3> hdims;
//     TinyVector<int,3> dims;
//     H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
//     dims = hdims;
//     val.resize(dims);
//     herr_t status = H5Dread(DatasetID, ComplexTypeID, MemSpaceID, DiskSpaceID,
// 			     H5P_DEFAULT, val.data());
//     return (status == 0); 
//   }
  
//   template<> bool
//   IOVarHDF5<complex<double>,4>::VarRead(Array<complex<double>,4> &val)
//   {
//     TinyVector<hsize_t,4> hdims;
//     TinyVector<int,4> dims;
//     H5Sget_simple_extent_dims(MemSpaceID, &(hdims[0]), NULL);
//     dims = hdims;
//     val.resize(dims);
//     herr_t status = H5Dread(DatasetID, ComplexTypeID, MemSpaceID, DiskSpaceID,
// 			     H5P_DEFAULT, val.data());
//     return (status == 0);     
//   }


  
  //////////////////////////////////////////
  /// String and bool VarWrite functions ///
  //////////////////////////////////////////

  template<> bool
  IOVarHDF5<string,0>::VarWrite(string val)
  {
    Array<char,1> charArray(val.size()+1);
    charArray=0;
    for (int index=0; index<val.length(); index++)
      charArray(index) = val[index];
    charArray(val.size()) = '\0';
    hid_t type = H5Dget_type(DatasetID);
    // Write the dataset to the file.
    herr_t status = H5Dwrite(DatasetID, type, MemSpaceID, DiskSpaceID, H5P_DEFAULT, 
			     charArray.data());
    H5Tclose(type);
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar.\n";
      return false;
    }
    else
      return true;
  }

  template<> bool
  IOVarHDF5<string,1>::VarWrite(const Array<string,1> &val)
  {
    Array<char,2> charArray(val.extent(0),MAX_HDF5_STRING_LENGTH);
    charArray=0;
    for (int i=0; i<val.extent(0); i++) {
      assert (val(i).length() < (MAX_HDF5_STRING_LENGTH-1));
      for (int s=0; s<val(i).length(); s++)
	charArray(i,s) = val(i)[s];
      // NULL terminate
      int n = val(i).length();
      charArray(i,n) = '\0';
    }
    hid_t type = H5Dget_type(DatasetID);
    // Write the dataset to the file.
    herr_t status = H5Dwrite(DatasetID, type, MemSpaceID, DiskSpaceID, H5P_DEFAULT, 
			     charArray.data());
    H5Tclose(type);
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar.\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<string,2>::VarWrite(const Array<string,2> &val)
  {
    Array<char,3> charArray(val.extent(0),
			    val.extent(1),
			    MAX_HDF5_STRING_LENGTH);
    charArray=0;
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++) {
	assert (val(i,j).length() < (MAX_HDF5_STRING_LENGTH-1));
	for (int s=0; s<val(i,j).length(); s++)
	  charArray(i,j,s) = val(i,j)[s];
      }
    hid_t type = H5Dget_type(DatasetID);
    // Write the dataset to the file.
    herr_t status = H5Dwrite(DatasetID, type, MemSpaceID, DiskSpaceID, H5P_DEFAULT, 
			     charArray.data());
    H5Tclose(type);
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar.\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<string,3>::VarWrite(const Array<string,3> &val)
  {
    Array<char,4> charArray
      (val.extent(0), val.extent(1), val.extent(2),
       MAX_HDF5_STRING_LENGTH);
    charArray=0;
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++)
	for (int k=0; k<val.extent(2); k++) {
	  assert (val(i,j,k).length() < (MAX_HDF5_STRING_LENGTH-1));
	  for (int s=0; s<val(i,j,k).length(); s++)
	    charArray(i,j,k,s) = val(i,j,k)[s];
	}
    hid_t type = H5Dget_type(DatasetID);
    // Write the dataset to the file.
    herr_t status = H5Dwrite(DatasetID, type, MemSpaceID, DiskSpaceID, H5P_DEFAULT, 
			     charArray.data());
    H5Tclose(type);
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar.\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<string,4>::VarWrite(const Array<string,4> &val)
  {
    Array<char,5> charArray
      (val.extent(0), val.extent(1), val.extent(2), val.extent(3),
       MAX_HDF5_STRING_LENGTH);
    charArray=0;
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++)
	for (int k=0; k<val.extent(2); k++)
	  for (int l=0; l<val.extent(3); l++) {
	    assert (val(i,j,k,l).length() < (MAX_HDF5_STRING_LENGTH-1));
	    for (int s=0; s<val(i,j,k,l).length(); s++)
	      charArray(i,j,k,l,s) = val(i,j,k,l)[s];
	  }
    hid_t type = H5Dget_type(DatasetID);
    // Write the dataset to the file.
    herr_t status = H5Dwrite(DatasetID, type, MemSpaceID, DiskSpaceID, H5P_DEFAULT, 
			     charArray.data());
    H5Tclose(type);
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar.\n";
      return false;
    }
    else
      return true;
  }


  template<> bool
  IOVarHDF5<bool,0>::VarWrite(bool val)
  {
    unsigned char cval;
    cval = val ? (unsigned char)1 : (unsigned char)0;
    herr_t status = H5Dwrite(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID, H5P_DEFAULT,
			     &cval);
        if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar(bool).\n";
      return false;
    }
    else
      return true;
  }  
  
  template<> bool
  IOVarHDF5<bool,1>::VarWrite(const Array<bool,1> &val)
  {
    Array<unsigned char,1> cval(val.shape());
    for (int i=0; i<val.extent(0); i++)
      cval(i) = val(i) ? (unsigned char)1 : (unsigned char)0;
    herr_t status = H5Dwrite(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID, H5P_DEFAULT,
			     cval.data());
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar(Array<bool,1>).\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<bool,2>::VarWrite(const Array<bool,2> &val)
  {
    Array<unsigned char,2> cval(val.shape());
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++)
	cval(i,j) = val(i,j) ? (unsigned char)1 : (unsigned char)0;
    herr_t status = H5Dwrite(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID, H5P_DEFAULT,
			     cval.data());
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar(Array<bool,1>).\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<bool,3>::VarWrite(const Array<bool,3> &val)
  {
    Array<unsigned char,3> cval(val.shape());
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++)
	for (int k=0; k<val.extent(2); k++)
	  cval(i,j,k) = val(i,j,k) ? (unsigned char)1 : (unsigned char)0;
    herr_t status = H5Dwrite(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID, H5P_DEFAULT,
			     cval.data());
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar(Array<bool,1>).\n";
      return false;
    }
    else
      return true;
  }
  
  template<> bool
  IOVarHDF5<bool,4>::VarWrite(const Array<bool,4> &val)
  {
    Array<unsigned char,4> cval(val.shape());
    for (int i=0; i<val.extent(0); i++)
      for (int j=0; j<val.extent(1); j++)
	for (int k=0; k<val.extent(2); k++)
	  for (int l=0; l<val.extent(3); l++)
	    cval(i,j,k,l) = val(i,j,k,l) ? (unsigned char)1 : (unsigned char)0;
    herr_t status = H5Dwrite(DatasetID, BoolTypeID, MemSpaceID, DiskSpaceID, H5P_DEFAULT,
			     cval.data());
    if (status < 0) {
      cerr << "Error writing string to HDF5 file in WriteVar(Array<bool,1>).\n";
      return false;
    }
    else
      return true;
    
  }
  


  
  IOVarBase *NewIOVarHDF5(hid_t dataSetID, string name, hid_t boolType, hid_t cmplxType)
  {
    /// First, figure out the rank
    hid_t diskSpaceID = H5Dget_space(dataSetID);
    hid_t memSpaceID  = H5Scopy(diskSpaceID);
    int rank = H5Sget_simple_extent_ndims(diskSpaceID);
    if (rank == 1) {
      hsize_t size, maxSize;
      H5Sget_simple_extent_dims(diskSpaceID, &size, &maxSize);
      if (maxSize == 1)
	rank = 0;
    }
    
    /// First, figure out what type it is.
    hid_t typeID = H5Dget_type(dataSetID);
    H5T_class_t classID = H5Tget_class(typeID);
    H5Tclose (typeID);
    if (classID == H5T_FLOAT) {
      if (rank == 0) 
	return new IOVarHDF5<double,0> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 1) 
	return new IOVarHDF5<double,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 2) 
	return new IOVarHDF5<double,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 3) 
	return new IOVarHDF5<double,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 4) 
	return new IOVarHDF5<double,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 5) 
	return new IOVarHDF5<double,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 6) 
	return new IOVarHDF5<double,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
    }
    if (classID == H5T_INTEGER) {
      if (rank == 0) 
	return new IOVarHDF5<int,0> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 1)
	return new IOVarHDF5<int,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 2)
	return new IOVarHDF5<int,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 3)
	return new IOVarHDF5<int,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 4)
	return new IOVarHDF5<int,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 5)
	return new IOVarHDF5<int,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 6)
	return new IOVarHDF5<int,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
    }
    if (classID == H5T_STRING) {
      if (rank == 0) 
	return new IOVarHDF5<string,0> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 1) 
	return new IOVarHDF5<string,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 2) 
	return new IOVarHDF5<string,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 3) 
	return new IOVarHDF5<string,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 4) 
	return new IOVarHDF5<string,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 5) 
	return new IOVarHDF5<string,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 6) 
	return new IOVarHDF5<string,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
    }
    if (classID == H5T_ENUM){
      if (rank == 0) 
	return new IOVarHDF5<bool,0> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 1) 
	return new IOVarHDF5<bool,1> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 2) 
	return new IOVarHDF5<bool,2> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 3) 
	return new IOVarHDF5<bool,3> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 4) 
	return new IOVarHDF5<bool,4> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 5) 
	return new IOVarHDF5<bool,5> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
      else if (rank == 6) 
	return new IOVarHDF5<bool,6> (name, dataSetID, diskSpaceID, memSpaceID, boolType, cmplxType, true);
    }
    if (classID == H5T_COMPOUND){
      if (rank == 0) 
	return new IOVarHDF5<complex<double>,0> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 1) 
	return new IOVarHDF5<complex<double>,1> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 2) 
	return new IOVarHDF5<complex<double>,2> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 3) 
	return new IOVarHDF5<complex<double>,3> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 4) 
	return new IOVarHDF5<complex<double>,4> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 5) 
	return new IOVarHDF5<complex<double>,5> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
      else if (rank == 6) 
	return new IOVarHDF5<complex<double>,6> (name, dataSetID, diskSpaceID, 
						 memSpaceID, boolType, cmplxType, true);
    }
  }
}
