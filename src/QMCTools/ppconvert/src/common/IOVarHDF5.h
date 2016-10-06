//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Paul R. C. Kent, kentpr@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#ifndef IO_VAR_HDF5_H
#define IO_VAR_HDF5_H

#include "IOVarBase.h"

namespace IO {
  using namespace blitz;
#define MAX_HDF5_STRING_LENGTH 200
  ///////////////////////////////////////////////////////////
  ///                HDF5 Specializations                 ///  
  ///////////////////////////////////////////////////////////
  template<typename T, int RANK> class IOVarHDF5;

//   template<typename T,  typename T0, typename T1, typename T2, typename T3, typename T4,  
// 	   typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
//   class HDF5SliceMaker
//   {
//   public:
//     static const int rank =      ArraySectionInfo<T0>::rank + ArraySectionInfo<T1>::rank + 
//     ArraySectionInfo<T2>::rank + ArraySectionInfo<T3>::rank + ArraySectionInfo<T4>::rank + 
//     ArraySectionInfo<T5>::rank + ArraySectionInfo<T6>::rank + ArraySectionInfo<T7>::rank + 
//     ArraySectionInfo<T8>::rank + ArraySectionInfo<T9>::rank + ArraySectionInfo<T10>::rank;

//     typedef IOVarHDF5<T,rank> SliceType;
//   };

  template<typename T,  int RANK, typename T0, typename T1, typename T2, typename T3, typename T4,  
	   typename T5, typename T6, typename T7, typename T8, typename T9, typename T10>
  class HDF5SliceMaker
  {
  public:
    static const int rank = RANK -
      (SliceCheck<T0>::isSlice + SliceCheck<T1>::isSlice + SliceCheck<T2>::isSlice +
       SliceCheck<T3>::isSlice + SliceCheck<T4>::isSlice + SliceCheck<T5>::isSlice +
       SliceCheck<T6>::isSlice + SliceCheck<T7>::isSlice + SliceCheck<T8>::isSlice +
       SliceCheck<T9>::isSlice + SliceCheck<T10>::isSlice);

    typedef IOVarHDF5<T,rank> SliceType;
  };

  template<typename T, int RANK>
  class IOVarHDF5 : public IOVarBase
  {
  protected:

    bool OwnDataset;
  public:
    hid_t DatasetID, DiskSpaceID, MemSpaceID, BoolTypeID, ComplexTypeID;

    template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	     typename T6, typename T7, typename T8, typename T9, typename T10>
    typename HDF5SliceMaker<T,RANK,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType     
    Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10);
    
    int GetRank();
    IODataType GetType();
    IOFileType GetFileType();

    int GetExtent (int dim);
    void Resize(int n);

    bool VarRead(T &val);
    bool VarRead(blitz::Array<T,RANK> &val);
    template<typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
		 T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10);


    bool VarWrite(T val);
    bool VarWrite(const blitz::Array<T,RANK> &val);
    template<typename TVAL, typename T0, typename T1, typename T2, typename T3, typename T4,
	     typename T5, typename T6, typename T7, typename T8, typename T9,
	     typename T10>
    bool VarWriteSlice(TVAL val,
		       T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
    {
      std::cerr << "val = " << val << std::endl;
      Slice(val, s0,s1,s2,s3,s4,s5,s6,s7,s8,s9,s10).VarWrite(val);
    }
    IOVarHDF5() : OwnDataset(false)
    {

    }
    
    IOVarHDF5( std::string name, hid_t datasetID, hid_t diskSpaceID, hid_t memSpaceID, 
	      hid_t boolTypeID, hid_t complexTypeID, bool ownDataset=false)
    {
      Name          = name;
      DatasetID     = datasetID;
      DiskSpaceID   = diskSpaceID;
      MemSpaceID    = memSpaceID;
      OwnDataset    = ownDataset;
      BoolTypeID    = boolTypeID;
      ComplexTypeID = complexTypeID;
    }
    ~IOVarHDF5()
    {
      H5Sclose(DiskSpaceID);
      H5Sclose(MemSpaceID);
      if (OwnDataset)
	H5Dclose(DatasetID);
    }

  };

  
  ////////////////////////////////////////////////////////////////
  /// Template specializations of reads and writes for std::string  ///
  /// and bools.                                               ///
  ////////////////////////////////////////////////////////////////
  template<> bool IOVarHDF5<std::string,0>::VarRead( std::string &val);
  template<> bool IOVarHDF5<std::string,1>::VarRead(blitz::Array<std::string,1> &val);
  template<> bool IOVarHDF5<std::string,2>::VarRead(blitz::Array<std::string,2> &val);
  template<> bool IOVarHDF5<std::string,3>::VarRead(blitz::Array<std::string,3> &val);
  template<> bool IOVarHDF5<std::string,4>::VarRead(blitz::Array<std::string,4> &val);
  template<> bool IOVarHDF5<bool,  0>::VarRead(bool &val);
  template<> bool IOVarHDF5<bool,  1>::VarRead(blitz::Array<bool,1> &val);
  template<> bool IOVarHDF5<bool,  2>::VarRead(blitz::Array<bool,2> &val);
  template<> bool IOVarHDF5<bool,  3>::VarRead(blitz::Array<bool,3> &val);
  template<> bool IOVarHDF5<bool,  4>::VarRead(blitz::Array<bool,4> &val);
  template<> bool IOVarHDF5<std::string,0>::VarWrite( std::string val);
  template<> bool IOVarHDF5<std::string,1>::VarWrite(const blitz::Array<std::string,1> &val);
  template<> bool IOVarHDF5<std::string,2>::VarWrite(const blitz::Array<std::string,2> &val);
  template<> bool IOVarHDF5<std::string,3>::VarWrite(const blitz::Array<std::string,3> &val);
  template<> bool IOVarHDF5<std::string,4>::VarWrite(const blitz::Array<std::string,4> &val);
  template<> bool IOVarHDF5<bool,  0>::VarWrite(bool val);
  template<> bool IOVarHDF5<bool,  1>::VarWrite(const blitz::Array<bool,1> &val);
  template<> bool IOVarHDF5<bool,  2>::VarWrite(const blitz::Array<bool,2> &val);
  template<> bool IOVarHDF5<bool,  3>::VarWrite(const blitz::Array<bool,3> &val);
  template<> bool IOVarHDF5<bool,  4>::VarWrite(const blitz::Array<bool,4> &val);





  template<typename T, int RANK> void
  IOVarHDF5<T,RANK>::Resize(int n) {
    /// The "+1" fixed compile error for RANK=0
    hsize_t dims[RANK+1];
    H5Sget_simple_extent_dims(DiskSpaceID, dims, NULL);
    if (n < dims[0]) {
      std::cerr << "Cannot resize and HDF5 dataset to a smaller size.\n";
      abort();
    }
    dims[0] = n;
    herr_t status = H5Dextend (DatasetID, dims);
    H5Sclose (DiskSpaceID);
    DiskSpaceID = H5Dget_space(DatasetID);
    H5Sclose(MemSpaceID);
    MemSpaceID = H5Scopy(DiskSpaceID);
  }

  template<typename T, int RANK> int
  IOVarHDF5<T,RANK>::GetExtent(int n) {
    assert (n < RANK);
    /// The "+1" fixed compile error for RANK=0
    hsize_t dims[RANK+1];
    H5Sget_simple_extent_dims(DiskSpaceID, dims, NULL);
    return dims[n];
  }


    
  ////////////////////////////////////////////////////////////////
  /// This function creates a new IOVarHDF5 from a dataset id. /// 
  ////////////////////////////////////////////////////////////////
  IOVarBase *NewIOVarHDF5(hid_t dataSetID, std::string name, hid_t boolType, hid_t cmplxType);
  

  ////////////////////////////////////////////////////////////////
  /// This function creates a new HDF5 dataset from a group id ///
  /// and a variable.  This is used for creating a new         ///
  /// variable for writing.                                    ///
  ////////////////////////////////////////////////////////////////
  template<typename T, int RANK>
  IOVarBase *NewIOVarHDF5(hid_t groupID, std::string name, const blitz::Array<T,RANK> &val,
			  hid_t boolType, hid_t cmplxType)
  {
    /// First, create a new DataSpace.
    TinyVector<hsize_t, RANK> h5Dims, h5MaxDims;
    h5Dims = val.shape();
    h5MaxDims = val.shape();
    h5MaxDims[0] = H5S_UNLIMITED;
    hid_t diskSpaceID = H5Screate_simple(RANK, &(h5Dims[0]), &(h5MaxDims[0]));
    hid_t typeID;
    bool mustCloseType = false;
    if (TypeConvert<T>::Type == DOUBLE_TYPE)
      typeID = H5T_NATIVE_DOUBLE;
    if (TypeConvert<T>::Type == INT_TYPE)
      typeID = H5T_NATIVE_INT;
    if (TypeConvert<T>::Type == STRING_TYPE) {
      typeID = H5Tcopy (H5T_C_S1);
      H5Tset_size (typeID, MAX_HDF5_STRING_LENGTH);
      mustCloseType = true;
    }
    else if (TypeConvert<T>::Type == BOOL_TYPE) {
      typeID = boolType;
    }
    else if (TypeConvert<T>::Type == COMPLEX_TYPE) {
      typeID = cmplxType;
    }
    
    hid_t chunkProp = H5Pcreate(H5P_DATASET_CREATE);
    /// Chunk_dims tells us how big of chunks to allocate in the file at
    /// a time.
    hsize_t chunk_dims[RANK];
    if (RANK < 3)
      chunk_dims[0]= std::min(val.extent(0), 32);
    else 
      chunk_dims[0]=1;
    for (int i=1;i<RANK;i++)
      chunk_dims[i]=val.extent(i);

    herr_t status = H5Pset_chunk(chunkProp, RANK, chunk_dims);  
    assert(status>=0);

    hid_t datasetID = 
      H5Dcreate (groupID, name.c_str(), typeID, diskSpaceID, chunkProp);

//     if (RANK == 1) {
//       hsize_t attrDim = 1;
//       hid_t attrSpace = H5Screate_simple(1, &attrDim, &attrDim);
//       hid_t attrID = H5Acreate(datasetID, "dim", H5T_NATIVE_INT,  attrSpace, H5P_DEFAULT);
//       int dim = 1;
//       herr_t status = H5Awrite(attrID, H5T_NATIVE_INT, &dim);
//       assert (status == 0);
//       H5Aclose (attrID);
//       H5Sclose (attrSpace);
//     }

    if (mustCloseType)
      H5Tclose (typeID);
    H5Sclose(diskSpaceID);

    IOVarHDF5<T,RANK> *newVar = 
      dynamic_cast<IOVarHDF5<T,RANK>*> (NewIOVarHDF5(datasetID, name, boolType, cmplxType));
    if (newVar == NULL) {
      std::cerr << "Error in dynamic_cast in NewIOVarHDF5 #1.\n";
      abort();
    }
    newVar->VarWrite(val);
    
    return newVar;
  }


  ////////////////////////////////////////////////////////////////
  /// This function creates a new HDF5 dataset from a group id ///
  /// and a variable.  This is used for creating a new         ///
  /// variable for writing.                                    ///
  ////////////////////////////////////////////////////////////////
  template<typename T> inline
  IOVarBase *NewIOVar0HDF5(hid_t groupID, std::string name, T val,
			   hid_t boolType, hid_t cmplxType)
  {
    /// First, create a new DataSpace.
    hsize_t h5Dims, h5MaxDims;
    h5Dims    = 1;
    h5MaxDims = 1;
    hid_t diskSpaceID = H5Screate_simple(1, &h5Dims, &h5MaxDims);
    hid_t typeID;
    bool mustCloseType = false;
    if (TypeConvert<T>::Type == DOUBLE_TYPE)
      typeID = H5T_NATIVE_DOUBLE;
    if (TypeConvert<T>::Type == INT_TYPE)
      typeID = H5T_NATIVE_INT;
    if (TypeConvert<T>::Type == STRING_TYPE) {
      typeID = H5Tcopy (H5T_C_S1);
      H5Tset_size (typeID, MAX_HDF5_STRING_LENGTH+1);
      mustCloseType = true;
    }
    else if (TypeConvert<T>::Type == BOOL_TYPE) {
      typeID = boolType;
    }
    else if (TypeConvert<T>::Type == COMPLEX_TYPE) {
      typeID = cmplxType;
    }
    
    hid_t datasetID = 
      H5Dcreate (groupID, name.c_str(), typeID, diskSpaceID, H5P_DEFAULT);

//     hsize_t attrDim = 1;
//     hid_t attrSpace = H5Screate_simple(1, &attrDim, &attrDim);
//     hid_t attrID = H5Acreate(datasetID, "dim", H5T_NATIVE_INT,  attrSpace, H5P_DEFAULT);
//     int dim = 0;
//     herr_t status = H5Awrite(attrID, H5T_NATIVE_INT, &dim);
//     assert (status == 0);
//     H5Aclose (attrID);
//     H5Sclose (attrSpace);

    if (mustCloseType)
      H5Tclose (typeID);
    H5Sclose(diskSpaceID);

    IOVarHDF5<T,0> *newVar =  dynamic_cast<IOVarHDF5<T,0>*> (NewIOVarHDF5(datasetID, name, boolType, cmplxType));
    if (newVar == NULL) {
      std::cerr << "Error in dynamic_cast in NewIOVarHDF5 #2.\n";
      abort();
    }
    newVar->VarWrite(val);
    
    return newVar;
  }


  ////////////////////////////////////////////////////////////////
  /// String specialization of above                           ///
  ////////////////////////////////////////////////////////////////
  template<> inline
  IOVarBase *NewIOVar0HDF5(hid_t groupID, std::string name, std::string val,
			   hid_t boolType, hid_t cmplxType)
  {
    /// First, create a new DataSpace.
    hsize_t h5Dims, h5MaxDims;
    h5Dims    = 1;
    h5MaxDims = 1;
    hid_t diskSpaceID = H5Screate_simple(1, &h5Dims, &h5MaxDims);
    hid_t typeID;
    typeID = H5Tcopy (H5T_C_S1);
    H5Tset_size (typeID, val.length()+1);
    
    hid_t datasetID = 
      H5Dcreate (groupID, name.c_str(), typeID, diskSpaceID, H5P_DEFAULT);

    H5Tclose (typeID);
    H5Sclose(diskSpaceID);

    IOVarHDF5<std::string,0> *newVar = dynamic_cast<IOVarHDF5<std::string,0>*>
      (NewIOVarHDF5(datasetID, name, boolType, cmplxType));
    if (newVar == NULL) {
      std::cerr << "Error in dynamic_cast in NewIOVarHDF5 #3.\n";
      abort();
    }
    newVar->VarWrite(val);
    
    return newVar;
  }



  template<typename T, int RANK> 
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10> bool
  IOVarHDF5<T,RANK>::VarRead(typename SliceInfo<T,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::T_slice &val,
			     T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
  {
    return Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarRead(val);
  }


//   template<typename T, int RANK> 
//   template<typename TVAL, typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
// 	   typename T6, typename T7, typename T8, typename T9, typename T10> bool
//   IOVarHDF5<T,RANK>::VarWriteSlice(TVAL val, T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T8 s9, T10 s10)
//   {
//     return Slice(s0, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10).VarWrite(val);
//   }

  template<typename T>
  class RangeConvert
  {
  public:
    static Range GetRange (T r) {
      return Range();
    }
  };
  
  template<>
  class RangeConvert<int>
  {
  public:
    static Range GetRange(int r)
    {
      return Range(r);
    }
  };

  template<>
  class RangeConvert<Range>
  {
  public:
    static Range GetRange(Range r)
    {
      return Range(r);
    }
  };

  template<class T, int RANK>
  template<typename T0, typename T1, typename T2, typename T3, typename T4, typename T5,
	   typename T6, typename T7, typename T8, typename T9, typename T10>
  typename HDF5SliceMaker<T,RANK,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType 
  IOVarHDF5<T,RANK>::Slice(T0 s0, T1 s1, T2 s2, T3 s3, T4 s4, T5 s5, T6 s6, T7 s7, T8 s8, T9 s9, T10 s10)
  {
    typedef typename HDF5SliceMaker<T,RANK,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::SliceType newSliceType;
    newSliceType newVar;

    newVar.DatasetID = DatasetID;
    newVar.DiskSpaceID = H5Scopy(H5Dget_space(DatasetID));
    newVar.BoolTypeID = BoolTypeID;
    newVar.ComplexTypeID = ComplexTypeID;
  
    hsize_t start[RANK], count[RANK], stride[RANK], dims[RANK], maxdims[RANK];
    hsize_t memDims[HDF5SliceMaker<T,RANK,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::rank];
    H5Sget_simple_extent_dims(newVar.DiskSpaceID, dims, maxdims);
  
    int memDimsIndex=0;
  
  
    /// Select the disk space hyperslab
    if (RANK > 0) {
      Range r0 = RangeConvert<T0>::GetRange(s0);
      start[0]  = r0.first(0);
      count[0]  = (r0.last(dims[0]-1)-start[0])/r0.stride() + 1;
      stride[0] = r0.stride();
      if (ArraySectionInfo<T0>::rank==1) {
	memDims[memDimsIndex]=count[0];
	memDimsIndex++;
      }
    }
    if (RANK > 1) {
      Range r1 = RangeConvert<T1>::GetRange(s1);
      start[1] = r1.first(0);
      count[1] = (r1.last(dims[1]-1)-start[1])/r1.stride() + 1;
      stride[1] = r1.stride();
      if (ArraySectionInfo<T1>::rank==1) {
	memDims[memDimsIndex]=count[1];
	memDimsIndex++;
      }
    }
    if (RANK > 2) {
      Range r2 = RangeConvert<T2>::GetRange(s2);
      start[2] = r2.first(0);
      count[2] = (r2.last(dims[2]-1)-start[2])/r2.stride() + 1;
      stride[2] = r2.stride();
      if (ArraySectionInfo<T2>::rank==1) {
	memDims[memDimsIndex]=count[2];
	memDimsIndex++;
      }
    }
    if (RANK > 3) {
      Range r3 = RangeConvert<T3>::GetRange(s3);
      start[3] = r3.first(0);
      count[3] = (r3.last(dims[3]-1)-start[3])/r3.stride() + 1;
      stride[3] = r3.stride();
      if (ArraySectionInfo<T3>::rank==1) {
	memDims[memDimsIndex]=count[3];
	memDimsIndex++;
      }
    }
    if (RANK > 4) {
      Range r4 = RangeConvert<T4>::GetRange(s4);
      start[4] = r4.first(0);
      count[4] = (r4.last(dims[4]-1)-start[4])/r4.stride() + 1;
      stride[4] = r4.stride();
      if (ArraySectionInfo<T4>::rank==1) {
	memDims[memDimsIndex]=count[4];
	memDimsIndex++;
      }
    }
    if (RANK > 5) {
      Range r5 = RangeConvert<T5>::GetRange(s5);
      start[5] = r5.first(0);
      count[5] = (r5.last(dims[5]-1)-start[5])/r5.stride() + 1;
      stride[5] = r5.stride();
      if (ArraySectionInfo<T5>::rank==1) {
	memDims[memDimsIndex]=count[5];
	memDimsIndex++;
      }
    }
    if (RANK > 6) {
      Range r6 = RangeConvert<T6>::GetRange(s6);
      start[6] = r6.first(0);
      count[6] = (r6.last(dims[6]-1)-start[6])/r6.stride() + 1;
      stride[6] = r6.stride();
      if (ArraySectionInfo<T6>::rank==1) {
	memDims[memDimsIndex]=count[6];
	memDimsIndex++;
      }
    }
    if (RANK > 7) {
      Range r7 = RangeConvert<T7>::GetRange(s7);
      start[7] = r7.first(0);
      count[7] = (r7.last(dims[7]-1)-start[7])/r7.stride() + 1;
      stride[7] = r7.stride();
      if (ArraySectionInfo<T7>::rank==1) {
	memDims[memDimsIndex]=count[7];
	memDimsIndex++;
      }
    }
    if (RANK > 8) {
      Range r8 = RangeConvert<T8>::GetRange(s8);
      start[8] = r8.first(0);
      count[8] = (r8.last(dims[8]-1)-start[8])/r8.stride() + 1;
      stride[8] = r8.stride();
      if (ArraySectionInfo<T8>::rank==1) {
	memDims[memDimsIndex]=count[8];
	memDimsIndex++;
      }
    }
    if (RANK > 9) {
      Range r9 = RangeConvert<T9>::GetRange(s9);
      start[9] = r9.first(0);
      count[9] = (r9.last(dims[9]-1)-start[9])/r9.stride() + 1;
      stride[9] = r9.stride();
      if (ArraySectionInfo<T9>::rank==1) {
	memDims[memDimsIndex]=count[9];
	memDimsIndex++;
      }
    }
    if (RANK > 10) {
      Range r10 = RangeConvert<T10>::GetRange(s10);
      start[10] = r10.first(0);
      count[10] = (r10.last(dims[10]-1)-start[10])/r10.stride() + 1;
      stride[10] = r10.stride();
      if (ArraySectionInfo<T10>::rank==1) {
	memDims[memDimsIndex]=count[10];
	memDimsIndex++;
      }
    }
    newVar.MemSpaceID = H5Screate_simple(HDF5SliceMaker<T,RANK,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10>::rank,
					 memDims, memDims);
    H5Sselect_hyperslab(newVar.DiskSpaceID, H5S_SELECT_SET, start, stride, count, NULL);
    return newVar;
  }

  /// This routine should cover double and int types.  Strings and bools
  /// need to be handled explicitly
  template<class T, int RANK> bool
  IOVarHDF5<T,RANK>::VarRead(T &val)
  {
    assert (RANK == 0);
    IODataType dataType = TypeConvert<T>::Type;
    hid_t memType;
    if      (dataType == DOUBLE_TYPE)  memType = H5T_NATIVE_DOUBLE;
    else if (dataType == INT_TYPE)     memType = H5T_NATIVE_INT;
    else if (dataType == COMPLEX_TYPE) memType = ComplexTypeID;
    else {
      T a;
      std::cerr << "Unknown data type in IOVarHDF5<" << TypeString(a) << ", " 
	   << RANK << ">" << std::endl;
    }

    /// Resize val to appropriate size
    hsize_t h5dim;
    H5Sget_simple_extent_dims(MemSpaceID, &h5dim, NULL);
    assert(h5dim == 1);
   
    /// Now, call HDF5 to do the actual reading.
    herr_t status = 
      H5Dread (DatasetID, memType, MemSpaceID, DiskSpaceID, H5P_DEFAULT, &val);
    
    return (status == 0);
  }


  /// This routine should cover double and int types.  Strings and bools
  /// need to be handled explicitly
  template<class T, int RANK> bool
  IOVarHDF5<T,RANK>::VarRead(blitz::Array<T,RANK> &val)
  {
    IODataType dataType = TypeConvert<T>::Type;
    hid_t memType;
    if      (dataType == DOUBLE_TYPE) memType = H5T_NATIVE_DOUBLE;
    else if (dataType == INT_TYPE)    memType = H5T_NATIVE_INT;
    else if (dataType == COMPLEX_TYPE) memType = ComplexTypeID;
    else {
      T a;
      std::cerr << "Unknown data type in IOVarHDF5<" << TypeString(a) << ", " 
	   << RANK << ">" << std::endl;
    }

    /// Resize val to appropriate size
    TinyVector<hsize_t,RANK> h5dims;
    TinyVector<int,RANK> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(h5dims[0]), NULL);
    dims = h5dims;
    //    std::cerr << "resizing val to " << dims << std::endl;
    val.resize(dims);
   
    /// Now, call HDF5 to do the actual reading.
    herr_t status = 
      H5Dread (DatasetID, memType, MemSpaceID, DiskSpaceID, H5P_DEFAULT, val.data());
    
    return (status == 0);
  }


  /// This routine should cover double and int types.  Strings and bools
  /// need to be handled explicitly
  template<class T, int RANK> bool
  IOVarHDF5<T,RANK>::VarWrite(T val)
  {
    /// Check to see if we have the write dimensions.
    assert (H5Sget_simple_extent_ndims(MemSpaceID) == 1);
    hsize_t dims;
    H5Sget_simple_extent_dims(MemSpaceID, &dims, NULL);
    assert (dims == 1);

    IODataType dataType = TypeConvert<T>::Type;
    hid_t memType;
    if      (dataType == DOUBLE_TYPE) memType = H5T_NATIVE_DOUBLE;
    else if (dataType == INT_TYPE)    memType = H5T_NATIVE_INT;
    else if (dataType == COMPLEX_TYPE) memType = ComplexTypeID;
    else {
      T a;
      std::cerr << "Unknown data type in IOVarHDF5<" << TypeString(a) << ", 0>" << std::endl;
    }
    /// Now, call HDF5 to do the actual writing.
    herr_t status = 
      H5Dwrite (DatasetID, memType, MemSpaceID, DiskSpaceID, H5P_DEFAULT, &val);
    
    return (status == 0);
  }



  /// This routine should cover double and int types.  Strings and bools
  /// need to be handled explicitly
  template<class T, int RANK> bool
  IOVarHDF5<T,RANK>::VarWrite(const blitz::Array<T,RANK> &val)
  {
    /// Check to see if we have the write dimensions.
    assert (H5Sget_simple_extent_ndims(MemSpaceID) == RANK);
    TinyVector<hsize_t,RANK> dims;
    H5Sget_simple_extent_dims(MemSpaceID, &(dims[0]), NULL);			      
    for (int i=0; i < RANK; i++)
      assert (dims[i] == val.extent(i));

    IODataType dataType = TypeConvert<T>::Type;
    hid_t memType;
    if      (dataType == DOUBLE_TYPE) memType = H5T_NATIVE_DOUBLE;
    else if (dataType == INT_TYPE)    memType = H5T_NATIVE_INT;
    else if (dataType == COMPLEX_TYPE) memType = ComplexTypeID;
    else {
      T a;
      std::cerr << "Unknown data type in IOVarHDF5<" << TypeString(a) << ", " 
	   << RANK << ">" << std::endl;
    }
    /// Now, call HDF5 to do the actual writing.
    herr_t status = 
      H5Dwrite (DatasetID, memType, MemSpaceID, DiskSpaceID, H5P_DEFAULT, val.data());
    
    return (status == 0);
  }


  template<typename T, int RANK> IOFileType
  IOVarHDF5<T,RANK>::GetFileType() 
  { return HDF5_TYPE; }
  
  template<typename T, int RANK> int 
  IOVarHDF5<T,RANK>::GetRank()
  { return RANK; }

  template<typename T, int RANK> IODataType 
  IOVarHDF5<T,RANK>::GetType()
  { 
    return TypeConvert<T>::Type;
  }

  IOVarBase *NewIOVarHDF5(hid_t dataSetID, hid_t boolType, hid_t cmplxType);



}      // namespace IO

#endif // ifndef IO_VAR_HDF5
