
#ifndef QMCPLUSPLUS_AFQMC_SMSPARSEMATRIX_H
#define QMCPLUSPLUS_AFQMC_SMSPARSEMATRIX_H

#include<tuple>
#include<algorithm>
#include"AFQMC/Utilities/tuple_iterator.hpp"
#include<iostream>
#include<vector>
#include<assert.h>
#include"AFQMC/config.0.h"
#include"Utilities/UtilityFunctions.h"
#include<mpi.h>
#include <sys/time.h>
#include <ctime>

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/exceptions.hpp>

#define ASSERT_SPARSEMATRIX 

#if defined(USE_EIGEN)
namespace qmcplusplus
{
}

#else  // In this case, use OhhmsPETE and your sparse matrix class

namespace qmcplusplus
{

// NOTE: In principle, it is possible to use different integer types for row/col and rowIndex.
// This allows to store very large matrices while keeping the storage associated with row/col
// to a minimum. In practice this doesn't work unless I have sparse blas routines that are capable
// of using mixed int precision for column and rowIndex arrays.
// As a result, I'm keeping only 1 type of integer. If it is set to long, then MKL must be compiled

// class that implements a sparse matrix in CSR format
template<class T>
class SMSparseMatrix
{
  public:

  template<typename spT> using ShmemAllocator = boost::interprocess::allocator<spT, boost::interprocess::managed_shared_memory::segment_manager>;
  template<typename spT> using SMVector = boost::interprocess::vector<spT, ShmemAllocator<spT>>;

  typedef T                Type_t;
  typedef T                value_type;
  typedef T*               pointer;
  typedef const T*         const_pointer;
  typedef const int*       const_intPtr;
  typedef int              intType;
  typedef int*             intPtr;
//  typedef const long*  const_intPtr;
//  typedef long    intType;
//  typedef long*   intPtr;
  typedef typename SMVector<T>::iterator iterator;
  typedef typename SMVector<T>::const_iterator const_iterator;
  typedef typename SMVector<intType>::iterator int_iterator;
  typedef typename SMVector<intType>::const_iterator const_int_iterator;
  typedef typename SMVector<intType>::iterator indx_iterator;
  typedef typename SMVector<intType>::const_iterator const_indx_iterator;
  typedef SMSparseMatrix<T>  This_t;

  SMSparseMatrix<T>():nr(0),nc(0),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  SMSparseMatrix<T>(int n):nr(n),nc(n),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL)
  {
    remover.ID="NULL";
    remover.head=false;
  }

  SMSparseMatrix<T>(int n,int m):nr(n),nc(m),compressed(false),zero_based(true),head(false),ID(""),SMallocated(false),storage_format(0),vals(NULL),rowIndex(NULL),myrows(NULL),colms(NULL),mutex(NULL),share_buff(NULL),segment(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  ~SMSparseMatrix<T>()
  {
    if(segment != NULL) {
     delete segment;
     boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  }

  // this should probably be disallowed
  SMSparseMatrix(const SMSparseMatrix<T> &rhs)
  {
    APP_ABORT(" Error: SMSparseMatrix copy constructor has been disabled.");
  }

  inline void setup(bool hd, std::string ii, MPI_Comm comm_) {
    head=hd;
    ID=ii;
    remover.ID=ii;
    remover.head=hd;
    comm=comm_;
  }

  inline void barrier() {
    MPI_Barrier(comm);
  }

  inline void reserve(std::size_t n, bool allow_reduce = false)
  {
    assert(n<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
    if(vals==NULL || (vals!=NULL && vals->capacity() < n) || (vals!=NULL && vals->capacity() > n && allow_reduce)) 
      allocate(n,allow_reduce);
    if(head) {
      assert(n < vals->max_size());
      assert(n < colms->max_size());
      assert((nr+1) < rowIndex->max_size());
      vals->reserve(n);
      myrows->reserve(n);
      colms->reserve(n); 
      rowIndex->resize(std::max(nr,nc)+1);
      share_buff->resize(1000*sizeof(std::complex<double>));
    }
    barrier();
  }

  // does not allow grow/shrink
  inline bool allocate_serial(std::size_t n)
  {
    assert(n<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
    if(!head) { SMallocated=true; return true; }
    if(vals!=NULL && vals->capacity() >= n) { SMallocated=true; return true; }

    memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+(2*n+2*std::max(nr,nc)+2)*sizeof(int)+1000*sizeof(std::complex<double>)+100000;

    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory); 
    } catch(boost::interprocess::interprocess_exception &ex) {
      std::cout<<"\n Found managed_shared_memory segment, removing. CAREFUL WITH PERSISTENT SHM MEMORY. !!! \n\n";
      boost::interprocess::shared_memory_object::remove(ID.c_str());
      segment=NULL;
    }

    if(segment==NULL) {
      try {
        segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory); 
      } catch(boost::interprocess::interprocess_exception &ex) {
        std::cerr<<"Problems setting up managed_shared_memory in SMSparseMatrix." <<std::endl;
        APP_ABORT(" Error: Problems setting up managed_shared_memory in SMSparseMatrix. \n");
        return false;
      }
    }

    try {
          
      alloc_ulong = new ShmemAllocator<intType>(segment->get_segment_manager());
      alloc_int = new ShmemAllocator<intType>(segment->get_segment_manager());
      alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());
      alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());
          
      mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();
          
      share_buff = segment->construct<SMVector<unsigned char>>("share_buff")(*alloc_uchar);
      share_buff->resize(1000*sizeof(std::complex<double>));

          
      rowIndex = segment->construct<SMVector<intType>>("rowIndex")(*alloc_ulong);
      assert((nr+1) < rowIndex->max_size());
      rowIndex->resize(std::max(nr,nc)+1);
          
      myrows = segment->construct<SMVector<int>>("myrows")(*alloc_int);
      assert(n < myrows->max_size());
      myrows->reserve(n);
          
      colms = segment->construct<SMVector<int>>("colms")(*alloc_int);
      assert(n < colms->max_size());
      colms->reserve(n);
          
      vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
      assert(n < vals->max_size());
      vals->reserve(n);

    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMSparseMatrix." <<std::endl;
      APP_ABORT(" Error: Problems setting up managed_shared_memory in SMSparseMatrix. \n");
      return false;
    }
    SMallocated=true;
    return true;
  }

  // all processes must call this routine
  inline bool allocate(std::size_t n, bool allow_reduce=false)
  {
    assert(n<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
    bool grow = false;
    std::size_t old_sz = (segment==NULL)?0:(segment->get_size());
    if(SMallocated) {
      if(vals!=NULL && vals->capacity() >= n && !allow_reduce) return true;
      grow = true; 
      if(!head) {  // delay delete call on head in case you need to shrink vector
        delete segment;
        segment=NULL;
      }
    }
    barrier();    
    if(head) { 
      memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+(2*n+2*std::max(nr,nc)+2)*sizeof(int)+1000*sizeof(std::complex<double>)+100000;

      if(grow) {
        if(memory > old_sz) {
          std::size_t extra = memory - old_sz;
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::grow(ID.c_str(), extra)) {
            std::cerr<<" Error growing shared memory in  SMSparseMatrix::allocate(). \n";
            APP_ABORT(" Error: growing shared memory in  SMSparseMatrix::allocate(). \n"); 
            return false;
          }
        } else {
          segment->destroy<SMVector<T>>("vals");
          segment->destroy<SMVector<int>>("myrows");
          segment->destroy<SMVector<int>>("colms");
          myrows = segment->construct<SMVector<int>>("myrows")(*alloc_T);
          assert(n < myrows->max_size());
          myrows->reserve(n);
          colms = segment->construct<SMVector<int>>("colms")(*alloc_T);
          assert(n < colms->max_size());
          colms->reserve(n);
          vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
          assert(n < vals->max_size());
          vals->reserve(n);
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::shrink_to_fit(ID.c_str())) {
            std::cerr<<" Error in shrink_to_fit shared memory in SMSparseMatrix::allocate(). \n";
            APP_ABORT(" Error: in shrink_to_fit shared memory in SMSparseMatrix::allocate(). \n"); 
            return false;
          }
        }

        try {
          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str());
          vals = segment->find<SMVector<T>>("vals").first;
          colms = segment->find<SMVector<int>>("colms").first;
          myrows = segment->find<SMVector<int>>("myrows").first;
          rowIndex = segment->find<SMVector<intType>>("rowIndex").first;
          share_buff = segment->find<SMVector<unsigned char>>("share_buff").first;
          mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
          assert(share_buff != 0);
          assert(mutex != 0);
          assert(vals != 0);
          assert(myrows != 0);
          assert(colms != 0);
          assert(rowIndex != 0);
          share_buff->resize(1000*sizeof(std::complex<double>)); 
          assert((std::max(nr,nc)+1) < rowIndex->max_size());
          rowIndex->resize(std::max(nr,nc)+1);
          assert(n < myrows->max_size());
          myrows->reserve(n);
          assert(n < colms->max_size());
          colms->reserve(n);
          assert(n < vals->max_size());
          vals->reserve(n);  
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems opening shared memory in SMDenseVector::allocate() ." <<std::endl;
          APP_ABORT(" Error: opening shared memory in SMDenseVector::allocate() ."); 
          return false;
        }

      } else { // grow

        try {
          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
        } catch(boost::interprocess::interprocess_exception &ex) {
          std::cout<<"\n Found managed_shared_memory segment, removing. CAREFUL WITH PERSISTENT SHM MEMORY. !!! \n\n";
          boost::interprocess::shared_memory_object::remove(ID.c_str());
          segment=NULL;
        }
    
        if(segment==NULL) {
          try {
            segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
          } catch(boost::interprocess::interprocess_exception &ex) {
            std::cerr<<"Problems setting up managed_shared_memory in SMSparseMatrix." <<std::endl;
            APP_ABORT(" Error: Problems setting up managed_shared_memory in SMSparseMatrix.");
            return false;
          }
        }

        try {

          alloc_ulong = new ShmemAllocator<intType>(segment->get_segment_manager());
          alloc_int = new ShmemAllocator<intType>(segment->get_segment_manager());
          alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());
          alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());

          mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();

          share_buff = segment->construct<SMVector<unsigned char>>("share_buff")(*alloc_uchar);
          share_buff->resize(1000*sizeof(std::complex<double>));

          rowIndex = segment->construct<SMVector<intType>>("rowIndex")(*alloc_ulong);
          assert((std::max(nr,nc)+1) < rowIndex->max_size());
          rowIndex->resize(std::max(nr,nc)+1);

          myrows = segment->construct<SMVector<int>>("myrows")(*alloc_int);
          assert(n < myrows->max_size());
          myrows->reserve(n);

          colms = segment->construct<SMVector<int>>("colms")(*alloc_int);
          assert(n < colms->max_size());
          colms->reserve(n);

          vals = segment->construct<SMVector<T>>("vals")(*alloc_T);
          assert(n < vals->max_size());
          vals->reserve(n);

        } catch(std::bad_alloc& ex) {
          std::cerr<<"Problems allocating shared memory in SMSparseMatrix." <<std::endl;
          APP_ABORT(" Error: Problems allocating shared memory in SMSparseMatrix."); 
          return false;
        }
      }
    }
    barrier();
    SMallocated=true;
    initializeChildren();
    return true;
  }

  // only call this when all arrays have been allocated and modified  
  inline bool initializeChildren()
  { 
    if(head) return true;
    // delete segment in case this routine is called multiple times.
    // SHM is not removed, just the mapping of the local process. 
    if(segment!=NULL) {
      delete segment;
      segment=NULL;
    }    
    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str()); 
      vals = segment->find<SMVector<T>>("vals").first;
      colms = segment->find<SMVector<int>>("colms").first;
      myrows = segment->find<SMVector<int>>("myrows").first;
      rowIndex = segment->find<SMVector<intType>>("rowIndex").first;
      share_buff = segment->find<SMVector<unsigned char>>("share_buff").first;
      mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
      assert(share_buff != 0);
      assert(mutex != 0);
      assert(vals != 0);
      assert(myrows != 0);
      assert(colms != 0);
      assert(rowIndex != 0);
    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMSparseMatrix: initializeChildren() ." <<std::endl;
      APP_ABORT("Problems allocating shared memory in SMSparseMatrix: initializeChildren() .");
      return false;
    }
    return true;
  }

  // does not allow grow/shrink, aborts if resizing beyond capacity
  inline void resize_serial(std::size_t nnz)
  {
    assert(nnz<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
    if(!head) return;
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz))
      APP_ABORT(" Error: Call to SMSparseMatrix::resize_serial without enough capacity. \n");
    if(head) {
      assert(nnz < vals->max_size());
      vals->resize(nnz);
      assert(nnz < myrows->max_size());
      myrows->resize(nnz);
      assert(nnz < colms->max_size());
      colms->resize(nnz);
      assert((nr+1) < rowIndex->max_size());
      rowIndex->resize(nr+1);
    }
  } 

  // this routine does not preserve information when allow_reduce=true  
  inline void resize(std::size_t nnz, bool allow_reduce=false)
  {
    assert(nnz<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) ) {
      allocate(nnz,allow_reduce);
    } else if(vals!=NULL && vals->capacity() > nnz && allow_reduce) {
      // not keeping information if reduction occurs 
      allocate(nnz,allow_reduce);
    }
    if(head) {
      assert(nnz < vals->max_size());
      assert(nnz < colms->max_size());
      assert(nnz < myrows->max_size());
      assert((nr+1) < rowIndex->max_size());
      vals->resize(nnz);
      myrows->resize(nnz);
      colms->resize(nnz);
      rowIndex->resize(nr+1);
    }
    barrier();
  }

  inline void clear() { 
    compressed=false;
    zero_based=true;
    if(!SMallocated) return; 
    if(!head) return;
    vals->clear();
    colms->clear();
    myrows->clear();
    rowIndex->clear();
  }

  inline void setDims(int n, int m)
  {
    nr=n;
    nc=m;
    compressed=false;
    zero_based=true;
    clear();
  }

  template<class T1>
  inline bool copyToBuffer(std::vector<T1>& v)
  {
    if( v.size()*sizeof(T1) > share_buff->capacity() )
      return false; 
    {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      std::memcpy( share_buff->data(), v.data(), v.size()*sizeof(T1) );
    }
  }
 
  template<class T1>
  inline bool copyFromBuffer(std::vector<T1>& v)
  {
    if( v.size()*sizeof(T1) > share_buff->capacity() )
      return false; 
    std::memcpy( v.data(), share_buff->data(), v.size()*sizeof(T1) );
  }  

  template<typename T1>
  inline void share(T1* x, int n, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    assert( sizeof(T1)*n < sizeof(unsigned char)*share_buff->size() );
    if(sender) {
      std::memcpy(&((*share_buff)[0]),x,sizeof(T1)*n);
      barrier();
    } else {
      barrier();
      std::memcpy(x,&((*share_buff)[0]),sizeof(T1)*n);
    }
    barrier();
  }

  template<typename T1>
  inline void share(T1& x, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    if(sender) {
      std::memcpy(&((*share_buff)[0]),&x,sizeof(T1));
      barrier();
    } else {
      barrier();
      std::memcpy(&x,&((*share_buff)[0]),sizeof(T1));
    }
    barrier();
  }

  template<typename T1>
  inline void share(std::vector<T1>& x, int n, bool sender) {
    if(!SMallocated)
      APP_ABORT("Error: Call to SMDenseVector::share with unallocated object. \n");
    assert( sizeof(T1)*n < sizeof(unsigned char)*share_buff->size() );
    assert( x.size() >= n);
    if(sender) {
      std::memcpy(&((*share_buff)[0]),x.data(),sizeof(T1)*n);
      barrier();
    } else {
      barrier();
      std::memcpy(x.data(),&((*share_buff)[0]),sizeof(T1)*n);
    }
    barrier();
  }

  inline void setCompressed() 
  {
    compressed=true;
  }

  inline bool isCompressed() const
  {
    return compressed;
  }

  inline std::size_t memoryUsage() { return memory; }

  inline std::size_t capacity() const
  {
    return (vals!=NULL)?(vals->capacity()):0;
  }
  inline std::size_t size() const
  {
    return (vals!=NULL)?(vals->size()):0;
  }
  inline int rows() const
  {
    return nr;
  }
  inline int cols() const
  {
    return nc;
  }

  inline const_pointer values(intType n=intType(0)) const 
  {
    return (vals!=NULL)?(&((*vals)[n])):NULL;
  }

  inline pointer values(intType n=intType(0)) 
  {
    return (vals!=NULL)?(&((*vals)[n])):NULL;
  }

  inline const_intPtr column_data(intType n=intType(0)) const 
  {
    return (colms!=NULL)?(&((*colms)[n])):NULL;
  }
  inline intPtr column_data(intType n=intType(0)) 
  {
    return (colms!=NULL)?(&((*colms)[n])):NULL;
  }

  inline const_intPtr row_data(intType n=intType(0)) const 
  {
    return (myrows!=NULL)?(&((*myrows)[n])):NULL;
  }
  inline intPtr row_data(intType n=intType(0)) 
  {
    return (myrows!=NULL)?(&((*myrows)[n])):NULL;
  }

  inline const_intPtr row_index(intType n=intType(0)) const 
  {
    return (rowIndex!=NULL)?(&((*rowIndex)[n])):NULL;
  }
  inline intPtr row_index(intType n=intType(0)) 
  {
    return  (rowIndex!=NULL)?(&((*rowIndex)[n])):NULL;
  }

  inline This_t& operator=(const SMSparseMatrix<T> &rhs) 
  { 
    APP_ABORT(" Error: SMSharedMatrix::operator= is disabled. \n");
    compressed=rhs.compressed;
    zero_based=rhs.zero_based;
    nr=rhs.nr;
    nc=rhs.nc;
    if(!head) return *this;
    (*vals)=*(rhs.vals);
    (*myrows)=*(rhs.myrows);
    (*colms)=*(rhs.colms);
    (*rowIndex)=*(rhs.rowIndex);
    return *this;
  }  

  // use binary search PLEASE!!! Not really used anyway
  inline int find_element(int i, int j) const {
    for (int k = (*rowIndex)[i]; k<(*rowIndex)[i+1]; k++) {
      if ((*colms)[k] == j) return k;
    }
    return -1;
  }

  inline Type_t operator()( int i, int j) const
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc && compressed); 
#endif
    int idx = find_element(i,j);
    if (idx == -1) return T(0);
    return (*vals)[idx]; 
  }

  inline void add(const int i, const int j, const T& v, bool needs_locks=false) 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<nr && j>=0 && j<nc);
#endif
    compressed=false;
    if(needs_locks) {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      myrows->push_back(i);
      colms->push_back(j);
      vals->push_back(v);
    } else {
      if(!head) return;
      myrows->push_back(i);
      colms->push_back(j);
      vals->push_back(v);
    }
    assert(vals->size()<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
  }

  inline void add(const std::vector<std::tuple<int,int,T>>& v, bool needs_locks=false)
  {
    compressed=false;
    if(needs_locks) {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      for(auto&& a: v) {
#ifdef ASSERT_SPARSEMATRIX
        assert(std::get<0>(a)>=0 && std::get<0>(a)<nr && std::get<1>(a)>=0 && std::get<1>(a)<nc);
#endif
        myrows->push_back(std::get<0>(a));
        colms->push_back(std::get<1>(a));
        vals->push_back(std::get<2>(a));
      }
    } else {
      if(!head) return;
      for(auto&& a: v) {
#ifdef ASSERT_SPARSEMATRIX
        assert(std::get<0>(a)>=0 && std::get<0>(a)<nr && std::get<1>(a)>=0 && std::get<1>(a)<nc);
#endif
        myrows->push_back(std::get<0>(a));
        colms->push_back(std::get<1>(a));
        vals->push_back(std::get<2>(a));
      }
    }
    assert(vals->size()<static_cast<std::size_t>(INT_MAX)); // right now limited to INT_MAX due to indexing problem.
  }

  inline bool remove_repeated_and_compress(MPI_Comm local_comm=MPI_COMM_SELF) 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
#endif

    if(myrows->size() <= 1) return true;
    compress(local_comm);

    if(head) {
      int_iterator first_r=myrows->begin(), last_r=myrows->end(); 
      int_iterator first_c=colms->begin(), last_c=colms->end(); 
      iterator first_v=vals->begin(), last_v = vals->end(); 
      int_iterator result_r = first_r;
      int_iterator result_c = first_c;
      iterator result_v = first_v;

      while ( ( (first_r+1) != last_r)  )
      {
        ++first_r; 
        ++first_c; 
        ++first_v; 
        if (!( (*result_r == *first_r) && (*result_c == *first_c) ) ) { 
          *(++result_r)=*first_r;
          *(++result_c)=*first_c;
          *(++result_v)=*first_v;
        } else {
          if( std::abs(*result_v - *first_v) > 1e-8 ) { //*result_v != *first_v) {
            app_error()<<" Error in call to SMSparseMatrix::remove_repeate_and_compressd. Same indexes with different values. \n";
            app_error()<<"ri, ci, vi: "
                       <<*result_r <<" "
                       <<*result_c <<" "
                       <<*result_v <<" "
                       <<"rj, cj, vj: "
                       <<*first_r <<" "
                       <<*first_c <<" "
                       <<*first_v <<std::endl;
            return false;
          }
        }
      }
      ++result_r;
      ++result_c;
      ++result_v;

      std::size_t sz1 = std::distance(myrows->begin(),result_r); 
      std::size_t sz2 = std::distance(colms->begin(),result_c); 
      std::size_t sz3 = std::distance(vals->begin(),result_v); 
      if(sz1 != sz2 || sz1 != sz2) {
        std::cerr<<"Error: Different number of erased elements in SMSparseMatrix::remove_repeate_and_compressed. \n" <<std::endl;  
        return false;
      }
  
      myrows->resize(sz1);
      colms->resize(sz1);
      vals->resize(sz1);

      // define rowIndex
      int curr=-1;
      for(intType n=0; n<myrows->size(); n++) {
        if( (*myrows)[n] != curr ) {
          int old = curr;
          curr = (*myrows)[n];
          for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
        }
      }
      for(int i=myrows->back()+1; i<rowIndex->size(); i++)
        (*rowIndex)[i] = static_cast<intType>(vals->size());
    }
    MPI_Barrier(local_comm);

    return true;
  }

  // implement parallel compression by using serial sorts on segments of the vector
  // and then either 1. using a parallel in_place_merge of the sorted segments
  // or using serial iterative (every 2 adjacent segments are joined) in_place_merge where every 
  // iteration half of the tasks do nothing (compared to the last iteration)  
  //
  // Possible algorithm:
  // 00. sort rows in parallel, generate rowIndex in serial, round robin for columns
  // 0. only processor numbers in powers of 2 work, so the largest power of 2 in the TG do work, 
  //    everybody else (from 2^lmax to nproc) do nothing 
  // 1. divide the array in Nseg = 2^lmax equalish segments. each working core performs serial in_place_sort on each segment. You might want to implement std::sort with boost::iterator_facade, which should be faster than your pure quick_sort
  // 2. now you have Nseg consecutive sorted subsequences. Now iteratively inplace_merge neighboring subsequences in parallel. Iterate from i=0:lmax and at each iteration 2^(i+1) processors fo work on a given inplace_merge. After each iteration, we adjust the boundaries of the new places to be merged and the number of processors involved.     
  //   - the parallel inplace_merge proceeds by performing iteratively block swap calls to break up the 2 sorted subsets into 2^(i+1) sorted consecutive subsets which can be inplace_merged serially by each working core. 
  //   While it is not ideal, you can use MPI_Barrier() to synchronize the work until you find something better. 
  //
  inline void compress(MPI_Comm local_comm=MPI_COMM_SELF)
  {

    auto comp = [](std::tuple<intType, intType, value_type> const& a, std::tuple<intType, intType, value_type> const& b){return std::get<0>(a) < std::get<0>(b) || (!(std::get<0>(b) < std::get<0>(a)) && std::get<1>(a) < std::get<1>(b));};

    if(local_comm==MPI_COMM_SELF) {
      std::sort(make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin(),colms->begin(),vals->begin()),
                make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->end(),colms->end(),vals->end()),
                comp);
      // define rowIndex
      rowIndex->resize(nr+1);
      int curr=-1; 
      for(intType n=0; n<myrows->size(); n++) {
        if( (*myrows)[n] != curr ) {
          int old = curr;
          curr = (*myrows)[n];
          for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
        }
      }
      for(int i=myrows->back()+1; i<rowIndex->size(); i++)
        (*rowIndex)[i] = static_cast<intType>(vals->size());
      compressed=true;   
    }

    double t1,t2,t3,t4;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    t1 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;

    int npr,rank;
    MPI_Comm_rank(local_comm,&rank); 
    MPI_Comm_size(local_comm,&npr); 

    MPI_Barrier(local_comm);

    gettimeofday(&tv, NULL);
    t2 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time waiting: " <<t2-t1 <<std::endl;

    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
    if(vals->size() == 0) return;

    int nlvl = static_cast<int>(  std::floor(std::log2(npr*1.0)) );
    int nblk = pow(2,nlvl);

    // sort equal segments in parallel
    std::vector<int> pos(nblk+1);     
    // a core processes elements from pos[rank]-pos[rank+1]
    FairDivide(vals->size(),nblk,pos); 
    
    // sort local segment
    if(rank < nblk)   
      std::sort(make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank],colms->begin()+pos[rank],vals->begin()+pos[rank]),
          make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank+1],colms->begin()+pos[rank+1],vals->begin()+pos[rank+1]),
          comp);
    MPI_Barrier(local_comm);

    gettimeofday(&tv, NULL);
    t1 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time sorting: " <<t1-t2 <<std::endl;

      gettimeofday(&tv, NULL);
      t4 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;

/*
      int mid=1;
      int end=2;
      for(int i=0; i<(npr-1); i++) {
        std::inplace_merge( make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[0],colms->begin()+pos[0],vals->begin()+pos[0]),
          make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[mid],colms->begin()+pos[mid],vals->begin()+pos[mid]),
          make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[end],colms->begin()+pos[end],vals->begin()+pos[end]),
          comp);
        mid++;
        end++;

        gettimeofday(&tv, NULL);
        t3 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
        app_log()<<" Time merging: " <<i <<" " <<t3-t4 <<std::endl;
        t4=t3;

      }
*/

    // 

    for(int i=0, k=rank, sz=1; i<nlvl; i++, sz*=2 ) {
      if(k%2==0 && rank<nblk) {
        //std::inplace_merge( vals->begin()+pos[rank], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz*2], comp);
        std::inplace_merge( 
           make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank],colms->begin()+pos[rank],vals->begin()+pos[rank]),
           make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank+sz],colms->begin()+pos[rank+sz],vals->begin()+pos[rank+sz]),
           make_tuple_iterator<int_iterator,int_iterator,iterator>(myrows->begin()+pos[rank+sz*2],colms->begin()+pos[rank+sz*2],vals->begin()+pos[rank+sz*2]),
           comp);
        k/=2;
      } else
        k=1;
      MPI_Barrier(local_comm);
    }

    if(head) {
      // define rowIndex
      rowIndex->resize(nr+1);
      int curr=-1; 
      for(intType n=0; n<myrows->size(); n++) {
        if( (*myrows)[n] != curr ) {
          int old = curr;
          curr = (*myrows)[n];
          for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
        }
      }
      for(int i=myrows->back()+1; i<rowIndex->size(); i++)
        (*rowIndex)[i] = static_cast<intType>(vals->size());
    
    } else {
      if(local_comm==MPI_COMM_SELF)
        APP_ABORT("Error in SMSparseMatrix::compress: Calling with MPI_COMM_SELF from node that is not head of node. \n\n\n");
    }
    MPI_Barrier(local_comm);
    compressed=true;   

    gettimeofday(&tv, NULL);
    t2 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time merging + indexing: " <<t2-t1 <<std::endl;
  }

/*
  inline void compress_old() 
  {
    if(!(myrows->size() == colms->size() && myrows->size() == vals->size()))
      std::cout<<"problem: " <<myrows->size() <<" " <<colms->size() <<" " <<vals->size() <<std::endl; 
#ifdef ASSERT_SPARSEMATRIX
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
#endif
    if(!head) { compressed=true; return; }
    if(vals->size() == 0) return;

    // order along myrows
    int n=myrows->size(); 
    sort_rows(0,n-1);
    if(!std::is_sorted(myrows->begin(),myrows->end())) 
      std::cout<<"ERROR: list is not sorted. \n" <<std::endl;

    // define rowIndex
    rowIndex->resize(nr+1);
    int curr=-1;
    for(int n=0; n<myrows->size(); n++) {
      if( (*myrows)[n] != curr ) {
        int old = curr;
        curr = (*myrows)[n];
        for(int i=old+1; i<=curr; i++) (*rowIndex)[i] = n;
      }
    }
    for(int i=myrows->back()+1; i<rowIndex->size(); i++)
      (*rowIndex)[i] = vals->size();
   
    // order within each rowIndex block
    for(int k=0; k<nr; k++) {
      if((*rowIndex)[k] == (*rowIndex)[k+1]) continue;       
      sort_colms((*rowIndex)[k],(*rowIndex)[k+1]-1);
    } 
   
    compressed=true;
  }

  void sort_rows(int left, int right) {
    int i = left, j = right;
    auto pivot = (*myrows)[(left + right) / 2];

    // partition 
    while (i <= j) {
      while ((*myrows)[i] < pivot)
        i++;
      while ((*myrows)[j] > pivot)
        j--;
      if (i <= j) {
        std::swap((*myrows)[i],(*myrows)[j]);
        std::swap((*colms)[i],(*colms)[j]);
        std::swap((*vals)[i++],(*vals)[j--]);
      }
    };

    // recursion 
    if (left < j)
      sort_rows(left, j);
    if (i < right)
      sort_rows(i, right);
  }

  void sort_colms(int left, int right) {
    int i = left, j = right;
    auto pivot = (*colms)[(left + right) / 2];

    // partition 
    while (i <= j) {
      while ((*colms)[i] < pivot)
        i++;
      while ((*colms)[j] > pivot)
        j--;
      if (i <= j) {
        std::swap((*colms)[i],(*colms)[j]);
        std::swap((*vals)[i++],(*vals)[j--]);
      }
    };

    // recursion 
    if (left < j)
      sort_colms(left, j);
    if (i < right)
      sort_colms(i, right);
  }
*/
  
  inline void transpose(MPI_Comm local_comm=MPI_COMM_SELF)   
  {
    if(!SMallocated) return;
    assert(myrows->size() == colms->size() && myrows->size() == vals->size());
    if(head) {
      // can parallelize this if you want
      for(int_iterator itR=myrows->begin(),itC=colms->begin(); itR!=myrows->end(); ++itR,++itC)
        std::swap(*itR,*itC);
    } else {
      if(local_comm == MPI_COMM_SELF)
        APP_ABORT(" Error in SMSparseMatrix::transpose: Calling with MPI_COMM_SELF from node that is not head of node. \n\n\n");
    }
    std::swap(nr,nc);
    compress(local_comm);
  }

  inline void check()
  {
    if(!head) return; 
    for(int i=0; i<rowIndex->size()-1; i++)
    {
      if((*rowIndex)[i+1] < (*rowIndex)[i]) std::cout<<"Error: SMSparseMatrix::check(): rowIndex-> \n" <<std::endl; 
  
    }

  }

  inline SMSparseMatrix<T>& operator*=(const RealType rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline SMSparseMatrix<T>& operator*=(const std::complex<RealType> rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline void toZeroBase() {
    if(!head) return; 
    if(zero_based) return;
    zero_based=true;
    for (intType& i : *colms ) i--; 
    for (intType& i : *myrows ) i--; 
    for (intType& i : *rowIndex ) i--; 
  }

  inline void toOneBase() {
    if(!head) return; 
    if(!zero_based) return;
    zero_based=false;
    for (intType& i : *colms ) i++; 
    for (intType& i : *myrows ) i++; 
    for (intType& i : *rowIndex ) i++; 
  }

  friend std::ostream& operator<<(std::ostream& out, const SMSparseMatrix<T>& rhs)
  {
    for(intType i=0; i<rhs.vals->size(); i++)
      out<<"(" <<(*(rhs.myrows))[i] <<"," <<(*(rhs.colms))[i] <<":" <<(*(rhs.vals))[i] <<")\n"; 
    return out;
  }

  friend std::istream& operator>>(std::istream& in, SMSparseMatrix<T>& rhs)
  {
    if(!rhs.head) return in;
    T v;
    int c,r;
    in>>r >>c >>v;
    (rhs.vals)->push_back(v); 
    (rhs.myrows)->push_back(r);
    (rhs.colms)->push_back(c);
    return in;
  }

  // this is ugly, but I need to code quickly 
  // so I'm doing this to avoid adding hdf5 support here 
  inline SMVector<T>* getVals() const { return vals; } 
  inline SMVector<intType>* getRows() const { return myrows; }
  inline SMVector<intType>* getCols() const { return colms; }
  inline SMVector<intType>* getRowIndex() const { return rowIndex; }

  inline iterator vals_begin() { assert(vals!=NULL); return vals->begin(); } 
  inline int_iterator rows_begin() { assert(myrows!=NULL); return myrows->begin(); }
  inline int_iterator cols_begin() { assert(colms!=NULL); return colms->begin(); }
  inline indx_iterator rowIndex_begin() { assert(rowIndex!=NULL); return rowIndex->begin(); }
  inline const_iterator vals_begin() const { return vals->begin(); } 
  inline const_int_iterator cols_begin() const { assert(colms!=NULL); return colms->begin(); }
  inline const_indx_iterator rowIndex_begin() const { assert(rowIndex!=NULL); return rowIndex->begin(); }
  inline const_iterator vals_end() const { assert(vals!=NULL); return vals->end(); } 
  inline const_int_iterator rows_end() const { assert(myrows!=NULL); return myrows->end(); }
  inline const_int_iterator cols_end() const { assert(colms!=NULL); return colms->end(); }
  inline const_indx_iterator rowIndex_end() const { assert(rowIndex!=NULL); return rowIndex->end(); }
  inline iterator vals_end() { assert(vals!=NULL); return vals->end(); } 
  inline int_iterator rows_end() { assert(myrows!=NULL); return myrows->end(); }
  inline int_iterator cols_end() { assert(colms!=NULL); return colms->end(); }
  inline indx_iterator rowIndex_end() { assert(rowIndex!=NULL); return rowIndex->end(); }

  inline bool isAllocated() {
    return (SMallocated)&&(vals!=NULL);
  }

  void setRowsFromRowIndex()
  {
    if(!head) return;
    intType shift = zero_based?0:1;
    myrows->resize(vals->size());
    for(int i=0; i<nr; i++)
     for(intType j=(*rowIndex)[i]; j<(*rowIndex)[i+1]; j++)
      (*myrows)[j]=i+shift;
  }

  bool zero_base() const { return zero_based; }
  int row_max() const { return max_in_row; }
  int format() const { return storage_format; } 

  private:

  boost::interprocess::interprocess_mutex *mutex;
  bool compressed;
  int nr,nc;
  SMVector<T> *vals;
  SMVector<intType> *colms,*myrows;
  SMVector<intType> *rowIndex;
  SMVector<unsigned char> *share_buff;
  bool head;
  std::string ID; 
  bool SMallocated;
  bool zero_based;
  int storage_format; // 0: CSR, 1: Compressed Matrix (ESSL) 
  int max_in_row; 
  std::size_t memory=0;

  //_mySort_snD_ my_sort;

  boost::interprocess::managed_shared_memory *segment;
  ShmemAllocator<T> *alloc_T;
  ShmemAllocator<boost::interprocess::interprocess_mutex> *alloc_mutex;
  ShmemAllocator<intType> *alloc_int;
  ShmemAllocator<intType> *alloc_ulong;
  ShmemAllocator<unsigned char> *alloc_uchar;

  // using MPI for barrier calls until I find solution
  MPI_Comm comm;
 
  struct shm_remove
  {
    bool head;
    std::string ID; 
    shm_remove() {
      if(head) boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
    ~shm_remove(){
      if(head) boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  } remover;

};


}

#endif
#endif
