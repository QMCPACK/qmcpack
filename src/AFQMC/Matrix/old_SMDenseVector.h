
#ifndef QMCPLUSPLUS_AFQMC_SMDENSEVECTOR_H
#define QMCPLUSPLUS_AFQMC_SMDENSEVECTOR_H

#include<iostream>
#include<vector>
#include<tuple>
#include<assert.h>
#include<algorithm>
#include<complex>
#include"../config.0.h"
//#include"AFQMC/config.0.h"

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/vector.hpp>
//#include <boost/interprocess/sync/interprocess_barrier.hpp>

#define ASSERT_VECTOR

namespace qmcplusplus
{

// wrapper for boost::interprocess::vector 
template<class T>
class SMDenseVector
{
  public:

  template<typename spT> using ShmemAllocator = boost::interprocess::allocator<spT, boost::interprocess::managed_shared_memory::segment_manager>;
  template<typename spT> using boost_SMVector = boost::interprocess::vector<spT, ShmemAllocator<spT>>;


  typedef T            Type_t;
  typedef T            value_type;
  typedef T*           pointer;
  typedef const T*     const_pointer;
  typedef const int*   const_indxPtr;
  typedef int*           indxPtr;
  typedef typename boost_SMVector<T>::iterator iterator;
  typedef typename boost_SMVector<T>::const_iterator const_iterator;
  typedef typename boost_SMVector<int>::iterator int_iterator;
  typedef typename boost_SMVector<int>::const_iterator const_int_iterator;
  typedef boost_SMVector<T>  This_t;

  SMDenseVector<T>():head(false),ID(""),SMallocated(false),vals(NULL),share_buff(NULL),mutex(NULL),npig(0) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  ~SMDenseVector<T>()
  {
    if(head && SMallocated) {
     segment->destroy<boost_SMVector<T>>("vals");
     segment->destroy<boost_SMVector<unsigned char>>("share_buff");
     segment->destroy<boost_SMVector<T>>("");
     delete segment;
     boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  }

  // this should probably be disallowed
  SMDenseVector(const SMDenseVector<T> &rhs)
  {
    ID = rhs.ID; // is this a good idea???
    head = rhs.head;
    return;
  }

  inline void setup(bool hd, std::string ii, int _npig, MPI_Comm comm_) {
    npig = _npig;
    head=hd;
    ID=ii;
    remover.ID=ii;
    remover.head=hd;
    comm=comm_;
  }

  inline void reserve(int n)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < n)) {
      allocate(n);
      barrier();
      initializeChildren();
      barrier();
    }
    barrier();
  }

  template<typename T1>
  inline void share(T1& x, bool sender) {
    if(sender) {
      std::memcpy(&((*share_buff)[0]),&x,sizeof(T1));
      barrier(); 
    } else {
      barrier(); 
      std::memcpy(&x,&((*share_buff)[0]),sizeof(T1));
    }
    barrier(); 
  }

  inline void barrier() {
    if(npig==1) return;
//    bool done = mybarrier->wait();
    MPI_Barrier(comm);   
  }

  inline bool deallocate()
  {
    SMallocated = false;
    if(head) {
      try{

cout<<"in deallocate() " <<std::endl;
//        segment->destroy<boost_SMVector<T>>("vals");
//        segment->destroy<boost_SMVector<unsigned char>>("share_buff");
//        segment->destroy<boost_SMVector<boost::interprocess::interprocess_mutex>>("mutex");
//        delete segment;
        segment=NULL;
cout<<"in deallocate(): calling remove " <<std::endl;
        boost::interprocess::shared_memory_object::remove(ID.c_str());
cout<<"in deallocate(): done calling remove " <<std::endl;

      } catch(std::bad_alloc&) {
        std::cerr<<"Problems de-allocating shared memory in SMDenseVector." <<std::endl;
        return false;
      }
    }
cout<<"done in deallocate() " <<std::endl;
    //barrier();
  } 

  // done like this to avoid dependencies on MPI
  // Barrier called outside this routine
  inline bool allocate(int n)
  {
    if(!head) { 
      SMallocated=true; 
      return true; 
    } 
    if(SMallocated) {
      if(vals!=NULL && vals->capacity() >= n) return true; 
      deallocate(); 
    }
    //uint64_t memory = sizeof(boost::interprocess::barrier)+sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+1000*sizeof(unsigned char)+8000; 
    uint64_t memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+1000*sizeof(unsigned char)+8000; 

    try {

      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
    } catch(boost::interprocess::interprocess_exception &ex) {
      std::cout<<" Found managed_shared_memory segment, removing. \n";
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str());
      boost::interprocess::shared_memory_object::remove(ID.c_str());
      segment=NULL;
    }

    if(segment==NULL) {
      try {
        segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
      } catch(boost::interprocess::interprocess_exception &ex) {
        std::cerr<<"Problems setting up managed_shared_memory in SMSparseMatrix." <<std::endl;
        return false;
      }
    }


    try {
      alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());
      alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());

      share_buff = segment->construct<boost_SMVector<unsigned char>>("share_buff")(*alloc_uchar);
      share_buff->resize(1000);

      mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();

      //mybarrier = segment->construct<boost::interprocess::barrier>("barrier")(npig);

      vals = segment->construct<boost_SMVector<T>>("vals")(*alloc_T);
      vals->reserve(n);

    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMDenseVector." <<std::endl;
      return false;
    }
    SMallocated=true;
    return true;
  }

  // only call this when all arrays have been allocated and modified  
  inline bool initializeChildren()
  { 
    if(head) return true;
    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str()); 
      vals = segment->find<boost_SMVector<T>>("vals").first;
      share_buff = segment->find<boost_SMVector<unsigned char>>("share_buff").first;
      mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
      //mybarrier = segment->find<boost::interprocess::barrier>("barrier").first;
    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMDenseVector: initializeChildren() ." <<std::endl;
      return false;
    }
    return true;
  }

  // resize is probably the best way to setup the vector 
  inline void resize(int nnz, bool wait=true)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz)) {
      allocate(nnz);
      if(wait) {
        barrier();
        initializeChildren();
        barrier();
      }
    }
    if(head) vals->resize(nnz);
    if(wait) barrier();
  }

  inline void clear() { 
    if(!head) return;
    if(!SMallocated) return;
    vals->clear();
  }

  inline int size() const
  {
    return (vals!=NULL)?(vals->size()):0;
  }

  inline const_pointer values() const 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

  inline pointer values() 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

  inline This_t& operator=(const SMDenseVector<T> &rhs) 
  { 
    resize(rhs.size());
    if(!head) return;
    (*vals)=*(rhs.vals);
  }  

  inline Type_t& operator()(unsigned int i)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<vals->size());
#endif
    return (*vals)[i]; 
  }

  inline Type_t& operator[](unsigned int i)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<vals->size());
#endif
    return (*vals)[i]; 
  }

  inline void add(const int i, const T& v, bool needs_locks=false) 
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<vals->size());
#endif
    if(needs_locks) {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
        (*vals)[i]=v;
    } else {
      if(!head) return;
      (*vals)[i]=v;
    }
  }

  inline int capacity() { return (vals==NULL)?0:vals->capacity(); }

  inline void push_back(const T& v, bool needs_locks=false)             
  {
    if(vals==NULL) return;
    if(vals->capacity() <= vals->size()+1) allocate(vals->size()+1000);
    if(needs_locks) {
        boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
        vals->push_back(v);
    } else {
      if(!head) return;
      vals->push_back(v);
    }
  }

  inline SMDenseVector<T>& operator*=(const RealType rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  inline SMDenseVector<T>& operator*=(const std::complex<RealType> rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  friend std::ostream& operator<<(std::ostream& out, const SMDenseVector<T>& rhs)
  {
    for(int i=0; i<rhs.vals->size(); i++)
      out<<"(" <<(*(rhs.myrows))[i] <<"," <<(*(rhs.colms))[i] <<":" <<(*(rhs.vals))[i] <<")\n"; 
    return out;
  }

  // this is ugly, but I need to code quickly 
  // so I'm doing this to avoid adding hdf5 support here 
  inline boost_SMVector<T>* getVector() const { return vals; } 

  inline iterator begin() { return vals->begin(); } 
  inline const_iterator begin() const { return vals->begin(); } 
  inline const_iterator end() const { return vals->end(); } 
  inline iterator end() { return vals->end(); } 

  boost::interprocess::interprocess_mutex* getMutex()
  {
    return mutex;
  } 

  private:

  boost::interprocess::interprocess_mutex *mutex;
//  boost::interprocess::barrier *mybarrier;
  boost_SMVector<T> *vals;
  boost_SMVector<unsigned char> *share_buff;
  bool head;
  std::string ID; 
  bool SMallocated;
  int npig;

  boost::interprocess::managed_shared_memory *segment;
  ShmemAllocator<T> *alloc_T;
  ShmemAllocator<boost::interprocess::interprocess_mutex> *alloc_mutex;
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
