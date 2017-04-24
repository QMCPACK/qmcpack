
#ifndef QMCPLUSPLUS_AFQMC_SMDENSEVECTOR_H
#define QMCPLUSPLUS_AFQMC_SMDENSEVECTOR_H

#include<iostream>
#include<vector>
#include<tuple>
#include <cassert>
#include<algorithm>
#include<complex>
#include<cmath>
#include"AFQMC/config.0.h"
#include"AFQMC/Utilities/Utils.h"

#include <boost/interprocess/managed_shared_memory.hpp>
#include <boost/interprocess/containers/vector.hpp>
#include <boost/interprocess/allocators/allocator.hpp>
#include <boost/interprocess/containers/vector.hpp>

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

  SMDenseVector<T>():head(false),ID(""),SMallocated(false),vals(NULL),share_buff(NULL),mutex(NULL),
                      segment(NULL),alloc_T(NULL),alloc_mutex(NULL),alloc_uchar(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  ~SMDenseVector<T>()
  {
    if(segment!=NULL) {
     delete segment;
     boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  }

  // this should probably be disallowed
  SMDenseVector(const SMDenseVector<T> &rhs)
  {
//    ID = rhs.ID; // is this a good idea???
//    head = rhs.head;
    APP_ABORT(" Error: SMDenseVector(SMDenseVector rhs) copy constructor has been disabled.");       
  }

  inline void setup(bool hd, std::string ii, MPI_Comm comm_) {
    head=hd;
    ID=ii;
    remover.ID=ii;
    remover.head=hd;
    comm=comm_;
  }

  inline void reserve(int nnz, bool allow_reduce = false)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) || (vals!=NULL && vals->capacity() > nnz && allow_reduce)) 
      allocate(nnz,allow_reduce);
    if(head) vals->reserve(nnz);
    barrier();
  }

  template<typename T1>
  void share(T1* x, int n, bool sender) {
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
  void share(std::vector<T1>& x, int n, bool sender) {
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

  inline void barrier() {
    MPI_Barrier(comm);   
  }

  inline bool deallocate()
  {
    SMallocated = false;
    barrier();
    if(!head) {
      try{
        delete segment;
        segment=NULL;
      } catch(std::bad_alloc&) {
        std::cerr<<"Problems deleting segment in SMDenseVector::deallocate()." <<std::endl;
        return false;
      }
    }
    barrier();
    if(head) {
      try{
        delete segment;
        segment=NULL;
        boost::interprocess::shared_memory_object::remove(ID.c_str());
      } catch(std::bad_alloc&) {
        std::cerr<<"Problems de-allocating shared memory in SMDenseVector." <<std::endl;
        return false;
      }
    }
    barrier();
  } 

  // this routine does not allow grow/shrink, meant in cases where only head can call it
  inline bool allocate_serial(int n)
  {
    if(!head) return false; /* XA: This was returning nothing, I assume false is the right thing to return here */
    if(vals!=NULL && vals->capacity() >= n) return true;

    memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+1000*sizeof(unsigned char)+8000;
 
    try {
      segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
    } catch(boost::interprocess::interprocess_exception &ex) {
      std::cout<<" Found managed_shared_memory segment, removing. Careful with persistent SHM segment. \n";
      boost::interprocess::shared_memory_object::remove(ID.c_str());
      segment=NULL;
    }

    if(segment==NULL) {
      try {
        segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
      } catch(boost::interprocess::interprocess_exception &ex) {
        std::cerr<<"Problems setting up managed_shared_memory in SMDenseVector." <<std::endl;
        return false;
      }
    }

    try {
      alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());
      alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());

      share_buff = segment->construct<boost_SMVector<unsigned char>>("share_buff")(*alloc_uchar);
      share_buff->resize(1000);

      mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();

      vals = segment->construct<boost_SMVector<T>>("vals")(*alloc_T);
      vals->reserve(n);

    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMDenseVector." <<std::endl;
      return false;
    }
    SMallocated=true;
    return true;

  }

  inline bool allocate(int n, bool allow_reduce=false)
  {
    bool grow = false;
    uint64_t old_sz = (segment==NULL)?0:(segment->get_size()); 
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
      memory = sizeof(boost::interprocess::interprocess_mutex)+n*sizeof(T)+1000*sizeof(unsigned char)+8000; 

      if(grow) {
        if(memory > old_sz) { 
          uint64_t extra = memory - old_sz;
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::grow(ID.c_str(), extra)) {
            std::cerr<<" Error growing shared memory in  SMDenseVector::allocate(). \n";
            return false;   
          } 
        } else {
          segment->destroy<boost_SMVector<T>>("vals");
          vals = segment->construct<boost_SMVector<T>>("vals")(*alloc_T);
          vals->reserve(n);          
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::shrink_to_fit(ID.c_str())) {
            std::cerr<<" Error in shrink_to_fit shared memory in  SMDenseVector::allocate(). \n";
            return false;   
          } 
        }

        try {
          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::open_only, ID.c_str());
          vals = segment->find<boost_SMVector<T>>("vals").first;
          share_buff = segment->find<boost_SMVector<unsigned char>>("share_buff").first;
          mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
          assert(vals != 0);
          assert(share_buff != 0);
          assert(mutex != 0);
          vals->reserve(n);
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems opening shared memory in SMDenseVector::allocate() ." <<std::endl;
          return false;
        }

      } else {

        try {

          segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
        } catch(boost::interprocess::interprocess_exception &ex) {
          std::cout<<" Found managed_shared_memory segment, removing. Careful with persistent SHM segment. \n";
          boost::interprocess::shared_memory_object::remove(ID.c_str());
          segment=NULL;
        }

        if(segment==NULL) {
          try {
            segment = new boost::interprocess::managed_shared_memory(boost::interprocess::create_only, ID.c_str(), memory);
          } catch(boost::interprocess::interprocess_exception &ex) {
            std::cerr<<"Problems setting up managed_shared_memory in SMDenseVector." <<std::endl;
            return false;
          }
        }

        try {
          alloc_T = new ShmemAllocator<T>(segment->get_segment_manager());
          alloc_uchar = new ShmemAllocator<unsigned char>(segment->get_segment_manager());

          share_buff = segment->construct<boost_SMVector<unsigned char>>("share_buff")(*alloc_uchar);
          share_buff->resize(1000);

          mutex = segment->construct<boost::interprocess::interprocess_mutex>("mutex")();

          vals = segment->construct<boost_SMVector<T>>("vals")(*alloc_T);
          vals->reserve(n);
 
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems allocating shared memory in SMDenseVector." <<std::endl;
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
      vals = segment->find<boost_SMVector<T>>("vals").first;
      share_buff = segment->find<boost_SMVector<unsigned char>>("share_buff").first;
      mutex = segment->find<boost::interprocess::interprocess_mutex>("mutex").first;
      assert(vals != 0);
      assert(share_buff != 0);
      assert(mutex != 0);
    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMDenseVector: initializeChildren() ." <<std::endl;
      return false;
    }
    return true;
  }

  // resize is probably the best way to setup the vector 
  inline void resize(int nnz, bool allow_reduce=false) 
  {
    if(vals==NULL) {
      allocate(nnz,allow_reduce);
    } else if(vals->capacity() < nnz) {
      std::vector<T> tmp;
      int sz = vals->size();
      if(head) {
        tmp.resize(sz);
        std::copy(vals->begin(),vals->begin()+sz,tmp.begin());
      }
      allocate(nnz,allow_reduce);
      if(head) {
        vals->resize(nnz);
        std::copy(tmp.begin(),tmp.begin()+sz,vals->begin());
      }
    } else if(vals->capacity() > nnz && allow_reduce) {
      std::vector<T> tmp;
      if(head) {
        tmp.resize(nnz);
        std::copy(vals->begin(),vals->begin()+nnz,tmp.begin());
      } 
      allocate(nnz,allow_reduce);
      if(head) {
        vals->resize(nnz);
        std::copy(tmp.begin(),tmp.begin()+nnz,vals->begin());
      }
    }
    if(head) vals->resize(nnz);
    barrier();
  }

  // does not allow grow/shrink
  inline void resize_serial(int nnz)
  {
    if(!head) return;
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) ) 
      APP_ABORT("Error: Calling SMDenseVector::resize_serial(n) without enough capacity. \n");
    vals->resize(nnz);
  }

  inline void clear() { 
    if(!head) return;
    if(!SMallocated) return;
    vals->clear();
  }

  inline unsigned long size() const
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

  inline bool isAllocated() {
    return (SMallocated)&&(vals!=NULL); 
  }

  inline This_t& operator=(const SMDenseVector<T> &rhs) 
  { 
    APP_ABORT(" Error: SMDenseVector(SMDenseVector rhs) operator= has been disabled.");       
    //resize(rhs.size());
    //if(!head) return *this;
    //(*vals)=*(rhs.vals);
    //return *this;
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

  inline uint64_t memoryUsage() { return memory; }

  inline int capacity() { return (vals==NULL)?0:vals->capacity(); }

  inline void push_back(const T& v, bool needs_locks=false)             
  {
    assert(vals != NULL);
    assert(vals->capacity() >= vals->size()+1 );
    if(needs_locks) {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      vals->push_back(v);
    } else {
      vals->push_back(v);
    }
  }

  inline void sort() {
    if(!head) return;
    if(vals==NULL) return;
    std::sort(vals->begin(),vals->end());
  }

  template<class Compare>
  inline void sort(Compare comp) {
    if(!head) return;
    if(vals==NULL) return;
    std::sort(vals->begin(),vals->end(),comp);
  }

  template<class Compare>
  inline void sort(Compare comp, MPI_Comm local_comm, bool inplace=true) {

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

    if(vals==NULL) return;
    if(vals->size() == 0) return;

    int nlvl = static_cast<int>(  std::floor(std::log2(npr*1.0)) );
    int nblk = pow(2,nlvl);

    // sort equal segments in parallel
    std::vector<int> pos(nblk+1);
    // a core processes elements from pos[rank]-pos[rank+1]
    FairDivide(vals->size(),nblk,pos);

    // sort local segment
    if(rank < nblk)
      std::sort(vals->begin()+pos[rank], vals->begin()+pos[rank+1], comp);

    MPI_Barrier(local_comm);

    gettimeofday(&tv, NULL);
    t1 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time sorting: " <<t1-t2 <<std::endl;

/*
    if(head) {

      gettimeofday(&tv, NULL);
      t4 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;

      for(int i=0, mid=1, end=2; i<(nblk-1); i++,mid++,end++) {
        std::inplace_merge( vals->begin()+pos[0], vals->begin()+pos[mid], vals->begin()+pos[end], comp); 

        gettimeofday(&tv, NULL);
        t3 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
        app_log()<<" Time merging: " <<i <<" " <<t3-t4 <<std::endl;
        t4=t3;
      }
    }
    MPI_Barrier(local_comm);
*/

    gettimeofday(&tv, NULL);
    t4 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;


    std::vector<T> temp;
    if (!inplace && rank<nblk)
      for(int i=0, sz=nblk; i<nlvl; i++,sz/=2) 
        if( rank%sz==0 ) {
          int nt=0;
          int nb = pow(2,i);
          for(int j=0; j<nb; j++)
            nt = std::max(nt, pos[(j+1)*sz]-pos[j*sz]);
          temp.resize(nt);  
          break;
        }
    
    for(int i=0, k=rank, sz=1; i<nlvl; i++, sz*=2 ) {
      if(k%2==0 && rank<nblk) {
        if(inplace) {
          std::inplace_merge( vals->begin()+pos[rank], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz*2], comp);  
        } else { 
          int nt = pos[rank+sz*2] - pos[rank];
          temp.resize(nt); // just in case
          std::merge( vals->begin()+pos[rank], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz*2], temp.begin(), comp);  
          std::copy(temp.begin(),temp.begin()+nt,vals->begin()+pos[rank]);
        }
        k/=2;
      } else
        k=1;  
      MPI_Barrier(local_comm);
    }



// FIX FIX FIX
// Why is this not working??????
/*
    for(int i=0; i<nlvl; i++) {
      int np = std::pow(2,i+1);  // number of processors on each group 
      int del = std::pow(2,i);   // number of "skips" on bounds array  
      int beg = (rank/np)*np;
      int mid = beg + del;
      int end = mid + del;
      int rk = rank%np;
      T* ptr_b = &((*vals)[0])+pos[beg];
      T* ptr_m = &((*vals)[0])+pos[mid];
      T* ptr_e = &((*vals)[0])+pos[end];
      if (rank < nblk)   
        parallel_inplace_merge(np,rk,ptr_b,ptr_m,ptr_e,local_comm,comp);
        //parallel_inplace_merge(np,rk,vals->begin()+pos[beg],vals->begin()+pos[mid],vals->begin()+pos[end],local_comm,comp);
      else {
        for(int k=0; k<3*(i+1); k++)
         MPI_Barrier(local_comm); 
      }
    }
*/

    gettimeofday(&tv, NULL);
    t2 =  double(tv.tv_sec)+double(tv.tv_usec)/1000000.0;
    app_log()<<" Time merging + indexing: " <<t2-t1 <<std::endl;

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

  inline iterator begin() { assert(vals!=NULL); return vals->begin(); } 
  inline const_iterator begin() const { assert(vals!=NULL); return vals->begin(); } 
  inline const_iterator end() const { assert(vals!=NULL); return vals->end(); } 
  inline iterator end() { assert(vals!=NULL); return vals->end(); } 
  inline T& back() { assert(vals!=NULL); return vals->back(); } 

  boost::interprocess::interprocess_mutex* getMutex()
  {
    return mutex;
  } 

  void setDims(int nr, int nc) {
    dims = {nr,nc};  
  }
  
  inline std::pair<int,int> getDims() { return dims; }
  inline int rows() { return dims.first; }
  inline int cols() { return dims.second; }

  private:

  boost::interprocess::interprocess_mutex *mutex;
  boost_SMVector<T> *vals;
  boost_SMVector<unsigned char> *share_buff;
  bool head;
  std::string ID; 
  bool SMallocated;
  uint64_t memory=0;
  std::pair<int, int> dims;   // kind of cheating, but allows me to use as a matrix when needed 

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
