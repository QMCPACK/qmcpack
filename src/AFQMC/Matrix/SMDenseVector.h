
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


  typedef T              Type_t;
  typedef T              value_type;
  typedef T*             pointer;
  typedef const T*       const_pointer;
  typedef std::size_t    indxType;
  typedef typename boost_SMVector<T>::iterator iterator;
  typedef typename boost_SMVector<T>::const_iterator const_iterator;
  typedef SMDenseVector<T>  This_t;

  SMDenseVector<T>():head(false),ID(""),SMallocated(false),vals(NULL),share_buff(NULL),mutex(NULL),
                      segment(NULL),alloc_T(NULL),alloc_mutex(NULL),alloc_uchar(NULL) 
  {
    remover.ID="NULL";
    remover.head=false;
  }

  SMDenseVector<T>(bool hd, std::string ii, MPI_Comm comm_):head(false),ID(""),SMallocated(false),vals(NULL),share_buff(NULL),mutex(NULL),
                      segment(NULL),alloc_T(NULL),alloc_mutex(NULL),alloc_uchar(NULL)
  {
    setup(hd,ii,comm_);
  }

  ~SMDenseVector<T>()
  {
    if(segment!=NULL) {
     delete segment;
     boost::interprocess::shared_memory_object::remove(ID.c_str());
    }
  }

  // this should probably be disallowed
  SMDenseVector(const SMDenseVector<T> &rhs) = delete;
  This_t& operator=(const SMDenseVector<T> &rhs) = delete; 

  SMDenseVector(SMDenseVector<T>&& other):head(false),ID(""),SMallocated(false),vals(NULL),
                    share_buff(NULL),mutex(NULL),segment(NULL),alloc_T(NULL),alloc_mutex(NULL),
                    alloc_uchar(NULL) {
    *this = std::move(other);
  }

  This_t& operator=(SMDenseVector<T>&& other)
  {
    if(this != &other) {

      if(segment!=NULL) {
        delete segment;
        boost::interprocess::shared_memory_object::remove(ID.c_str());
      }

      remover.ID=other.remover.ID;
      remover.head=other.remover.head;

      mutex = other.mutex;
      vals = other.vals;
      share_buff = other.share_buff;
      segment = other.segment;
      alloc_T = other.alloc_T;
      alloc_mutex = other.alloc_mutex;
      alloc_uchar = other.alloc_uchar;
      comm = other.comm;
      head = other.head;
      ID = other.ID;
      SMallocated = other.SMallocated;
      memory = other.memory;
      dims = other.dims;

      other.setNull();
    }
    return *this;
  }

   void setup(bool hd, std::string ii, MPI_Comm comm_) {
    head=hd;
    ID=ii;
    remover.ID=ii;
    remover.head=hd;
    comm=comm_;
  }

   void reserve(std::size_t nnz, bool allow_reduce = false)
  {
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) || (vals!=NULL && vals->capacity() > nnz && allow_reduce)) 
      allocate(nnz,allow_reduce);
    if(head) {
      assert(nnz < vals->max_size());
      vals->reserve(nnz);
    } 
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

   void barrier() {
    MPI_Barrier(comm);   
  }

   bool deallocate()
  {
    SMallocated = false;
    barrier();
    if(!head) {
      try{
        delete segment;
        segment=NULL;
      } catch(std::bad_alloc&) {
        std::cerr<<"Problems deleting segment in SMDenseVector::deallocate()." <<std::endl;
        APP_ABORT("Problems deleting segment in SMDenseVector::deallocate().\n");
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
        APP_ABORT("Problems de-allocating shared memory in SMDenseVector.\n");
        return false;
      }
    }
    barrier();
  } 

  // this routine does not allow grow/shrink, meant in cases where only head can call it
   bool allocate_serial(std::size_t n)
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
        APP_ABORT("Problems setting up managed_shared_memory in SMDenseVector.\n");
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
      assert(n < vals->max_size());
      vals->reserve(n);

    } catch(std::bad_alloc&) {
      std::cerr<<"Problems allocating shared memory in SMDenseVector." <<std::endl;
      APP_ABORT("Problems allocating shared memory in SMDenseVector.\n");
      return false;
    }
    SMallocated=true;
    return true;

  }

   bool allocate(std::size_t n, bool allow_reduce=false)
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
            APP_ABORT(" Error growing shared memory in  SMDenseVector::allocate(). \n");
            return false;   
          } 
        } else {
          segment->destroy<boost_SMVector<T>>("vals");
          vals = segment->construct<boost_SMVector<T>>("vals")(*alloc_T);
          assert(n < vals->max_size());
          vals->reserve(n);          
          delete segment;
          segment=NULL;
          if(!boost::interprocess::managed_shared_memory::shrink_to_fit(ID.c_str())) {
            std::cerr<<" Error in shrink_to_fit shared memory in  SMDenseVector::allocate(). \n";
            APP_ABORT(" Error in shrink_to_fit shared memory in  SMDenseVector::allocate(). \n");
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
          assert(n < vals->max_size());
          vals->reserve(n);
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems opening shared memory in SMDenseVector::allocate() ." <<std::endl;
          APP_ABORT("Problems opening shared memory in SMDenseVector::allocate() .\n");
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
            APP_ABORT("Problems setting up managed_shared_memory in SMDenseVector.\n");
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
          assert(n < vals->max_size());
          vals->reserve(n);
 
        } catch(std::bad_alloc&) {
          std::cerr<<"Problems allocating shared memory in SMDenseVector." <<std::endl;
          APP_ABORT("Problems allocating shared memory in SMDenseVector.\n");
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
   bool initializeChildren()
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
      APP_ABORT("Problems allocating shared memory in SMDenseVector: initializeChildren() .\n");
      return false;
    }
    return true;
  }

  // resize is probably the best way to setup the vector 
   void resize(std::size_t nnz, bool allow_reduce=false) 
  {
    if(vals==NULL) {
      allocate(nnz,allow_reduce);
    } else if(vals->capacity() < nnz) {
      std::vector<T> tmp;
      indxType sz = vals->size();
      if(head) {
        tmp.resize(sz);
        std::copy(vals->begin(),vals->begin()+sz,tmp.begin());
      }
      allocate(nnz,allow_reduce);
      if(head) {
        assert(nnz < vals->max_size());
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
        assert(nnz < vals->max_size());
        vals->resize(nnz);
        std::copy(tmp.begin(),tmp.begin()+nnz,vals->begin());
      }
    }
    if(head) {
      assert(nnz < vals->max_size());
      vals->resize(nnz);
    }
    barrier();
  }

  // does not allow grow/shrink
   void resize_serial(std::size_t nnz)
  {
    if(!head) return;
    if(vals==NULL || (vals!=NULL && vals->capacity() < nnz) ) 
      APP_ABORT("Error: Calling SMDenseVector::resize_serial(n) without enough capacity. \n");
    assert(nnz < vals->max_size());
    vals->resize(nnz);
  }

   void clear() { 
    if(!head) return;
    if(!SMallocated) return;
    vals->clear();
  }

   std::size_t size() const
  {
    return (vals!=NULL)?(vals->size()):0;
  }

   const_pointer values() const 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

   pointer values() 
  {
    return (vals!=NULL)?(&((*vals)[0])):NULL;
  }

   bool isAllocated() {
    return (SMallocated)&&(vals!=NULL); 
  }

   Type_t& operator()(std::size_t i)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<vals->size());
#endif
    return (*vals)[i]; 
  }

   Type_t& operator[](std::size_t i)
  {
#ifdef ASSERT_SPARSEMATRIX
    assert(i>=0 && i<vals->size());
#endif
    return (*vals)[i]; 
  }

  template<typename IType>
   void add(const IType i, const T& v, bool needs_locks=false) 
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

   std::size_t memoryUsage() { return memory; }

   std::size_t capacity() { return (vals==NULL)?0:vals->capacity(); }

   void push_back(const T& v, bool needs_locks=false)             
  {
    assert(vals != NULL);
    if(needs_locks) {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      assert(vals->capacity() >= vals->size()+1 );
      vals->push_back(v);
    } else {
      assert(vals->capacity() >= vals->size()+1 );
      vals->push_back(v);
    }
  }

   void push_back(const std::vector<T>& v, bool needs_locks=false)
  {
    assert(vals != NULL);
    if(needs_locks) {
      boost::interprocess::scoped_lock<boost::interprocess::interprocess_mutex> lock(*mutex);
      assert(vals->capacity() >= vals->size()+v.size() );
      for(auto&& i: v) 
        vals->push_back(i);
    } else {
      assert(vals->capacity() >= vals->size()+v.size() );
      for(auto&& i: v) 
        vals->push_back(i);
    }
  }

   void sort() {
    if(!head) return;
    if(vals==NULL) return;
    std::sort(vals->begin(),vals->end());
  }

  template<class Compare>
   void sort(Compare comp, MPI_Comm local_comm=MPI_COMM_SELF, bool inplace=true) {

    if(vals==NULL) return;
    if(local_comm==MPI_COMM_SELF) {
      std::sort(vals->begin(),vals->end(),comp);
    }

    int npr,rank;
    MPI_Comm_rank(local_comm,&rank);
    MPI_Comm_size(local_comm,&npr);

    MPI_Barrier(local_comm);

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

    std::vector<T> temp;
    if (!inplace && rank<nblk) {
      int nt=0;
      for(int i=0, k=rank, sz=1; i<nlvl; i++, sz*=2 ) {
        if(k%2==0 && rank<nblk) {
          nt = std::max(nt,pos[rank+sz*2] - pos[rank]);
          k/=2;
        } else
          k=1;
      }
      temp.resize(nt);
    }
    
    for(int i=0, k=rank, sz=1; i<nlvl; i++, sz*=2 ) {
      if(k%2==0 && rank<nblk) {
        if(inplace) {
          std::inplace_merge( vals->begin()+pos[rank], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz*2], comp);  
        } else { 
          std::size_t nt = pos[rank+sz*2] - pos[rank];
           assert( temp.size() >= nt );  
          std::merge( vals->begin()+pos[rank], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz], vals->begin()+pos[rank+sz*2], temp.begin(), comp);  
          std::copy(temp.begin(),temp.begin()+nt,vals->begin()+pos[rank]);
        }
        k/=2;
      } else
        k=1;  
      MPI_Barrier(local_comm);
    }

  }

   SMDenseVector<T>& operator*=(const RealType rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

   SMDenseVector<T>& operator*=(const std::complex<RealType> rhs ) 
  {
    if(!head) return *this; 
    for(iterator it=vals->begin(); it!=vals->end(); it++)
      (*it) *= rhs;
    return *this; 
  }

  friend std::ostream& operator<<(std::ostream& out, const SMDenseVector<T>& rhs)
  {
    for(std::size_t i=0; i<rhs.vals->size(); i++)
      out<<"(" <<(*(rhs.myrows))[i] <<"," <<(*(rhs.colms))[i] <<":" <<(*(rhs.vals))[i] <<")\n"; 
    return out;
  }

  // this is ugly, but I need to code quickly 
  // so I'm doing this to avoid adding hdf5 support here 
   boost_SMVector<T>* getVector() const { return vals; } 

   iterator begin() { assert(vals!=NULL); return vals->begin(); } 
   const_iterator begin() const { assert(vals!=NULL); return vals->begin(); } 
   const_iterator end() const { assert(vals!=NULL); return vals->end(); } 
   iterator end() { assert(vals!=NULL); return vals->end(); } 
   T& back() { assert(vals!=NULL); return vals->back(); } 

  boost::interprocess::interprocess_mutex* getMutex()
  {
    return mutex;
  } 

  void setDims(int nr, int nc) {
    dims = {nr,nc};  
  }
  
   std::pair<int,int> getDims() { return dims; }
   int rows() { return dims.first; }
   int cols() { return dims.second; }

  private:

  void setNull() {
    mutex=NULL;
    vals=NULL;
    share_buff=NULL;
    head=false;
    SMallocated=false;
    memory=0;
    segment=NULL;
    alloc_T=NULL;
    alloc_mutex=NULL;
    alloc_uchar=NULL;
    remover.head=false;
    remover.ID="";
  }

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
