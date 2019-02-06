#ifndef AFQMC_CONFIG_H 
#define AFQMC_CONFIG_H 

#include <string>
#include <algorithm>
#include<cstdlib>
#include<ctype.h>
#include <vector>
#include <map>
#include <complex>
#include <tuple>
#include <fstream>
#include "Configuration.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"

#include "AFQMC/config.0.h"

#include "AFQMC/Memory/custom_pointers.hpp"

#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/coo_matrix.hpp"

#include "mpi3/shared_window.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"

#include "Utilities/NewTimer.h"
#include "AFQMC/Utilities/myTimer.h"
extern myTimer Timer; 

namespace qmcplusplus
{

  extern TimerList_t AFQMCTimers;
  enum AFQMCTimerIDs {    
    block_timer,
    pseudo_energy_timer,
    energy_timer,
    vHS_timer,
    assemble_X_timer,
    vbias_timer,
    G_for_vbias_timer,
    propagate_timer,
    E_comm_overhead_timer,
    vHS_comm_overhead_timer
  };
  extern TimerNameList_t<AFQMCTimerIDs> AFQMCTimerNames;  

namespace afqmc
{

  // ultil we switch to c++17, to reduce extra lines 
  using tp_ul_ul = std::tuple<std::size_t,std::size_t>;

  enum WALKER_TYPES {UNDEFINED_WALKER_TYPE, CLOSED, COLLINEAR, NONCOLLINEAR};
  // when QMC_CUDA is not set, DEVICE and TG_LOCAL are the same
  enum ALLOCATOR_TYPES {STD,NODE,STD_DEVICE,SHARED_LOCAL_DEVICE,SHARED_DEVICE};

  template<typename T> using s1D = std::tuple<IndexType,T>;
  template<typename T> using s2D = std::tuple<IndexType,IndexType,T>;
  template<typename T> using s3D = std::tuple<IndexType,IndexType,IndexType,T>;
  template<typename T> using s4D = std::tuple<IndexType,IndexType,IndexType,IndexType,T>;

  enum SpinTypes {Alpha,Beta};  

  // allocators
  template<class T>
  using shared_allocator = boost::mpi3::intranode::allocator<T>;
  template<class T>
  using shm_pointer = typename shared_allocator<T>::pointer; 

#if defined(QMC_CUDA)
  template<class T>  using device_allocator = qmc_cuda::cuda_gpu_allocator<T>;
  template<class T>  using device_ptr = qmc_cuda::cuda_gpu_ptr<T>;
  template<class T>  using localTG_allocator = device_allocator<T>;
  template<class T>  using node_allocator = device_allocator<T>;
  template<class T, class TG> 
  localTG_allocator<T> make_localTG_allocator(TG&) {return localTG_allocator<T>{};}
  template<class T, class TG> 
  node_allocator<T> make_node_allocator(TG&) {return node_allocator<T>{};}
/*   Temporary fix for the conflict problem between cpu and gpu pointers. Find proper fix */
  template<class T> device_ptr<T> make_device_ptr(device_ptr<T> p) { return p; }
  template<class T> device_ptr<T> make_device_ptr(T* p) 
  {
    print_stacktrace;
    throw std::runtime_error(" Invalid pointer conversion: cuda_gpu_ptr<T> to T*.");
  }   
  template<class T> device_ptr<T> make_device_ptr(boost::mpi3::intranode::array_ptr<T> p) 
  { 
    print_stacktrace;
    throw std::runtime_error(" Invalid pointer conversion: cuda_gpu_ptr<T> to T*.");
  }   
#else
  template<class T>  using device_allocator = std::allocator<T>;
  template<class T>  using device_ptr = T*;
  template<class T>  using localTG_allocator = shared_allocator<T>; 
  template<class T>  using node_allocator = shared_allocator<T>;
  template<class T, class TG> 
  localTG_allocator<T> make_localTG_allocator(TG& t_) {return localTG_allocator<T>{t_.TG_local()};}
  template<class T, class TG> 
  node_allocator<T> make_node_allocator(TG& t_) {return node_allocator<T>{t_.Node()};}
/*   Temporary fix for the conflict problem between cpu and gpu pointers. Find proper fix */
  template<class T> device_ptr<T> make_device_ptr(T* p) { return p; }
  template<class T> device_ptr<T> make_device_ptr(boost::mpi3::intranode::array_ptr<T> p) 
    { return device_ptr<T>{to_address(p)}; }
#endif

  // new types
  using SpCType_shm_csr_matrix = ma::sparse::csr_matrix<SPComplexType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPComplexType>,
                                ma::sparse::is_root>;
  using SpVType_shm_csr_matrix = ma::sparse::csr_matrix<SPValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPValueType>,
                                ma::sparse::is_root>;
  using CType_shm_csr_matrix = ma::sparse::csr_matrix<ComplexType,int,std::size_t,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                ma::sparse::is_root>;
  using VType_shm_csr_matrix = ma::sparse::csr_matrix<ValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<ValueType>,
                                ma::sparse::is_root>;

//#ifdef PsiT_IN_SHM
  using PsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                ma::sparse::is_root>;
#ifdef QMC_CUDA
  using devPsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int,
                                device_allocator<ComplexType>>; 
#else
  using devPsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                ma::sparse::is_root>;
#endif
//#else
//  using PsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int>;
//  using devPsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int>;
//#endif


#if defined(QMC_CUDA)
  using P1Type = ma::sparse::csr_matrix<ComplexType,int,int,
                                 localTG_allocator<ComplexType>>;
#else
  using P1Type = ma::sparse::csr_matrix<ComplexType,int,int,
                                localTG_allocator<ComplexType>,
                                ma::sparse::is_root>;
#endif

  enum HamiltonianTypes {Factorized,THC,KPTHC,KPFactorized,UNKNOWN};

  template<std::ptrdiff_t D> 
  using iextensions = typename boost::multi::iextensions<D>;
  //using extensions = typename boost::multi::layout_t<D>::extensions_type;  

  // general matrix definitions
  template< class Alloc = std::allocator<int> >
  using IntegerVector =  boost::multi::array<int,1,Alloc>;
  template< class Alloc = std::allocator<ValueType> >
  using ValueVector =  boost::multi::array<ValueType,1,Alloc>;
  template< class Alloc = std::allocator<ComplexType> >
  using ComplexVector =  boost::multi::array<ComplexType,1,Alloc>;
  template< class Alloc = std::allocator<SPComplexType> >
  using SPComplexVector =  boost::multi::array<SPComplexType,1,Alloc>;
  template< class Ptr = ComplexType* >
  using ComplexVector_ref =  boost::multi::array_ref<ComplexType,1,Ptr>;
  template< class Ptr = SPComplexType* >
  using SPComplexVector_ref =  boost::multi::array_ref<SPComplexType,1,Ptr>;

  template< class Alloc = std::allocator<int> >
  using IntegerMatrix =  boost::multi::array<int,2,Alloc>;
  template< class Alloc = std::allocator<ValueType> >
  using ValueMatrix =  boost::multi::array<ValueType,2,Alloc>;
  template< class Alloc = std::allocator<ComplexType> >
  using ComplexMatrix =  boost::multi::array<ComplexType,2,Alloc>;
  template< class Alloc = std::allocator<SPComplexType> >
  using SPComplexMatrix =  boost::multi::array<SPComplexType,2,Alloc>;
  template< class Ptr = ComplexType* >
  using ComplexMatrix_ref =  boost::multi::array_ref<ComplexType,2,Ptr>;
  template< class Ptr = SPComplexType* >
  using SPComplexMatrix_ref =  boost::multi::array_ref<SPComplexType,2,Ptr>;

  template< class Alloc = std::allocator<ComplexType> >
  using Complex3Tensor =  boost::multi::array<ComplexType,3,Alloc>;
  template< class Alloc = std::allocator<SPComplexType> >
  using SPComplex3Tensor =  boost::multi::array<SPComplexType,3,Alloc>;
  template< class Ptr = ComplexType* >
  using Complex3Tensor_ref =  boost::multi::array_ref<ComplexType,3,Ptr>;
  template< class Ptr = SPComplexType* >
  using SPComplex3Tensor_ref =  boost::multi::array_ref<SPComplexType,3,Ptr>;

  template<std::ptrdiff_t D, class Alloc = std::allocator<ComplexType> >
  using ComplexArray =  boost::multi::array<ComplexType,D,Alloc>;
  template<std::ptrdiff_t D, class Alloc = std::allocator<SPComplexType> >
  using SPComplexArray =  boost::multi::array<SPComplexType,D,Alloc>;
  template<std::ptrdiff_t D, class Ptr = ComplexType* >
  using ComplexArray_ref =  boost::multi::array_ref<ComplexType,D,Ptr>;
  template<std::ptrdiff_t D, class Ptr = SPComplexType* >
  using SPComplexArray_ref =  boost::multi::array_ref<SPComplexType,D,Ptr>;


  struct AFQMCInfo 
  {
    public:

    // default constructor
    AFQMCInfo():
        name(""),NMO(-1),NMO_FULL(-1),NAEA(-1),NAEB(-1),NCA(0),NCB(0),NETOT(-1),
        MS2(-99),spinRestricted(true),ISYM(-1)
    {
    }

    AFQMCInfo(std::string nm, int nmo_, int naea_, int naeb_):
        name(nm),NMO(nmo_),NMO_FULL(nmo_),NAEA(naea_),NAEB(naeb_),NCA(0),NCB(0),
        NETOT(-1),MS2(-99),spinRestricted(true),ISYM(-1)
    {
    }

    AFQMCInfo( const AFQMCInfo& other) = default;
    AFQMCInfo& operator=( const AFQMCInfo& other) = default;

    // destructor
    ~AFQMCInfo() {}

    // identifier
    std::string name;

    // number of orbitals
    int NMO_FULL;

    // number of active orbitals
    int NMO;

    // number of active electrons alpha/beta 
    int NAEA, NAEB;

    // number of core electrons alpha/beta
    int NCA,NCB;

    // total number of electrons
    int NETOT; 

    // ms2
    int MS2; 

    // isym
    int ISYM;  

    // if true then RHF calculation, otherwise it is UHF 
    bool spinRestricted;

    // copies values from object
    void copyInfo(const AFQMCInfo& a) {
      name = a.name;
      NMO_FULL=a.NMO_FULL;
      NMO=a.NMO;
      NAEA=a.NAEA;
      NAEB=a.NAEB;
      NCA=a.NCA;
      NCB=a.NCB;
      NETOT=a.NETOT;
      MS2=a.MS2;
      ISYM=a.ISYM;
      spinRestricted=a.spinRestricted;
    }

    // no fully spin polarized yet, not sure what it will break 
    bool checkAFQMCInfoState() {
      if(NMO_FULL<1 || NAEA<1 || NAEB<1 || NCA<0 || NCB<0 ) //|| NETOT!= NCA+NCB+NAEA+NAEB ) //|| MS2<0 )
        return false;
      return true; 
    } 

    void printAFQMCInfoState(std::ostream& out) {
      out<<"AFQMC info: \n"
         <<"name: " <<name <<"\n"
         <<"NMO_FULL: " <<NMO_FULL <<"\n"
         <<"NAEA: " <<NAEA  <<"\n"
         <<"NAEB: " <<NAEB  <<"\n"
         <<"NCA: " <<NCA  <<"\n"
         <<"NCB: " <<NCB  <<"\n"
         <<"NETOT: " <<NETOT  <<"\n"
         <<"MS2: " <<MS2  <<"\n"
         <<"spinRestricted: " <<spinRestricted <<std::endl;
    }

    bool parse(xmlNodePtr cur)
    {

      if(cur == NULL)
        return false;

      xmlNodePtr curRoot=cur;
      OhmmsAttributeSet oAttrib;
      oAttrib.add(name,"name");
      oAttrib.put(cur);

      std::string sR("yes");
      ParameterSet m_param;
      m_param.add(NMO_FULL,"NMO_FULL","int");
      m_param.add(NMO_FULL,"NMO","int");
      m_param.add(NAEA,"NAEA","int");
      m_param.add(NAEB,"NAEB","int");
      m_param.add(NCA,"NCA","int");
      m_param.add(NCB,"NCB","int");
      m_param.add(NETOT,"NETOT","int");
      m_param.add(MS2,"MS2","int");
      m_param.add(sR,"spinRestricted","string");
      m_param.put(cur);

      spinRestricted=false;
      std::string sR0(sR);
      std::transform(sR0.begin(),sR0.end(),sR.begin(),(int (*)(int))tolower);
      if(sR == "yes" || sR == "true") spinRestricted = true;

      NMO = NMO_FULL-NCA;
      if(NETOT==-1) NETOT = NCA+NCB+NAEA+NAEB;

      return true;
    }

  };

}
}

#endif
