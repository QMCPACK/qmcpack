#ifndef AFQMC_CONFIG_H
#define AFQMC_CONFIG_H

#include <string>
#include <algorithm>
#include <cstdlib>
#include <ctype.h>
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

//#include "mpi3/shared_window.hpp"
#include "AFQMC/Memory/SharedMemory/shm_ptr_with_raw_ptr_dispatch.hpp"
#include "multi/array.hpp"
#include "multi/array_ref.hpp"
#include "multi/memory/fallback.hpp"

#include "Utilities/TimerManager.h"

namespace qmcplusplus
{
extern TimerList_t AFQMCTimers;
enum AFQMCTimerIDs
{
  block_timer,
  pseudo_energy_timer,
  energy_timer,
  vHS_timer,
  assemble_X_timer,
  vbias_timer,
  G_for_vbias_timer,
  propagate_timer,
  back_propagate_timer,
  E_comm_overhead_timer,
  vHS_comm_overhead_timer,
  popcont_timer,
  ortho_timer,
  setup_timer,
  extra_timer,
  T1_t,
  T2_t,
  T3_t,
  T4_t,
  T5_t,
  T6_t,
  T7_t,
  T8_t
};
extern const TimerNameList_t<AFQMCTimerIDs> AFQMCTimerNames;

namespace afqmc
{
// ultil we switch to c++17, to reduce extra lines
using tp_ul_ul = std::tuple<std::size_t, std::size_t>;

enum WALKER_TYPES
{
  UNDEFINED_WALKER_TYPE,
  CLOSED,
  COLLINEAR,
  NONCOLLINEAR
};
// when ENABLE_CUDA is not set, DEVICE and TG_LOCAL are the same
enum ALLOCATOR_TYPES
{
  STD,
  NODE,
  STD_DEVICE,
  SHARED_LOCAL_DEVICE,
  SHARED_DEVICE
};

inline WALKER_TYPES initWALKER_TYPES(int i)
{
  if (i == 0)
    return UNDEFINED_WALKER_TYPE;
  else if (i == 1)
    return CLOSED;
  else if (i == 2)
    return COLLINEAR;
  else if (i == 3)
    return NONCOLLINEAR;
  return UNDEFINED_WALKER_TYPE;
}

template<typename T>
using s1D = std::tuple<IndexType, T>;
template<typename T>
using s2D = std::tuple<IndexType, IndexType, T>;
template<typename T>
using s3D = std::tuple<IndexType, IndexType, IndexType, T>;
template<typename T>
using s4D = std::tuple<IndexType, IndexType, IndexType, IndexType, T>;

enum SpinTypes
{
  Alpha,
  Beta
};

// allocators
template<class T>
using shared_allocator = shm::allocator_shm_ptr_with_raw_ptr_dispatch<T>;
template<class T>
using shm_pointer = typename shared_allocator<T>::pointer;

#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
template<class T>
using device_allocator = device::device_allocator<T>;
template<class T>
using device_ptr = device::device_pointer<T>;
template<class T>
using localTG_allocator = device_allocator<T>;
template<class T>
using node_allocator = device_allocator<T>;
template<class T, class TG>
localTG_allocator<T> make_localTG_allocator(TG&)
{
  return localTG_allocator<T>{};
}
template<class T, class TG>
node_allocator<T> make_node_allocator(TG&)
{
  return node_allocator<T>{};
}
/*   Temporary fix for the conflict problem between cpu and gpu pointers. Find proper fix */
template<class T>
device_ptr<T> make_device_ptr(device_ptr<T> p)
{
  return p;
}
template<class T>
device_ptr<T> make_device_ptr(T* p)
{
  print_stacktrace;
  throw std::runtime_error(" Invalid pointer conversion: device_pointer<T> to T*.");
}
//template<class T>
//device_ptr<T> make_device_ptr(boost::mpi3::intranode::array_ptr<T> p)
//{
//  print_stacktrace;
//  throw std::runtime_error(" Invalid pointer conversion: device_pointer<T> to T*.");
//}
template<class T>
device_ptr<T> make_device_ptr(shm::shm_ptr_with_raw_ptr_dispatch<T> p)
{
  print_stacktrace;
  throw std::runtime_error(" Invalid pointer conversion: device_pointer<T> to T*.");
}

using device_memory_resource = device::memory_resource;
using shm_memory_resource    = device::memory_resource;
template<class T>
using device_constructor = device::constructor<T>;
template<class T>
using shm_constructor = device::constructor<T>;

#else
template<class T>
using device_allocator = std::allocator<T>;
template<class T>
using device_ptr = T*;
template<class T>
using localTG_allocator = shared_allocator<T>;
template<class T>
using node_allocator = shared_allocator<T>;
template<class T, class TG>
localTG_allocator<T> make_localTG_allocator(TG& t_)
{
  return localTG_allocator<T>{t_.TG_local()};
}
template<class T, class TG>
node_allocator<T> make_node_allocator(TG& t_)
{
  return node_allocator<T>{t_.Node()};
}
/*   Temporary fix for the conflict problem between cpu and gpu pointers. Find proper fix */
template<class T>
device_ptr<T> make_device_ptr(T* p)
{
  return p;
}
//template<class T>
//device_ptr<T> make_device_ptr(boost::mpi3::intranode::array_ptr<T> p) = delete;
//{ //return device_ptr<T>{to_address(p)}; }*/
//  print_stacktrace;*/
//  throw std::runtime_error(" Invalid pointer conversion: device_pointer<T> to T*.");*/
//}
template<class T>
device_ptr<T> make_device_ptr(shm::shm_ptr_with_raw_ptr_dispatch<T> p)
{
  return device_ptr<T>{to_address(p)};
}

using device_memory_resource = boost::multi::memory::resource<>;
using shm_memory_resource    = shm::memory_resource_shm_ptr_with_raw_ptr_dispatch;
template<class T>
using device_constructor = device_allocator<T>;
template<class T>
using shm_constructor = shared_allocator<T>;

#endif

template<class T>
using host_constructor     = std::allocator<T>;
using host_memory_resource = boost::multi::memory::resource<>;

// new types
using SpCType_shm_csr_matrix =
    ma::sparse::csr_matrix<SPComplexType, int, std::size_t, shared_allocator<SPComplexType>, ma::sparse::is_root>;
using SpVType_shm_csr_matrix =
    ma::sparse::csr_matrix<SPValueType, int, std::size_t, shared_allocator<SPValueType>, ma::sparse::is_root>;
using CType_shm_csr_matrix =
    ma::sparse::csr_matrix<ComplexType, int, std::size_t, shared_allocator<ComplexType>, ma::sparse::is_root>;
using VType_shm_csr_matrix =
    ma::sparse::csr_matrix<ValueType, int, std::size_t, shared_allocator<ValueType>, ma::sparse::is_root>;

//#ifdef PsiT_IN_SHM
template<typename T>
using PsiT_Matrix_t = ma::sparse::csr_matrix<T, int, int, shared_allocator<T>, ma::sparse::is_root>;
using PsiT_Matrix   = PsiT_Matrix_t<ComplexType>;
#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
using devcsr_Matrix = ma::sparse::csr_matrix<ComplexType, int, int, device_allocator<ComplexType>>;
#else
using devcsr_Matrix = ma::sparse::csr_matrix<ComplexType, int, int, shared_allocator<ComplexType>, ma::sparse::is_root>;
#endif
//#else
//  using PsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int>;
//  using devPsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int>;
//#endif


#if defined(ENABLE_CUDA) || defined(ENABLE_HIP)
using P1Type = ma::sparse::csr_matrix<ComplexType, int, int, localTG_allocator<ComplexType>>;
#else
using P1Type        = ma::sparse::csr_matrix<ComplexType, int, int, localTG_allocator<ComplexType>, ma::sparse::is_root>;
#endif

enum HamiltonianTypes
{
  Factorized,
  THC,
  KPTHC,
  KPFactorized,
  RealDenseFactorized,
  UNKNOWN
};

template<std::ptrdiff_t D>
using iextensions = typename boost::multi::iextensions<D>;
//using extensions = typename boost::multi::layout_t<D>::extensions_type;

// general matrix definitions
template<class Alloc = std::allocator<int>>
using IntegerVector = boost::multi::array<int, 1, Alloc>;
template<class Alloc = std::allocator<ValueType>>
using ValueVector = boost::multi::array<ValueType, 1, Alloc>;
template<class Alloc = std::allocator<ComplexType>>
using ComplexVector = boost::multi::array<ComplexType, 1, Alloc>;
template<class Alloc = std::allocator<SPComplexType>>
using SPComplexVector = boost::multi::array<SPComplexType, 1, Alloc>;
template<class Ptr = ComplexType*>
using ComplexVector_ref = boost::multi::array_ref<ComplexType, 1, Ptr>;
template<class Ptr = SPComplexType*>
using SPComplexVector_ref = boost::multi::array_ref<SPComplexType, 1, Ptr>;

template<class Alloc = std::allocator<int>>
using IntegerMatrix = boost::multi::array<int, 2, Alloc>;
template<class Alloc = std::allocator<ValueType>>
using ValueMatrix = boost::multi::array<ValueType, 2, Alloc>;
template<class Alloc = std::allocator<ComplexType>>
using ComplexMatrix = boost::multi::array<ComplexType, 2, Alloc>;
template<class Alloc = std::allocator<SPComplexType>>
using SPComplexMatrix = boost::multi::array<SPComplexType, 2, Alloc>;
template<class Ptr = ComplexType*>
using ComplexMatrix_ref = boost::multi::array_ref<ComplexType, 2, Ptr>;
template<class Ptr = SPComplexType*>
using SPComplexMatrix_ref = boost::multi::array_ref<SPComplexType, 2, Ptr>;

template<class Alloc = std::allocator<ComplexType>>
using Complex3Tensor = boost::multi::array<ComplexType, 3, Alloc>;
template<class Alloc = std::allocator<SPComplexType>>
using SPComplex3Tensor = boost::multi::array<SPComplexType, 3, Alloc>;
template<class Ptr = ComplexType*>
using Complex3Tensor_ref = boost::multi::array_ref<ComplexType, 3, Ptr>;
template<class Ptr = SPComplexType*>
using SPComplex3Tensor_ref = boost::multi::array_ref<SPComplexType, 3, Ptr>;

template<std::ptrdiff_t D, class Alloc = std::allocator<ComplexType>>
using ComplexArray = boost::multi::array<ComplexType, D, Alloc>;
template<std::ptrdiff_t D, class Alloc = std::allocator<SPComplexType>>
using SPComplexArray = boost::multi::array<SPComplexType, D, Alloc>;
template<std::ptrdiff_t D, class Ptr = ComplexType*>
using ComplexArray_ref = boost::multi::array_ref<ComplexType, D, Ptr>;
template<std::ptrdiff_t D, class Ptr = SPComplexType*>
using SPComplexArray_ref = boost::multi::array_ref<SPComplexType, D, Ptr>;


struct AFQMCInfo
{
public:
  // default constructor
  AFQMCInfo()
      : name(""),
        NMO(-1),
        NMO_FULL(-1),
        NAEA(-1),
        NAEB(-1),
        NCA(0),
        NCB(0),
        NETOT(-1),
        MS2(-99),
        ISYM(-1),
        spinRestricted(true)
  {}

  AFQMCInfo(std::string nm, int nmo_, int naea_, int naeb_)
      : name(nm),
        NMO(nmo_),
        NMO_FULL(nmo_),
        NAEA(naea_),
        NAEB(naeb_),
        NCA(0),
        NCB(0),
        NETOT(-1),
        MS2(-99),
        ISYM(-1),
        spinRestricted(true)
  {}

  AFQMCInfo(const AFQMCInfo& other) = default;
  AFQMCInfo& operator=(const AFQMCInfo& other) = default;

  // destructor
  ~AFQMCInfo() {}

  // identifier
  std::string name;

  // number of active orbitals
  int NMO;

  // number of orbitals
  int NMO_FULL;

  // number of active electrons alpha/beta
  int NAEA, NAEB;

  // number of core electrons alpha/beta
  int NCA, NCB;

  // total number of electrons
  int NETOT;

  // ms2
  int MS2;

  // isym
  int ISYM;

  // if true then RHF calculation, otherwise it is UHF
  bool spinRestricted;

  // copies values from object
  void copyInfo(const AFQMCInfo& a)
  {
    name           = a.name;
    NMO_FULL       = a.NMO_FULL;
    NMO            = a.NMO;
    NAEA           = a.NAEA;
    NAEB           = a.NAEB;
    NCA            = a.NCA;
    NCB            = a.NCB;
    NETOT          = a.NETOT;
    MS2            = a.MS2;
    ISYM           = a.ISYM;
    spinRestricted = a.spinRestricted;
  }

  // no fully spin polarized yet, not sure what it will break
  bool checkAFQMCInfoState()
  {
    if (NMO_FULL < 1 || NAEA < 1 || NAEB < 1 || NCA < 0 || NCB < 0) //|| NETOT!= NCA+NCB+NAEA+NAEB ) //|| MS2<0 )
      return false;
    return true;
  }

  void printAFQMCInfoState(std::ostream& out)
  {
    out << "AFQMC info: \n"
        << "name: " << name << "\n"
        << "NMO_FULL: " << NMO_FULL << "\n"
        << "NAEA: " << NAEA << "\n"
        << "NAEB: " << NAEB << "\n"
        << "NCA: " << NCA << "\n"
        << "NCB: " << NCB << "\n"
        << "NETOT: " << NETOT << "\n"
        << "MS2: " << MS2 << "\n"
        << "spinRestricted: " << spinRestricted << std::endl;
  }

  bool parse(xmlNodePtr cur)
  {
    if (cur == NULL)
      return false;

    OhmmsAttributeSet oAttrib;
    oAttrib.add(name, "name");
    oAttrib.put(cur);

    std::string sR("yes");
    ParameterSet m_param;
    m_param.add(NMO_FULL, "NMO_FULL");
    m_param.add(NMO_FULL, "NMO");
    m_param.add(NAEA, "NAEA");
    m_param.add(NAEB, "NAEB");
    m_param.add(NCA, "NCA");
    m_param.add(NCB, "NCB");
    m_param.add(NETOT, "NETOT");
    m_param.add(MS2, "MS2");
    m_param.add(sR, "spinRestricted");
    m_param.put(cur);

    spinRestricted = false;
    std::string sR0(sR);
    std::transform(sR0.begin(), sR0.end(), sR.begin(), (int (*)(int))tolower);
    if (sR == "yes" || sR == "true")
      spinRestricted = true;

    NMO = NMO_FULL - NCA;
    if (NETOT == -1)
      NETOT = NCA + NCB + NAEA + NAEB;

    return true;
  }
};

} // namespace afqmc
} // namespace qmcplusplus

#endif
