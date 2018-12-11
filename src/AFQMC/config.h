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
#include <OhmmsPETE/OhmmsMatrix.h>
#include <OhmmsPETE/Tensor.h>
#include <OhmmsPETE/OhmmsVector.h>
#include <OhmmsPETE/TinyVector.h>

#include "AFQMC/config.0.h"
#include "AFQMC/Matrix/SparseMatrix.h"
#include "AFQMC/Matrix/SMSparseMatrix.h"
#include "AFQMC/Matrix/SMDenseVector.h"

#include "AFQMC/Matrix/csr_matrix.hpp"
#include "AFQMC/Matrix/coo_matrix.hpp"
#include<boost/multi_array.hpp>

#include "mpi3/shared_window.hpp"

#include "Utilities/NewTimer.h"
#include "AFQMC/Utilities/myTimer.h"
extern myTimer Timer; 

#define AFQMC_DEBUG 3 
#define AFQMC_TIMER 

#define MAXIMUM_EMPLACE_BUFFER_SIZE 102400 

// maximum size in Bytes for a dataset with walker data on WalkerIO
#define WALKER_HDF_BLOCK_SIZE 100000000 

// maximum size in Bytes for a block of data in CSR matrix HDF IO 
#define CSR_HDF_BLOCK_SIZE 2000000 

// careful here that RealType is consistent with this!!!
#define MKL_INT         int
#define MKL_Complex8    std::complex<float> 
#define MKL_Complex16   std::complex<double> 

#define byRows   999
#define byCols   111

#define PsiT_IN_SHM

namespace qmcplusplus
{

  enum WALKER_TYPES {UNDEFINED_WALKER_TYPE, CLOSED, COLLINEAR, NONCOLLINEAR};

  //using boost::multi_array_types::index_gen;
  using boost::extents;
  using boost::indices;
  using range_t = boost::multi_array_types::index_range;

  enum SpinTypes {Alpha,Beta};  

  // new types
  using SpCType_shm_csr_matrix = ma::sparse::csr_matrix<SPComplexType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPComplexType>,
                                boost::mpi3::intranode::is_root>;
  using SpVType_shm_csr_matrix = ma::sparse::csr_matrix<SPValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<SPValueType>,
                                boost::mpi3::intranode::is_root>;
  using CType_shm_csr_matrix = ma::sparse::csr_matrix<ComplexType,int,std::size_t,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                boost::mpi3::intranode::is_root>;
  using VType_shm_csr_matrix = ma::sparse::csr_matrix<ValueType,int,std::size_t,
                                boost::mpi3::intranode::allocator<ValueType>,
                                boost::mpi3::intranode::is_root>;

#ifdef PsiT_IN_SHM
  using PsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                boost::mpi3::intranode::is_root>;
#else
  using PsiT_Matrix = ma::sparse::csr_matrix<ComplexType,int,int>;
#endif


  using P1Type = ma::sparse::csr_matrix<ComplexType,int,int,
                                boost::mpi3::intranode::allocator<ComplexType>,
                                boost::mpi3::intranode::is_root>;

  // old types

  typedef Vector<IndexType>     IndexVector;
  typedef Vector<RealType>      RealVector;
  typedef Vector<ValueType>     ValueVector;
  typedef Vector<SPValueType>   SPValueVector;
  typedef Vector<ComplexType>   ComplexVector;
  typedef Vector<SPComplexType>   SPComplexVector;

  typedef SMDenseVector<IndexType>     IndexSMVector;
  typedef SMDenseVector<RealType>      RealSMVector;
  typedef SMDenseVector<ValueType>     ValueSMVector;
  typedef SMDenseVector<SPValueType>   SPValueSMVector;
  typedef SMDenseVector<ComplexType>   ComplexSMVector;
  typedef SMDenseVector<SPComplexType>   SPComplexSMVector;

  typedef Matrix<IndexType>     IndexMatrix;
  typedef Matrix<RealType>      RealMatrix;
  typedef Matrix<ValueType>     ValueMatrix;
  typedef Matrix<SPValueType>     SPValueMatrix;
  typedef Matrix<ComplexType>   ComplexMatrix;
  typedef Matrix<SPComplexType>   SPComplexMatrix;

  typedef SparseMatrix<IndexType>     IndexSpMat;
  typedef SparseMatrix<RealType>      RealSpMat;
  typedef SparseMatrix<ValueType>     ValueSpMat;
  typedef SparseMatrix<SPValueType>   SPValueSpMat;
  typedef SparseMatrix<ComplexType>   ComplexSpMat;

  typedef SMSparseMatrix<IndexType>     IndexSMSpMat;
  typedef SMSparseMatrix<RealType>      RealSMSpMat;
  typedef SMSparseMatrix<ValueType>     ValueSMSpMat;
  typedef SMSparseMatrix<SPValueType>   SPValueSMSpMat;
  typedef SMSparseMatrix<ComplexType>   ComplexSMSpMat;
  typedef SMSparseMatrix<SPComplexType>   SPComplexSMSpMat;

  enum HamiltonianTypes {Factorized,SymmetricFactorized,s4DInts,THC};

  extern TimerList_t AFQMCTimers;
  enum AFQMCTimerIDs {    
    block_timer,
    pseudo_energy_timer,
    energy_timer,
    vHS_timer,
    vbias_timer,
    G_for_vbias_timer,
    propagate_timer,
    E_comm_overhead_timer,
    vHS_comm_overhead_timer
  };
  extern TimerNameList_t<AFQMCTimerIDs> AFQMCTimerNames;  

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

#endif
