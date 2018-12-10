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

#include "AFQMC/Utilities/myTimer.h"
extern myTimer Timer; 

#define AFQMC_DEBUG 3 
#define AFQMC_TIMER 

// careful here that RealType is consistent with this!!!
#define MKL_INT         int
#define MKL_Complex8    std::complex<float> 
#define MKL_Complex16   std::complex<double> 

namespace qmcplusplus
{
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

struct AFQMCInfo 
{
  public:

  // default constructor
  AFQMCInfo():name(""),NMO(-1),NMO_FULL(-1),NAEA(-1),NAEB(-1),NCA(0),NCB(0),NETOT(-1),MS2(-99)
            ,spinRestricted(-1),ISYM(-1)
  {}

  // destructor
  ~AFQMCInfo() {}

  // name
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
    name=a.name;
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
