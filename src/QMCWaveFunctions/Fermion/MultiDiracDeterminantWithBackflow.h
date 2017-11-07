//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
/**@file DiracDeterminantBaseBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_MULTIDIRACDETERMINANTWITHBASE_H
#define QMCPLUSPLUS_MULTIDIRACDETERMINANTWITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase_help.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "Message/Communicate.h"
#include "Numerics/DeterminantOperators.h"
//#include "Numerics/OhmmsBlas.h"

namespace qmcplusplus
{

class MultiDiracDeterminantBase: public OrbitalBase
{

public:
  bool Optimizable;
  void registerTimers();
  NewTimer UpdateTimer, RatioTimer, InverseTimer, buildTableTimer, readMatTimer, evalWTimer, evalOrbTimer,evalOrb1Timer;
  NewTimer readMatGradTimer,buildTableGradTimer,ExtraStuffTimer;
  // Optimizable parameters
  opt_variables_type myVars;

  typedef SPOSetBase::IndexVector_t IndexVector_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
  typedef SPOSetBase::HessMatrix_t  HessMatrix_t;
  typedef SPOSetBase::HessType      HessType;
  typedef OrbitalSetTraits<ValueType>::HessVector_t  HessVector_t;
  typedef TinyVector<HessType, 3>   GGGType;
  typedef Vector<GGGType>           GGGVector_t;
  typedef Matrix<GGGType>           GGGMatrix_t;
  typedef Array<HessType,3>         HessArray_t;

  /** constructor
   *@param spos the single-particle orbital set
   *@param first index of the first particle
   */
  MultiDiracDeterminantBase(SPOSetBasePtr const &spos, int first=0);

  ///default destructor
  ~MultiDiracDeterminantBase();

  /**copy constructor
   * @param s existing DiracDeterminantBase
   *
   * This constructor makes a shallow copy of Phi.
   * Other data members are allocated properly.
   */
  MultiDiracDeterminantBase(const MultiDiracDeterminantBase& s);

  MultiDiracDeterminantBase& operator=(const MultiDiracDeterminantBase& s);

  /** return a clone of Phi
   */
  SPOSetBasePtr clonePhi() const;

  SPOSetBasePtr getPhi()
  {
    return Phi;
  };

  inline IndexType rows() const
  {
    return NumPtcls;
  }

  inline IndexType cols() const
  {
    return NumOrbitals;
  }

  /** set the index of the first particle in the determinant and reset the size of the determinant
   *@param first index of first particle
   *@param nel number of particles in the determinant
   *@param norb total number of orbitals (including unoccupied)
   */
  void set(int first, int nel, int norb);

  void set(int first, int nel);

  ///set BF pointers
  void setBF(BackflowTransformation* bf)
  {
    usingBF=true;
    BFTrans=bf;
  }

  ///optimizations  are disabled
  inline void checkInVariables(opt_variables_type& active)
  {
    //Phi->checkInVariables(active);
  }

  inline void checkOutVariables(const opt_variables_type& active)
  {
    //Phi->checkOutVariables(active);
  }

  void resetParameters(const opt_variables_type& active)
  {
    //Phi->resetParameters(active);
  }

  inline void reportStatus(std::ostream& os)
  {
  }
  void resetTargetParticleSet(ParticleSet& P)
  {
    Phi->resetTargetParticleSet(P);
  }

  virtual RealType getAlternatePhaseDiff()
  {
    return 0.0;
  }
  virtual RealType getAlternatePhaseDiff(int iat)
  {
    return 0.0;
  }

  ///reset the size: with the number of particles and number of orbtials
  virtual void resize(int nel, int morb);

  void registerData(ParticleSet& P, WFBufferType& buf);

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, WFBufferType& buf);

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat);

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  /****************************************************************************
   * These functions should not be called.
   ***************************************************************************/

  ValueType ratio(ParticleSet& P, int iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }

  GradType evalGrad(ParticleSet& P, int iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }

  RealType evaluate(ParticleSet& P,
                    ParticleSet::ParticleGradient_t& G,
                    ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }


  /****************************************************************************
   * END END END
   ***************************************************************************/

  // create necessary structures used in the evaluation of the determinants
  // this works with confgList, which shouldn't change during a simulation
  void createDetData(ci_configuration2& ref, std::vector<int>& data,
                     std::vector<std::pair<int,int> >& pairs, std::vector<double>& sign);
  /*
      void testDets() {

         std::vector<int> indx(16);
         ValueMatrix_t Mat(8,8);
         NewTimer dummy("dummy");
         dummy.reset();

         std::ifstream in("matrix.txt");
         for(int i=0; i<8; i++)
         for(int j=0; j<8; j++)
            in>>Mat(i,j);
         in.close();

         std::vector<double> dets(9);
         dets[2] =  -0.085167; // 2x2
         dets[3] = 0.051971; // 3x3
         dets[4] = 0.029865; // 4x4
         dets[5] = -0.055246; // 5x5
         dets[6] = -0.0087414; // 6x6
         dets[7] = 0.025257; // 7x7
         dets[8] = 0.067119; // 8x8

         for(int q=2; q<=8; q++) {
           int k=0;
           for(int i=0; i<q; i++)
            indx[k++] = i;
           for(int i=0; i<q; i++)
            indx[k++] = i;
           std::cout <<"determinant of " <<q <<": " <<dets[q] <<"  -  " <<CalculateRatioFromMatrixElements(q,Mat,indx.begin()) << std::endl;
         }
         int k=0;
         for(int i=0; i<4; i++)
          indx[k++] = i;
         for(int i=0; i<4; i++)
          indx[k++] = i;

         std::cout <<"Timing n=4 \n";
         ValueType res;
         dummy.start();
         for(int i=0; i<1000; i++)
           res=CalculateRatioFromMatrixElements(4,Mat,indx.begin());
         dummy.stop();
         double direct = dummy.get_total();
         dummy.reset();
         dummy.start();
         for(int i=0; i<1000; i++)
           res=NewCalculateRatioFromMatrixElements(4,Mat,indx.begin());
         dummy.stop();
         double func = dummy.get_total();
         std::cout <<"Direct, function: " <<direct <<"  " <<func << std::endl << std::endl;

         k=0;
         for(int i=0; i<5; i++)
          indx[k++] = i;
         for(int i=0; i<5; i++)
          indx[k++] = i;

         std::cout <<"Timing n=5 \n";
         dummy.start();
         for(int i=0; i<1000; i++)
           res=CalculateRatioFromMatrixElements(5,Mat,indx.begin());
         dummy.stop();
         direct = dummy.get_total();
         dummy.reset();
         dummy.start();
         for(int i=0; i<1000; i++)
           res=NewCalculateRatioFromMatrixElements(5,Mat,indx.begin());
         dummy.stop();
         func = dummy.get_total();
         std::cout <<"Direct, function: " <<direct <<"  " <<func << std::endl << std::endl;

         k=0;
         for(int i=0; i<6; i++)
          indx[k++] = i;
         for(int i=0; i<6; i++)
          indx[k++] = i;

         std::cout <<"Timing n=6 \n";
         dummy.start();
         for(int i=0; i<1000; i++)
           res=CalculateRatioFromMatrixElements(6,Mat,indx.begin());
         dummy.stop();
         direct = dummy.get_total();
         dummy.reset();
         dummy.start();
         for(int i=0; i<1000; i++)
           res=NewCalculateRatioFromMatrixElements(6,Mat,indx.begin());
         dummy.stop();
         func = dummy.get_total();
         std::cout <<"Direct, function: " <<direct <<"  " <<func << std::endl << std::endl;

         exit(1);
      }
  */

  inline ValueType CalculateRatioFromMatrixElements(int n, ValueMatrix_t& dotProducts, std::vector<int>::iterator it)
  {
    switch(n)
    {
    case 0:
      return 1.0;
    case 1:
      return dotProducts(*it,*(it+1));
      break;
    case 2:
    {
      register int i=*it;
      register int j=*(it+1);
      register int a=*(it+2);
      register int b=*(it+3);
      return dotProducts(i,a)*dotProducts(j,b)-dotProducts(i,b)*dotProducts(j,a);
      break;
    }
    case 3:
    {
      register int i1=*it;
      register int i2=*(it+1);
      register int i3=*(it+2);
      register int a1=*(it+3);
      register int a2=*(it+4);
      register int a3=*(it+5);
      return DetCalculator.evaluate(
               dotProducts(i1,a1),dotProducts(i1,a2),dotProducts(i1,a3),
               dotProducts(i2,a1),dotProducts(i2,a2),dotProducts(i2,a3),
               dotProducts(i3,a1),dotProducts(i3,a2),dotProducts(i3,a3));
      break;
    }
    case 4:
    {
      register int i1=*it;
      register int i2=*(it+1);
      register int i3=*(it+2);
      register int i4=*(it+3);
      register int a1=*(it+4);
      register int a2=*(it+5);
      register int a3=*(it+6);
      register int a4=*(it+7);
      return DetCalculator.evaluate(
               dotProducts(i1,a1),dotProducts(i1,a2),dotProducts(i1,a3),dotProducts(i1,a4),
               dotProducts(i2,a1),dotProducts(i2,a2),dotProducts(i2,a3),dotProducts(i2,a4),
               dotProducts(i3,a1),dotProducts(i3,a2),dotProducts(i3,a3),dotProducts(i3,a4),
               dotProducts(i4,a1),dotProducts(i4,a2),dotProducts(i4,a3),dotProducts(i4,a4));
      break;
    }
    case 5:
    {
      register int i1=*it;
      register int i2=*(it+1);
      register int i3=*(it+2);
      register int i4=*(it+3);
      register int i5=*(it+4);
      register int a1=*(it+5);
      register int a2=*(it+6);
      register int a3=*(it+7);
      register int a4=*(it+8);
      register int a5=*(it+9);
      return DetCalculator.evaluate(
               dotProducts(i1,a1),dotProducts(i1,a2),dotProducts(i1,a3),dotProducts(i1,a4),dotProducts(i1,a5),
               dotProducts(i2,a1),dotProducts(i2,a2),dotProducts(i2,a3),dotProducts(i2,a4),dotProducts(i2,a5),
               dotProducts(i3,a1),dotProducts(i3,a2),dotProducts(i3,a3),dotProducts(i3,a4),dotProducts(i3,a5),
               dotProducts(i4,a1),dotProducts(i4,a2),dotProducts(i4,a3),dotProducts(i4,a4),dotProducts(i4,a5),
               dotProducts(i5,a1),dotProducts(i5,a2),dotProducts(i5,a3),dotProducts(i5,a4),dotProducts(i5,a5));
      break;
    }
    /*
              case 6:
              {
                register int i1=*it;
                register int i2=*(it+1);
                register int i3=*(it+2);
                register int i4=*(it+3);
                register int i5=*(it+4);
                register int i6=*(it+5);
                register int a1=*(it+6);
                register int a2=*(it+7);
                register int a3=*(it+8);
                register int a4=*(it+9);
                register int a5=*(it+10);
                register int a6=*(it+11);
                return DetCalculator.evaluate(
                        dotProducts(i1,a1),dotProducts(i1,a2),dotProducts(i1,a3),dotProducts(i1,a4),dotProducts(i1,a5),dotProducts(i1,a6),
                        dotProducts(i2,a1),dotProducts(i2,a2),dotProducts(i2,a3),dotProducts(i2,a4),dotProducts(i2,a5),dotProducts(i2,a6),
                        dotProducts(i3,a1),dotProducts(i3,a2),dotProducts(i3,a3),dotProducts(i3,a4),dotProducts(i3,a5),dotProducts(i3,a6),
                        dotProducts(i4,a1),dotProducts(i4,a2),dotProducts(i4,a3),dotProducts(i4,a4),dotProducts(i4,a5),dotProducts(i4,a6),
                        dotProducts(i5,a1),dotProducts(i5,a2),dotProducts(i5,a3),dotProducts(i5,a4),dotProducts(i5,a5),dotProducts(i5,a6),
                        dotProducts(i6,a1),dotProducts(i6,a2),dotProducts(i6,a3),dotProducts(i6,a4),dotProducts(i6,a5),dotProducts(i6,a6));
                break;
              }
    */
    default:
    {
      return DetCalculator.evaluate(dotProducts,it,n);
    }
    }
    return 0.0;
  }
  /*
      inline ValueType NewCalculateRatioFromMatrixElements(int n, ValueMatrix_t& dotProducts, std::vector<int>::iterator it)
      {
         switch(n)
         {
            case 0:
              return 1.0;
            case 1:
              return dotProducts(*it,*(it+1));
              break;
            case 2:
            {
              register int i=*it;
              register int j=*(it+1);
              register int a=*(it+2);
              register int b=*(it+3);
              return dotProducts(i,a)*dotProducts(j,b)-dotProducts(i,b)*dotProducts(j,a);
              break;
            }
            default:
            {
              return DetCalculator.evaluate(dotProducts,it,n);
            }
         }
      }
  */

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, ValueVector_t& ratios, ValueMatrix_t &psiinv, ValueMatrix_t &psi, ValueMatrix_t& dotProducts, std::vector<int>& data, std::vector<std::pair<int,int> >& pairs, std::vector<double>& sign)
  {
    ValueType det0 = ratios[ref];
    buildTableTimer.start();
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    buildTableTimer.stop();
    readMatTimer.start();
    std::vector<int>::iterator it2 = data.begin();
    //ValueVector_t::iterator itr = ratios.begin();
    //vector<double>::iterator its = sign.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        //itr++;
        //its++;
        count++;
        continue;
      }
      ratios[count] = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      //*(itr++) = *(its++)*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
    readMatTimer.stop();
  }

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, GradMatrix_t& ratios, ValueMatrix_t& psiinv, ValueMatrix_t& psi, ValueMatrix_t& dotProducts, std::vector<int>& data, std::vector<std::pair<int,int> >& pairs, std::vector<double>& sign, int dx)
  {
    ValueType det0 = ratios(ref,iat)[dx];
    buildTableGradTimer.start();
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    buildTableGradTimer.stop();
    readMatGradTimer.start();
    std::vector<int>::iterator it2 = data.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        count++;
        continue;
      }
      ratios(count,iat)[dx] = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
    readMatGradTimer.stop();
  }

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, ValueMatrix_t& ratios, ValueMatrix_t& psiinv, ValueMatrix_t& psi, ValueMatrix_t& dotProducts, std::vector<int>& data, std::vector<std::pair<int,int> >& pairs, std::vector<double>& sign)
  {
    ValueType det0 = ratios(ref,iat);
    int num=psi.extent(1);
    std::vector<std::pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    std::vector<int>::iterator it2 = data.begin();
    int count= 0;  // number of determinants processed
    while(it2 != data.end())
    {
      int n = *it2; // number of excitations
      if(count == ref)
      {
        it2+=3*n+1;  // number of integers used to encode the current excitation
        count++;
        continue;
      }
      ratios(count,iat) = sign[count]*det0*CalculateRatioFromMatrixElements(n,dotProducts,it2+1);
      count++;
      it2+=3*n+1;
    }
  }

//   Finish this at some point
  inline void InverseUpdateByColumn_GRAD(ValueMatrix_t& Minv
                                         , GradVector_t& newcol, ValueVector_t& rvec
                                         , ValueVector_t& rvecinv, int colchanged
                                         , ValueType c_ratio, int dx)
  {
    c_ratio=1.0/c_ratio;
    int m = Minv.rows();
    BLAS::gemv('N', m, m, c_ratio, Minv.data(), m, newcol[0].data()+dx, 3, 0.0, rvec.data(), 1);
    rvec[colchanged]=1.0-c_ratio;
    BLAS::copy(m,Minv.data()+colchanged,m,rvecinv.data(),1);
    BLAS::ger(m,m,-1.0,rvec.data(),1,rvecinv.data(),1,Minv.data(),m);
  }
  /*
      inline void InverseUpdateByColumn(ValueMatrix_t& Minv
            , GradVector_t& dM, ValueVector_t& rvec
            , ValueVector_t& rvecinv, int colchanged
            , ValueType c_ratio, std::vector<int>::iterator& it)
      {
            ValueType c_ratio=1.0/ratioLapl;
            BLAS::gemv('N', NumPtcls, NumPtcls, c_ratio, dpsiMinv.data(), NumPtcls, tv, 1, T(), workV1, 1);
            workV1[colchanged]=1.0-c_ratio;
            BLAS::copy(m,pinv+colchanged,m,workV2,1);
            BLAS::ger(m,m,-1.0,workV1,1,workV2,1,dpsiMinv.data(),m);
      }
  */

  void setDetInfo(int ref, std::vector<ci_configuration2> list);

  void evaluateDetsForPtclMove(ParticleSet& P, int iat);

  void evaluateDetsAndGradsForPtclMove(ParticleSet& P, int iat);

  void evaluateGrads(ParticleSet& P, int iat);

  void evaluateAllForPtclMove(ParticleSet& P, int iat);

  // full evaluation of all the structures from scratch, used in evaluateLog for example
  void evaluateForWalkerMove(ParticleSet& P, bool fromScratch=true);

  ///total number of particles
  int NP;
  ///number of single-particle orbitals which belong to this Dirac determinant
  int NumOrbitals;
  ///number of particles which belong to this Dirac determinant
  int NumPtcls;
  ///index of the first particle with respect to the particle set
  int FirstIndex;
  ///index of the last particle with respect to the particle set
  int LastIndex;
  ///index of the particle (or row)
  int WorkingIndex;
  ///a set of single-particle orbitals used to fill in the  values of the matrix
  SPOSetBasePtr Phi;
  /// number of determinants handled by this object
  int NumDets;
  ///
  std::vector<ci_configuration2> confgList;
// the reference determinant never changes, so there is no need to store it.
// if its value is zero, then use a data from backup, but always use this one
// by default
  int ReferenceDeterminant;

  /// store determinants (old and new)
  ValueVector_t detValues, new_detValues;
  GradMatrix_t grads,new_grads;
  GradMatrix_t grads_bf;
  ValueMatrix_t lapls,new_lapls;

  /// psiM(i,j) \f$= \psi_j({\bf r}_i)\f$
  /// TpsiM(i,j) \f$= psiM(j,i) \f$
  ValueMatrix_t psiM,psiM_temp, TpsiM, psiMinv, psiMinv_temp;
  /// dpsiM(i,j) \f$= \nabla_i \psi_j({\bf r}_i)\f$
  GradMatrix_t  dpsiM,dpsiM_temp;
  // temporaty storage
  ValueMatrix_t dpsiMinv;
  /// d2psiM(i,j) \f$= \nabla_i^2 \psi_j({\bf r}_i)\f$
  ValueMatrix_t d2psiM,d2psiM_temp;

  HessMatrix_t grad_grad_psiM;
  HessVector_t grad_gradV;
  HessMatrix_t grad_grad_psiM_temp;
  GGGMatrix_t  grad_grad_grad_psiM;
  ParticleSet::ParticleGradient_t Gtemp;
  HessMatrix_t Ajk_sum,Qmat;
  GradMatrix_t Fmat;
  GradVector_t Fmatdiag;
  GradVector_t Fmatdiag_temp;
  GradMatrix_t dFa;

  /// value of single-particle orbital for particle-by-particle update
  ValueVector_t psiV, psiV_temp;
  GradVector_t dpsiV;
  ValueVector_t d2psiV;
  ValueVector_t workV1, workV2;

  ValueMatrix_t dotProducts;

  Vector<ValueType> WorkSpace;
  Vector<IndexType> Pivot;

  ValueType curRatio,cumRatio;
  ValueType *FirstAddressOfGrads;
  ValueType *LastAddressOfGrads;
  ValueType *FirstAddressOfdpsiM;
  ValueType *LastAddressOfdpsiM;

  /* mmorales:
   *  i decided to store the excitation information of all determinants in the following
   *  compact form: (all these integers are found consecutively in the array)
   *    For each determinant:
   *     -n : number of excitations
   *     -i1,i2,...,in : occupied orbital to be replaced (these must be numbers from 0:Nptcl-1)
   *     -a1,a2,...,an : excited states that replace the orbitals (these can be anything)
   */
  std::vector<int> detData;
  std::vector<std::pair<int,int> > uniquePairs;
  std::vector<double> DetSigns;

  int backup_reference;
  std::vector<int> backup_detData;
  std::vector<std::pair<int,int> > backup_uniquePairs;
  std::vector<double> backup_DetSigns;

  // transformation
  BackflowTransformation *BFTrans;
  bool usingBF;

  MyDeterminant<ValueType> DetCalculator;

};

//  #include "MultiDiracDeterminantBase_help.h"

}
#endif
