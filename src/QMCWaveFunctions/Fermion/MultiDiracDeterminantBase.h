//////////////////////////////////////////////////////////////////
// (c) Copyright 1998-2002,2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
/**@file DiracDeterminantBaseBase.h
 * @brief Declaration of DiracDeterminantBase with a S(ingle)P(article)O(rbital)SetBase
 */
#ifndef QMCPLUSPLUS_MULTIDIRACDETERMINANTWITHBASE_H
#define QMCPLUSPLUS_MULTIDIRACDETERMINANTWITHBASE_H
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/Fermion/ci_configuration2.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminantBase_help.h"
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

  void setBF(BackflowTransformation* bf) {}

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

  inline void reportStatus(ostream& os)
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

  RealType registerData(ParticleSet& P, PooledData<RealType>& buf);

  void registerDataForDerivatives(ParticleSet& P, PooledData<RealType>& buf, int storageType=0);

  virtual void memoryUsage_DataForDerivatives(ParticleSet& P,long& orbs_only, long& orbs, long& invs, long& dets)
  {
    orbs_only += NumPtcls*NumOrbitals;
    orbs += NumPtcls*NumOrbitals*5;
    invs += NumPtcls*NumPtcls;
    dets += NumDets*(1 + 4*NumPtcls);
  }

  void copyToDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);

  void copyFromDerivativeBuffer(ParticleSet& P, PooledData<RealType>& buf);

  RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf, bool fromscratch=false);

  void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf);

  /** dump the inverse to the buffer
   */
  void dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    APP_ABORT("  Need to implement MSDFast::dumpToBuffer. \n");
  }

  /** copy the inverse from the buffer
   */
  void dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    APP_ABORT("  Need to implement MSDFast::dumpFromBuffer. \n");
  }

  /** move was accepted, update the real container
   */
  void acceptMove(ParticleSet& P, int iat);

  /** move was rejected. copy the real container to the temporary to move on
   */
  void restore(int iat);

  void update(ParticleSet& P,
              ParticleSet::ParticleGradient_t& dG,
              ParticleSet::ParticleLaplacian_t& dL,
              int iat);

  RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf);

  OrbitalBasePtr makeClone(ParticleSet& tqp) const;

  /****************************************************************************
   * These functions should not be called.
   ***************************************************************************/

  ValueType ratio(ParticleSet& P, int iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return ValueType();
  }

  GradType evalGrad(ParticleSet& P, int iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return GradType();
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return ValueType();
  }

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return ValueType();
  }

  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return 0.0;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    APP_ABORT("  MultiDiracDeterminantBase: This should not be called. \n");
    return ValueType();
  }


  /****************************************************************************
   * END END END
   ***************************************************************************/

  // create necessary structures used in the evaluation of the determinants
  // this works with confgList, which shouldn't change during a simulation
  void createDetData(ci_configuration2& ref, vector<int>& data,
                     vector<pair<int,int> >& pairs, vector<double>& sign);
  /*
      void testDets() {

         vector<int> indx(16);
         ValueMatrix_t Mat(8,8);
         NewTimer dummy("dummy");
         dummy.reset();

         ifstream in("matrix.txt");
         for(int i=0; i<8; i++)
         for(int j=0; j<8; j++)
            in>>Mat(i,j);
         in.close();

         vector<double> dets(9);
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
           cout<<"determinant of " <<q <<": " <<dets[q] <<"  -  " <<CalculateRatioFromMatrixElements(q,Mat,indx.begin()) <<endl;
         }
         int k=0;
         for(int i=0; i<4; i++)
          indx[k++] = i;
         for(int i=0; i<4; i++)
          indx[k++] = i;

         cout<<"Timing n=4 \n";
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
         cout<<"Direct, function: " <<direct <<"  " <<func <<endl <<endl;

         k=0;
         for(int i=0; i<5; i++)
          indx[k++] = i;
         for(int i=0; i<5; i++)
          indx[k++] = i;

         cout<<"Timing n=5 \n";
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
         cout<<"Direct, function: " <<direct <<"  " <<func <<endl <<endl;

         k=0;
         for(int i=0; i<6; i++)
          indx[k++] = i;
         for(int i=0; i<6; i++)
          indx[k++] = i;

         cout<<"Timing n=6 \n";
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
         cout<<"Direct, function: " <<direct <<"  " <<func <<endl <<endl;

         exit(1);
      }
  */

  inline ValueType CalculateRatioFromMatrixElements(int n, ValueMatrix_t& dotProducts, vector<int>::iterator it)
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
      inline ValueType NewCalculateRatioFromMatrixElements(int n, ValueMatrix_t& dotProducts, vector<int>::iterator it)
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

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, ValueVector_t& ratios, ValueMatrix_t &psiinv, ValueMatrix_t &psi, ValueMatrix_t& dotProducts, vector<int>& data, vector<pair<int,int> >& pairs, vector<double>& sign)
  {
    ValueType det0 = ratios[ref];
    buildTableTimer.start();
    int num=psi.extent(1);
    vector<pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    buildTableTimer.stop();
    readMatTimer.start();
    vector<int>::iterator it2 = data.begin();
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

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, GradMatrix_t& ratios, ValueMatrix_t& psiinv, ValueMatrix_t& psi, ValueMatrix_t& dotProducts, vector<int>& data, vector<pair<int,int> >& pairs, vector<double>& sign, int dx)
  {
    ValueType det0 = ratios(ref,iat)[dx];
    buildTableGradTimer.start();
    int num=psi.extent(1);
    vector<pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    buildTableGradTimer.stop();
    readMatGradTimer.start();
    vector<int>::iterator it2 = data.begin();
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

  inline void BuildDotProductsAndCalculateRatios(int ref, int iat, ValueMatrix_t& ratios, ValueMatrix_t& psiinv, ValueMatrix_t& psi, ValueMatrix_t& dotProducts, vector<int>& data, vector<pair<int,int> >& pairs, vector<double>& sign)
  {
    ValueType det0 = ratios(ref,iat);
    int num=psi.extent(1);
    vector<pair<int,int> >::iterator it(pairs.begin()), last(pairs.end());
    while(it != last)
    {
      dotProducts((*it).first,(*it).second) = simd::dot(psiinv[(*it).first],psi[(*it).second],num);
      it++;
    }
    vector<int>::iterator it2 = data.begin();
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
            , ValueType c_ratio, vector<int>::iterator& it)
      {
            ValueType c_ratio=1.0/ratioLapl;
            BLAS::gemv('N', NumPtcls, NumPtcls, c_ratio, dpsiMinv.data(), NumPtcls, tv, 1, T(), workV1, 1);
            workV1[colchanged]=1.0-c_ratio;
            BLAS::copy(m,pinv+colchanged,m,workV2,1);
            BLAS::ger(m,m,-1.0,workV1,1,workV2,1,dpsiMinv.data(),m);
      }
  */

  void setDetInfo(int ref, vector<ci_configuration2> list);

  inline void
  evaluateDetsForPtclMove(ParticleSet& P, int iat)
  {
    UpdateMode=ORB_PBYP_RATIO;
    RatioTimer.start();
    evalOrbTimer.start();
    Phi->evaluate(P,iat,psiV);
    evalOrbTimer.stop();
    WorkingIndex = iat-FirstIndex;
    if(NumPtcls==1)
    {
      vector<ci_configuration2>::iterator it(confgList.begin());
      vector<ci_configuration2>::iterator last(confgList.end());
      ValueVector_t::iterator det(new_detValues.begin());
      while(it != last)
      {
        int orb = (it++)->occup[0];
        *(det++) = psiV(orb);
      }
    }
    else
    {
      vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
// mmorales: the only reason this is here is because
// NonlocalECP do not necessarily call rejectMove after
// calling ratio(), and even if the move is rejected
// this matrix needs to be restored
// If we always restore after ratio, then this is not needed
// For efficiency reasons, I don't do this for ratioGrad or ratio(P,dG,dL)
      ExtraStuffTimer.start();
      psiMinv_temp = psiMinv;
      for(int i=0; i<NumPtcls; i++)
        psiV_temp[i] = psiV(*(it++));
      ValueType ratioRef = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
      new_detValues[ReferenceDeterminant] = ratioRef * detValues[ReferenceDeterminant];
      InverseUpdateByColumn(psiMinv_temp,psiV_temp,workV1,workV2,WorkingIndex,ratioRef);
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiV(i);
      ExtraStuffTimer.stop();
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant,WorkingIndex,new_detValues,psiMinv_temp,TpsiM,dotProducts,detData,uniquePairs,DetSigns);
// check comment above
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
    }
    RatioTimer.stop();
  }

  inline void
  evaluateDetsAndGradsForPtclMove(ParticleSet& P, int iat)
  {
    UpdateMode=ORB_PBYP_PARTIAL;
    evalOrb1Timer.start();
    Phi->evaluate(P,iat,psiV,dpsiV,d2psiV);
    evalOrb1Timer.stop();
    WorkingIndex = iat-FirstIndex;
    if(NumPtcls==1)
    {
      vector<ci_configuration2>::iterator it(confgList.begin());
      vector<ci_configuration2>::iterator last(confgList.end());
      ValueVector_t::iterator det(new_detValues.begin());
      GradMatrix_t::iterator grad(new_grads.begin());
      while(it != last)
      {
        int orb = (it++)->occup[0];
        *(det++) = psiV[orb];
        *(grad++) = dpsiV[orb];
      }
    }
    else
    {
      ExtraStuffTimer.start();
//mmorales: check comment above
      psiMinv_temp = psiMinv;
      vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
      GradType ratioGradRef;
      for(int i=0; i<NumPtcls; i++)
      {
        psiV_temp[i] = psiV(*it);
        ratioGradRef += psiMinv_temp(i,WorkingIndex)*dpsiV(*it);
        it++;
      }
      ValueType ratioRef = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
      new_grads(ReferenceDeterminant,WorkingIndex) = ratioGradRef*detValues[ReferenceDeterminant];
      new_detValues[ReferenceDeterminant] = ratioRef * detValues[ReferenceDeterminant];
      InverseUpdateByColumn(psiMinv_temp,psiV_temp,workV1,workV2,WorkingIndex,ratioRef);
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiV[i];
      ExtraStuffTimer.stop();
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant,WorkingIndex,new_detValues,psiMinv_temp,TpsiM,dotProducts,detData,uniquePairs,DetSigns);
      for(int idim=0; idim<OHMMS_DIM; idim++)
      {
        ExtraStuffTimer.start();
        //dpsiMinv = psiMinv_temp;
        dpsiMinv = psiMinv;
        it = confgList[ReferenceDeterminant].occup.begin();
        for(int i=0; i<NumPtcls; i++)
          psiV_temp[i] = dpsiV(*(it++))[idim];
        InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,WorkingIndex,ratioGradRef[idim]);
        for(int i=0; i<NumOrbitals; i++)
          TpsiM(i,WorkingIndex) = dpsiV[i][idim];
        ExtraStuffTimer.stop();
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant,WorkingIndex,new_grads,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns,idim);
      }
// check comment above
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
    }
  }

  inline void
  evaluateGrads(ParticleSet& P, int iat)
  {
    WorkingIndex = iat-FirstIndex;
    vector<int>::iterator it;
    if(NumPtcls==1)
    {
      vector<ci_configuration2>::iterator it(confgList.begin());
      vector<ci_configuration2>::iterator last(confgList.end());
      GradMatrix_t::iterator grad(grads.begin());
      while(it != last)
      {
        int orb = (it++)->occup[0];
        *(grad++) = dpsiM(0,orb);
      }
    }
    else
    {
      for(int idim=0; idim<OHMMS_DIM; idim++)
      {
        //dpsiMinv = psiMinv_temp;
        dpsiMinv = psiMinv;
        it = confgList[ReferenceDeterminant].occup.begin();
        ValueType ratioG = 0.0;
        for(int i=0; i<NumPtcls; i++)
        {
          psiV_temp[i] = dpsiM(WorkingIndex,*it)[idim];
          ratioG += psiMinv(i,WorkingIndex)*dpsiM(WorkingIndex,*it)[idim];
          it++;
        }
        grads(ReferenceDeterminant,WorkingIndex)[idim] = ratioG*detValues[ReferenceDeterminant];
        InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,WorkingIndex,ratioG);
        for(int i=0; i<NumOrbitals; i++)
          TpsiM(i,WorkingIndex) = dpsiM(WorkingIndex,i)[idim];
        BuildDotProductsAndCalculateRatios(ReferenceDeterminant,WorkingIndex,grads,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns,idim);
      }
// check comment above
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
    }
  }

  inline void
  evaluateAllForPtclMove(ParticleSet& P, int iat)
  {
    UpdateMode=ORB_PBYP_ALL;
    Phi->evaluate(P,iat,psiV,dpsiV,d2psiV);
    WorkingIndex = iat-FirstIndex;
    if(NumPtcls==1)
    {
      vector<ci_configuration2>::iterator it(confgList.begin());
      vector<ci_configuration2>::iterator last(confgList.end());
      ValueVector_t::iterator det(new_detValues.begin());
      ValueMatrix_t::iterator lap(new_lapls.begin());
      GradMatrix_t::iterator grad(new_grads.begin());
      while(it != last)
      {
        int orb = (it++)->occup[0];
        *(det++) = psiV(orb);
        *(lap++) = d2psiV(orb);
        *(grad++) = dpsiV(orb);
      }
    }
    else
    {
//mmorales: check comment above
      psiMinv_temp = psiMinv;
      vector<int>::iterator it(confgList[ReferenceDeterminant].occup.begin());
      //GradType ratioGradRef;
      //ValueType ratioLapl = 0.0;
      for(int i=0; i<NumPtcls; i++)
      {
        psiV_temp[i] = psiV(*it);
        //ratioGradRef += psiMinv_temp(i,WorkingIndex)*dpsiV(*it);
        //ratioLapl += psiMinv_temp(i,WorkingIndex)*d2psiV(*it);
        it++;
      }
      ValueType det0, ratioRef = DetRatioByColumn(psiMinv_temp, psiV_temp, WorkingIndex);
      //new_lapls(ReferenceDeterminant,WorkingIndex) = ratioLapl*detValues[ReferenceDeterminant];
      //new_grads(ReferenceDeterminant,WorkingIndex) = ratioGradRef*detValues[ReferenceDeterminant];
      det0 = ratioRef * detValues[ReferenceDeterminant];
      new_detValues[ReferenceDeterminant] = det0;
      InverseUpdateByColumn(psiMinv_temp,psiV_temp,workV1,workV2,WorkingIndex,ratioRef);
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiV[i];
      BuildDotProductsAndCalculateRatios(ReferenceDeterminant,WorkingIndex,new_detValues,psiMinv_temp,TpsiM,dotProducts,detData,uniquePairs,DetSigns);
      for(int jat=0; jat<NumPtcls; jat++)
      {
        it = confgList[ReferenceDeterminant].occup.begin();
        GradType gradRatio;// = 0.0;
        ValueType ratioLapl = 0.0;
        if(jat==WorkingIndex)
        {
          for(int i=0; i<NumPtcls; i++)
          {
            gradRatio += psiMinv_temp(i,jat)*dpsiV(*it);
            ratioLapl += psiMinv_temp(i,jat)*d2psiV(*it);
            it++;
          }
          new_grads(ReferenceDeterminant,jat) = det0*gradRatio;
          new_lapls(ReferenceDeterminant,jat) = det0*ratioLapl;
          for(int idim=0; idim<OHMMS_DIM; idim++)
          {
            dpsiMinv = psiMinv_temp;
            it = confgList[ReferenceDeterminant].occup.begin();
            for(int i=0; i<NumPtcls; i++)
              psiV_temp[i] = dpsiV(*(it++))[idim];
            InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,jat,gradRatio[idim]);
            for(int i=0; i<NumOrbitals; i++)
              TpsiM(i,jat) = dpsiV[i][idim];
            BuildDotProductsAndCalculateRatios(ReferenceDeterminant,jat,new_grads,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns,idim);
          }
          dpsiMinv = psiMinv_temp;
          it = confgList[ReferenceDeterminant].occup.begin();
          for(int i=0; i<NumPtcls; i++)
            psiV_temp[i] = d2psiV(*(it++));
          InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,jat,ratioLapl);
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,jat) = d2psiV[i];
          BuildDotProductsAndCalculateRatios(ReferenceDeterminant,jat,new_lapls,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns);
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,jat) = psiV[i];
        }
        else
        {
          for(int i=0; i<NumPtcls; i++)
          {
            gradRatio += psiMinv_temp(i,jat)*dpsiM(jat,*it);
            ratioLapl += psiMinv_temp(i,jat)*d2psiM(jat,*it);
            it++;
          }
          new_grads(ReferenceDeterminant,jat) = det0*gradRatio;
          new_lapls(ReferenceDeterminant,jat) = det0*ratioLapl;
          for(int idim=0; idim<OHMMS_DIM; idim++)
          {
            dpsiMinv = psiMinv_temp;
            it = confgList[ReferenceDeterminant].occup.begin();
            for(int i=0; i<NumPtcls; i++)
              psiV_temp[i] = dpsiM(jat,*(it++))[idim];
            InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,jat,gradRatio[idim]);
            for(int i=0; i<NumOrbitals; i++)
              TpsiM(i,jat) = dpsiM(jat,i)[idim];
            BuildDotProductsAndCalculateRatios(ReferenceDeterminant,jat,new_grads,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns,idim);
          }
          dpsiMinv = psiMinv_temp;
          it = confgList[ReferenceDeterminant].occup.begin();
          for(int i=0; i<NumPtcls; i++)
            psiV_temp[i] = d2psiM(jat,*(it++));
          InverseUpdateByColumn(dpsiMinv,psiV_temp,workV1,workV2,jat,ratioLapl);
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,jat) = d2psiM(jat,i);
          BuildDotProductsAndCalculateRatios(ReferenceDeterminant,jat,new_lapls,dpsiMinv,TpsiM,dotProducts,detData,uniquePairs,DetSigns);
          for(int i=0; i<NumOrbitals; i++)
            TpsiM(i,jat) = psiM(jat,i);
        }
      } // jat
// check comment above
      for(int i=0; i<NumOrbitals; i++)
        TpsiM(i,WorkingIndex) = psiM(WorkingIndex,i);
    }
  }


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
  vector<ci_configuration2> confgList;
// the reference determinant never changes, so there is no need to store it.
// if its value is zero, then use a data from backup, but always use this one
// by default
  int ReferenceDeterminant;

  /// store determinants (old and new)
  ValueVector_t detValues, new_detValues;
  GradMatrix_t grads,new_grads;
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
   *  i decided to stored the excitation information of all determinants in the following
   *  compact form: (all these integers are found consecutively in the array)
   *    For each determinant:
   *     -n : number of excitations
   *     -i1,i2,...,in : occupied orbital to be replaced (these must be numbers from 0:Nptcl-1)
   *     -a1,a2,...,an : excited states that replace the orbitals (these can be anything)
   */
  vector<int> detData;
  vector<pair<int,int> > uniquePairs;
  vector<double> DetSigns;

  int backup_reference;
  vector<int> backup_detData;
  vector<pair<int,int> > backup_uniquePairs;
  vector<double> backup_DetSigns;

  MyDeterminant<ValueType> DetCalculator;

};

//  #include "MultiDiracDeterminantBase_help.h"

}
#endif
/***************************************************************************
 * $RCSfile$   $Author: kesler $
 * $Revision: 3574 $   $Date: 2009-02-19 17:11:24 -0600 (Thu, 19 Feb 2009) $
 * $Id: DiracDeterminantBase.h 3574 2009-02-19 23:11:24Z kesler $
 ***************************************************************************/
