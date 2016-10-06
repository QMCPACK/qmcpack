//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Bryan Clark, bclark@Princeton.edu, Princeton University
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"



namespace qmcplusplus
{

class ExcitationClass
{

  typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
  typedef SPOSetBase::ValueVector_t ValueVector_t;
  typedef SPOSetBase::GradMatrix_t  GradMatrix_t;
  typedef SPOSetBase::GradVector_t  GradVector_t;
  typedef double ValueType;
  typedef TinyVector<ValueType,3> GradType;

public:
  std::vector<int> unoccupied_orbitals_to_use;
  std::vector<int> orbitals_to_replace;
  ValueMatrix_t dotProducts; // orbitals_to_replace.size() x unoccupied_orbitals_to_use.size()
  GradMatrix_t dotProducts_grad; // orbitals_to_replace.size() x unoccupied_orbitals_to_use.size()
  TinyVector<int,2> ii;
  TinyVector<int,2> jj;

  TinyVector<int,2> second_replaces_first;
  std::list<TinyVector<int,2> > replaceList;
  int orbitals_to_replace_index;
  int unoccupied_orbitals_to_use_index;
  void testMe()
  {
    int i=5;
  }




  TinyVector<int,2>  begin()
  {
    orbitals_to_replace_index=0;
    unoccupied_orbitals_to_use_index=0;
    if (orbitals_to_replace.size()==0 || unoccupied_orbitals_to_use.size()==0)
    {
      second_replaces_first[0]=-1;
      second_replaces_first[1]=-1;
      return second_replaces_first;
    }
    second_replaces_first[0]=orbitals_to_replace[orbitals_to_replace_index];
    second_replaces_first[1]=unoccupied_orbitals_to_use[unoccupied_orbitals_to_use_index];
    return second_replaces_first;
  }

  TinyVector<int,2> next()
  {
    unoccupied_orbitals_to_use_index++;
    if (unoccupied_orbitals_to_use_index==unoccupied_orbitals_to_use.size())
    {
      unoccupied_orbitals_to_use_index=0;
      orbitals_to_replace_index++;
      if (orbitals_to_replace_index==orbitals_to_replace.size())
      {
        second_replaces_first[0]=-1;
        second_replaces_first[1]=-1;
        return second_replaces_first;
      }
    }
    second_replaces_first[0]=orbitals_to_replace[orbitals_to_replace_index];
    second_replaces_first[1]=unoccupied_orbitals_to_use[unoccupied_orbitals_to_use_index];
    return second_replaces_first;
  }


  //gs_inverse needs to be orbital x ptcl
  //psiM_actual needs to be orbtial x ptcl
  void BuildDotProducts(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual)
  {
    int NumPtcls=psiM_actual.extent(1);
    for (int ii=0; ii<orbitals_to_replace.size(); ii++)
      for (int jj=0; jj<unoccupied_orbitals_to_use.size(); jj++)
      {
        int i=orbitals_to_replace[ii];
        int j=unoccupied_orbitals_to_use[jj];
        dotProducts(ii,jj)=0.0;
        for (int ptcl=0; ptcl<NumPtcls; ptcl++)
          dotProducts(ii,jj)+=gs_inverse(i,ptcl)*psiM_actual(j,ptcl);
      }
  }


//     void BuildDotProducts_grad_slow(GradMatrix_t &gs_inverse, GradMatrix_t &dpsiM_actual)
//     {
//       int NumPtcls=dpsiM_actual.extent(1);
//       for (int ii=0;ii<orbitals_to_replace.size();ii++)
// 	for (int jj=0;jj<unoccupied_orbitals_to_use.size();jj++){
// 	  int i=orbitals_to_replace[ii];
// 	  int j=unoccupied_orbitals_to_use[jj];
// 	  dotProducts_grad(ii,jj)=0.0;
// 	  for (int ptcl=0;ptcl<NumPtcls;ptcl++)
// 	    dotProducts(ii,jj)+=gs_inverse(i,ptcl)*dpsiM_actual(ptcl,j);
// 	}


//     }


///Assumes: gs_inverse is ptcl x orbital
///Assumes: dpsiM_actual is ptcl x orbital

  void BuildDotProducts_grad(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual,GradMatrix_t &dpsiM_actual, GradVector_t &op2,GradType one_over_ratio,int WorkingIndex)
  {
    //note: op1(j) should be gs_inverse(WorkingIndex,j);
    GradMatrix_t gs_inverse_update;
    gs_inverse_update.resize(gs_inverse.extent(0),gs_inverse.extent(1));
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int i=0; i<gs_inverse.extent(0); i++)
      {
        for (int j=0; j<gs_inverse.extent(1); j++)
        {
          gs_inverse_update(i,j)[dim]=gs_inverse(i,j)+one_over_ratio[dim]*op2(i)[dim]*gs_inverse(WorkingIndex,j);
        }
      }
    }
    for (int i=0; i<gs_inverse_update.extent(0); i++)
      for (int j=0; j<gs_inverse_update.extent(1); j++)
        std::cerr <<i<<" "<<j<<" "<<gs_inverse_update(i,j)[0]<< std::endl;
    ////BE CAREFUL IF YOU CHANGE THE TRANPSOSE OF dpsiM!!
    int NumPtcls=dpsiM_actual.extent(0);
    assert(NumPtcls==psiM_actual.extent(1));
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int ii=0; ii<orbitals_to_replace.size(); ii++)
        for (int jj=0; jj<unoccupied_orbitals_to_use.size(); jj++)
        {
          int i=orbitals_to_replace[ii];
          int j=unoccupied_orbitals_to_use[jj];
          dotProducts_grad(ii,jj)[dim]=0.0;
          std::cerr <<"enter-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
          for (int ptcl=0; ptcl<NumPtcls; ptcl++)
          {
            //	    GradType toAdd=one_over_ratio*op2(i)*op1(j);
            std::cerr <<"on-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
            dotProducts_grad(ii,jj)[dim]+=gs_inverse_update(ptcl,i)[dim]*psiM_actual(j,ptcl);
            std::cerr <<"on-fix af "<<ii<<" "<<jj<<" "<<i<<" "<<j<<" "<<ptcl<<" "<<dotProducts_grad(ii,jj)<<" "<<gs_inverse_update(ptcl,i)[dim]<<" "<<psiM_actual(j,ptcl)<< std::endl;
          }
          std::cerr <<"pre-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
          dotProducts_grad(ii,jj)[dim]-=gs_inverse_update(WorkingIndex,i)[dim]*psiM_actual(j,WorkingIndex);
          dotProducts_grad(ii,jj)[dim]+=gs_inverse_update(WorkingIndex,i)[dim]*dpsiM_actual(WorkingIndex,j)[dim];
          std::cerr <<"post-fix "<<dotProducts_grad(ii,jj)<< std::endl;
        }
    }
  }


  void CalcSingleExcitations_grad(ValueVector_t &coefs, GradType &val,int &coefIndex)
  {
    for (int i=0; i<orbitals_to_replace.size(); i++)
      for (int j=0; j<unoccupied_orbitals_to_use.size(); j++)
      {
        for (int dim=0; dim<OHMMS_DIM; dim++)
        {
          val[dim]+=dotProducts_grad(i,j)[dim];// *coefs(coefIndex);
        }
        coefIndex++;
      }
  }


  //coefIndex should be set to start the double excitations
  void CalcDoubleExcitations_grad(ValueVector_t &coefs,GradType &val,int &coefIndex)
  {
    for (int ii1=0; ii1<orbitals_to_replace.size(); ii1++)
    {
      for (int ii2=ii1+1; ii2<orbitals_to_replace.size(); ii2++)
      {
        for (int jj1=0; jj1<unoccupied_orbitals_to_use.size(); jj1++)
        {
          for (int jj2=jj1+1; jj2<unoccupied_orbitals_to_use.size(); jj2++)
          {
            for (int dim=0; dim<OHMMS_DIM; dim++)
            {
              ValueType a=dotProducts_grad(ii1,jj1)[dim];
              ValueType d=dotProducts_grad(ii2,jj2)[dim];
              ValueType b=dotProducts_grad(ii1,jj2)[dim];
              ValueType c=dotProducts_grad(ii2,jj1)[dim];
              ValueType det=a*d-b*c;
              val[dim]+=coefs(coefIndex)*det;
            }
            coefIndex++;
          }
        }
      }
    }
  }


  //////////////////////////////////////////////////////////


  void BuildDotProducts_lapl(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual,ValueMatrix_t &dpsiM_actual, ValueVector_t &op2,GradType one_over_ratio,int WorkingIndex)
  {
    //note: op1(j) should be gs_inverse(WorkingIndex,j);
    GradMatrix_t gs_inverse_update;
    gs_inverse_update.resize(gs_inverse.extent(0),gs_inverse.extent(1));
    int dim=0;
    {
      for (int i=0; i<gs_inverse.extent(0); i++)
      {
        for (int j=0; j<gs_inverse.extent(1); j++)
        {
          gs_inverse_update(i,j)[dim]=gs_inverse(i,j)+one_over_ratio[dim]*op2(i)*gs_inverse(WorkingIndex,j);
        }
      }
    }
    for (int i=0; i<gs_inverse_update.extent(0); i++)
      for (int j=0; j<gs_inverse_update.extent(1); j++)
        std::cerr <<i<<" "<<j<<" "<<gs_inverse_update(i,j)[0]<< std::endl;
    ////BE CAREFUL IF YOU CHANGE THE TRANPSOSE OF dpsiM!!
    int NumPtcls=dpsiM_actual.extent(0);
    assert(NumPtcls==psiM_actual.extent(1));
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int ii=0; ii<orbitals_to_replace.size(); ii++)
        for (int jj=0; jj<unoccupied_orbitals_to_use.size(); jj++)
        {
          int i=orbitals_to_replace[ii];
          int j=unoccupied_orbitals_to_use[jj];
          dotProducts_grad(ii,jj)[dim]=0.0;
          std::cerr <<"enter-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
          for (int ptcl=0; ptcl<NumPtcls; ptcl++)
          {
            //	    GradType toAdd=one_over_ratio*op2(i)*op1(j);
            std::cerr <<"on-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
            dotProducts_grad(ii,jj)[dim]+=gs_inverse_update(ptcl,i)[dim]*psiM_actual(j,ptcl);
            std::cerr <<"on-fix af "<<ii<<" "<<jj<<" "<<i<<" "<<j<<" "<<ptcl<<" "<<dotProducts_grad(ii,jj)<<" "<<gs_inverse_update(ptcl,i)[dim]<<" "<<psiM_actual(j,ptcl)<< std::endl;
          }
          std::cerr <<"pre-fix "<<ii<<" "<<jj<<" "<<i<<" "<<j<<dotProducts_grad(ii,jj)<< std::endl;
          dotProducts_grad(ii,jj)[dim]-=gs_inverse_update(WorkingIndex,i)[dim]*psiM_actual(j,WorkingIndex);
          dotProducts_grad(ii,jj)[dim]+=gs_inverse_update(WorkingIndex,i)[dim]*dpsiM_actual(WorkingIndex,j);
          std::cerr <<"post-fix "<<dotProducts_grad(ii,jj)<< std::endl;
        }
    }
  }


  void CalcSingleExcitations_lapl(ValueVector_t &coefs, GradType &val,int &coefIndex)
  {
    for (int i=0; i<orbitals_to_replace.size(); i++)
      for (int j=0; j<unoccupied_orbitals_to_use.size(); j++)
      {
        int dim=0;
        {
          val[dim]+=dotProducts_grad(i,j)[dim];// *coefs(coefIndex);
        }
        coefIndex++;
      }
  }


  //coefIndex should be set to start the double excitations
  void CalcDoubleExcitations_lapl(ValueVector_t &coefs,GradType &val,int &coefIndex)
  {
    for (int ii1=0; ii1<orbitals_to_replace.size(); ii1++)
    {
      for (int ii2=ii1+1; ii2<orbitals_to_replace.size(); ii2++)
      {
        for (int jj1=0; jj1<unoccupied_orbitals_to_use.size(); jj1++)
        {
          for (int jj2=jj1+1; jj2<unoccupied_orbitals_to_use.size(); jj2++)
          {
            int dim=0;
            {
              ValueType a=dotProducts_grad(ii1,jj1)[dim];
              ValueType d=dotProducts_grad(ii2,jj2)[dim];
              ValueType b=dotProducts_grad(ii1,jj2)[dim];
              ValueType c=dotProducts_grad(ii2,jj1)[dim];
              ValueType det=a*d-b*c;
              val[dim]+=coefs(coefIndex)*det;
            }
            coefIndex++;
          }
        }
      }
    }
  }

  ////////////////////////////////////






  //val should be sent with 1.0 likely upon entering this
  void CalcSingleExcitations(ValueVector_t &coefs,ValueType &val,int &coefIndex)
  {
    for (int i=0; i<orbitals_to_replace.size(); i++)
      for (int j=0; j<unoccupied_orbitals_to_use.size(); j++)
      {
        val+=dotProducts(i,j);// *coefs(coefIndex);
        coefIndex++;
      }
  }

  //assumes inclusive
  void set_excitations(int start_uo,int end_uo,int start_rep,int end_rep,ValueVector_t &coefs)
  {
    unoccupied_orbitals_to_use.resize(end_uo-start_uo+1);
    for (int i=start_uo; i<=end_uo; i++)
    {
      unoccupied_orbitals_to_use[i-start_uo]=i;
    }
    orbitals_to_replace.resize(end_rep-start_rep+1);
    for (int i=start_rep; i<=end_rep; i++)
    {
      orbitals_to_replace[i-start_rep]=i;
    }
    dotProducts.resize(orbitals_to_replace.size(),unoccupied_orbitals_to_use.size());
    dotProducts_grad.resize(orbitals_to_replace.size(),unoccupied_orbitals_to_use.size());
    int m=orbitals_to_replace.size();
    int n=unoccupied_orbitals_to_use.size();
    coefs.resize(m*n+m*(m-1)*n*(n-1)/4);
    //for now let's set al teh ocefs to 1;
    coefs=1.0;
  }


  //coefIndex should be set to start the double excitations
  void CalcDoubleExcitations(ValueVector_t &coefs,ValueType &val,int &coefIndex)
  {
    for (int ii1=0; ii1<orbitals_to_replace.size(); ii1++)
    {
      for (int ii2=ii1+1; ii2<orbitals_to_replace.size(); ii2++)
      {
        for (int jj1=0; jj1<unoccupied_orbitals_to_use.size(); jj1++)
        {
          for (int jj2=jj1+1; jj2<unoccupied_orbitals_to_use.size(); jj2++)
          {
            ValueType a=dotProducts(ii1,jj1);
            ValueType d=dotProducts(ii2,jj2);
            ValueType b=dotProducts(ii1,jj2);
            ValueType c=dotProducts(ii2,jj1);
            ValueType det=a*d-b*c;
            val+=coefs(coefIndex)*det;
            coefIndex++;
          }
        }
      }
    }
  }

//     //gs_inverse needs to be orbital x ptcl
//     //psiM_actual needs to be orbtial x ptcl
//     //val should be typically sent in as 1.0
//     void BuildDotProducts(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual,ValueVector_t &coefs,ValueType &val)
//     {

//       for (int i=0;i<orbitals_to_replace.size();i++)
// 	for (int j=0;j<unoccupied_orbitals_to_use.size();j++){
// 	  dotProducts(i,j)=0.0;
// 	  for (int ptcl=0;ptcl<NumPtcls;ptcl++)
// 	    dotProducts(i,j)+=gs_inverse(i,ptcl)*psiM_actual(j,ptcl);

// 	}

//     }


};

}
