
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/SPOSetBase.h"



namespace qmcplusplus {

  class ExcitationClass {

    typedef SPOSetBase::ValueMatrix_t ValueMatrix_t;
    typedef SPOSetBase::ValueVector_t ValueVector_t;
    typedef double ValueType ;

  public:
     vector<int> unoccupied_orbitals_to_use;
    vector<int> orbitals_to_replace;
    ValueMatrix_t dotProducts; // orbitals_to_replace.size() x unoccupied_orbitals_to_use.size()
    void testMe()
    {
      int i=5;
    }

    //gs_inverse needs to be orbital x ptcl
    //psiM_actual needs to be orbtial x ptcl
    void BuildDotProducts(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual)
    {
      int NumPtcls=psiM_actual.extent(1);
      for (int ii=0;ii<orbitals_to_replace.size();ii++)
	for (int jj=0;jj<unoccupied_orbitals_to_use.size();jj++){
	  int i=orbitals_to_replace[ii];
	  int j=unoccupied_orbitals_to_use[jj];
	  dotProducts(i,j)=0.0;
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++)
	    dotProducts(i,j)+=gs_inverse(i,ptcl)*psiM_actual(j,ptcl);
	}
    }
    


//     //gs_inverse needs to be orbital x ptcl
//     //psiM_actual needs to be orbtial x ptcl
//     void UpdateDotProducts(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual,
// 			   ValueVector_t &gs_inverse_new_WorkingIndex,ValueVector_t &psiM_actual_new_WorkingIndex,
// 			   int WorkingIndex)
//     {
//       for (int ii=0;ii<orbitals_to_replace.size();ii++)
// 	for (int jj=0;jj<unoccupied_orbitals_to_use.size();jj++){
// 	  int i=orbitals_to_replace[ii];
// 	  int j=unoccupied_orbitals_to_use[jj];
// 	  dotProducts(i,j)-=gs_inverse(i,WorkingIndex)*psiM_actual(j,WorkingIndex);	
// 	  dotProducts(i,j)+=gs_inverse_new_WorkingIndex(i)*psiM_actual_new_WorkingIndex(j);
// 	}
//     }

    //val should be sent with 1.0 likely upon entering this
    void CalcSingleExcitations(ValueVector_t &coefs,ValueType &val,int &coefIndex)
    {
      for (int i=0;i<orbitals_to_replace.size();i++)
	for (int j=0;j<unoccupied_orbitals_to_use.size();j++){
	  val+=coefs(coefIndex)*dotProducts(i,j);
	  coefIndex++;
	}
	  


    }

    //assumes inclusive
    void set_excitations(int start_uo,int end_uo,int start_rep,int end_rep,ValueVector_t &coefs)
    {
      unoccupied_orbitals_to_use.resize(end_uo-start_uo+1);
      for (int i=start_uo;i<=end_uo;i++)
	unoccupied_orbitals_to_use[i-start_uo]=start_uo;
      orbitals_to_replace.resize(end_rep-start_rep+1);
      for (int i=start_rep;i<=end_rep;i++)
	orbitals_to_replace[i-start_rep]=start_rep;
      dotProducts.resize(orbitals_to_replace.size(),unoccupied_orbitals_to_use.size());
      int m=orbitals_to_replace.size();
      int n=unoccupied_orbitals_to_use.size();
      coefs.resize(m*(m+1)*n*(n+1)/4);
      //for now let's set al teh ocefs to 1;
      coefs=1.0;
    }
    

    //coefIndex should be set to start the double excitations
    void CalcDoubleExcitations(ValueVector_t &coefs,ValueType &val,int &coefIndex)
    {
      for (int ii1=0;ii1<orbitals_to_replace.size();ii1++){
	for (int ii2=ii1+1;ii2<orbitals_to_replace.size();ii2++){
	  for (int jj1=0;jj1<unoccupied_orbitals_to_use.size();jj1++){
	    for (int jj2=jj1+1;jj2<unoccupied_orbitals_to_use.size();jj2++){
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
