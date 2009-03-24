
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
      TinyVector<int,2> ii;
      TinyVector<int,2> jj;

    TinyVector<int,2> second_replaces_first;
    list<TinyVector<int,2> > replaceList;
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
      if (orbitals_to_replace.size()==0 || unoccupied_orbitals_to_use.size()==0){
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
      if (unoccupied_orbitals_to_use_index==unoccupied_orbitals_to_use.size()){
	unoccupied_orbitals_to_use_index=0;
	orbitals_to_replace_index++;
	if (orbitals_to_replace_index==orbitals_to_replace.size()){
	  second_replaces_first[0]=-1;
	  second_replaces_first[1]=-1;
	  return second_replaces_first;
	}
      }
      second_replaces_first[0]=orbitals_to_replace[orbitals_to_replace_index];
      second_replaces_first[1]=unoccupied_orbitals_to_use[unoccupied_orbitals_to_use_index];
      return second_replaces_first;
    }
    


//     list<TinyVector<int,2> > begin_doubleExcitation()
//     {
//       ii[0]=0;
//       ii[1]=ii[0]+1;
//       jj[0]=0;
//       jj[1]=jj[0]+1;
//       list<TinyVector<int,2> > myList;
//       myList.push_back(ii);
//       myList.push_back(jj)
      
//       for (int ii1=0;ii1<orbitals_to_replace.size();ii1++){
// 	for (int ii2=ii1+1;ii2<orbitals_to_replace.size();ii2++){
// 	  for (int jj1=0;jj1<unoccupied_orbitals_to_use.size();jj1++){
// 	    for (int jj2=jj1+1;jj2<unoccupied_orbitals_to_use.size();jj2++){




      
//     }

//     list<TinyVector<int,2> > next_doubleExcitation()
//     {


//     }

    //gs_inverse needs to be orbital x ptcl
    //psiM_actual needs to be orbtial x ptcl
    void BuildDotProducts(ValueMatrix_t &gs_inverse,ValueMatrix_t &psiM_actual)
    {
      int NumPtcls=psiM_actual.extent(1);

      for (int ii=0;ii<orbitals_to_replace.size();ii++)
	for (int jj=0;jj<unoccupied_orbitals_to_use.size();jj++){
	  int i=orbitals_to_replace[ii];
	  int j=unoccupied_orbitals_to_use[jj];
	  dotProducts(ii,jj)=0.0;
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++)
	    dotProducts(ii,jj)+=gs_inverse(i,ptcl)*psiM_actual(j,ptcl);
	}
    }
    //dpsiM is currently ptcl x orbital which is bad!
    void BuildDotProducts(ValueMatrix_t &gs_inverse,ValueVector_t &op1,ValueVector_t &op2,ValueMatrix_t &dpsiM_actual)
    {
      int NumPtcls=dpsiM_actual.extent(0);
      for (int ii=0;ii<orbitals_to_replace.size();ii++)
	for (int jj=0;jj<unoccupied_orbitals_to_use.size();jj++){
	  int i=orbitals_to_replace[ii];
	  int j=unoccupied_orbitals_to_use[jj];
	  dotProducts(ii,jj)=0.0;
	  for (int ptcl=0;ptcl<NumPtcls;ptcl++)
	    dotProducts(ii,jj)+=(gs_inverse(i,ptcl)+op1(i)*op2(ptcl))*psiM_actual(j,ptcl);
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

	  val+=dotProducts(i,j);// *coefs(coefIndex);
	  coefIndex++;
	}



    }

    //assumes inclusive
    void set_excitations(int start_uo,int end_uo,int start_rep,int end_rep,ValueVector_t &coefs)
    {
      unoccupied_orbitals_to_use.resize(end_uo-start_uo+1);
      for (int i=start_uo;i<=end_uo;i++){
	unoccupied_orbitals_to_use[i-start_uo]=i;

      }
      orbitals_to_replace.resize(end_rep-start_rep+1);
      for (int i=start_rep;i<=end_rep;i++){
	orbitals_to_replace[i-start_rep]=i;

      }
      dotProducts.resize(orbitals_to_replace.size(),unoccupied_orbitals_to_use.size());
      int m=orbitals_to_replace.size();
      int n=unoccupied_orbitals_to_use.size();
      coefs.resize(m*n+m*(m-1)*n*(n-1)/4);
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
