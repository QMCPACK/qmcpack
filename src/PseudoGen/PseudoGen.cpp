//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include "SQD/HFConfiguration.h"
#include "SQD/SphericalPotential/RegularLinearTransform.h"
#include "SQD/SphericalPotential/RegularLogTransform.h"
#include "SQD/SphericalPotential/NuclearLinearTransform.h"
#include "SQD/SphericalPotential/NuclearLogTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "PseudoGen/PseudoGen.h"
#include "Optimize/Minimize.h"

namespace ohmmshf {

  /**
     @param pot the potential
     @param psi the wavefunction
     @brief The contructor.
  */

  PseudoGen::PseudoGen(RadialPotentialSet& pot, 
		       SphericalOrbitalTraits::BasisSetType& psi):
    Pot(pot), Psi(psi), num_closed_shells(0), maxiter(1000), 
    eig_tol(1e-12), scf_tol(1e-8), ratio(0.35), 
    GridType("none"), rmatch(1.0), cg_stepsize(0.01), 
    cg_tolerance(1.0e-12), cg_epsilon(1.0e-6), weight_eig(0.5), 
    weight_norm(0.5), grid_ptr(NULL), orb_ptr(NULL), 
    pot_ptr(NULL), opt_ptr(NULL), myGrid(NULL), NumCostCalls(0){ 

    IDtag.resize(2);
    IDtag[0] = "SJ_lambda";
    IDtag[1] = "r_core";

  }

  /**
     @param aroot the root for all output files
     @brief Sets the root for all output files.
  */

  void PseudoGen::setRoot(const string& aroot) {
    RootFileName = aroot;
    LogFileName = RootFileName + ".log";
    log_stream.open(LogFileName.c_str());
  }


  bool PseudoGen::optimize() {
    putOptParams();
    ConjugateGradient CG;
    CG.Tolerance = cg_tolerance;
    CG.StepSize = cg_stepsize;
    CG.epsilon = cg_epsilon;
    CG.Minimize(*this);
    plot_ascii();
    plot_siesta_grid();
    Cost();
    plot_ascii();
    plot_siesta_grid();
    return true;
  }

 bool PseudoGen::run() {
   putOptParams();
   Cost();
   return true;
 }

  bool PseudoGen::putOptParams(){
    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++){
      //locate the id in the variable list for Psi
      int id = vreg.find(IDtag[i]);
      if(id>= 0) {
	//find the number of variables represented by the id
	int size = vreg.Sizes[id];
	value_type* temp = vreg.Pointers[id];
	//loop over the number of variables for the id
	//assign the value to the optimization parameter set
	for(int j=0; j<size; j++){
	  OptParams.push_back(temp[j]);
	}
      } else {
	ERRORMSG("Could not find parameter " << IDtag[i])
	  return false;
	OptParams.push_back(0.0);
      }
    }
    XMLReport("Initial variables")
      copy(OptParams.begin(),OptParams.end(),ostream_iterator<value_type>(cout,"\n"));
  
    return true;
  }

  bool PseudoGen::getOptParams(){
      
    XMLReport("Updated variables")
      copy(OptParams.begin(),OptParams.end(),ostream_iterator<value_type>(cout,"\n"));
      
    int offset = 0;
    //loop over all the unique id's
    for(int i=0; i<IDtag.size(); i++){
      //locate the id in the variable list for Psi
      int id = vreg.find(IDtag[i]);
      if(id>= 0) {
	//find the number of variables represented by the id
	int size = vreg.Sizes[id];
	//loop over the number of variables for the id
	//assign the variable its new value
	value_type* temp = vreg.Pointers[id];
	for(int j=0; j<size; j++,temp++){
	  *(temp) = OptParams[j+offset];
	}
	//increment offset to account for the fact that id
	//may represent more than one variable
	offset += size;
      } 
    }
      
    return true;
  }
    
  double PseudoGen::Cost(){

    getOptParams();

    double TotE = 0.0;
    int norb = Psi.size();
    if(GridType == "linear"){
      RegularLinearTransform<RadialOrbital_t> *afake=NULL;
      TotE = runHF(afake,norb);
    } else if(GridType == "log"){
      RegularLogTransform<RadialOrbital_t> *afake=NULL;
      TotE = runHF(afake,norb);
    }

   
    double cost = 0.0;
    string label("spdf"); 
    //  std::ofstream log_stream(LogFileName.c_str());
    log_stream.precision(8); 

    log_stream << "Iteration = " << NumCostCalls++ 
	       << " and Total Energy = " << TotE << endl;

    value_type sum_norm = 0.0;
    value_type sum_eig = 0.0;
    log_stream << "orb" << '\t' << "PPeigVal" 
	       << '\t' << "AEeigVal" << endl;
    for(int ob=0; ob < norb; ob++){
      log_stream << Psi.N[ob]<< label[Psi.L[ob]] << '\t' 
		 << PPeigVal[ob] << '\t' << AEeigVal[ob] 
		 << endl;
      RadialOrbital_t psi_norm(Psi(ob));
      RadialOrbital_t psi_sq(Psi(ob));
      for(int j=0; j<Psi(ob).size(); j++)
	psi_sq(j) = Psi(ob,j)*Psi(ob,j);
      integrate_RK2_forward(psi_sq,psi_norm);

      //	AEorbitals[ob].setgrid(rpn);
      //	AEorbitals_norm[ob].setgrid(rpn);

      //	psi_norm.setgrid(rpn);
      //	value_type rpninv = 1.0/rpn;
      //	AEorbitals[ob].evaluate(rpn,rpninv);
      //	AEorbitals_norm[ob].evaluate(rpn,rpninv);
      //	Psi(ob).evaluate(rpn,rpninv);
      //	psi_norm.evaluate(rpn,rpninv);

      //	cout << "Log derivatives = "
      //             << rpn*AEorbitals[ob].dY/AEorbitals[ob].Y - 1.0 
      // << " " << rpn*Psi(ob).dY/Psi(ob).Y - 1.0 <<  " " << endl;

      int x = Psi.m_grid->index(rmatch);
      sum_norm += fabs(1.0-psi_norm(x)/AEorbitals_norm[ob](x));
      sum_eig += fabs(1.0-PPeigVal[ob]/AEeigVal[ob]);
    }
    sum_norm /= static_cast<value_type>(norb);
    sum_eig /= static_cast<value_type>(norb);
    cost = (sum_norm*weight_norm+sum_eig*weight_eig);

    log_stream << "Differential in eigenvalues:   " << sum_eig << endl;
    log_stream << "Differential in partial norms: " << sum_norm << endl;
    log_stream << "Cost = " << cost << endl;
    log_stream << endl;

    if(cost > 100.0) return 100;
    
    if (Params(0) < 0.0) return 100;
    else if (Params(1) > rmatch) return 100;
    else return cost;
      
  }
    
  /**
     @param fake Transformation object
     @param norb the number of eigen vectors to be obtained
     @brief Perform self-consistent Hartree-Fock calculations.
     *
     *The first argument is used to tell the compiler which transform
     *object is used but is not meaningful for the calculations. 
     *This is to force the compiler to inline everything possible. 
     *A better solution could be implemented using traits.
     */
  template<class Transform_t> 
  inline double
  PseudoGen::runHF(Transform_t* fake, int norb) {
    typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
    value_type Vtotal,KEnew, KEold,E; 
    value_type lowerbound, upperbound;
    vector<value_type> energy(Pot.size());
    Pot.reset(Psi);
    int iter = 0;
    Vtotal = Pot.evaluate(Psi,energy,norb);
    Pot.mix(0.0);
    KEnew = Pot.calcKE(Psi,0,norb);

    string label("spdf");  
    //   std::ofstream log_stream(LogFileName.c_str());
    // log_stream.precision(8);
      
    do {
      KEold = KEnew;
      value_type eigsum = 0.0;
      //loop over the orbitals
      for(int ob=0; ob < norb; ob++){ 

        //set the number of nodes of the eigen vector 
        Pot.V[ob].setNumOfNodes(Pot.getNumOfNodes(Psi.N[ob],Psi.L[ob]));

	//calculate the lower and upper bounds for the eigenvalue
        Pot.EnergyBound(lowerbound,upperbound);

	//set up the transformer
	Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob],Psi.CuspParam, 
		       Pot.getMass());

	//initialize the numerov solver
	Numerov_t numerov(es,Psi(ob));

	//calculate the eigenvalue and the corresponding orbital
	eigsum += (PPeigVal[ob] = 
		   numerov.solve(lowerbound, upperbound, eig_tol));

	//	LOGMSG(Psi.N[ob]<< label[Psi.L[ob]] << '\t' << PPeigVal[ob]);
      }
	
      //normalize the orbitals
      Psi.normalize(norb);
      //restrict the orbitals
      Psi.applyRestriction(norb);
      //calculate the new kinetic energy	
      KEnew = Pot.calcKE(Psi,eigsum,norb);
      //the total energy
      E = KEnew + Vtotal;
      //for the new orbitals Psi, calculate the new SCF potentials
      Vtotal = Pot.evaluate(Psi,energy,norb);
      //restrict the potential
      Pot.applyRestriction(Psi);
      //mix the new SCF potential with the old
      Pot.mix(ratio);
      //  cout.precision(10);
      //       cout << "Iteration #" << iter+1 << endl;
      //       cout << "KE    = " << setw(15) << KEnew 
      // 		 << "  PE     = " << setw(15) << Vtotal << endl; 
      //       cout << "PE/KE = " << setw(15) << Vtotal/KEnew 
      // 		 << "  Energy = " << setw(15) << E << endl;
      //       cout << endl;
      iter++;
      //continue the loop until the kinetic energy converges
    }while(fabs(KEnew-KEold)>scf_tol && iter<maxiter);
    log_stream.precision(10);
    log_stream << "Total Hartree-Fock iterations = " << iter << endl;
    log_stream << "KE    = " << setw(15) << KEnew 
	 << "  PE     = " << setw(15) << Vtotal << endl; 
    log_stream << "PE/KE = " << setw(15) << Vtotal/KEnew 
	 << "  Energy = " << setw(15) << E << endl;
    log_stream << endl;
    log_stream << "V_External = " << energy[0] << endl;
    log_stream << "V_Hartree = "  << energy[1] << endl;
    log_stream << "V_Exchange = " << energy[2] << endl;
    log_stream << endl;

    return E;

  }

}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

  
  
