#ifndef OHMMS_PSEUDOGEN_HARTREEFOCK_H
#define OHMMS_PSEUDOGEN_HARTREEFOCK_H
#include "Utilities/SimpleParser.h"
#include "AtomicHF/HFConfiguration.h"
#include "AtomicHF/HFAtomicOrbitals.h"
#include "AtomicHF/RadialPotentialSet.h"
#include "Numerics/PPTransform.h"
#include "Numerics/Numerov.h"
#include "Numerics/RadialFunctorUtility.h"
#include "Numerics/Transform2GridFunctor.h"
#include "OhmmsData/libxmldefs.h"
#include "Optimize/Minimize.h"
#include "Numerics/HDFNumericAttrib.h"
#include <libxml++/libxml++.h>
#include <fstream>
#include <iostream>

namespace ohmmshf { 

  struct PseudoHF: public MinimizeFunction {
    
    typedef HFAtomicOrbitals::value_type value_type;
    typedef HFAtomicOrbitals::RadialOrbital_t RadialOrbital_t;
    typedef HFAtomicOrbitals::RadialGrid_t RadialGrid_t;
    
    int maxiter;
    value_type eig_tol;
    value_type ratio;
    value_type scf_tol;
    value_type cg_tol, cg_stepsize, cg_epsilon;
    RadialPotentialSet& Pot;
    HFAtomicOrbitals& Psi;
    vector<scalar> OptParams;
    vector<string> IDtag;
    vector<value_type> AEepsilon;
    vector<value_type> PPepsilon;

    vector<RadialOrbital_t> AEorbitals;
    vector<RadialOrbital_t> AEorbitals_norm;
    VarRegistry<value_type> vreg;
    int norb;
    int Counter;
    int pflag;

    value_type Zeff;
    value_type rmatch;
    value_type weight_eig, weight_norm;
    scalar& Params(int i) { return OptParams[i]; }
    scalar Params(int i) const { return OptParams[i]; }
    int NumParams() { return OptParams.size(); }
    void WriteStuff() { cout << "calling WriteStuff" << endl; }
    
    PseudoHF(RadialPotentialSet& pot, HFAtomicOrbitals& psi,
	     const xmlpp::Node* root):
      Pot(pot), Psi(psi), cg_stepsize(0.01), cg_tol(1.0e-8), 
      rmatch(2.0), cg_epsilon(1.0e-6), eig_tol(1e-12), ratio(0.25), 
      norb(4), maxiter(1000), scf_tol(1e-5), Counter(0),
      weight_eig(0.5), weight_norm(0.5), Zeff(4.0), pflag(0) {
  
      for(int aeorbs=0; aeorbs<norb; aeorbs++){
	AEorbitals.push_back(RadialOrbital_t(Psi.m_grid));
	AEorbitals_norm.push_back(RadialOrbital_t(Psi.m_grid));
      }
      AEepsilon.resize(norb);
      PPepsilon.resize(norb);
      IDtag.resize(2);
 
      AEepsilon[0] = 0.0;
      AEepsilon[1] = 0.0;
      AEepsilon[2] = 0.0;
      AEepsilon[3] = 0.0;

      IDtag[0] = "SJ_lambda";
      IDtag[1] = "r_core";
      
      string afilename = "Ge.AE.h5";
      string species = "Ge";
      hid_t afile = H5Fopen(afilename.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      XMLReport("Opening file: " << afilename)
	hid_t group_id = H5Gopen(afile,species.c_str());
      if(Psi.Restriction == "spin+space") {
	int index = 0;
	Counter = 6;
	for(int i=0; i<2; i++){
	  char grpname[128];	
	  sprintf(grpname,"orbital%04d",Counter++);
	  cout << "Opening group " << grpname << endl;
	  hid_t group_id_orb = H5Gopen(group_id,grpname);
	  Vector<int> NLMS;
	  HDFAttribIO<Vector<int> > QuantumNoHDFIn(NLMS);
	  QuantumNoHDFIn.read(group_id_orb,"nlms");
	  XMLReport("A sherical orbital (n,l,m) " << NLMS[0] << " " << NLMS[1] 
		    << " " << NLMS[2])  
      
	    Vector<double> u_r_temp;
	  RadialOrbital_t temp(Psi(0));
	  HDFAttribIO<Vector<double> > U_RIn(u_r_temp);
	  U_RIn.read(group_id_orb,"u_r");

	  for(int j=0; j<AEorbitals[index].size(); j++){
	    AEorbitals[index](j) = u_r_temp[j];
	    AEorbitals[index+1](j) = u_r_temp[j];
	    temp(j) = u_r_temp[j]*u_r_temp[j];
	  }
	  integrate_RK2_forward(temp,AEorbitals_norm[index]);
	  integrate_RK2_forward(temp,AEorbitals_norm[index+1]);
	  index+=2;
	  H5Gclose(group_id_orb);   
	} 
      }
      if(Psi.Restriction == "none") {
	Counter = 28;
	for(int i=0; i<4; i++){
	  char grpname[128];	
	  sprintf(grpname,"orbital%04d",Counter++);
	  cout << "Opening group " << grpname << endl;
	  hid_t group_id_orb = H5Gopen(group_id,grpname);
	  Vector<int> NLMS;
	  HDFAttribIO<Vector<int> > QuantumNoHDFIn(NLMS);
	  QuantumNoHDFIn.read(group_id_orb,"nlms");
	  XMLReport("A sherical orbital (n,l,m) " << NLMS[0] << " " << NLMS[1] 
		    << " " << NLMS[2])  
      
	    Vector<double> u_r_temp;
	  RadialOrbital_t temp(Psi(0));
	  HDFAttribIO<Vector<double> > U_RIn(u_r_temp);
	  U_RIn.read(group_id_orb,"u_r");

	  for(int j=0; j<AEorbitals[i].size(); j++){
	    AEorbitals[i](j) = u_r_temp[j];
	    temp(j) = u_r_temp[j]*u_r_temp[j];
	  }
	  integrate_RK2_forward(temp,AEorbitals_norm[i]);
	  H5Gclose(group_id_orb);   
	}
      }
      H5Gclose(group_id);
      H5Fclose(afile);
      

      XMLReport("Closing external file")

	if(put(root))
	  XMLReport("Sucessfully parsed the XML file.")
	else ERRORMSG("XML file not sucessfully parsed.")
	       }
    
    bool run() {
      putOptParams();
      ConjugateGradient CG;
      CG.Tolerance = cg_tol;
      CG.StepSize = cg_stepsize;
      CG.epsilon = cg_epsilon;
      CG.Minimize(*this);
      if(pflag){
	plot_ascii();
	plot_mygrid();
	plot_siesta_grid();
      }
      return true;
    }

    bool putOptParams(){
      
      for(int i=0; i<IDtag.size(); i++){
	int id = vreg.find(IDtag[i]);
	if(id>= 0) {
	  int size = vreg.Sizes[id];
	  value_type* temp = vreg.Pointers[id];
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

    bool getOptParams(){
      
      XMLReport("Updated variables")
	copy(OptParams.begin(),OptParams.end(),ostream_iterator<value_type>(cout,"\n"));
      
      int offset = 0;
      for(int i=0; i<IDtag.size(); i++){
	int id = vreg.find(IDtag[i]);
	if(id>= 0) {
	  int size = vreg.Sizes[id];
	  value_type* temp = vreg.Pointers[id];
	  for(int j=0; j<size; j++,temp++){
	    *(temp) = OptParams[j+offset];
	  }
	  offset += size;
	} 
      }
      
      return true;
    }
    


    bool put(const xmlpp::Node* q){
      using namespace xmlpp;
      value_type rc = 1.0;
      value_type lambda = 1.0;
      NodeSet mset = q->find("./PseudoSolve");
      if(mset.empty()){
	ERRORMSG("No PseudoPotential generation information preset.");
	return false;
      } else {
	NodeSet pset = mset[0]->find("./Parameter");
	for(int i=0; i<pset.size(); i++){
	  Element* cur = dynamic_cast<Element*>(pset[i]);
	  Attribute* att = cur->get_attribute("name");
	  if(att) {
	    const string& aname = att->get_value();
	    xmlNode* ccur = cur->cobj();
	    if(aname == "cg_tol"){
	      putContent(cg_tol,ccur);
	      XMLReport("conjugate gradient tolerance = " << cg_tol)
		} else if(aname == "eig_tol"){
		  putContent(eig_tol,ccur);
		  XMLReport("eigen value tolerance = " << eig_tol)
		    } else if(aname == "en_tol"){
		      putContent(scf_tol,ccur);
		      XMLReport("scf tolerance = " << scf_tol)
			} else if(aname == "mix_ratio"){
			  putContent(ratio,ccur); 
			  XMLReport("mixing ratio = " << ratio)
			    } else if(aname == "r_match"){
			      putContent(rmatch,ccur);
			      XMLReport("matching radius = " << rmatch)
				} else if(aname == "weight_norm"){
				  putContent(weight_norm,ccur);
				  XMLReport("weight for partial norms = " << weight_norm)
				    } else if(aname == "weight_eig"){
				      putContent(weight_eig,ccur);
				      XMLReport("weight for eigenvalues = " << weight_eig)
					} else if(aname == "cg_stepsize") {
					  putContent(cg_stepsize,ccur);
					  XMLReport("cg step size = " << cg_stepsize)
					    } else if(aname == "r_core") {
					      putContent(rc,ccur);
					      XMLReport("core radius = " << rc)
						} else if(aname == "lambda") {
						  putContent(lambda,ccur);
						  XMLReport("lambda = " << lambda)
						    } else if(aname == "eigenvalue_400-1") {
						      putContent(AEepsilon[0],ccur);
						      XMLReport("Eigen Value for 400-1 = " << AEepsilon[0])
							} else if(aname == "eigenvalue_4001") {
							  putContent(AEepsilon[1],ccur);
							  XMLReport("Eigen Value for 4001 = " << AEepsilon[1])
							    } else if(aname == "eigenvalue_4111") {
							      putContent(AEepsilon[2],ccur);
							      XMLReport("Eigen Value for 4111 = " << AEepsilon[2])
								} else if(aname == "eigenvalue_4101") {
								  putContent(AEepsilon[3],ccur);
								  XMLReport("Eigen Value for 4101 = " << AEepsilon[3])
								    } else if(aname == "plot") {
								      putContent(pflag,ccur);
								      XMLReport("Ploting flag is set to = " << pflag)
								    }
	  }
	}
      }
      Pot.add(new PseudoPotential(vreg,Zeff,rc,lambda));
      XMLReport("Adding SJ-PseudoPotential.")
	
	return true; 
    }


    void plot_ascii(){
      char* fname = "Ge.pp.ASCII";
      ofstream gplot(fname);
      cout << "Writing Pseudopotential to file " << fname << endl;     
 
      for(int i=0; i<Psi.m_grid->size(); i++){
      	value_type r =  Psi.m_grid->r(i);
	value_type SJ_num = 1.0-exp(-Params(0)*r);
	value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
	gplot << setw(15) << r << setw(15) <<  (-1.0*Zeff/r)*(SJ_num/SJ_den) << endl;
      }
      gplot.close();
    }
    
    void plot_siesta_grid(){
      vector<string> vlist;
      int npts;   
      char* infname = "Ge.tm2.psf";
      ifstream fin(infname,ios_base::in);
      if(!fin){
	cout << "Could not open file " << infname << endl;
	exit(-1);
      }
      
      cout << "Reading siesta pseudopotential file " << infname 
	   << " for grid information." << endl;



      getwords(vlist,fin);
      getwords(vlist,fin);
      getwords(vlist,fin);
      getwords(vlist,fin);
      npts = atoi(vlist[2].c_str());
      getwords(vlist,fin);
      
      value_type rin;
      vector<value_type> grid_temp;
      for(int ig=0; ig<npts; ig++) {
	fin >> rin;
	grid_temp.push_back(rin);
      }   
      RadialGrid_t* red_grid = new NumericalGrid<value_type>(grid_temp);
      
      cout << endl;
  
      RadialOrbital_t red_vc(red_grid);
    
      RadialOrbital_t VCharge(Psi(0));
      for(int ob=0; ob<norb; ob++){
	for(int j=0; j<VCharge.size(); j++){
	  VCharge(j)+= Psi(ob,j)*Psi(ob,j);
	  // VCharge(j)+= AEorbitals[ob](j)*AEorbitals[ob](j);
	}
      }

      int imin = 0;
      value_type deriv = (VCharge(imin+1)-VCharge(imin))/VCharge.dr(imin);
      VCharge.spline(imin,deriv,VCharge.size()-1,0.0);
      
      Transform2GridFunctor<OneDimGridFunctor<value_type>,OneDimGridFunctor<value_type> > transform(VCharge,red_vc);
      
      transform.generate();   
      cout << endl; 
      cout << "Output on a reduced grid with rmin = " << red_grid->rmin() 
	   << " and rmax = " << red_grid->rmax() << " and size = "
	   << npts << endl;
      cout << "Total valence charge = " << integrate_RK2(VCharge) << endl;
      cout << "Total valence charge = " << integrate_RK2(red_vc) << endl;
      RadialOrbital_t PP(red_grid);
    
      value_type prefactor = -2.0*Zeff;
      for(int i=0; i<npts; i++){
      	value_type r =  red_grid->r(i);
	value_type SJ_num = 1.0-exp(-Params(0)*r);
	value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
	PP(i) = prefactor*(SJ_num/SJ_den);
      }
      
      char* fname = "Ge.siesta_grid.psf";
      ofstream siesta(fname);
      cout << "Writing Pseudopential to file " << fname << endl;

      siesta << " Ge  ca nrl nc" << endl;
      siesta << " ATM3       no_date   Troullier-Martins" << endl;
      siesta << " 4s 2.00  r=" << setw(3) << rmatch
	     << "/4p 2.00  r=" << setw(3) << rmatch 
	    << "/3d  2.00 r=" << setw(3) << rmatch << endl;
      siesta << "  3  0 " << npts << endl;
      siesta << " Radial grid follows" << endl;
      int nlines = npts/4;
      int remainder = npts-nlines*4;
      siesta.precision(12);
      siesta.setf(ios::scientific,ios::floatfield);
      siesta << setiosflags(ios::uppercase);
      int j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << red_grid->r(j)
	       << setw(20) << red_grid->r(j+1) 
	       << setw(20) << red_grid->r(j+2)
	       << setw(20) << red_grid->r(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
      for(int i=j; i<npts; i++)
	siesta << setw(20) << red_grid->r(i);
      siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  0" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  1" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
      for(int i=j; i<npts; i++)
	siesta << setw(20) << PP(i);
      siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  2" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  3" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      double CC = 0.0;
      siesta << " Core charge follows" << endl;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << CC
	       << setw(20) << CC 
	       << setw(20) << CC
	       << setw(20) << CC 
	       << endl;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << CC;
	siesta << endl;
      }
      siesta << " Valence charge follows" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << red_vc(j)
	       << setw(20) << red_vc(j+1)
	       << setw(20) << red_vc(j+2)
	       << setw(20) << red_vc(j+3)
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << red_vc(i);
	siesta << endl;
      }
      siesta.close();

 
    }
    
  void plot_mygrid(){
    
    char* fname = "Ge.mygrid.psf";
    cout << "Writing Pseudopotential to file " << fname << endl;
    
    RadialOrbital_t VCharge(Psi(0));
    
    for(int ob=0; ob<norb; ob++){
      for(int j=0; j<VCharge.size(); j++){
	VCharge(j)+= Psi(ob,j)*Psi(ob,j);
      }
    }
    
    int imin = 0;
    value_type deriv = (VCharge(imin+1)-VCharge(imin))/VCharge.dr(imin);
    VCharge.spline(imin,deriv,VCharge.size()-1,0.0);
    cout << endl;
    cout << "Output on a grid with rmin = " << Psi.m_grid->rmin() 
	 << " and rmax = " << Psi.m_grid->rmax() << " and size = "
	 << Psi.m_grid->size() << endl;
    cout << "Total valence charge = " << integrate_RK2(VCharge) << endl;
    RadialOrbital_t PP(Psi(0));
    
    int npts = Psi.m_grid->size();
    
    value_type prefactor = -2.0*Zeff;
    for(int i=0; i<npts; i++){
      value_type r =  Psi.m_grid->r(i);
      value_type SJ_num = 1.0-exp(-Params(0)*r);
      value_type SJ_den = 1.0+exp(-Params(0)*(r-Params(1)));
      PP(i) = prefactor*(SJ_num/SJ_den);
    }

    ofstream siesta(fname);
    siesta << " Ge  ca nrl nc" << endl;
    siesta << " ATM3       no_date   Troullier-Martins" << endl;
    siesta << " 4s 2.00  r=" << setw(3) << rmatch
	   << "/4p 2.00  r=" << setw(3) << rmatch 
	   << "/3d  2.00 r=" << setw(3) << rmatch << endl;
    siesta << "  3  0 " << npts << endl;
    siesta << " Radial grid follows" << endl;
    int nlines = npts/4;
    int remainder = npts-nlines*4;
    siesta.precision(12);
    siesta.setf(ios::scientific,ios::floatfield);
    siesta << setiosflags(ios::uppercase);
    int j=0;
    for(int i=0; i<nlines; i++){
      siesta << setw(20) << Psi.m_grid->r(j)
	     << setw(20) << Psi.m_grid->r(j+1) 
	     << setw(20) << Psi.m_grid->r(j+2)
	     << setw(20) << Psi.m_grid->r(j+3) 
	     << endl;
	j+=4;
    }
    if(remainder){
      for(int i=j; i<npts; i++)
	siesta << setw(20) << Psi.m_grid->r(i);
      siesta << endl;
    }
    siesta << " Down Pseudopotential follows (l on next line)" << endl;
    siesta << "  0" << endl;
    j=0;
    for(int i=0; i<nlines; i++){
      siesta << setw(20) << PP(j)
	     << setw(20) << PP(j+1) 
	     << setw(20) << PP(j+2)
	     << setw(20) << PP(j+3) 
	     << endl;
      j+=4;
    }
    if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  1" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  2" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      siesta << " Down Pseudopotential follows (l on next line)" << endl;
      siesta << "  3" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << PP(j)
	       << setw(20) << PP(j+1) 
	       << setw(20) << PP(j+2)
	       << setw(20) << PP(j+3) 
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << PP(i);
	siesta << endl;
      }
      double CC = 0.0;
      siesta << " Core charge follows" << endl;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << CC
	       << setw(20) << CC 
	       << setw(20) << CC
	       << setw(20) << CC 
	       << endl;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << CC;
	siesta << endl;
      }
      siesta << " Valence charge follows" << endl;
      j=0;
      for(int i=0; i<nlines; i++){
	siesta << setw(20) << VCharge(j)
	       << setw(20) << VCharge(j+1)
	       << setw(20) << VCharge(j+2)
	       << setw(20) << VCharge(j+3)
	       << endl;
	j+=4;
      }
      if(remainder){
	for(int i=j; i<npts; i++)
	  siesta << setw(20) << VCharge(i);
	siesta << endl;
      }
      siesta.close();
 
    }


    value_type Cost(){
      typedef PPLogTransform<RadialOrbital_t> Transform_t;
      typedef Numerov<Transform_t, RadialOrbital_t> Numerov_t;
      value_type Vtotal, KEnew, KEold, E;
      vector<value_type> energy(Pot.size());
      value_type lowerbound, upperbound;
      Pot.reset(Psi);
      int iter = 0;
      Vtotal = Pot.evaluate(Psi,energy,norb);
      Pot.mix(0.0);
      KEnew = Pot.calcKE(Psi,0,norb);
      string label("spdf");  
      cout.precision(10);
      scalar cost = 0.0;
      getOptParams();
      Psi.EigenBoundParam = -2.0*Zeff/Params(1);
   
      do {
	KEold = KEnew;
	value_type eigsum = 0.0;
	for(int ob=0; ob < norb; ob++){ 
	  Transform_t es(Pot.V[ob], Psi.N[ob], Psi.L[ob], 
			 Psi.CuspParam, Psi.EigenBoundParam);
	  Numerov_t numerov(es,Psi(ob));
	  es.EnergyBound(lowerbound,upperbound);
	  eigsum += (PPepsilon[ob] = 
		     numerov.solve(lowerbound, upperbound, eig_tol));
	}


	Psi.normalize(norb);
	Psi.applyRestriction(norb);
	KEnew = Pot.calcKE(Psi,eigsum,norb);
	E = KEnew + Vtotal;
	Vtotal = Pot.evaluate(Psi,energy,norb);	
	Pot.applyRestriction(Psi);
	Pot.mix(ratio);
	iter++;
      }while(fabs(KEnew-KEold)>scf_tol && iter<maxiter);

      cout << "Iteration = " << iter << " and Energy = " << E << endl;
      value_type sum_norm = 0.0;
      value_type sum_eig = 0.0;
      for(int ob=0; ob < norb; ob++){
	cout << Psi.N[ob]<< label[Psi.L[ob]] << '\t' << PPepsilon[ob] << endl;
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

	//	cout << "Log derivatives = " << rpn*AEorbitals[ob].dY/AEorbitals[ob].Y - 1.0 << " " << rpn*Psi(ob).dY/Psi(ob).Y - 1.0 <<  " " << endl;

	int x = Psi.m_grid->index(rmatch);
	sum_norm += fabs(1.0-psi_norm(x)/AEorbitals_norm[ob](x));
	sum_eig += fabs(1.0-PPepsilon[ob]/AEepsilon[ob]);
      }
      cost = (sum_norm*weight_norm+sum_eig*weight_eig)/static_cast<value_type>(norb);
      cout << endl;
      cout << "Differential in eigenvalues:   " << sum_eig << endl;
      cout << "Differential in partial norms: " << sum_norm << endl;
      cout << "Cost = " << cost << endl;

      if(cost > 100.0) return 100;

      if (Params(0) < 0.0) return 100;
      else if (Params(1) > rmatch) return 100;
      else return cost;
  
    }
    
  };
}

#endif
