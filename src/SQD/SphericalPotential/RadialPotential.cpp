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
#include <math.h>
#include <fstream>
#include "SQD/SphericalPotential/RadialPotential.h"
#include "Numerics/Clebsch_Gordan.h"
#include "Numerics/RadialFunctorUtility.h"
namespace ohmmshf {

  /**
   *@param cg The Clebsch-Gordan matrix elements
   *@param norb the number of orbitals
   *@brief The Constructor for the HartreePotential
   */
  
  HartreePotential::HartreePotential(Clebsch_Gordan* cg, int norb):
    CG_coeff(cg) { 
    storage.resize(norb*(norb+1)/2); 
    for(int nn=0; nn<storage.size(); nn++) storage[nn] = 0.0;
  }
  
  /**
   *@param psi the wavefunction
   *@param V increment the potential
   *@parma norb the number of orbitals
   *@return The Hartree energy
   \f[
   E_{Hartree} = 
   \frac{1}{2} \sum_{ij} \sum_{k=0}^{\min(2l_i,2l_j)} (-1)^{-m_i-m_j}
   \frac{(2l_j+1)(2l_i+1)}{(2k+1)^2} \\
   \langle l_j m_j l_j (-m_j)| k 0 \rangle \langle l_j 0 l_j 0| k 0 \rangle
   \langle l_i m_i l_i (-m_i)| k 0 \rangle \langle l_i 0 l_i 0| k 0 \rangle
   {\cal R}^k(ij;ij) 
   \f]
   *
   *@brief Calculates the Hartree potential for each orbital and the
   *Hartree energy.  
   *
   The Hartree potential for the \textit{ith} orbital
   \f[
   \hat{V}_{H} u_{n_i l_i}(r) =  
   \sum_j \sum_{k=0}^{\min(2l_i,2l_j)} (-1)^{-m_i-m_j}
   \frac{(2l_j+1)(2l_i+1)}{(2k+1)^2} \langle l_j m_j l_j (-m_j)| k 0 \rangle 
   \langle l_j 0 l_j 0| k 0 \rangle 
   \langle l_i m_i l_i (-m_i)| k 0 \rangle \langle l_i 0 l_i 0| k 0 \rangle 
   \frac{{\cal Y}_k(n_jl_j;n_jl_j/r)}{r} u_{n_i l_i}(r)
   \f]
   *
   see RadialFunctorUtility.h for the implementation of 
   \f$ {\cal Y}_k(n_jl_j;n_jl_j/r) /r \f$ and 
   \f$ {\cal R}^k(ij;ij). \f$  
   *
   *@note The sumnation over the orbitals {ij} includes \f$ i=j. \f$
   */
  
  RadialPotentialBase::value_type
  HartreePotential::evaluate(const BasisSetType& psi, 
			     RadialOrbitalSet_t& V, int norb) {
  
    int kmax, k;
    int pt;
    value_type coeff, ith_orb_coeff, jth_orb_coeff, energy_coeff = 0;
    value_type Ehartree=0;
  
    RadialOrbital_t Ykii_r(psi(0));
    RadialOrbital_t Ykjj_r(psi(0));
    RadialOrbital_t Psisq_x_Yk(psi(0));
    int npts = psi.m_grid->size();
    int nn=0;

    for(int i=0; i < norb; i++) {

      int mi = psi.M[i];
      int li = psi.L[i];
      int two_li_plus_one = 2*li + 1;
    
      for(int j=i; j < norb; j++, nn++) {
      
	int mj = psi.M[j];
	int lj = psi.L[j];
	int two_lj_plus_one = 2*lj + 1;
      
	int kmax = (li > lj) ? 2*lj : 2*li;

	for(int k=kmax; k >= 0; k -= 2) {
	  int two_k_plus_one = 2*k+1;

	  coeff = static_cast<value_type>(two_li_plus_one*two_lj_plus_one)/
	    static_cast<value_type>(two_k_plus_one)/static_cast<value_type>(two_k_plus_one)
	    * CG_coeff->cg(li,li,k,0,0) * CG_coeff->cg(lj,lj,k,0,0)
	    * CG_coeff->cg(li,li,k,mi,-mi) * CG_coeff->cg(lj,lj,k,mj,-mj)
	    * pow(-1.0, mi+mj);

	  if(i == j) coeff /= 2.0; //double counting

	  ith_orb_coeff = psi.Occ[i] * coeff;
	  jth_orb_coeff = psi.Occ[j] * coeff;
	  energy_coeff = psi.Occ[j] * psi.Occ[i] * coeff;

	  Ykofr(Ykii_r, psi(i), psi(i), k);
	  Ykofr(Ykjj_r, psi(j), psi(j), k);
	  for(int gp=0; gp<npts; gp++){
	    V[i](gp) += jth_orb_coeff*Ykjj_r(gp);
	    V[j](gp) += ith_orb_coeff*Ykii_r(gp);
	  }
	
	  storage[nn] = 0.0;
	  storage[nn] += Phisq_x_Yk(Ykjj_r, psi(i), psi(i), 0.5*energy_coeff);
	  storage[nn] += Phisq_x_Yk(Ykii_r, psi(j), psi(j), 0.5*energy_coeff);
	  Ehartree += storage[nn];
	}
      }
    }

    return Ehartree;

  }

  void HartreePotential::getStorage(const BasisSetType& psi,
				    const std::string& RootFileName){
    
    string fileforoutput = RootFileName + ".hartree";
    ofstream fout(fileforoutput.c_str());
    string llabel("spdf");
    string slabel("d0u");  
    int norb = psi.size();
    int nn=0;
    value_type temp;
    fout << "orb#1" << '\t' << "orb#2" << '\t' << "Hartree" << endl;
    fout.precision(12);
    fout.setf(ios::scientific,ios::floatfield);
    for(int i=0; i<norb; i++) {
      for(int j=i; j<norb; j++, nn++) {
	temp = storage[nn];
	if(fabs(temp) < 1e-12) temp = 0.0;
	fout << psi.N[i] << llabel[psi.L[i]] << slabel[psi.S[i]+1]
	     << '\t' << psi.N[j] << llabel[psi.L[j]] << slabel[psi.S[j]+1]
	     << '\t' << temp << endl;
      }
    }
  }

  /**
   *@param cg The Clebsch-Gordan matrix elements
   *@param norb the number of orbitals
   *@brief The Constructor for the ExchangePotential
   */
  
  ExchangePotential::ExchangePotential(Clebsch_Gordan* cg, int norb):
    CG_coeff(cg) {
    storage.resize(norb*(norb+1)/2); 
    for(int nn=0; nn<storage.size(); nn++) storage[nn] = 0.0; 
  }
  
  /**
   *@param psi the wavefunction
   *@param V increment the potential
   *@parma norb the number of orbitals
   *@return The Exchange energy
   \f[
   E_{Exchange} =
   -\sum_j \delta_{\sigma_i,\sigma_j} \sum_{k=|l_j-l_i|}^{l_i+l_j}
   \frac{(2l_i+1)(2l_j+1)}{(2k+1)^2}
   \langle l_i (-m_i) l_j m_j| k (m_i-m_j) \rangle^2 
   \langle l_i 0 l_j 0| k 0 \rangle^2 {\cal R}^k(ij;ji)
   \f]
   *
   *@brief Calculates and the Exchange potential for each orbital and
   the Exchange energy.  
   *
   The Exchange potential for the \textit{ith} orbital
   \f[
   \hat{V}_{E} u_{n_j l_j}(r) =  
   -\sum_j \delta_{\sigma_i,\sigma_j} \sum_{k=|l_j-l_i|}^{l_i+l_j} 
   \frac{(2l_i+1)(2l_j+1)}{(2k+1)^2} 
   \langle l_i (-m_i) l_j m_j| k (m_i-m_j) \rangle^2 
   \langle l_i 0 l_j 0| k 0 \rangle^2 
   \frac{{\cal Y}_k(n_jl_j;n_il_i/r)}{r} u_{n_j l_j}(r)
   \f]
   *
   see RadialFunctorUtility.h for the implementation of 
   \f$ {\cal Y}_k(n_jl_j;n_il_i/r) /r \f$ and 
   \f$ {\cal R}^k(ij;ji). \f$  Also includes a function to 
   make the exchange a local potential.
   *
   *@note The sumnation over the orbitals {ij} includes \f$ i=j. \f$
   */

  RadialPotentialBase::value_type
  ExchangePotential::evaluate(const BasisSetType& psi, 
			      RadialOrbitalSet_t& V, int norb) {
    
    value_type ith_orb_coeff, jth_orb_coeff, coeff;
    value_type energy_coeff=0;
    value_type Eexchange=0;  
    //zero_all_orbitals(); V is reset before entering
    RadialOrbital_t Ykij_r(psi(0));
    int nn=0;

    //Loop over all pairs of electrons once
    for(int i=0; i < norb; i++) {
      int si = psi.S[i];
      int mi = psi.M[i];
      int li = psi.L[i];
      int two_li_plus_one = 2*li + 1;

      for(int j=i; j < norb; j++, nn++) {
	int sj = psi.S[j];
	int mj = psi.M[j];
	int lj = psi.L[j];
	int two_lj_plus_one = 2*lj + 1;

	int kmax = li + lj;
	int kmin = abs(li - lj);

	if( si == sj ) {
	  for(int k=kmax; k >= kmin; k-=2) {
	    int two_k_plus_one = 2*k + 1;

	    coeff = static_cast<value_type>(two_li_plus_one * two_lj_plus_one) / 
	      static_cast<value_type>(two_k_plus_one*two_k_plus_one)
	      * CG_coeff->cg(li,lj,k,0,0) * CG_coeff->cg(li,lj,k,-mi,mj) 
	      * CG_coeff->cg(li,lj,k,0,0) * CG_coeff->cg(li,lj,k,-mi,mj);

	    if(i == j) coeff /= 2.0; //double counting

	    ith_orb_coeff = psi.Occ[i] * coeff;
	    jth_orb_coeff = psi.Occ[j] * coeff;
	    energy_coeff = psi.Occ[j] * psi.Occ[i] * coeff;

	    Ykofr(Ykij_r, psi(i), psi(j), k); // Ykofr_phi1_phi2

	    Make_Loc_Pot(V[i], Ykij_r, psi(i), psi(j),jth_orb_coeff);           
	    Make_Loc_Pot(V[j], Ykij_r, psi(j), psi(i),ith_orb_coeff);           
      
	    storage[nn] = -1.0*Phisq_x_Yk(Ykij_r, psi(i), psi(j), energy_coeff);
	    Eexchange += storage[nn];
	  }
	}
      }
    }
  
    return Eexchange;
  }

  void ExchangePotential::getStorage(const BasisSetType& psi,
				     const std::string& RootFileName){
    
    string fileforoutput = RootFileName + ".exchange";
    ofstream fout(fileforoutput.c_str());
    string llabel("spdf");  
    string slabel("d0u");
    int norb = psi.size();
    int nn=0;
    value_type temp;
    fout << "orb#1" << '\t' << "orb#2" << '\t' << "Exchange" << endl;
    fout.precision(12);
    fout.setf(ios::scientific,ios::floatfield);
    for(int i=0; i<norb; i++) {
      for(int j=i; j<norb; j++, nn++) {
	temp = storage[nn];
	if(fabs(temp) < 1e-12) temp = 0.0;
	fout << psi.N[i] << llabel[psi.L[i]] << slabel[psi.S[i]+1]
	     << '\t' << psi.N[j] << llabel[psi.L[j]] << slabel[psi.S[j]+1]
	     << '\t' << temp << endl;
      }
    }
  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

