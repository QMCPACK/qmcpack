//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim 
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
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
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#ifndef OHMMS_QMC_LOCALENERGYESTIMATOR_H
#define OHMMS_QMC_LOCALENERGYESTIMATOR_H
#include <fstream>
#include "OhmmsData/libxmldefs.h"
#include "Estimators/ScalarEstimatorBase.h"

namespace ohmmsqmc {

  /*** A class to evaluate the local energy 
   *
   *The LocalEnergyEstimator evaluates 
   <ul>
   <li> LocalEnergy
   <li> Variance of the LocalEnergy
   </ul>
   and the values of each QHCHamiltonianBase elements added by an application.
   Typical local energies are
   <li> Kinetic Energy
   <li> External Potential Energy
   <li> Electron-Electron Energy
   <li> Ion-Ion Energy (for molecules only)
   <li> Conserved Quantity (VMC only)
   </ul>
   The method of evaluating the estimators is as follows
   \f[
   \langle O \rangle = \frac{\sum_{i=1}^N w_i O_i}{\sum_{i=1}^N w_i}
   \f]
   where the sum runs over the total number of accumulated walkers 
   and \f$ w_i \f$ is the weight of the \f$ ith \f$ walker.
   *
   *The formula for the LocalEnergy 
   \f[ E_L({\bf R}) = \frac{\hat{H} \Psi({\bf R})}{ \Psi({\bf R})} \f]
  */
  template<class T>
  class LocalEnergyEstimator: public ScalarEstimatorBase<T> {

    ///locator of the first data this object handles
    int LocalEnergyIndex;

    ///local data
    T energy_sum, energy_sq_sum;
    ///vector to contain the names of all the constituents of the local energy
    std::vector<string> elocal_name;
    ///vector to contain all the constituents of the local energy
    std::vector<T>  elocal;

  public:

    typedef typename ScalarEstimatorBase<T>::Walker_t Walker_t;
  
    LocalEnergyEstimator(QMCHamiltonian& h):energy_sum(0.0), 
					    energy_sq_sum(0.0) { 
      int nterms=h.size();
      elocal.resize(nterms);
      elocal_name.resize(nterms+2);
      elocal_name[0] = "LocalEnergy";
      elocal_name[1] = "Variance";
      int ii=2;
      for(int i=0; i < nterms; i++, ii++) elocal_name[ii] = h.getName(i);
    }

    /**
       @param record a container class for storing scalar records (name,value)
       @brief add the local energy, variance and all the constituents 
       of the local energy to the scalar record container
    */
    void add2Record(RecordNamedProperty<T>& record) {
      LocalEnergyIndex = record.add(elocal_name[0].c_str());
      for(int i=1; i<elocal_name.size(); i++)
	record.add(elocal_name[i].c_str());
    }

    inline void accumulate(const Walker_t& awalker, T wgt) {
      T e = awalker.Properties(LocalEnergy);
      energy_sum += wgt*e; energy_sq_sum += wgt*e*e;
      for(int i=0, ii=0; i<elocal.size(); i++) {
	elocal[i] += wgt*(awalker.E[i]);
      }
    }

    ///reset all the cumulative sums to zero
    inline void reset() { 
      energy_sum = T(); energy_sq_sum = T(); 
      for(int i=0; i<elocal.size(); i++) {
	elocal[i] = T();
      }
    }

    /*!
     *\param record a container class for storing scalar records (name,value)
     *\param wgtinv the inverse weight
     *\brief calculate the averages and reset to zero
     */
    inline void report(RecordNamedProperty<T>& record, T wgtinv) {
      b_average =  energy_sum*wgtinv;
      b_variance = energy_sq_sum*wgtinv-b_average*b_average;
      register int ir=LocalEnergyIndex;
      record[ir++] = b_average;
      record[ir++] = b_variance;
      for(int i=0; i<elocal.size(); i++) {
	record[ir++] = elocal[i]*wgtinv;
	elocal[i] = T();
      }
      energy_sum = T(); energy_sq_sum = T(); 
    }
  };

  //   template<class T>
  //   class LocalEnergyEstimator: public EstimatorBase<T> {

  //     std::string RootName;
  //     std::ofstream fout;
  //     std::vector<std::string>  elocal_name;
  //     std::vector<T>  elocal;
  //   public:

  //     LocalEnergyEstimator(QMCHamiltonian& h) { 
  //       int nterms=h.size();
  //       elocal.resize(nterms);
  //       elocal_name.resize(nterms+2);
  //       elocal_name[0] = "LocalEnergy";
  //       elocal_name[1] = "Variance";
  //       int ii=2;
  //       for(int i=0; i < nterms; i++, ii++) elocal_name[ii] = h.getName(i);
  //     }

  //     inline void resetReportSettings(const string& aname) {
  //       if(RootName != aname) {
  //         finalize();
  //       }
  //       RootName = aname;
  //       string fname(aname);
  //       fname.append(".e.dat");
  //       fout.open(fname.c_str());
  //       fout.setf(ios::right,ios::scientific);
  //       fout << "#      index ";
  //       for(int i=0; i < elocal_name.size(); i++) 
  //        fout << setw(14) << elocal_name[i];
  //       fout << endl;
  //       fout<<setprecision(8);
  //     }

  //     inline void finalize() {
  //       if(fout.is_open()) fout.close();
  //     }

  //     inline void accumulate(const MCWalkerConfiguration& W) {

  //       for(MCWalkerConfiguration::const_iterator it = W.begin(); 
  // 	  it != W.end(); it++){
  // 	T w = (*it)->Properties(Weight);
  // 	for(int i=0; i<elocal.size(); i++) {
  // 	  elocal[i] += w*((*it)->E[i]);
  // 	}
  // 	add(w,(*it)->Properties(LocalEnergy));
  //       }
  //     }

  //     inline void report(int iter){
  //       if(iter%stride == 0) { 
  // 	register T winv = 1.0/d_v[WEIGHT];
  //         register T eavg = d_v[VALUE]*winv;
  //         register T evar = d_v[VALUE_SQ]*winv-eavg*eavg;

  // 	fout << setw(10) << iter;
  //         fout << setw(15) << eavg << setw(15) << evar;
  // 	for(int i=0; i<elocal.size(); i++) {
  //            fout << setw(15) << elocal[i]*winv; elocal[i]=0.0;
  //         }
  //         fout << endl;

  //         d_v =0.0;
  //         d_v[AVERAGE]=eavg;
  //         d_v[VARIANCE]=evar;
  //       }
  //     }

  //     bool get(ostream& os) const {return true;}
  //     bool put(istream& is) {return true;}
  //     bool put(xmlNodePtr cur) {return true;}

  //   };
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
