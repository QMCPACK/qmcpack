//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim and Jordan Vincent
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim and Jordan Vincent
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
#ifndef OHMMS_QMC_FIXEDNODE_BRANCHER_H
#define OHMMS_QMC_FIXEDNODE_BRANCHER_H

#include <deque>
#include <algorithm>
#include <numeric>
#include "OhmmsData/ParameterSet.h"

namespace ohmmsqmc {

  /** For use in Fixed-Node Diffusion Monte Carlo. 
   *
   *@brief Provides routines for the fixed-node diffusion
   Monte Carlo algorithm.  
   *
   Calculates the Branching Green's function
   \f[
   G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])
   \f]
   which is used for the weight and mulitplicity of each walker
   \f[ Weight =  G_{branch} \f]
   \f[ Mulitplicity =  G_{branch} + \nu \f]
   and to update the energy offset
   \f[ E_T = <E_G> - feed \log \left(\frac{P(t)}{P_0}\right) \f]
   where \f$P(t)\f$ is the current population, \f$P_0\f$ is the 
   ideal population, and \f$<E_G>\f$ is an estimate of the
   local energy. 
  */
  template<class T>
  class MolecuFixedNodeBranch {

  public:

    ///the timestep
    T Tau;
    ///feedback parameter to control the population
    T Feed;
    ///energy offset to control branching
    T E_T;

    ///Feed*log(N)
    T logN;

    ///ideal population
    int Nideal;
    ///maximum population
    int Nmax;
    ///minumum population
    int Nmin;

    ///maximum copies of a walker
    int MaxCopy;

    ///size of container Eg
    int EgBufferSize;

    ///counts the number of times update has been called
    int Counter;

    int InFeed;

    ///container to store old values of e_ref
    std::deque<T> Eg;

    ///Constructor
    MolecuFixedNodeBranch(T tau, int nideal): 
      Tau(tau), E_T(0.0), Nideal(nideal), 
      EgBufferSize(100), Counter(0), MaxCopy(10) {
      //feed = 1.0/(50.0*tau); 
      InFeed = 50;
      Eg.resize(EgBufferSize);
      for(int i=0; i<EgBufferSize; i++) Eg[i] = 0.0; 
    }
    
    ///return true if the nodal surface is crossed
    inline bool operator()(T psi0, T psi1) const { return psi0*psi1 < 0;}
    
    /**
     *@param tau effective time step
     *@param emixed mixed energy \f$(E_L(R)+E_L(R'))/2\f$
     *@param reject rejection probability
     *@return \f$G_{branch}\f$
     *@brief Calculates the Branching Green's function
     \f[G_{branch} = \exp(-\tau \left[(E_L(R)+E_L(R'))/2-E_T\right])\f]
     *@note Use the rejection probability \f$q\f$ to limit \f$G_{branch}\f$
     \f[ G_{branch} = \min \left(\frac{1}{2q},G_{branch}\right). \f]
     */
    inline T branchGF(T tau, T emixed, T reject) const { 
      return exp(-tau*(emixed-E_T));
      //return min(0.5/(reject+1e-12),exp(-tau*(emix-E_T)));
    }
    
    ///set \f$ <E_G> = eg \f$
    inline void setEguess(T eg){
      E_T = eg;
      for(int i=0; i<EgBufferSize; i++) Eg[i] = eg;
    } 


    /** call MCWalkerConfiguration::branch
     *@param iter the iteration
     *@param w the walker ensemble
     *@return the number of walkers after branching
     */
    inline int branch(int iter, MCWalkerConfiguration& w) {
      return w.branch(10,Nmax,Nmin);
    }

    /**
     *@param pop the current population \f$ P(t) \f$
     *@param eavg the average value of the local energy
     *@return the energy offset \f$E_T\f$
     *@brief Update the energy offset
     \f[ E_T = <E_G> - feed \log \(\frac{P(t)}{P_0}\) \f]
    */
    inline T update(T pop_now, T eavg) {
      Counter++;
      //pop off the last value of Eg
      Eg.pop_back();
      //insert a new value at the beggining of the deque
      Eg.push_front(eavg);
      int mlimit = std::min(Counter/2+1,EgBufferSize);
      T Esum = std::accumulate(Eg.begin(),Eg.begin()+mlimit,T());
      //calculate the average
      //average over the last half of the simulation
      //for(int i=0; i<limit; i++) Esum = Eg[i];
      E_T = Esum/static_cast<T>(mlimit)-Feed*log(static_cast<T>(pop_now))+logN;
      return E_T;
    }


    /**
     *@param q the current xmlNode 
     *@brief Parse the xml file for parameters
     <ul>
     <li> en_ref: a reference energy
     <li> num_gen: number of generations \f$N_G\f$ to reach 
     equilibrium, used in the feedback parameter
     \f$ feed = \frac{1}{N_G \tau} \f$ 
     <\ul>
    */
    bool put(xmlNodePtr cur, OhmmsInform *LogOut){
      ParameterSet m_param;
      m_param.add(E_T,"en_ref","AU");
      m_param.add(InFeed,"num_gen","int");
      m_param.add(MaxCopy,"max_copy","int");
      m_param.add(Nideal,"target_walkers","int");
      m_param.put(cur);
      for(int i=0; i<EgBufferSize; i++) Eg[i] = E_T;
      reset();
      LogOut->getStream() << "reference energy = " << Eg[0] << endl;
      LogOut->getStream() << "number of generations = " << Feed << endl;
      LogOut->getStream() << "feedback = " << Feed << endl;
      return true;
    }

    void reset() {
      Nmax = 2*Nideal;
      Nmin = static_cast<int>(Nideal/2);
      Feed = 1.0/(static_cast<T>(InFeed)*Tau);
      logN = Feed*log(static_cast<T>(Nideal));

      XMLReport("Target walkers = " << Nideal)
      XMLReport("Branching: Referece energy = " << Eg[0])
      XMLReport("MaxCopy for branching = " << MaxCopy)
      XMLReport("Branching: Number of generations = " << InFeed)
      XMLReport("Branching: Feedback parameter = " << Feed)
    }

    void write(hid_t grp) {
      hsize_t dim = static_cast<hsize_t>(EgBufferSize);
      //vector<T> etrial(Eg.begin(),Eg.end());
      hid_t dataspace  = H5Screate_simple(1, &dim, NULL);
      hid_t dataset =  
        H5Dcreate(grp, "TrialEnergies", H5T_NATIVE_DOUBLE, dataspace, H5P_DEFAULT);
      hid_t ret = 
        H5Dwrite(dataset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,&Eg[0]);
      H5Dclose(dataset);
      H5Sclose(dataspace);
      //make a char array which contains the species names separated by space
      std::stringstream s;
      s << Nideal<< ' ' << InFeed << ' ' << Counter << ' ' << E_T << '\0';
      dim = s.str().size();
      dataspace = H5Screate_simple(1, &dim, NULL);
      dataset = H5Dcreate(grp, "BranchParameters", H5T_NATIVE_CHAR, dataspace, H5P_DEFAULT);
      ret =  H5Dwrite(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT,s.str().c_str());
      H5Dclose(dataset);
      H5Sclose(dataspace);

    }

    void read(const string& aroot) {
      string h5file(aroot);
      h5file.append(".config.h5");
      hid_t h_file =  H5Fopen(h5file.c_str(),H5F_ACC_RDWR,H5P_DEFAULT);
      hid_t h_config = H5Gopen(h_file,"config_collection");
      hsize_t dataset = H5Dopen(h_config, "TrialEnergies");
      if(dataset) {
        hsize_t ret = H5Dread(dataset, H5T_NATIVE_DOUBLE, 
            H5S_ALL, H5S_ALL, H5P_DEFAULT, &Eg[0]);
        H5Dclose(dataset);
        char temp[256];
        dataset = H5Dopen(h_config, "BranchParameters");
        ret = H5Dread(dataset, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, temp);
        H5Dclose(dataset);
        std::istringstream s(temp);
        s >> Nideal>>InFeed >>Counter >> E_T;
        reset();
        Counter = std::min(Counter,EgBufferSize);
      }
      H5Gclose(h_config);
      H5Fclose(h_file);

    }

  private:
    ///default constructor (disabled)
    MolecuFixedNodeBranch(){}
  };
  
}
#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/

