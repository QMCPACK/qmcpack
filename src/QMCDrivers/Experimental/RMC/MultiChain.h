//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#ifndef QMCPLUSPLUS_MULTICHAIN_H
#define QMCPLUSPLUS_MULTICHAIN_H
#include "Particle/MCWalkerConfiguration.h"
#include "OhmmsPETE/OhmmsVector.h"
#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/IteratorUtility.h"
#include "Numerics/HDFNumericAttrib.h"
#include "QMCDrivers/SpaceWarp.h"
#include "QMCDrivers/DriftOperators.h"
#include <deque>
namespace qmcplusplus
{

struct Bead: public MCWalkerConfiguration::Walker_t
{

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::RealType RealType;
  typedef MCWalkerConfiguration::ParticlePos_t ParticlePos_t;
  typedef MCWalkerConfiguration::ParticleLaplacian_t ParticleLaplacian_t;

  //Drift is moved  from Walker
  ParticlePos_t Drift;
  Vector<int> BeadSignWgt;
  std::vector<ParticlePos_t*> Gradients;
  std::vector<ParticleLaplacian_t*> Laplacians;
  std::vector<ParticlePos_t*> DriftVectors;
  Matrix<RealType> Action;
  RealType TransProb[2];
  bool ScaleDrift;
  std::vector<RealType> Tau_eff;

  //Extra Observable Spots for Pressure, etc.
  Vector<RealType> deltaRSquared;
  int stepmade,timesTouched;
  std::vector<int> ptclAge;

  inline Bead(const Bead& a)
  {
    makeCopyBead(a);
  }

  inline ~Bead()
  {
    delete_iter(Gradients.begin(),Gradients.end());
    delete_iter(Laplacians.begin(),Laplacians.end());
    delete_iter(DriftVectors.begin(),DriftVectors.end());
  }

  inline Bead(const Walker_t& a):Walker_t(a)
  {
    resize_bead(R.size());
  }

  //       deltaRSquared.resize(3);
  inline Bead(const Walker_t& a, std::string scaleBeadDrift): Walker_t(a)
  {
    ScaleDrift = ((scaleBeadDrift=="true")||(scaleBeadDrift=="yes"));
    resize_bead(R.size());
  }

  inline void resize_bead(int n)
  {
    Drift.resize(R.size());
    int rows=Properties.rows();
    Resize_Grad_and_Action(rows,R.size());
    BeadSignWgt.resize(rows);
    Tau_eff.resize(rows);
  }

  inline Bead& operator=(const Bead& a)
  {
    makeCopyBead(a);
    return *this;
  }

  inline void makeCopyBead(const Bead& a)
  {
    makeCopy(a);
    Drift.resize(R.size());
    int rows=a.Gradients.size();
    Resize_Grad_and_Action(rows,a.size());
    Action=a.Action;
    Tau_eff.resize(rows);
    for(int i=0; i<rows; i++)
      *Gradients[i] = *(a.Gradients[i]);
    for(int i=0; i<rows; i++)
      *Laplacians[i] = *(a.Laplacians[i]);
    for(int i=0; i<rows; i++)
      *DriftVectors[i] = *(a.DriftVectors[i]);
    Tau_eff = a.Tau_eff ;
    BeadSignWgt.resize(rows);
    BeadSignWgt=a.BeadSignWgt;
    TransProb[0]=a.TransProb[0];
    TransProb[1]=a.TransProb[1];
    stepmade=a.stepmade;
    timesTouched=a.timesTouched;
    deltaRSquared=a.deltaRSquared;
  }

  inline void Resize_Grad_and_Action(int n, int m)
  {
    int curg=Gradients.size();
    while(curg<n)
    {
      Gradients.push_back(new ParticlePos_t(m));
      Laplacians.push_back(new ParticleLaplacian_t(m));
      ++curg;
    }
    curg=DriftVectors.size();
    while(curg<n)
    {
      DriftVectors.push_back(new ParticlePos_t(m));
      ++curg;
    }
    Tau_eff.resize(n);
    deltaRSquared.resize(3);
    deltaRSquared=0.0;
    Action.resize(n,3);
    Action=0.0;
    ptclAge.resize(R.size());
    for (int i=0; i<ptclAge.size(); i++)
      ptclAge[i]=0;
  }

  void registerData(Buffer_t& buf);
  void copyFromBuffer(Buffer_t& buf);
  void copyToBuffer(Buffer_t& buf);

  inline void getDrift(std::vector<RealType>& LogNorm)
  {
    int npsi(Properties.rows());
    //compute Drift
    RealType denom(0.e0),wgtpsi;
    Drift=0.e0;
    for(int ipsi=0; ipsi<npsi; ipsi++)
    {
      wgtpsi=BeadSignWgt[ipsi]*std::exp(2.0*( Properties(ipsi,LOGPSI)- LogNorm[ipsi]
                                              -Properties(0,LOGPSI)   + LogNorm[0]));
      denom += wgtpsi;
      Drift += (wgtpsi*(*DriftVectors[ipsi]));
    }
    denom=1.0/denom;
    Drift *= denom;
  }

  inline void getScaledDrift(std::vector<RealType>& LogNorm, RealType Tau)
  {
    int npsi(Properties.rows());
    //compute Drift
    RealType denom(0.e0),wgtpsi;
    Drift=0.e0;
    ParticlePos_t TMPgrad(Drift);
    for(int ipsi=0; ipsi<npsi; ipsi++)
    {
      wgtpsi=BeadSignWgt[ipsi]*std::exp(2.0*( Properties(ipsi,LOGPSI)- LogNorm[ipsi]
                                              -Properties(0,LOGPSI)   + LogNorm[0]));
      denom += wgtpsi;
      if (ScaleDrift)
      {
        //setScaledDriftPbyP(Tau,*Gradients[ipsi],TMPgrad);
        Tau_eff[ipsi] = setLargestScaledDriftPbyP(Tau,*Gradients[ipsi],TMPgrad);
        (*DriftVectors[ipsi]) = TMPgrad;
        Drift += (wgtpsi*TMPgrad);
      }
      else
      {
        Tau_eff[ipsi] = Tau;
        (*DriftVectors[ipsi]) = Tau*(*Gradients[ipsi]);
        Drift += Tau*(wgtpsi*(*Gradients[ipsi]));
      }
    }
    denom=1.0/denom;
    Drift *= denom;
  }

  inline void getScaledDriftSingle(std::vector<RealType>& LogNorm, RealType Tau, int ipsi)
  {
    if (ScaleDrift)
    {
//        setScaledDriftPbyP(Tau,*Gradients[ipsi],(*DriftVectors[ipsi]));
      Tau_eff[ipsi] = setLargestScaledDriftPbyP(Tau,*Gradients[ipsi],(*DriftVectors[ipsi]));
    }
    else
    {
      Tau_eff[ipsi] = Tau;
      (*DriftVectors[ipsi]) = Tau* (*Gradients[ipsi]);
    }
  }

  inline void getDrift(std::vector<RealType>& Jacobian, SpaceWarp& PtclWarp,std::vector<RealType>& LogNorm)
  {
    QMCTraits::PosType WarpDrift;
    int npsi(Properties.rows());
    int nptcl(Gradients[0]->size());
    //compute Drift
    RealType denom(0.e0),wgtpsi;
    Drift=0.e0;
    for(int ipsi=0; ipsi<npsi; ipsi++)
    {
      wgtpsi=BeadSignWgt[ipsi]*Jacobian[ipsi]*std::exp(2.0*( Properties(ipsi,LOGPSI)-LogNorm[ipsi]
             -Properties(0,LOGPSI)   +LogNorm[0]));
      denom += wgtpsi;
      for(int iptcl=0; iptcl< nptcl; iptcl++)
      {
        WarpDrift=dot(  (*Gradients[ipsi])[iptcl],PtclWarp.get_Jacob_matrix(iptcl,ipsi)  )
                  +5.0e-1*PtclWarp.get_grad_ln_Jacob(iptcl,ipsi) ;
        Drift[iptcl] += (wgtpsi*WarpDrift);
      }
    }
    denom=1.0/denom;
    //Drift = denom*Drift;
    Drift *= denom;
  }
};

//struct MultiChain: public std::deque<Bead*> {
struct MultiChain
{

  typedef std::deque<Bead*>               Container_t;
  typedef Container_t::iterator           iterator;
  typedef Container_t::const_iterator     const_iterator;
  typedef Container_t::reference          reference;
  typedef Container_t::const_reference    const_reference;

  typedef MCWalkerConfiguration::Walker_t Walker_t;
  typedef MCWalkerConfiguration::RealType RealType;
  typedef Bead::Buffer_t                  Buffer_t;

  /// Direction of growth
  int GrowthDirection;

  ///The index of the middle
  int Middle;

  ///The index of the Last
  int Last;

  ///The number of H/Psi pairs
  int nPsi;

  //age of entire reptile, number of steps taken.
  int Age;

  ///size of bead buffe for each Bead
  int BeadBufferSize;
  ///size of MultiChain buffer for the collective data
  int MyBufferSize;

  ///hdf5 handle to multichain
  hid_t h_config;

  RealType GlobalWgt;

  Vector<RealType> GlobalAction,UmbrellaWeight;
  Vector<int>      GlobalSignWgt,RefSign;

  Container_t      Beads;
  Vector< Vector<RealType> > OldObs;

  /**  constructor
   * @param abead Bead used to generate len Beads to form a chain
   * @param len size of the chain
   * @param direction initial growth direction
   * @param npsi number of Psi/H pairs
   */
  //MultiChain(Walker_t* awalker,int len, int direction, int npsi):
  MultiChain(Bead* abead,int len, int direction, int npsi);

  /** destructor
   *
   * Need to clean up the walkers in the repository and the polymer chain
   */
  ~MultiChain();

  inline reference operator[](int n)
  {
    return Beads[n];
  }
  inline const_reference operator[](int n) const
  {
    return Beads[n];
  }
  inline iterator begin()
  {
    return Beads.begin();
  }
  inline iterator end()
  {
    return Beads.end();
  }
  inline reference front()
  {
    return Beads.front();
  }
  inline reference back()
  {
    return Beads.back();
  }

  //An easy way to grab the center bead for observables calculated in center.
  inline iterator middle()
  {
    iterator Bit = Beads.begin()+Middle;
    return Bit;
  }
  inline reference center()
  {
    return *(this->middle());
  }

  INLINE_ALL void push_front(Bead* abead)
  {
    Beads.push_front(abead);
  }
  INLINE_ALL void push_back(Bead* abead)
  {
    Beads.push_back(abead);
  }
  INLINE_ALL void pop_back()
  {
    Beads.pop_back();
  }
  INLINE_ALL void pop_front()
  {
    Beads.pop_front();
  }
  inline void flip()
  {
    GrowthDirection = std::abs(GrowthDirection-1); //flip the direction
  }

  inline void setRefSign(int ipsi, RealType phase)
  {
    RefSign[ipsi]=(phase>M_PI_2)?-1:1;
  }

  inline int getSign(RealType phase)
  {
    return (phase>M_PI_2)?-1:1;
  }

  inline void Reset_Ages()
  {
    for(Container_t::iterator Pit=Beads.begin(); Pit!=Beads.end(); Pit++)
    {
      (*Pit)->stepmade=0;
      (*Pit)->timesTouched=0;
    };
  };

  /** copy the restart data from buf
   * @param buf buffer to read from
   */
  void copyFromBuffer(Buffer_t& buf);

  /** add the restart data to buf
   * @param buf buffer to write
   *
   * add takes care of memory allocation and assignment
   */
  void copyToBuffer(Buffer_t& buf);

  /** read multi-chain configuration from a file
   * @param aroot root name
   * @return true, if the input file contains valid data
   */
  bool read(const std::string& aroot);

  /** read MultiChain tree a group *
   * @param grp hdf5 group
   * @return true, if the input file contains valid data
   *
   * Typically, grp is the file id
   */
  bool read(hid_t grp);

  /** open hdf5 and initialize the dataspace
   * @param aroot root name of hdf5 output
   */
  void open(const std::string& aroot);
  /** record the data for restart */
  void record();
  /** close the file */
  void close();

};
}
#endif
