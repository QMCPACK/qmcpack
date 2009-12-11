//////////////////////////////////////////////////////////////////
// (c) Copyright 2003  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
#ifndef QMCPLUSPLUS_WALKER_H
#define QMCPLUSPLUS_WALKER_H

#include "OhmmsPETE/OhmmsMatrix.h"
#include "Utilities/PooledData.h"
#include <assert.h>
#include <deque>
namespace qmcplusplus {

  /** an enum denoting index of physical properties 
   *
   * LOCALPOTENTIAL should be always the last enumeation 
   */
  enum {LOGPSI=0,       /*!< log(fabs(psi)) instead of square of the many-body wavefunction \f$|\Psi|^2\f$ */
	SIGN,           /*!< value of the many-body wavefunction \f$\Psi(\{R\})\f$ */
        UMBRELLAWEIGHT, /*!< sum of wavefunction ratios for multiple H and Psi */
        R2ACCEPTED,     /*!< r^2 for accepted moves */
        R2PROPOSED,     /*!< r^2 for proposed moves */
        DRIFTSCALE,     /*!< scaling value for the drift */
	LOCALENERGY,    /*!< local energy, the sum of all the components */
	LOCALPOTENTIAL, /*!< local potential energy = local energy - kinetic energy */
        NUMPROPERTIES   /*!< the number of properties */
       };
  
  /** A container class to represent a walker.
   *
   * A walker stores the particle configurations {R}  and a property container.
   * RealTypehe template (P)articleSet(A)ttribute is a generic container  of position types.
   * RealTypehe template (G)radient(A)ttribute is a generic container of gradients types.
   * Data members for each walker
   * - ID : identity for a walker. default is 0. 
   * - Age : generation after a move is accepted.
   * - Weight : weight to take the ensemble averages
   * - Multiplicity : multiplicity for branching. Probably can be removed.
   * - Properties  : 2D container. RealTypehe first index corresponds to the H/Psi index and second index >=NUMPROPERTIES.
   * - DataSet : anonymous container. 
   */
  template<typename t_traits, typename p_traits>
  struct Walker 
  {
    enum {DIM=t_traits::DIM};
    /** typedef for real data type */
    typedef typename t_traits::RealType RealType;
    /** typedef for value data type. */
    typedef typename t_traits::ValueType ValueType;
    /** array of particles */
    typedef typename p_traits::ParticlePos_t ParticlePos_t;
    /** array of gradients */
    typedef typename p_traits::ParticleGradient_t ParticleGradient_t;
    /** array of laplacians */
    typedef typename p_traits::ParticleLaplacian_t ParticleLaplacian_t;

    ///typedef for the property container, fixed size
    typedef Matrix<RealType>      PropertyContainer_t;
    typedef PooledData<RealType>  Buffer_t;

    ///id reserved for forward walking
    long ID;
    ///id reserved for forward walking
    long ParentID;
    ///DMCgeneration
    int Generation;
    ///Age of this walker age is incremented when a walker is not moved after a sweep
    int Age;
    ///Weight of the walker
    RealType Weight;
    /** Number of copies for branching
     *
     * When Multiplicity = 0, this walker will be destroyed.
     */
    RealType Multiplicity;

    /**the configuration vector (3N-dimensional vector to store
       the positions of all the particles for a single walker)*/
    ParticlePos_t R;
    /** \f$ \nabla_i d\log \Psi for the i-th particle */
    ParticleGradient_t G;
    /** \f$ \nabla^2_i d\log \Psi for the i-th particle */
    ParticleLaplacian_t L;
    /////drift of the walker \f$ Drift({\bf R}) = \tau v_{drift}({\bf R}) \f$
    //ParticlePos_t Drift;

    ///scalar properties of a walker
    PropertyContainer_t  Properties;
    
    ///Property history vector
    vector<vector<RealType> >  PropertyHistory;
    vector<int> PHindex;

    ///buffer for the data for particle-by-particle update
    Buffer_t DataSet;

    ///default constructor
    inline Walker() : ID(0),ParentID(0), Generation(0),Age(0),
                  Weight(1.0e0),Multiplicity(1.0e0) 
    { 
      Properties.resize(1,NUMPROPERTIES);
      reset();
    }

    ///create a walker for n-particles
    inline explicit Walker(int nptcl) : ID(0),ParentID(0), Generation(0),Age(0),
                           Weight(1.0e0),Multiplicity(1.0e0)
    {
      Properties.resize(1,NUMPROPERTIES);
      resize(nptcl); 
      reset();
    }

    inline int addPropertyHistory(int leng)
    {
      int newL = PropertyHistory.size();
      vector<RealType> newVecHistory=vector<RealType>(leng,0.0);
      PropertyHistory.push_back(newVecHistory);
      PHindex.push_back(0);
      return newL;
    }
    
    inline void resetPropertyHistory( )
    {
      for(int i=0;i<PropertyHistory.size();i++)
      {
        PHindex[i]=0;
	for(int k=0;k<PropertyHistory[i].size();k++)
	{
	  PropertyHistory[i][k]=0.0;
	}
      }
    }

    inline void addPropertyHistoryPoint(int index, RealType data)
    {
      PropertyHistory[index][PHindex[index]]=(data);
      PHindex[index]++;
      if (PHindex[index]==PropertyHistory[index].size()) PHindex[index]=0;
//       PropertyHistory[index].pop_back();
    }
   
    inline RealType getPropertyHistorySum(int index, int endN)
    {
      RealType mean=0.0; 
      typename vector<RealType>::const_iterator phStart;
      phStart=PropertyHistory[index].begin()+PHindex[index];
      for(int i=0;i<endN;phStart++,i++){
        if (phStart>=PropertyHistory[index].end()) phStart -= PropertyHistory[index].size();
        mean+= (*phStart);
      }
      return mean ;
    }

    inline ~Walker() { }

    ///assignment operator
    inline Walker& operator=(const Walker& a) {
      if(this != &a) makeCopy(a);
      return *this;
    }

    ///return the number of particles per walker
    inline int size() const { return R.size(); }

    ///resize for n particles
    inline void resize(int nptcl) 
    {
      R.resize(nptcl); 
      G.resize(nptcl);
      L.resize(nptcl);
      //Drift.resize(nptcl); 
    }

    ///copy the content of a walker
    inline void makeCopy(const Walker& a) 
    {
      ID=a.ID;
      ParentID=a.ParentID;
      Generation=a.Generation;
      Age=a.Age;
      Weight=a.Weight;
      Multiplicity=a.Multiplicity;
      if(R.size()!=a.R.size()) resize(a.R.size());
      R = a.R;
      G = a.G;
      L = a.L;
      //Drift = a.Drift;
      Properties.copy(a.Properties);
      DataSet=a.DataSet;
      
      if(PropertyHistory.size()!=a.PropertyHistory.size()) PropertyHistory.resize(a.PropertyHistory.size());
      for(int i=0;i<PropertyHistory.size();i++) PropertyHistory[i]=a.PropertyHistory[i];
      PHindex=a.PHindex;
    }

    //return the address of the values of Hamiltonian terms
    inline RealType* restrict getPropertyBase() {
      return Properties.data();
    }

    //return the address of the values of Hamiltonian terms
    inline const RealType* restrict getPropertyBase() const {
      return Properties.data();
    }

    ///return the address of the i-th properties
    inline RealType* restrict getPropertyBase(int i) {
      return Properties[i];
    }

    ///return the address of the i-th properties
    inline const RealType* restrict getPropertyBase(int i) const {
      return Properties[i];
    }


    /** reset the property of a walker
     *@param logpsi \f$\log |\Psi|\f$
     *@param sigN  sign of the trial wavefunction
     *@param ene the local energy
     *
     *Assign the values and reset the age
     * but leave the weight and multiplicity 
     */
    inline void resetProperty(RealType logpsi, RealType sigN, RealType ene) 
    {
      Age=0;
      //Weight=1.0;
      Properties(LOGPSI)=logpsi;
      Properties(SIGN)=sigN;
      Properties(LOCALENERGY) = ene;
    }

    /** reset the property of a walker
     * @param logpsi \f$\log |\Psi|\f$
     * @param sigN  sign of the trial wavefunction
     * @param ene the local energy
     * @param r2a \f$r^2\f$ for the accepted moves
     * @param r2p \f$r^2\f$ for the proposed moves
     * @param vq \f$\bar{V}/V\f$ scaling to control node divergency in JCP 93
     *
     *Assign the values and reset the age
     * but leave the weight and multiplicity 
     */
    inline void resetProperty(RealType logpsi, RealType sigN, RealType ene, RealType r2a, RealType r2p, RealType vq) 
    {
      Age=0;
      Properties(LOGPSI)=logpsi;
      Properties(SIGN)=sigN;
      Properties(LOCALENERGY) = ene;
      Properties(R2ACCEPTED) = r2a;
      Properties(R2PROPOSED) = r2p;
      Properties(DRIFTSCALE) = vq;
    }

/** marked to die
     *
     * Multiplicity and weight are set to zero.
     */
    inline void willDie() {
      Multiplicity=0;
      Weight=0.0;
    }

    /** reset the walker weight, multiplicity and age */
    inline void reset() {
      Age=0;
      Multiplicity=1.0e0;
      Weight=1.0e0;
    }

    inline void resizeProperty(int n, int m) {
      Properties.resize(n,m);
    }


    /** byte size for a packed message 
     *
     * ID, Age, Properties, R, Drift, DataSet is packed
     */
    inline int byteSize() 
    {
      int numPH(0);
      for(int iat=0; iat<PropertyHistory.size();iat++) numPH += PropertyHistory[iat].size();
      return 2*sizeof(long)+2*sizeof(int)+ PHindex.size()*sizeof(int)
        +(Properties.size()+DataSet.size()+ numPH )*sizeof(RealType)
        +R.size()*(DIM*sizeof(RealType)+(DIM+1)*sizeof(ValueType));//R+G+L
        //+R.size()*(DIM*2*sizeof(RealType)+(DIM+1)*sizeof(ValueType));//R+Drift+G+L
    }

    template<class Msg>
    inline Msg& putMessage(Msg& m) {
      const int nat=R.size();
      m << ID << ParentID << Generation << Age;
      m.Pack(&(R[0][0]),nat*OHMMS_DIM);
#if defined(QMC_COMPLEX)
      m.Pack(reinterpret_cast<RealType*>(&(G[0][0])),nat*OHMMS_DIM*2);
      m.Pack(reinterpret_cast<RealType*>(L.first_address()),nat*2);
#else
      m.Pack(&(G[0][0]),nat*OHMMS_DIM);
      m.Pack(L.first_address(),nat);
#endif
      m.Pack(Properties.data(),Properties.size());
      m.Pack(DataSet.data(),DataSet.size());
      //Properties.putMessage(m);
      //DataSet.putMessage(m);
      for(int iat=0; iat<PropertyHistory.size();iat++)m.Pack(&(PropertyHistory[iat][0]),PropertyHistory[iat].size());
      m.Pack(&(PHindex[0]),PHindex.size());
      return m;
    }

    template<class Msg>
    inline Msg& getMessage(Msg& m) {
      const int nat=R.size();
      m>>ID >> ParentID >> Generation >> Age;
      m.Unpack(&(R[0][0]),nat*OHMMS_DIM);
#if defined(QMC_COMPLEX)
      m.Unpack(reinterpret_cast<RealType*>(&(G[0][0])),nat*OHMMS_DIM*2);
      m.Unpack(reinterpret_cast<RealType*>(L.first_address()),nat*2);
#else
      m.Unpack(&(G[0][0]),nat*OHMMS_DIM);
      m.Unpack(L.first_address(),nat);
#endif
      m.Unpack(Properties.data(),Properties.size());
      m.Unpack(DataSet.data(),DataSet.size());
      //Properties.getMessage(m);
      //DataSet.getMessage(m);
      for(int iat=0; iat<PropertyHistory.size();iat++)m.Unpack(&(PropertyHistory[iat][0]),PropertyHistory[iat].size());
      m.Unpack(&(PHindex[0]),PHindex.size());
      return m;
    }

  };

  template<class RealType, class PA>
    ostream& operator<<(ostream& out, const Walker<RealType,PA>& rhs)
    {
      copy(rhs.Properties.begin(), rhs.Properties.end(), 
      	   ostream_iterator<double>(out," "));
      out << endl;
      out << rhs.R;
      return out;
    }
}

#endif
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
