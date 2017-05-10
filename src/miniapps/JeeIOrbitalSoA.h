//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_EEIJASTROW_OPTIMIZED_SOA_H
#define QMCPLUSPLUS_EEIJASTROW_OPTIMIZED_SOA_H
#include "Configuration.h"
#if QMC_BUILD_LEVEL<5
#include "QMCWaveFunctions/OrbitalBase.h"
#endif
#include "Particle/DistanceTableData.h"
#include <simd/allocator.hpp>
#include <simd/algorithm.hpp>
#include <map>
#include <numeric>

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 *  @brief Specialization for three-body Jastrow function using multiple functors
 *
 *Each pair-type can have distinct function \f$u(r_{ij})\f$.
 *For electrons, distinct pair correlation functions are used
 *for spins up-up/down-down and up-down/down-up.
 */
template<class FT>
class JeeIOrbitalSoA: public OrbitalBase
{
  ///type of each component U, dU, d2U;
  using valT=typename FT::real_type;
  ///element position type
  using posT=TinyVector<valT,OHMMS_DIM>;
  ///use the same container 
  using RowContainer=DistanceTableData::RowContainer;
  ///table index for i-el, el-el is always zero
  int myTableID;
  //nuber of particles
  int Nelec, Nion;
  //number of groups of the target particleset
  int eGroups, iGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;
  ///diff value
  RealType DiffVal;

  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  ParticleAttrib<valT> Uat;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  ParticleAttrib<posT> dUat;
  valT *FirstAddressOfdU, *LastAddressOfdU;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  ParticleAttrib<valT> d2Uat;
  valT cur_Uat;
  aligned_vector<valT> cur_u, cur_du, cur_d2u;
  aligned_vector<valT> old_u, old_du, old_d2u;
  ///container for the Jastrow functions
  Array<FT*,3> F;

  std::map<std::string,FT*> J3Unique;
  //YYYY
  std::map<FT*,int> J3UniqueIndex;

  struct IonDataCompact
  {
    std::vector<int> elecs_inside;
    aligned_vector<RealType> elecs_dist;
    RealType cutoff_radius;
    IonDataCompact() : cutoff_radius(0.0) { }
  };

  std::vector<IonDataCompact> IonDataList;

  std::vector<valT> Ion_cutoff;

  // Used for evaluating derivatives with respect to the parameters
  int NumVars;
  Array<std::pair<int,int>,3> VarOffset;
  Vector<RealType> dLogPsi;
  Array<PosType,2> gradLogPsi;
  Array<RealType,2> lapLogPsi;

  // Temporary store for parameter derivatives of functor
  // The first index is the functor index in J3Unique.  The second is the parameter index w.r.t. to that
  // functor
  std::vector<std::vector<RealType> >               du_dalpha;
  std::vector<std::vector<PosType> >             dgrad_dalpha;
  std::vector<std::vector<Tensor<RealType,3> > > dhess_dalpha;

public:

  ///alias FuncType
  using FuncType=FT;

  JeeIOrbitalSoA(ParticleSet& ions, ParticleSet& elecs): Ions(ions), NumVars(0)
  {
    myTableID=elecs.addTable(Ions,DT_SOA);
    init(elecs);
  }

  ~JeeIOrbitalSoA() { }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    JeeIOrbitalSoA<FT>* eeIcopy=new JeeIOrbitalSoA<FT>(*this);
    return eeIcopy;
  }

  void init(ParticleSet& p)
  {
    Nelec=p.getTotalNum();
    Nion = Ions.getTotalNum();
    iGroups=Ions.getSpeciesSet().getTotalNum();
    eGroups=p.groups();

    Uat.resize(Nelec);
    dUat.resize(Nelec);
    FirstAddressOfdU = &(dUat[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dUat.size()*OHMMS_DIM;
    d2Uat.resize(Nelec);
    cur_u.resize(Nelec);
    cur_du.resize(Nelec);
    cur_d2u.resize(Nelec);
    old_u.resize(Nelec);
    old_du.resize(Nelec);
    old_d2u.resize(Nelec);

    F.resize(iGroups,eGroups,eGroups);
    F=nullptr;
    IonDataList.resize(Nion);
  }

  void initUnique()
  {
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    du_dalpha.resize(J3Unique.size());
    dgrad_dalpha.resize(J3Unique.size());
    dhess_dalpha.resize(J3Unique.size());
    int ifunc=0;
    while(it != it_end)
    {
      J3UniqueIndex[it->second]=ifunc;
      FT &functor = *(it->second);
      int numParams = functor.getNumParameters();
      du_dalpha[ifunc].resize(numParams);
      dgrad_dalpha[ifunc].resize(numParams);
      dhess_dalpha[ifunc].resize(numParams);
      ++it;
      ifunc++;
    }
  }

  void addFunc(int iSpecies, int eSpecies1, int eSpecies2, FT* j)
  {
    if(eSpecies1==eSpecies2)
    {
      //if only up-up is specified, assume spin-unpolarized correlations
      if(eSpecies1==0)
        for (int eG1=0; eG1<eGroups; eG1++)
          for (int eG2=0; eG2<eGroups; eG2++)
          {
            if(F(iSpecies,eG1,eG2)==0)
              F(iSpecies,eG1,eG2)=j;
          }
    }
    else
    {
      F(iSpecies,eSpecies1,eSpecies2) = j;
      F(iSpecies,eSpecies2,eSpecies1) = j;
    }
    if(j)
    {
      RealType rcut = 0.5 * j->cutoff_radius;
      for (int i=0; i<Nion; i++)
        if (Ions.GroupID[i] == iSpecies)
          IonDataList[i].cutoff_radius = rcut;
    }
    else
    {
      APP_ABORT("JeeIOrbitalSoA::addFunc  Jastrow function pointer is NULL");
    }
    std::strstream aname;
    aname << iSpecies << "_" << eSpecies1 << "_" << eSpecies2;
    J3Unique[aname.str()]=j;
    initUnique();
  }


  /** check that correlation information is complete
   */
  void check_complete()
  {
    //check that correlation pointers are either all 0 or all assigned
    bool complete = true;
    for(int i=0; i<iGroups; ++i)
    {
      int nfilled = 0;
      bool partial;
      for(int e1=0; e1<eGroups; ++e1)
        for(int e2=0; e2<eGroups; ++e2)
          if(F(i,e1,e2)!=0)
            nfilled++;
      partial = nfilled>0 && nfilled<eGroups*eGroups;
      if(partial)
        app_log() << "J3 eeI is missing correlation for ion "<<i<< std::endl;
      complete = complete && !partial;
    }
    if(!complete)
    {
      APP_ABORT("JeeIOrbitalSoA::check_complete  J3 eeI is missing correlation components\n  see preceding messages for details");
    }
    //first set radii
    for(int i=0; i<Nion; ++i)
    {
      FT* f = F(Ions.GroupID[i],0,0);
      if(f!=0)
        IonDataList[i].cutoff_radius = .5*f->cutoff_radius;
    }
    //then check radii
    bool all_radii_match = true;
    for(int i=0; i<iGroups; ++i)
    {
      if(F(i,0,0)!=0)
      {
        bool radii_match = true;
        RealType rcut = F(i,0,0)->cutoff_radius;
        for(int e1=0; e1<eGroups; ++e1)
          for(int e2=0; e2<eGroups; ++e2)
            radii_match = radii_match && F(i,e1,e2)->cutoff_radius==rcut;
        if(!radii_match)
          app_log() << "eeI functors for ion species " << i
                    << " have different radii"<< std::endl;
        all_radii_match = all_radii_match && radii_match;
      }
    }
    if(!all_radii_match)
    {
      APP_ABORT("JeeIOrbitalSoA::check_radii  J3 eeI are inconsistent for some ion species\n  see preceding messages for details");
    }
  }



  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P) {}

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while(it != it_end)
    {
      (*it).second->checkInVariables(active);
      (*it).second->checkInVariables(myVars);
      ++it;
    }
//       reportStatus(app_log());
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while (it != it_end)
    {
      (*it).second->myVars.getIndex(active);
      myVars.insertFrom((*it).second->myVars);
      ++it;
    }
    myVars.getIndex(active);
    NumVars=myVars.size();
    if (NumVars)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars,Nelec);
      lapLogPsi.resize(NumVars,Nelec);
      VarOffset.resize(Nion, Nelec, Nelec);
      int varoffset=myVars.Index[0];
      for (int i=0; i<Nion; i++)
      {
        if(IonDataList[i].cutoff_radius>0.0)
        {
          for (int j=0; j<Nelec; j++)
            for (int k=0; k<Nelec; k++)
            {
              FT &func_ijk = *F(i, j, k);
              VarOffset(i,j,k).first  = func_ijk.myVars.Index.front()-varoffset;
              VarOffset(i,j,k).second = func_ijk.myVars.Index.size()+VarOffset(i,j,k).first;
            }
        }
      }
    }
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    if(!Optimizable)
      return;
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while(it != it_end)
    {
      (*it++).second->resetParameters(active);
    }
    for(int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if(ii>=0)
        myVars[i]= active[ii];
    }
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while(it != it_end)
    {
      (*it).second->myVars.print(os);
      ++it;
    }
  }


  /**
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@param G returns the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
   *@param L returns the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$
   *@return \f$exp(-J({\bf R}))\f$
   *@brief While evaluating the value of the Jastrow for a set of
   *particles add the gradient and laplacian contribution of the
   *Jastrow to G(radient) and L(aplacian) for local energy calculations
   *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f]
   *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
   *@note The DistanceTableData contains only distinct pairs of the
   *particles belonging to one set, e.g., SymmetricDTD.
   */
  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    evaluateGL(P,G,L,true);
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  ValueType ratio(ParticleSet& P, int iat)
  {
    // construct nearest ion list
    // computeU3(u, gradu, hessu)
    // accumulateG
#if 0
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableID];
    curVal=0.0;
    RealType newval = 0.0;
    RealType oldval = 0.0;
    int ee0 = ee_table->M[iat] - (iat+1);
    for (int i=0; i<Nion; i++)
    {
      IonData &ion = IonDataList[i];
      RealType r_Ii = eI_table->Temp[i].r1;
      int nn0 = eI_table->M[i];
      if (r_Ii < ion.cutoff_radius)
      {
        for (int j=0; j<ion.elecs_inside.size(); j++)
        {
          int jat = ion.elecs_inside[j];
          if (jat != iat)
          {
            RealType r_ij = ee_table->Temp[jat].r1;
            RealType r_Ij = eI_table->r(nn0+jat);
            FT &func = *F(i, iat, jat);
            RealType u = func.evaluate(r_ij, r_Ii, r_Ij);
            curVal[jat] += u;
            newval -= u;
          }
        }
      }
      //if (Fs[i]) curVal += Fs[i]->evaluate(d_table->Temp[i].r1);
    }
    for (int jat=0; jat<Nelec; jat++)
      oldval -= U[iat*Nelec+jat];
    DiffVal = newval - oldval;
    return std::exp(DiffVal);
    //return std::exp(U[iat]-curVal);
    // DiffVal=0.0;
    // const int* pairid(PairID[iat]);
    // for(int jat=0, ij=iat*N; jat<N; jat++,ij++) {
    //   if(jat == iat) {
    //     curVal[jat]=0.0;
    //   } else {
    //     curVal[jat]=F[pairid[jat]]->evaluate(ee_table->Temp[jat].r1);
    //     DiffVal += U[ij]-curVal[jat];
    //     //DiffVal += U[ij]-F[pairid[jat]]->evaluate(ee_table->Temp[jat].r1);
    //   }
    // }
    // return std::exp(DiffVal);
#endif
  }

  //to be removed from QMCPACK: these are not used anymore with PbyPFast
  inline void update(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& dG,
                     ParticleSet::ParticleLaplacian_t& dL,
                     int iat) {}

  ValueType ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
  {return ValueType(1);}

  GradType evalGrad(ParticleSet& P, int iat)
  {
    return GradType(dUat[iat]);
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    // construct nearest ion list
    // computeU3(u, gradu, hessu)
    // accumulateG
#if 0
    curVal  = 0.0;
    curGrad_i = PosType();
    curLap_i  = 0.0;
    curGrad_j = PosType();
    curLap_j = 0.0;
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableID];
    int ee0 = ee_table->M[iat] - (iat+1);
    DiffVal = 0.0;
    for (int i=0; i<Nion; i++)
    {
      IonData &ion = IonDataList[i];
      RealType r_Ii     = eI_table->Temp[i].r1;
      RealType r_Ii_inv = 1.0/r_Ii;
      int nn0 = eI_table->M[i];
      if (r_Ii < ion.cutoff_radius)
      {
        for (int j=0; j<ion.elecs_inside.size(); j++)
        {
          int jat = ion.elecs_inside[j];
          if (jat != iat)
          {
            RealType r_ij = ee_table->Temp[jat].r1;
            RealType r_ij_inv = 1.0/r_ij;
            RealType r_Ij = eI_table->r(nn0+jat);
            RealType r_Ij_inv = 1.0/r_Ij;
            FT &func = *F(i, iat, jat);
            PosType gradF;
            Tensor<RealType,OHMMS_DIM> hessF;
            RealType u = func.evaluate(r_ij, r_Ii, r_Ij, gradF, hessF);
            PosType gr_ee =   -gradF[0]*r_ij_inv * ee_table->Temp[jat].dr1;
            PosType du_i, du_j;
            RealType d2u_i, d2u_j;
            du_i = gradF[1]*r_Ii_inv * eI_table->Temp[i].dr1 - gr_ee;
            du_j = gradF[2]*r_Ij_inv * eI_table->dr(nn0+jat) + gr_ee;
            d2u_i = (hessF(0,0) + 2.0*r_ij_inv*gradF[0] + 2.0*hessF(0,1) *
                     dot(ee_table->Temp[jat].dr1,
                         eI_table->Temp[i].dr1)*r_ij_inv*r_Ii_inv
                     + hessF(1,1) + 2.0*r_Ii_inv*gradF[1]);
            d2u_j = (hessF(0,0) + 2.0*r_ij_inv*gradF[0] - 2.0*hessF(0,2) *
                     dot(ee_table->Temp[jat].dr1,
                         eI_table->dr(nn0+jat))*r_ij_inv*r_Ij_inv
                     + hessF(2,2) + 2.0*r_Ij_inv*gradF[2]);
            curVal   [jat] += u;
            curGrad_j[jat] += du_j;
            curLap_j [jat] += d2u_j;
            curGrad_i[jat] += du_i;
            curLap_i [jat] += d2u_i;
            DiffVal -=   u;
          }
        }
      }
    }
    for (int jat=0; jat<Nelec; jat++)
    {
      if (iat != jat)
      {
        int ij = iat*Nelec+jat;
        DiffVal +=   U[ij];
        grad_iat -= curGrad_i[jat];
      }
    }
    return std::exp(DiffVal);
#endif
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
#if 0
    const DistanceTableData* eI_table=P.DistTables[myTableID];
    //      std::cerr << "acceptMove called.\n";
    DiffValSum += DiffVal;
    for(int jat=0,ij=iat*Nelec,ji=iat; jat<Nelec; jat++,ij++,ji+=Nelec)
      if (jat != iat)
      {
        dU[ij]  = curGrad_i[jat];
        dU[ji]  = curGrad_j[jat];
        d2U[ij] = curLap_i[jat];
        d2U[ji] = curLap_j[jat];
        U[ij] =  U[ji] = curVal[jat];
      }
    // Now, update elecs_inside for each ion
    for (int i=0; i < IonDataList.size(); i++)
    {
      IonData &ion = IonDataList[i];
      bool inside = eI_table->Temp[i].r1 < ion.cutoff_radius;
      IonData::eListType::iterator iter;
      iter = find(ion.elecs_inside.begin(),
                  ion.elecs_inside.end(), iat);
      if (inside && iter == ion.elecs_inside.end())
        ion.elecs_inside.push_back(iat);
      else
        if (!inside && iter != ion.elecs_inside.end())
          ion.elecs_inside.erase(iter);
    }
#endif
  }

  /** intenal function to compute \f$\sum_j u(r_j), du/dr, d2u/dr2\f$ */
  inline void computeU3(ParticleSet& P, int jel, const RealType* restrict dist,
    RealType* restrict u, RealType* restrict du, RealType* restrict d2u)
  {
  }

  inline void recompute(ParticleSet& P)
  {
    constexpr RealType cone(1);
    constexpr RealType cminus(-1);
    constexpr RealType ctwo(2);
    const DistanceTableData& ee_table=(*P.DistTables[0]);
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);

    // First build a neighbour list for each ion
    for (int iat=0; iat<Nion; iat++)
    {
      IonDataCompact &ion = IonDataList[iat];
      ion.elecs_inside.clear();
      ion.elecs_dist.clear();
      if (ion.cutoff_radius>0)
        for (int jel=0; jel<Nelec; jel++)
        {
          const RealType rij=eI_table.Distances[jel][iat];
          if (rij < ion.cutoff_radius)
          {
            ion.elecs_inside.push_back(jel);
            ion.elecs_dist.push_back(rij);
          }
        }
    }

    for(int jel=0; jel<Nelec; ++jel)
    {
      const int jg=P.GroupID[jel];
      Uat[jel]   = 0;
      dUat[jel]  = PosType();
      d2Uat[jel] = 0;
      for(int iat=0; iat<Nion; ++iat)
      {
        if(eI_table.Distances[jel][iat]<IonDataList[iat].cutoff_radius)
        {
          const int ig=Ions.GroupID[iat];
          IonDataCompact &ion = IonDataList[iat];
          const RealType r_Ij     = eI_table.Distances[jel][iat];
          const PosType disp_Ij   = cminus*eI_table.Displacements[jel][iat];
          const RealType r_Ij_inv = cone/r_Ij;

          for(int kelid=0; kelid<ion.elecs_inside.size(); kelid++)
          {
            const int kel=ion.elecs_inside[kelid];
            if(jel==kel) continue;
            const int kg=P.GroupID[kel];
            const FT& feeI(*F(ig,jg,kg));
            const RealType r_Ik     = ion.elecs_dist[kelid];
            const RealType r_jk     = ee_table.Distances[jel][kel];
            const PosType disp_jk   = ee_table.Displacements[jel][kel];
            const RealType r_jk_inv = cone/r_jk;

            // compute the contribution to jel
            TinyVector<valT,3> gradF;
            Tensor<valT,3> hessF;

            valT u = feeI.evaluate(r_jk, r_Ij, r_Ik, gradF, hessF);
            // sign is flipped in gradient
            PosType du_j = gradF[0]*r_jk_inv * disp_jk - gradF[1]*r_Ij_inv * disp_Ij;
            RealType d2u_j = hessF(0,0) + hessF(1,1)
                           + ctwo*r_jk_inv*gradF[0] + ctwo*r_Ij_inv*gradF[1]
                           - ctwo*hessF(0,1)*dot(disp_jk,disp_Ij)*r_jk_inv*r_Ij_inv;

            Uat[jel] += u;
            dUat[jel] += du_j;
            d2Uat[jel] -= d2u_j;
          }
        }
      }
    }
  }

  inline RealType registerData(ParticleSet& P, PooledData<RealType>& buf)
  {
    evaluateLog(P,P.G,P.L);
    buf.add(Uat.begin(), Uat.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
    buf.add(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  inline RealType updateBuffer(ParticleSet& P, PooledData<RealType>& buf,
                               bool fromscratch=false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.put(Uat.begin(), Uat.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    buf.put(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.get(Uat.begin(), Uat.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    buf.get(d2Uat.begin(), d2Uat.end());
  }

  inline RealType evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
  {
    buf.put(Uat.begin(), Uat.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    buf.put(d2Uat.begin(), d2Uat.end());
    return LogValue;
  }

  void evaluateGL(ParticleSet& P,
             ParticleSet::ParticleGradient_t& G,
             ParticleSet::ParticleLaplacian_t& L,
             bool fromscratch)
  {
    if(fromscratch) recompute(P);
    LogValue=valT(0);
    for(int iat=0; iat<Nelec; ++iat)
    {
      LogValue += Uat[iat];
      G[iat] += dUat[iat];
      L[iat] += d2Uat[iat];
    }

    constexpr valT mhalf(-0.5);
    LogValue=mhalf*LogValue;
  }

  void evaluateDerivatives(ParticleSet& P,
                           const opt_variables_type& optvars,
                           std::vector<RealType>& dlogpsi,
                           std::vector<RealType>& dhpsioverpsi)
  {
#if 0
    bool recalculate(false);
    std::vector<bool> rcsingles(myVars.size(),false);
    for (int k=0; k<myVars.size(); ++k)
    {
      int kk=myVars.where(k);
      if (kk<0)
        continue;
      if (optvars.recompute(kk))
        recalculate=true;
      rcsingles[k]=true;
    }
    if (recalculate)
    {
      const DistanceTableData* ee_table=P.DistTables[0];
      const DistanceTableData* eI_table=P.DistTables[myTableID];
      // First, create lists of electrons within the sphere of each ion
      for (int i=0; i<Nion; i++)
      {
        IonData &ion = IonDataList[i];
        ion.elecs_inside.clear();
        int iel=0;
        if (ion.cutoff_radius > 0.0)
          for (int nn=eI_table->M[i]; nn<eI_table->M[i+1]; nn++, iel++)
            if (eI_table->r(nn) < ion.cutoff_radius)
              ion.elecs_inside.push_back(iel);
      }
      dLogPsi=0.0;
      gradLogPsi = PosType();
      lapLogPsi = 0.0;
      RealType u;
      PosType gradF;
      Tensor<RealType,3> hessF;
      // Now, evaluate three-body term for each ion
      for (int i=0; i<Nion; i++)
      {
        IonData &ion = IonDataList[i];
        int nn0 = eI_table->M[i];
        for (int j=0; j<ion.elecs_inside.size(); j++)
        {
          int jel = ion.elecs_inside[j];
          RealType r_Ij     = eI_table->r(nn0+jel);
          RealType r_Ij_inv = eI_table->rinv(nn0+jel);
          int ee0 = ee_table->M[jel]-(jel+1);
          for (int k=j+1; k<ion.elecs_inside.size(); k++)
          {
            int kel = ion.elecs_inside[k];
            RealType r_Ik     = eI_table->r(nn0+kel);
            RealType r_Ik_inv = eI_table->rinv(nn0+kel);
            RealType r_jk     = ee_table->r(ee0+kel);
            RealType r_jk_inv = ee_table->rinv(ee0+kel);
            FT &func = *F(i, jel, kel);
            int idx = J3UniqueIndex[F(i, jel, kel)];
            func.evaluateDerivatives(r_jk, r_Ij, r_Ik, du_dalpha[idx],
                                     dgrad_dalpha[idx], dhess_dalpha[idx]);
            u = func.evaluate (r_jk, r_Ij, r_Ik, gradF, hessF);
            LogValue -= u;
            // Save for ratio
            U[jel*Nelec+kel] += u;
            U[kel*Nelec+jel] += u;
            int first = VarOffset(i,jel,kel).first;
            int last  = VarOffset(i,jel,kel).second;
            std::vector<RealType> &dlog = du_dalpha[idx];
            std::vector<PosType>  &dgrad = dgrad_dalpha[idx];
            std::vector<Tensor<RealType,3> > &dhess = dhess_dalpha[idx];
            for (int p=first,ip=0; p<last; p++,ip++)
            {
              RealType dval =  dlog[ip];
              PosType dg  = dgrad[ip];
              Tensor<RealType,3> dh  = dhess[ip];
              PosType gr_ee =    dg[0]*r_jk_inv * ee_table->dr(ee0+kel);
              PosType du_j, du_k;
              RealType d2u_j, d2u_k;
              du_j = dg[1]*r_Ij_inv * eI_table->dr(nn0+jel) - gr_ee;
              du_k = dg[2]*r_Ik_inv * eI_table->dr(nn0+kel) + gr_ee;
              d2u_j = (dh(0,0) + 2.0*r_jk_inv*dg[0] -
                       2.0*dh(0,1)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+jel))*r_jk_inv*r_Ij_inv
                       + dh(1,1) + 2.0*r_Ij_inv*dg[1]);
              d2u_k = (dh(0,0) + 2.0*r_jk_inv*dg[0] +
                       2.0*dh(0,2)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+kel))*r_jk_inv*r_Ik_inv
                       + dh(2,2) + 2.0*r_Ik_inv*dg[2]);
              dLogPsi[p] -= dval;
              gradLogPsi(p,jel) -= du_j;
              gradLogPsi(p,kel) -= du_k;
              lapLogPsi(p,jel)  -= d2u_j;
              lapLogPsi(p,kel)  -= d2u_k;
            }
          }
        }
      }
      for (int k=0; k<myVars.size(); ++k)
      {
        int kk=myVars.where(k);
        if (kk<0)
          continue;
        dlogpsi[kk]=dLogPsi[k];
        RealType sum = 0.0;
        for (int i=0; i<Nelec; i++)
        {
#if defined(QMC_COMPLEX)
          sum -= 0.5*lapLogPsi(k,i);
          for(int jdim=0; jdim<OHMMS_DIM; ++jdim)
            sum -= P.G[i][jdim].real()*gradLogPsi(k,i)[jdim];
#else
          sum -= 0.5*lapLogPsi(k,i) + dot(P.G[i], gradLogPsi(k,i));
#endif
        }
        dhpsioverpsi[kk] = sum;
      }
    }
#endif
  }
};

}
#endif
