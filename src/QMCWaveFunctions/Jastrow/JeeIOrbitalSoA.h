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
  ///number of particles + padded
  size_t Nelec_padded;
  //number of groups of the target particleset
  int eGroups, iGroups;
  ///reference to the sources (ions)
  const ParticleSet& Ions;
  ///diff value
  RealType DiffVal;

  ///\f$Uat[i] = sum_(j) u_{i,j}\f$
  Vector<valT> Uat,oldUk,newUk;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  using gContainer_type=VectorSoaContainer<valT,OHMMS_DIM>;
  gContainer_type dUat,olddUk,newdUk;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  Vector<valT> d2Uat,oldd2Uk,newd2Uk;
  /// current values during PbyP
  valT cur_Uat,cur_d2Uat;
  posT cur_dUat, dUat_temp;
  ///container for the Jastrow functions
  Array<FT*,3> F;

  std::map<std::string,FT*> J3Unique;
  //YYYY
  std::map<FT*,int> J3UniqueIndex;

  /// the cutoff for e-I pairs
  std::vector<valT> Ion_cutoff;
  /// the electrons around ions within the cutoff radius, grouped by species
  Array<std::vector<int>,2> elecs_inside;
  Array<std::vector<valT>,2> elecs_inside_dist;

  /// compressed distances
  aligned_vector<valT> Distjk_Compressed, DistkI_Compressed;
  std::vector<int> DistIndice;

  VectorSoaContainer<valT,9> mVGL;

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

  JeeIOrbitalSoA(const ParticleSet& ions, ParticleSet& elecs, bool is_master=false)
    : Ions(ions), NumVars(0)
  {
    OrbitalName = "JeeIOrbitalSoA";
    myTableID=elecs.addTable(Ions,DT_SOA);
    elecs.DistTables[myTableID]->Need_full_table_loadWalker=true;
    init(elecs);
  }

  ~JeeIOrbitalSoA() { }

  OrbitalBasePtr makeClone(ParticleSet& elecs) const
  {
    JeeIOrbitalSoA<FT>* eeIcopy= new JeeIOrbitalSoA<FT>(Ions, elecs, false);
    std::map<const FT*,FT*> fcmap;
    for (int iG=0; iG<iGroups; iG++)
      for (int eG1=0; eG1<eGroups; eG1++)
        for (int eG2=0; eG2<eGroups; eG2++)
        {
          if(F(iG,eG1,eG2)==0)
            continue;
          typename std::map<const FT*,FT*>::iterator fit=fcmap.find(F(iG,eG1,eG2));
          if(fit == fcmap.end())
          {
            FT* fc=new FT(*F(iG,eG1,eG2));
            eeIcopy->addFunc(iG, eG1, eG2, fc);
            fcmap[F(iG,eG1,eG2)]=fc;
          }
        }
    // Ye: I don't like the following memory allocated by default.
    eeIcopy->myVars.clear();
    eeIcopy->myVars.insertFrom(myVars);
    eeIcopy->NumVars=NumVars;
    eeIcopy->dLogPsi.resize(NumVars);
    eeIcopy->gradLogPsi.resize(NumVars,Nelec);
    eeIcopy->lapLogPsi.resize(NumVars,Nelec);
    eeIcopy->VarOffset=VarOffset;
    eeIcopy->Optimizable = Optimizable;
    return eeIcopy;
  }

  void init(ParticleSet& p)
  {
    Nelec=p.getTotalNum();
    Nelec_padded=getAlignedSize<valT>(Nelec);
    Nion = Ions.getTotalNum();
    iGroups=Ions.getSpeciesSet().getTotalNum();
    eGroups=p.groups();

    Uat.resize(Nelec);
    dUat.resize(Nelec);
    d2Uat.resize(Nelec);

    oldUk.resize(Nelec);
    olddUk.resize(Nelec);
    oldd2Uk.resize(Nelec);
    newUk.resize(Nelec);
    newdUk.resize(Nelec);
    newd2Uk.resize(Nelec);

    F.resize(iGroups,eGroups,eGroups);
    F=nullptr;
    elecs_inside.resize(Nion,eGroups);
    elecs_inside_dist.resize(Nion,eGroups);
    Ion_cutoff.resize(Nion, 0.0);

    mVGL.resize(Nelec);
    DistkI_Compressed.resize(Nelec);
    Distjk_Compressed.resize(Nelec);
    DistIndice.resize(Nelec);
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
          Ion_cutoff[i] = rcut;
    }
    else
    {
      APP_ABORT("JeeIOrbitalSoA::addFunc  Jastrow function pointer is NULL");
    }
    std::stringstream aname;
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
        Ion_cutoff[i] = .5*f->cutoff_radius;
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
      VarOffset.resize(iGroups, eGroups, eGroups);
      int varoffset=myVars.Index[0];
      for (int ig=0; ig<iGroups; ig++)
        for (int jg=0; jg<eGroups; jg++)
          for (int kg=0; kg<eGroups; kg++)
          {
            FT *func_ijk = F(ig, jg, kg);
            if(func_ijk==nullptr) continue;
            VarOffset(ig,jg,kg).first  = func_ijk->myVars.Index.front()-varoffset;
            VarOffset(ig,jg,kg).second = func_ijk->myVars.Index.size()+VarOffset(ig,jg,kg).first;
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

  void build_compact_list(ParticleSet& P)
  {
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);

    for(int iat=0; iat<Nion; ++iat)
      for(int jg=0; jg<eGroups; ++jg)
      {
        elecs_inside(iat,jg).clear();
        elecs_inside_dist(iat,jg).clear();
      }

    for(int jg=0; jg<eGroups; ++jg)
      for(int jel=P.first(jg); jel<P.last(jg); jel++)
        for(int iat=0; iat<Nion; ++iat)
          if(eI_table.Distances[jel][iat]<Ion_cutoff[iat])
          {
            elecs_inside(iat,jg).push_back(jel);
            elecs_inside_dist(iat,jg).push_back(eI_table.Distances[jel][iat]);
          }
  }

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
    UpdateMode=ORB_PBYP_RATIO;

    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);
    cur_Uat=computeU(P, iat, P.GroupID[iat], eI_table.Temp_r.data(), ee_table.Temp_r.data());
    DiffVal=Uat[iat]-cur_Uat;
    return std::exp(DiffVal);
  }

  void evaluateRatiosAlltoOne(ParticleSet& P, std::vector<ValueType>& ratios)
  {
    const DistanceTableData* d_table=P.DistTables[0];
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);

    for(int jg=0; jg<eGroups; ++jg)
    {
      const valT sumU=computeU(P, -1, jg, eI_table.Temp_r.data(), ee_table.Temp_r.data());

      for(int j=P.first(jg); j<P.last(jg); ++j)
      {
        // remove self-interaction
        valT Uself(0);
        for(int iat=0; iat<Nion; ++iat)
        {
          const valT &r_Ij = eI_table.Temp_r[iat];
          const valT &r_Ik = eI_table.Distances[j][iat];
          if(r_Ij<Ion_cutoff[iat]&&r_Ik<Ion_cutoff[iat])
          {
            const int ig=Ions.GroupID[iat];
            Uself+=F(ig,jg,jg)->evaluate(ee_table.Temp_r[j],r_Ij,r_Ik);
          }
        }
        ratios[j]=std::exp(Uat[j]+Uself-sumU);
      }
    }
  }

  GradType evalGrad(ParticleSet& P, int iat)
  {
    return GradType(dUat[iat]);
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    UpdateMode=ORB_PBYP_PARTIAL;

    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);
    computeU3(P, iat, eI_table.Temp_r.data(), eI_table.Temp_dr, ee_table.Temp_r.data(), ee_table.Temp_dr,
              cur_Uat, cur_dUat, cur_d2Uat, newUk, newdUk, newd2Uk);
    DiffVal=Uat[iat]-cur_Uat;
    grad_iat+=cur_dUat;
    return std::exp(DiffVal);
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);
    // get the old value, grad, lapl
    computeU3(P, iat, eI_table.Distances[iat], eI_table.Displacements[iat], ee_table.Distances[iat], ee_table.Displacements[iat],
              Uat[iat], dUat_temp, d2Uat[iat], oldUk, olddUk, oldd2Uk);
    if(UpdateMode == ORB_PBYP_RATIO)
    {//ratio-only during the move; need to compute derivatives
      computeU3(P, iat, eI_table.Temp_r.data(), eI_table.Temp_dr, ee_table.Temp_r.data(), ee_table.Temp_dr,
                cur_Uat, cur_dUat, cur_d2Uat, newUk, newdUk, newd2Uk);
    }

    #pragma omp simd
    for(int jel=0; jel<Nelec; jel++)
    {
      Uat[jel]   += newUk[jel]-oldUk[jel];
      d2Uat[jel] += newd2Uk[jel]-oldd2Uk[jel];
    }
    for(int idim=0; idim<OHMMS_DIM; ++idim)
    {
      valT* restrict save_g=dUat.data(idim);
      const valT* restrict new_g=newdUk.data(idim);
      const valT* restrict old_g=olddUk.data(idim);
      #pragma omp simd aligned(save_g,new_g,old_g)
      for(int jel=0; jel<Nelec; jel++)
        save_g[jel]+=new_g[jel]-old_g[jel];
    }

    Uat[iat]   = cur_Uat;
    dUat(iat)  = cur_dUat;
    d2Uat[iat] = cur_d2Uat;

    const int ig = P.GroupID[iat];
    // update compact list elecs_inside
    for (int jat=0; jat < Nion; jat++)
    {
      bool inside = eI_table.Temp_r[jat] < Ion_cutoff[jat];
      std::vector<int>::iterator iter = find(elecs_inside(jat,ig).begin(), elecs_inside(jat,ig).end(), iat);
      std::vector<RealType>::iterator iter_dist = elecs_inside_dist(jat,ig).begin()+std::distance(elecs_inside(jat,ig).begin(),iter);
      if(inside)
      {
        if(iter==elecs_inside(jat,ig).end())
        {
          elecs_inside(jat,ig).push_back(iat);
          elecs_inside_dist(jat,ig).push_back(eI_table.Temp_r[jat]);
        }
        else
        {
          *iter_dist = eI_table.Temp_r[jat];
        }
      }
      else
      {
        if(iter!=elecs_inside(jat,ig).end())
        {
          elecs_inside(jat,ig).erase(iter);
          elecs_inside_dist(jat,ig).erase(iter_dist);
        }
      }
    }
  }

  inline void recompute(ParticleSet& P)
  {
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);

    build_compact_list(P);

    for(int jel=0; jel<Nelec; ++jel)
    {
      computeU3(P, jel, eI_table.Distances[jel], eI_table.Displacements[jel], ee_table.Distances[jel], ee_table.Displacements[jel],
                Uat[jel], dUat_temp, d2Uat[jel], newUk, newdUk, newd2Uk, true);
      dUat(jel) = dUat_temp;
      // add the contribution from the upper triangle
      #pragma omp simd
      for(int kel=0; kel<jel; kel++)
      {
        Uat[kel] += newUk[kel];
        d2Uat[kel] += newd2Uk[kel];
      }
      for(int idim=0; idim<OHMMS_DIM; ++idim)
      {
        valT* restrict save_g=dUat.data(idim);
        const valT* restrict new_g=newdUk.data(idim);
        #pragma omp simd aligned(save_g,new_g)
        for(int kel=0; kel<jel; kel++)
          save_g[kel]+=new_g[kel];
      }
    }
  }

  inline valT computeU(ParticleSet& P, int jel, int jg,
                        const RealType* distjI, const RealType* distjk)
  {
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);

    valT Uj = valT(0);
    for(int iat=0; iat<Nion; ++iat)
      if(distjI[iat]<Ion_cutoff[iat])
      {
        const int ig=Ions.GroupID[iat];
        const valT r_Ij     = distjI[iat];

        for(int kg=0; kg<eGroups; ++kg)
        {
          const FT& feeI(*F(ig,jg,kg));
          int kel_counter=0;
          for(int kind=0; kind<elecs_inside(iat,kg).size(); kind++)
          {
            const int kel=elecs_inside(iat,kg)[kind];
            if(kel!=jel)
            {
              DistkI_Compressed[kel_counter]=elecs_inside_dist(iat,kg)[kind];
              Distjk_Compressed[kel_counter]=distjk[kel];
              kel_counter++;
            }
          }
          Uj += feeI.evaluateV(kel_counter, Distjk_Compressed.data(), r_Ij, DistkI_Compressed.data());
        }
      }
    return Uj;
  }

  inline void computeU3(ParticleSet& P, int jel,
                        const RealType* distjI, const RowContainer& displjI,
                        const RealType* distjk, const RowContainer& displjk,
                        valT& Uj, posT& dUj, valT& d2Uj,
                        Vector<valT>& Uk, gContainer_type& dUk, Vector<valT>& d2Uk, bool triangle=false)
  {
    constexpr valT czero(0);
    constexpr valT cone(1);
    constexpr valT cminus(-1);
    constexpr valT ctwo(2);
    constexpr valT lapfac=OHMMS_DIM-cone;
    Uj = czero;
    dUj = posT();
    d2Uj = czero;

    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const int jg=P.GroupID[jel];

    const int kelmax=triangle?jel:Nelec;
    std::fill_n(Uk.data(),kelmax,czero);
    std::fill_n(d2Uk.data(),kelmax,czero);
    for(int idim=0; idim<OHMMS_DIM; ++idim)
      std::fill_n(dUk.data(idim),kelmax,czero);

    valT* restrict val=mVGL.data(0);
    valT* restrict gradF0=mVGL.data(1);
    valT* restrict gradF1=mVGL.data(2);
    valT* restrict gradF2=mVGL.data(3);
    valT* restrict hessF00=mVGL.data(4);
    valT* restrict hessF11=mVGL.data(5);
    valT* restrict hessF22=mVGL.data(6);
    valT* restrict hessF01=mVGL.data(7);
    valT* restrict hessF02=mVGL.data(8);

    for(int iat=0; iat<Nion; ++iat)
      if(distjI[iat]<Ion_cutoff[iat])
      {
        const int ig=Ions.GroupID[iat];
        const valT r_Ij     = distjI[iat];
        const posT disp_Ij  = cminus*displjI[iat];

        for(int kg=0; kg<eGroups; ++kg)
        {
          const FT& feeI(*F(ig,jg,kg));
          int kel_counter=0;
          for(int kind=0; kind<elecs_inside(iat,kg).size(); kind++)
          {
            const int kel=elecs_inside(iat,kg)[kind];
            if(kel<kelmax && kel!=jel)
            {
              DistkI_Compressed[kel_counter]=elecs_inside_dist(iat,kg)[kind];
              Distjk_Compressed[kel_counter]=distjk[kel];
              DistIndice[kel_counter]=kel;
              kel_counter++;
            }
          }

          feeI.evaluateVGL(kel_counter, Distjk_Compressed.data(), r_Ij, DistkI_Compressed.data(),
                           val, gradF0, gradF1, gradF2, hessF00, hessF11, hessF22, hessF01, hessF02);

          for(int kel_index=0; kel_index<kel_counter; kel_index++)
          {
            int kel=DistIndice[kel_index];
            const posT disp_Ik  = cminus*eI_table.Displacements[kel][iat];
            const posT disp_jk  = displjk[kel];

            // compute the contribution to jel
            Uj += val[kel_index];
            dUj += gradF0[kel_index] * disp_jk - gradF1[kel_index] * disp_Ij;
            d2Uj -= hessF00[kel_index] + hessF11[kel_index]
                    + lapfac*(gradF0[kel_index] + gradF1[kel_index])
                    - ctwo*hessF01[kel_index]*dot(disp_jk,disp_Ij);

            // compute the contribution to kel
            Uk[kel] += val[kel_index];
            dUk(kel) = dUk[kel] - gradF0[kel_index] * disp_jk - gradF2[kel_index] * disp_Ik;
            d2Uk[kel] -= hessF00[kel_index] + hessF22[kel_index]
                         + lapfac*(gradF0[kel_index] + gradF2[kel_index])
                         + ctwo*hessF02[kel_index]*dot(disp_jk,disp_Ik);
          }
        }
      }
  }

  inline void registerData(ParticleSet& P, WFBufferType& buf)
  {
    if ( Bytes_in_WFBuffer == 0 )
    {
      Bytes_in_WFBuffer = buf.current();
      buf.add(Uat.begin(), Uat.end());
      buf.add(dUat.data(), dUat.end());
      buf.add(d2Uat.begin(), d2Uat.end());
      Bytes_in_WFBuffer = buf.current()-Bytes_in_WFBuffer;
      // free local space
      Uat.free();
      dUat.free();
      d2Uat.free();
    }
    else
    {
      buf.forward(Bytes_in_WFBuffer);
    }
  }

  inline RealType updateBuffer(ParticleSet& P, WFBufferType& buf,
                               bool fromscratch=false)
  {
    evaluateGL(P, P.G, P.L, false);
    buf.forward(Bytes_in_WFBuffer);
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    Uat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    dUat.attachReference(Nelec, Nelec_padded, buf.lendReference<valT>(Nelec_padded*OHMMS_DIM));
    d2Uat.attachReference(buf.lendReference<valT>(Nelec), Nelec);
    build_compact_list(P);
  }

  void evaluateGL(ParticleSet& P,
             ParticleSet::ParticleGradient_t& G,
             ParticleSet::ParticleLaplacian_t& L,
             bool fromscratch=false)
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
      constexpr valT czero(0);
      constexpr valT cone(1);
      constexpr valT cminus(-1);
      constexpr valT ctwo(2);
      constexpr valT lapfac=OHMMS_DIM-cone;

      const DistanceTableData& ee_table=(*P.DistTables[0]);
      const DistanceTableData& eI_table=(*P.DistTables[myTableID]);

      build_compact_list(P);

      dLogPsi=czero;
      gradLogPsi = PosType();
      lapLogPsi = czero;

      for(int iat=0; iat<Nion; ++iat)
      {
        const int ig=Ions.GroupID[iat];
        for(int jg=0; jg<eGroups; ++jg)
          for(int jind=0; jind<elecs_inside(iat,jg).size(); jind++)
          {
            const int jel=elecs_inside(iat,jg)[jind];
            const valT r_Ij     = elecs_inside_dist(iat,jg)[jind];
            const posT disp_Ij  = cminus*eI_table.Displacements[jel][iat];
            const valT r_Ij_inv = cone/r_Ij;

            for(int kg=0; kg<eGroups; ++kg)
              for(int kind=0; kind<elecs_inside(iat,kg).size(); kind++)
              {
                const int kel=elecs_inside(iat,kg)[kind];
                if(kel<jel)
                {
                  const FT& feeI(*F(ig,jg,kg));

                  const valT r_Ik     = elecs_inside_dist(iat,kg)[kind];
                  const posT disp_Ik  = cminus*eI_table.Displacements[kel][iat];
                  const valT r_Ik_inv = cone/r_Ik;

                  const valT r_jk     = ee_table.Distances[jel][kel];
                  const posT disp_jk  = ee_table.Displacements[jel][kel];
                  const valT r_jk_inv = cone/r_jk;

                  FT &func = *F(ig, jg, kg);
                  int idx = J3UniqueIndex[F(ig, jg, kg)];
                  func.evaluateDerivatives(r_jk, r_Ij, r_Ik, du_dalpha[idx],
                                       dgrad_dalpha[idx], dhess_dalpha[idx]);
                  int first = VarOffset(ig,jg,kg).first;
                  int last  = VarOffset(ig,jg,kg).second;
                  std::vector<RealType> &dlog = du_dalpha[idx];
                  std::vector<PosType>  &dgrad = dgrad_dalpha[idx];
                  std::vector<Tensor<RealType,3> > &dhess = dhess_dalpha[idx];

                  for (int p=first,ip=0; p<last; p++,ip++)
                  {
                    RealType& dval = dlog[ip];
                    PosType& dg = dgrad[ip];
                    Tensor<RealType,3>& dh = dhess[ip];

                    dg[0]*=r_jk_inv;
                    dg[1]*=r_Ij_inv;
                    dg[2]*=r_Ik_inv;

                    PosType gr_ee = dg[0] * disp_jk;

                    gradLogPsi(p,jel) -= dg[1] * disp_Ij - gr_ee;
                    lapLogPsi(p,jel)  -= (dh(0,0) + lapfac*dg[0] -
                         ctwo*dh(0,1)*dot(disp_jk,disp_Ij)*r_jk_inv*r_Ij_inv
                         + dh(1,1) + lapfac*dg[1]);

                    gradLogPsi(p,kel) -= dg[2] * disp_Ik + gr_ee;
                    lapLogPsi(p,kel)  -= (dh(0,0) + lapfac*dg[0] +
                         ctwo*dh(0,2)*dot(disp_jk,disp_Ik)*r_jk_inv*r_Ik_inv
                         + dh(2,2) + lapfac*dg[2]);

                    dLogPsi[p] -= dval;
                  }
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
  }
};

}
#endif
