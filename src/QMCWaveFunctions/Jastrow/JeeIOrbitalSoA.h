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
  ParticleAttrib<valT> Uat,oldUk,newUk;
  ///\f$dUat[i] = sum_(j) du_{i,j}\f$
  ParticleAttrib<posT> dUat,olddUk,newdUk;
  valT *FirstAddressOfdU, *LastAddressOfdU;
  ///\f$d2Uat[i] = sum_(j) d2u_{i,j}\f$
  ParticleAttrib<valT> d2Uat,oldd2Uk,newd2Uk;
  valT cur_Uat,cur_d2Uat;
  posT cur_dUat;
  ///container for the Jastrow functions
  Array<FT*,3> F;

  std::map<std::string,FT*> J3Unique;
  //YYYY
  std::map<FT*,int> J3UniqueIndex;

  std::vector<valT> Ion_cutoff;

  /// compressed distances
  aligned_vector<valT> Distjk_Compressed, DistkI_Compressed;
  std::vector<int> DistIndice;

  using VGL_type=VectorSoaContainer<valT,9>;
  VGL_type mVGL;

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
#if 0
    eeIcopy->myVars.clear();
    eeIcopy->myVars.insertFrom(myVars);
    eeIcopy->NumVars=NumVars;
    eeIcopy->dLogPsi.resize(NumVars);
    eeIcopy->gradLogPsi.resize(NumVars,Nelec);
    eeIcopy->lapLogPsi.resize(NumVars,Nelec);
    eeIcopy->VarOffset=VarOffset;
    eeIcopy->Optimizable = Optimizable;
#endif
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

    oldUk.resize(Nelec);
    olddUk.resize(Nelec);
    oldd2Uk.resize(Nelec);
    newUk.resize(Nelec);
    newdUk.resize(Nelec);
    newd2Uk.resize(Nelec);

    F.resize(iGroups,eGroups,eGroups);
    F=nullptr;
    Ion_cutoff.resize(Nion);

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
#if 0
    myVars.clear();
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while(it != it_end)
    {
      (*it).second->checkInVariables(active);
      (*it).second->checkInVariables(myVars);
      ++it;
    }
#endif
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
#if 0
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
        if(Ion_cutoff[i]>0.0)
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
#endif
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
#if 0
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
#endif
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
#if 0
    typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    while(it != it_end)
    {
      (*it).second->myVars.print(os);
      ++it;
    }
#endif
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
    cur_Uat=computeU(P, iat, eI_table.Temp_r.data(), ee_table.Temp_r.data());
    DiffVal=Uat[iat]-cur_Uat;
    return std::exp(DiffVal);
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
              Uat[iat], dUat[iat], d2Uat[iat], oldUk, olddUk, oldd2Uk);
    if(UpdateMode == ORB_PBYP_RATIO)
    {//ratio-only during the move; need to compute derivatives
      computeU3(P, iat, eI_table.Temp_r.data(), eI_table.Temp_dr, ee_table.Temp_r.data(), ee_table.Temp_dr,
                cur_Uat, cur_dUat, cur_d2Uat, newUk, newdUk, newd2Uk);
    }

    for(int jel=0; jel<Nelec; jel++)
    {
      Uat[jel]   += newUk[jel]-oldUk[jel];
      dUat[jel]  += newdUk[jel]-olddUk[jel];
      d2Uat[jel] += newd2Uk[jel]-oldd2Uk[jel];
    }

    Uat[iat]   = cur_Uat;
    dUat[iat]  = cur_dUat;
    d2Uat[iat] = cur_d2Uat;
  }

  inline void recompute(ParticleSet& P)
  {
    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const DistanceTableData& ee_table=(*P.DistTables[0]);
    for(int jel=0; jel<Nelec; ++jel)
    {
      computeU3(P, jel, eI_table.Distances[jel], eI_table.Displacements[jel], ee_table.Distances[jel], ee_table.Displacements[jel],
                Uat[jel], dUat[jel], d2Uat[jel], newUk, newdUk, newd2Uk, true);
      // add the contribution from the upper triangle
      for(int kel=0; kel<jel; kel++)
      {
        Uat[kel] += newUk[kel];
        dUat[kel] += newdUk[kel];
        d2Uat[kel] += newd2Uk[kel];
      }
    }
  }

  inline valT computeU(ParticleSet& P, int jel,
                        const RealType* distjI, const RealType* distjk)
  {
    valT Uj = valT(0);

    const DistanceTableData& eI_table=(*P.DistTables[myTableID]);
    const int jg=P.GroupID[jel];

    for(int iat=0; iat<Nion; ++iat)
      if(distjI[iat]<Ion_cutoff[iat])
      {
        const int ig=Ions.GroupID[iat];
        const valT r_Ij     = distjI[iat];

        for(int kg=0; kg<eGroups; ++kg)
        {
          const FT& feeI(*F(ig,jg,kg));
          int kel_counter=0;
          for(int kel=P.first(kg); kel<P.last(kg); kel++)
            if(eI_table.Distances[kel][iat]<Ion_cutoff[iat] && kel!=jel)
            {
              DistkI_Compressed[kel_counter]=eI_table.Distances[kel][iat];
              Distjk_Compressed[kel_counter]=distjk[kel];
              kel_counter++;
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
                        ParticleAttrib<valT>& Uk, ParticleAttrib<posT>& dUk, ParticleAttrib<valT>& d2Uk, bool triangle=false)
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
    for(int kel=0; kel<kelmax; kel++)
    {
      Uk[kel] = czero;
      dUk[kel] = posT();
      d2Uk[kel] = czero;
    }

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
        const valT r_Ij_inv = cone/r_Ij;

        for(int kg=0; kg<eGroups; ++kg)
        {
          const FT& feeI(*F(ig,jg,kg));
          int kel_counter=0;
          for(int kel=P.first(kg); kel<std::min(P.last(kg),kelmax); kel++)
            if(eI_table.Distances[kel][iat]<Ion_cutoff[iat] && kel!=jel)
            {
              DistkI_Compressed[kel_counter]=eI_table.Distances[kel][iat];
              Distjk_Compressed[kel_counter]=distjk[kel];
              DistIndice[kel_counter]=kel;
              kel_counter++;
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
            dUk[kel] += - gradF0[kel_index] * disp_jk - gradF2[kel_index] * disp_Ik;
            d2Uk[kel] -= hessF00[kel_index] + hessF22[kel_index]
                         + lapfac*(gradF0[kel_index] + gradF2[kel_index])
                         + ctwo*hessF02[kel_index]*dot(disp_jk,disp_Ik);
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
    APP_ABORT("JeeIOrbitalSoA::evaluateDerivatives not implemented yet!");
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
