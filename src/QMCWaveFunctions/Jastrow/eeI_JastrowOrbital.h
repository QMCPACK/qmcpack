//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_COMMON_EEN_JASTROW_H
#define QMCPLUSPLUS_COMMON_EEN_JASTROW_H
#include "Configuration.h"
#include  <map>
#include  <numeric>
#include "QMCWaveFunctions/OrbitalBase.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "LongRange/StructFact.h"
#include "ParticleBase/ParticleAttribOps.h"
#include <cmath>

namespace qmcplusplus
{

struct IonData
{
  typedef std::vector<int> eListType;
  OrbitalBase::RealType cutoff_radius;
  eListType elecs_inside;
  IonData() : cutoff_radius(0.0) { }
};


/** @ingroup OrbitalComponent
 *  @brief Specialization for three-body Jastrow function using multiple functors
 *
 *Each pair-type can have distinct function \f$u(r_{ij})\f$.
 *For electrons, distinct pair correlation functions are used
 *for spins up-up/down-down and up-down/down-up.
 */
template<class FT>
class eeI_JastrowOrbital: public OrbitalBase
{

  //flag to prevent parallel output
  bool Write_Chiesa_Correction;
  ///table index for i-el, el-el is always zero
  int myTableIndex;
  //nuber of particles
  int Nelec, Nion;
  //N*N
  int NN;
  //number of groups of the target particleset
  int eGroups, iGroups;
  RealType DiffVal, DiffValSum;
  RealType *FirstAddressOfdU, *LastAddressOfdU;
  ParticleAttrib<RealType> U,d2U,curLap_i, curLap_j, curVal;
  ParticleAttrib<PosType> dU,curGrad_i, curGrad_j;
  ParticleAttrib<PosType> curGrad0;
  ParticleAttrib<RealType> refVal;
  ParticleSet *eRef, *IRef;
  // The first index is the ion, the second two are
  // the electrons
  int TripletID(int i, int j, int k) {
    return (IRef->GroupID[i]*eRef->groups() + eRef->GroupID[j])*eRef->groups() + eRef->GroupID[k];
  }

  std::map<std::string,FT*> J3Unique;
  std::map<FT*,int> J3UniqueIndex;
  bool FirstTime;
  RealType KEcorr;

  std::vector<IonData> IonDataList;

  // Temporary store for parameter derivatives of functor
  // The first index is the functor index in J3Unique.  The second is the parameter index w.r.t. to that
  // functor
  std::vector<std::vector<RealType> >               du_dalpha;
  std::vector<std::vector<PosType> >             dgrad_dalpha;
  std::vector<std::vector<Tensor<RealType,3> > > dhess_dalpha;

  // Used for evaluating derivatives with respect to the parameters
  int NumVars;
  Array<std::pair<int,int>,3> VarOffset;
  Vector<RealType> dLogPsi;
  Array<PosType,2> gradLogPsi;
  Array<RealType,2> lapLogPsi;

public:

  typedef FT FuncType;

  ///container for the Jastrow functions
  Array<FT*,3> F;

  RealType ChiesaKEcorrection()
  {
    return 0.0;
  }

  eeI_JastrowOrbital(ParticleSet& ions, ParticleSet& elecs, bool is_master)
    : Write_Chiesa_Correction(is_master), KEcorr(0.0)
  {
    OrbitalName = "eeI_JastrowOrbital";
    eRef = &elecs;
    IRef = &ions;
    myTableIndex=elecs.addTable(ions,DT_AOS);
    init(elecs);
    FirstTime = true;
    NumVars=0;
  }

  ~eeI_JastrowOrbital() { }

  void init(ParticleSet& p)
  {
    Nelec=p.getTotalNum();
    NN = Nelec*Nelec;
    Nion = IRef->getTotalNum();
    U.resize(Nelec*Nelec+1);
    d2U.resize(Nelec*Nelec);
    dU.resize(Nelec*Nelec);
    curGrad_i.resize(Nelec);
    curGrad_j.resize(Nelec);
    curGrad0.resize(Nelec);
    curLap_i.resize(Nelec);
    curLap_j.resize(Nelec);
    curVal.resize(Nelec);
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    int nisp=iGroups=IRef->getSpeciesSet().getTotalNum();
    int nesp=eGroups=p.groups();
    F.resize(nisp,nesp,nesp);
    F = 0;
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

  void addFunc(int iSpecies,
               int eSpecies1, int eSpecies2, FT* j)
  {
    if(eSpecies1==eSpecies2)
    {
      //if only up-up is specified, assume spin-unpolarized correlations
      if(eSpecies1==0)
      {
        int ijk = iSpecies * eGroups*eGroups;
        for (int eG1=0; eG1<eGroups; eG1++)
          for (int eG2=0; eG2<eGroups; eG2++, ijk++)
            if(F(iSpecies, eG1, eG2)==0)
              F(iSpecies, eG1, eG2) = j;
      }
    }
    else
    {
      F(iSpecies,eSpecies1,eSpecies2) = j;
      F(iSpecies, eSpecies2, eSpecies1) = j;
    }
    if(j)
    {
      RealType rcut = 0.5 * j->cutoff_radius;
      for (int i=0; i<Nion; i++)
        if (IRef->GroupID[i] == iSpecies)
          IonDataList[i].cutoff_radius = rcut;
    }
    else
    {
      APP_ABORT("eeI_JastrowOrbital::addFunc  Jastrow function pointer is NULL");
    }
    std::stringstream aname;
    aname << iSpecies << "_" << eSpecies1 << "_" << eSpecies2;
    J3Unique[aname.str()]=j;
    initUnique();
    FirstTime = false;
  }


  /** check that correlation information is complete
   */
  void check_complete()
  {
    //check that correlation pointers are either all 0 or all assigned
    bool complete = true;
    int ni = F.size(0);
    int ne = F.size(1);
    int ne2 = ne*ne;
    for(int i=0; i<ni; ++i)
    {
      int nfilled = 0;
      bool partial;
      for(int e1=0; e1<ne; ++e1)
        for(int e2=0; e2<ne; ++e2)
          if(F(i,e1,e2)!=0)
            nfilled++;
      partial = nfilled>0 && nfilled<ne2;
      if(partial)
        app_log() << "J3 eeI is missing correlation for ion "<<i<< std::endl;
      complete = complete && !partial;
    }
    if(!complete)
    {
      APP_ABORT("eeI_JastrowOrbital::check_complete  J3 eeI is missing correlation components\n  see preceding messages for details");
    }
    //first set radii
    for(int i=0; i<Nion; ++i)
    {
      FT* f = F(IRef->GroupID[i],0,0);
      if(f!=0)
        IonDataList[i].cutoff_radius = .5*f->cutoff_radius;
    }
    //then check radii
    bool all_radii_match = true;
    for(int i=0; i<ni; ++i)
    {
      if(F(i,0,0)!=0)
      {
        bool radii_match = true;
        RealType rcut = F(i,0,0)->cutoff_radius;
        for(int e1=0; e1<ne; ++e1)
          for(int e2=0; e2<ne; ++e2)
            radii_match = radii_match && F(i,e1,e2)->cutoff_radius==rcut;
        if(!radii_match)
          app_log() << "eeI functors for ion species " << i
                    << " have different radii"<< std::endl;
        all_radii_match = all_radii_match && radii_match;
      }
    }
    if(!all_radii_match)
    {
      APP_ABORT("eeI_JastrowOrbital::check_radii  J3 eeI are inconsistent for some ion species\n  see preceding messages for details");
    }
  }



  //evaluate the distance table with els
  void resetTargetParticleSet(ParticleSet& P)
  {
    eRef = &P;
    //      if(dPsi) dPsi->resetTargetParticleSet(P);
  }

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
    // myVars.getIndex(active);
    // Optimizable=myVars.is_optimizable();
    // typename std::map<std::string,FT*>::iterator it(J3Unique.begin()),it_end(J3Unique.end());
    // while(it != it_end)
    // {
    //   (*it++).second->checkOutVariables(active);
    // }
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
    //myVars.print(std::cout);
    if (NumVars)
    {
      dLogPsi.resize(NumVars);
      gradLogPsi.resize(NumVars,Nelec);
      lapLogPsi.resize(NumVars,Nelec);
      // for (int i=0; i<NumVars; ++i) {
      //   gradLogPsi[i].resize(Nelec);//=new GradVectorType(Nelec);
      //   lapLogPsi[i].resize(Nelec);//=new ValueVectorType(Nelec);
      // }
      int nisp=iGroups=IRef->getSpeciesSet().getTotalNum();
      int nesp=eGroups=eRef->groups();
      VarOffset.resize(Nion, Nelec, Nelec);
      int varoffset=myVars.Index[0];
      for (int i=0; i<Nion; i++)
      {
        if(IonDataList[i].cutoff_radius>0.0)
        {
          for (int j=0; j<Nelec; j++)
            for (int k=0; k<Nelec; k++)
            {
              FT &func_ijk = *F.data()[TripletID(i, j, k)];
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
    //if (FirstTime) {
    // if(!IsOptimizing)
    // {
    //   app_log() << "  Chiesa kinetic energy correction = "
    //     << ChiesaKEcorrection() << std::endl;
    //   //FirstTime = false;
    // }
    //      if(dPsi) dPsi->resetParameters( active );
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
    ChiesaKEcorrection();
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
    // HACK HACK HACK
    // evaluateLogAndStore(P,G,L);
    // return LogValue;
    LogValue=0.0;
    for (int jk=0; jk<Nelec*Nelec; jk++)
      U[jk] = 0.0;
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
        // std::cerr << "jel = " << jel << " dtable j = " << eI_table->J[nn0+jel] << std::endl;
        RealType r_Ij     = eI_table->r(nn0+jel);
        RealType r_Ij_inv = eI_table->rinv(nn0+jel);
        int ee0 = ee_table->M[jel]-(jel+1);
        for (int k=j+1; k<ion.elecs_inside.size(); k++)
        {
          int kel = ion.elecs_inside[k];
          // std::cerr << "kel = " << kel << " dtable k = " << ee_table->J[ee0+kel] << std::endl;
          // std::cerr << "jel,kel = " << jel << ", " << kel << std::endl;
          RealType r_Ik     = eI_table->r(nn0+kel);
          RealType r_Ik_inv = eI_table->rinv(nn0+kel);
          RealType r_jk     = ee_table->r(ee0+kel);
          RealType r_jk_inv = ee_table->rinv(ee0+kel);
          FT &func = *F.data()[TripletID(i, jel, kel)];
          u = func.evaluate (r_jk, r_Ij, r_Ik, gradF, hessF);
          LogValue -= u;
          // Save for ratio
          U[jel*Nelec+kel] += u;
          U[kel*Nelec+jel] += u;
          PosType gr_ee =    gradF[0]*r_jk_inv * ee_table->dr(ee0+kel);
          PosType du_j, du_k;
          RealType d2u_j, d2u_k;
          du_j = gradF[1]*r_Ij_inv * eI_table->dr(nn0+jel) - gr_ee;
          du_k = gradF[2]*r_Ik_inv * eI_table->dr(nn0+kel) + gr_ee;
          d2u_j = (hessF(0,0) + 2.0*r_jk_inv*gradF[0] -
                   2.0*hessF(0,1)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+jel))*r_jk_inv*r_Ij_inv
                   + hessF(1,1) + 2.0*r_Ij_inv*gradF[1]);
          d2u_k = (hessF(0,0) + 2.0*r_jk_inv*gradF[0] +
                   2.0*hessF(0,2)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+kel))*r_jk_inv*r_Ik_inv
                   + hessF(2,2) + 2.0*r_Ik_inv*gradF[2]);
          G[jel] -= du_j;
          G[kel] -= du_k;
          L[jel] -= d2u_j;
          L[kel] -= d2u_k;
          // PosType gr_ee =    gradF[0]*r_jk_inv * ee_table->dr(ee0+kel);
          // G[jel] +=  gr_ee - gradF[1]*r_Ij_inv * eI_table->dr(nn0+jel);
          // G[kel] -=  gr_ee + gradF[2]*r_Ik_inv * eI_table->dr(nn0+kel);
          // L[jel] -= (hessF(0,0) + 2.0*r_jk_inv*gradF[0] -
          // 	       2.0*hessF(0,1)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+jel))*r_jk_inv*r_Ij_inv
          // 	       + hessF(1,1) + 2.0*r_Ij_inv*gradF[1]);
          // L[kel] -= (hessF(0,0) + 2.0*r_jk_inv*gradF[0] +
          // 	       2.0*hessF(0,2)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+kel))*r_jk_inv*r_Ik_inv
          // 	       + hessF(2,2) + 2.0*r_Ik_inv*gradF[2]);
        }
      }
    }
    return LogValue;
    // RealType dudr, d2udr2;
    // PosType gr;
    // for(int i=0; i<ee_table->size(SourceIndex); i++) {
    // 	for(int nn=ee_table->M[i]; nn<ee_table->M[i+1]; nn++) {
    // 	  int j = ee_table->J[nn];
    // 	  //LogValue -= F[ee_table->PairID[nn]]->evaluate(ee_table->r(nn), dudr, d2udr2);
    // 	  RealType uij = F[ee_table->PairID[nn]]->evaluate(ee_table->r(nn), dudr, d2udr2);
    //     LogValue -= uij;
    // 	  U[i*Nelec+j]=uij; U[j*Nelec+i]=uij; //save for the ratio
    // 	  //multiply 1/r
    // 	  dudr *= ee_table->rinv(nn);
    // 	  gr = dudr*ee_table->dr(nn);
    // 	  //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)
    // 	  RealType lap(d2udr2+2.0*dudr);
    // 	  //multiply -1
    // 	  G[i] += gr;
    // 	  G[j] -= gr;
    // 	  L[i] -= lap;
    // 	  L[j] -= lap;
    // 	}
    // }
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }

  inline GradType evalGradSourceFD(ParticleSet& P,
                                   ParticleSet& source, int isrc)
  {
    GradType G;
    const RealType eps=1.0e-6;
    int N = P.G.size();
    ParticleSet::ParticleGradient_t  Gt(N);
    ParticleSet::ParticleLaplacian_t Lt(N);
    PosType itmp = source.R[isrc];
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      PosType delta;
      // delta[dim] = eps;
      // source.makeMove(isrc, delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim] + eps;
      P.update();
      G[dim] = evaluateLog(P, Gt, Lt);
      // source.makeMove(isrc, -2.0*delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim] - eps;
      P.update();
      G[dim] -= evaluateLog (P, Gt, Lt);
      // source.makeMove(isrc, delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim];
      P.update();
      evaluateLog (P, Gt, Lt);
    }
    G = 0.5/eps * G;
    return G;
  }

  inline GradType evalGradSource(ParticleSet& P,
                                 ParticleSet& source, int isrc)
  {
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
    IonData &ion = IonDataList[isrc];
    ion.elecs_inside.clear();
    int iel=0;
    if (ion.cutoff_radius > 0.0)
      for (int nn=eI_table->M[isrc]; nn<eI_table->M[isrc+1]; nn++, iel++)
        if (eI_table->r(nn) < ion.cutoff_radius)
          ion.elecs_inside.push_back(iel);
    GradType G;
    int nn0 = eI_table->M[isrc];
    RealType u;
    PosType gradF;
    Tensor<RealType,3> hessF;
    for (int j=0; j<ion.elecs_inside.size(); j++)
    {
      int jel = ion.elecs_inside[j];
      RealType r_Ij     = eI_table->r(nn0+jel);
      PosType dr_Ij     = eI_table->dr(nn0+jel);
      RealType r_Ij_inv = eI_table->rinv(nn0+jel);
      int ee0 = ee_table->M[jel]-(jel+1);
      for (int k=j+1; k<ion.elecs_inside.size(); k++)
      {
        int kel = ion.elecs_inside[k];
        RealType r_Ik     = eI_table->r(nn0+kel);
        PosType dr_Ik     = eI_table->dr(nn0+kel);
        RealType r_Ik_inv = eI_table->rinv(nn0+kel);
        RealType r_jk     = ee_table->r(ee0+kel);
        RealType r_jk_inv = ee_table->rinv(ee0+kel);
        FT &func = *F.data()[TripletID(isrc, jel, kel)];
        u = func.evaluate (r_jk, r_Ij, r_Ik, gradF, hessF);
        G += (gradF[1] * r_Ij_inv * dr_Ij +
              gradF[2] * r_Ik_inv * dr_Ik);
      }
    }
    return G;
  }


  inline GradType
  evalGradSourceFD(ParticleSet& P, ParticleSet& source, int isrc,
                   TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
                   TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    GradType G;
    const RealType eps=1.0e-6;
    int N = P.G.size();
    ParticleSet::ParticleGradient_t  grad_plus(N), grad_minus(N);
    ParticleSet::ParticleLaplacian_t lapl_plus(N), lapl_minus(N);
    PosType itmp = source.R[isrc];
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      grad_plus  = GradType();
      grad_minus = GradType();
      lapl_plus  = 0.0;
      lapl_minus = 0.0;
      PosType delta;
      //delta[dim] = eps;
      // source.makeMove(isrc, delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim] + eps;
      P.update();
      G[dim] = evaluateLog(P, grad_plus, lapl_plus);
      // source.makeMove(isrc, -2.0*delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim] - eps;
      P.update();
      G[dim] -= evaluateLog (P, grad_minus, lapl_minus);
      // source.makeMove(isrc, delta);
      // source.acceptMove(isrc);
      source.R[isrc][dim] = itmp[dim];
      P.update();
      //evaluateLog (P, P.G, P.L);
      for (int i=0; i<N; i++)
      {
        grad_grad[dim][i] = (grad_plus[i] - grad_minus[i])/(2.0*eps);
        lapl_grad[dim][i] = (lapl_plus[i] - lapl_minus[i])/(2.0*eps);
      }
    }
    source.R[isrc] = itmp;
    G = 0.5/eps * G;
    return G;
  }


  inline GradType
  evalGradSource(ParticleSet& P, ParticleSet& source, int isrc,
                 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
                 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    //return (evalGradSourceFD(P, source, isrc, grad_grad, lapl_grad));
    // for (int dim=0; dim<OHMMS_DIM; dim++) {
    // 	grad_grad[dim] = GradType();
    // 	lapl_grad[dim] = RealType();
    // }
    IonData &ion = IonDataList[isrc];
    // ion.elecs_inside.clear();
    // int iel=0;
    // if (ion.cutoff_radius > 0.0)
    // 	for (int nn=eI_table->M[isrc]; nn<eI_table->M[isrc+1]; nn++, iel++)
    // 	  if (eI_table->r(nn) < ion.cutoff_radius)
    // 	    ion.elecs_inside.push_back(iel);
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
    GradType G;
    int nn0 = eI_table->M[isrc];
    RealType u;
    PosType gradF;
    Tensor<RealType,3> hessF;
    TinyVector<Tensor<RealType,3>,3> d3F;
    for (int j=0; j<ion.elecs_inside.size(); j++)
    {
      int jel = ion.elecs_inside[j];
      RealType r_Ij     = eI_table->r(nn0+jel);
      PosType dr_Ij     = eI_table->dr(nn0+jel);
      RealType r_Ij_inv = eI_table->rinv(nn0+jel);
      PosType dr_Ij_hat = r_Ij_inv * dr_Ij;
      int ee0 = ee_table->M[jel]-(jel+1);
      for (int k=j+1; k<ion.elecs_inside.size(); k++)
      {
        int kel = ion.elecs_inside[k];
        RealType r_Ik     = eI_table->r(nn0+kel);
        PosType dr_Ik     = eI_table->dr(nn0+kel);
        RealType r_Ik_inv = eI_table->rinv(nn0+kel);
        PosType dr_Ik_hat = r_Ik_inv * dr_Ik;
        RealType r_jk     = ee_table->r(ee0+kel);
        RealType r_jk_inv = ee_table->rinv(ee0+kel);
        PosType dr_jk_hat = r_jk_inv * ee_table->dr(ee0+kel);
        FT &func = *F.data()[TripletID(isrc, jel, kel)];
        u = func.evaluate (r_jk, r_Ij, r_Ik, gradF, hessF, d3F);
        if (j < k)
          G += (gradF[1] * r_Ij_inv * dr_Ij +
                gradF[2] * r_Ik_inv * dr_Ik);
        for (int dim_ion=0; dim_ion < OHMMS_DIM; dim_ion++)
        {
          for (int dim_el=0; dim_el < OHMMS_DIM; dim_el++)
          {
            // Should be 6 terms total.
            grad_grad[dim_ion][jel][dim_el] -=
              hessF(0,1)*dr_jk_hat[dim_el]*dr_Ij_hat[dim_ion] -
              (hessF(1,1) - r_Ij_inv*gradF[1])*dr_Ij_hat[dim_el]*dr_Ij_hat[dim_ion] +
              hessF(0,2)*dr_jk_hat[dim_el]*dr_Ik_hat[dim_ion] -
              hessF(1,2)*dr_Ij_hat[dim_el]*dr_Ik_hat[dim_ion];
            grad_grad[dim_ion][kel][dim_el] -=
              -hessF(0,2)*dr_jk_hat[dim_el]*dr_Ik_hat[dim_ion] -
              (hessF(2,2) - r_Ik_inv*gradF[2])*dr_Ik_hat[dim_el]*dr_Ik_hat[dim_ion] +
              -hessF(0,1)*dr_jk_hat[dim_el]*dr_Ij_hat[dim_ion] -
              hessF(1,2)*dr_Ik_hat[dim_el]*dr_Ij_hat[dim_ion];
          }
          grad_grad[dim_ion][jel][dim_ion] += r_Ij_inv * gradF[1];
          grad_grad[dim_ion][kel][dim_ion] += r_Ik_inv * gradF[2];
          lapl_grad[dim_ion][jel] +=
            (d3F[0](0,1) + 2.0*r_jk_inv*hessF(0,1) +
             -2.0*(d3F[0](1,1) - r_Ij_inv*hessF(0,1)) * dot(dr_jk_hat, dr_Ij_hat) +
             d3F[1](1,1) - 2.0*r_Ij_inv*r_Ij_inv*gradF[1] + 2.0*r_Ij_inv*hessF(1,1)) * dr_Ij_hat[dim_ion];
          lapl_grad[dim_ion][jel] +=
            (d3F[0](0,2) + 2.0*r_jk_inv*hessF(0,2) - 2.0 * d3F[0](1,2)*dot(dr_jk_hat, dr_Ij_hat) +
             d3F[1](1,2) + 2.0*r_Ij_inv*hessF(1,2))*dr_Ik_hat[dim_ion];
          lapl_grad[dim_ion][jel] -= 2.0*r_Ij_inv*hessF(0,1)*dr_jk_hat[dim_ion];
          lapl_grad[dim_ion][kel] +=
            (d3F[0](0,2) + 2.0*r_jk_inv*hessF(0,2) +
             +2.0*(d3F[0](2,2) - r_Ik_inv*hessF(0,2)) * dot(dr_jk_hat, dr_Ik_hat) +
             d3F[2](2,2) - 2.0*r_Ik_inv*r_Ik_inv*gradF[2] + 2.0*r_Ik_inv*hessF(2,2)) * dr_Ik_hat[dim_ion];
          lapl_grad[dim_ion][kel] +=
            (d3F[0](0,1) + 2.0*r_jk_inv*hessF(0,1) + 2.0 * d3F[0](1,2)*dot(dr_jk_hat, dr_Ik_hat) +
             d3F[1](2,2) + 2.0*r_Ik_inv*hessF(1,2))*dr_Ij_hat[dim_ion];
          lapl_grad[dim_ion][kel] += 2.0*r_Ik_inv*hessF(0,2)*dr_jk_hat[dim_ion];
        }
      }
    }
    return G;
  }

  inline void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    const int iat=VP.activePtcl;
    const int nk=ratios.size();
    int nat=iat*Nelec;
    RealType x=std::accumulate(&(U[nat]),&(U[nat+Nelec]),0.0);
    std::vector<RealType> newval(nk,x);
    const DistanceTableData* ee_table=VP.DistTables[0];
    const DistanceTableData* eI_table=VP.DistTables[myTableIndex];
    const DistanceTableData* eI_0=VP.refPtcl.DistTables[myTableIndex];

    for (int i=0,nn=0; i<Nion; i++)
    {
      IonData &ion = IonDataList[i];
      int nn0=eI_0->M[i];
      for(int k=0; k<nk; ++k, ++nn)
      {
        RealType r_Ii = eI_table->r(nn);
        if (r_Ii < ion.cutoff_radius)
        {
          for (int j=0; j<ion.elecs_inside.size(); j++)
          {
            int jat = ion.elecs_inside[j];
            if (jat != iat)
            {
              RealType r_ij = ee_table->r(jat*nk+k);
              RealType r_Ij = eI_0->r(nn0+jat);
              FT &func = *F.data()[TripletID(i, iat, jat)];
              newval[k] -= func.evaluate(r_ij, r_Ii, r_Ij);
            }
          }
        }
      }
    }
    for(int k=0; k<ratios.size(); ++k)
      ratios[k]=std::exp(newval[k]);
  }

  ValueType ratio(ParticleSet& P, int iat)
  {
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
            FT &func = *F.data()[TripletID(i, iat, jat)];
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
  }

  GradType evalGrad(ParticleSet& P, int iat)
  {
    GradType gr;
    for(int jat=0,ij=iat*Nelec; jat<Nelec; ++jat,++ij)
      gr -= dU[ij];
    return gr;
  }

  ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    curVal  = 0.0;
    curGrad_i = PosType();
    curLap_i  = 0.0;
    curGrad_j = PosType();
    curLap_j = 0.0;
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
            FT &func = *F.data()[TripletID(i, iat, jat)];
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
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
  }


  inline void evaluateLogAndStore(ParticleSet& P,
                                  ParticleSet::ParticleGradient_t& G,
                                  ParticleSet::ParticleLaplacian_t& L)
  {
    //      std::cerr << "evaluateLogAndStore called.\n";
    LogValue=0.0;
    const DistanceTableData* ee_table=P.DistTables[0];
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
    RealType u;
    PosType gradF;
    Tensor<RealType,3> hessF;
    // Zero out cached data
    for (int jk=0; jk<Nelec*Nelec; jk++)
    {
      U[jk] = 0.0;
      dU[jk] = PosType();
      d2U[jk] = 0.0;
    }
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
          FT &func = *F.data()[TripletID(i, jel, kel)];
          u = func.evaluate (r_jk, r_Ij, r_Ik, gradF, hessF);
          LogValue -= u;
          PosType gr_ee =    gradF[0]*r_jk_inv * ee_table->dr(ee0+kel);
          PosType du_j, du_k;
          RealType d2u_j, d2u_k;
          du_j = gradF[1]*r_Ij_inv * eI_table->dr(nn0+jel) - gr_ee;
          du_k = gradF[2]*r_Ik_inv * eI_table->dr(nn0+kel) + gr_ee;
          d2u_j = (hessF(0,0) + 2.0*r_jk_inv*gradF[0] -
                   2.0*hessF(0,1)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+jel))*r_jk_inv*r_Ij_inv
                   + hessF(1,1) + 2.0*r_Ij_inv*gradF[1]);
          d2u_k = (hessF(0,0) + 2.0*r_jk_inv*gradF[0] +
                   2.0*hessF(0,2)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+kel))*r_jk_inv*r_Ik_inv
                   + hessF(2,2) + 2.0*r_Ik_inv*gradF[2]);
          G[jel] -= du_j;
          G[kel] -= du_k;
          L[jel] -= d2u_j;
          L[kel] -= d2u_k;
          int jk = jel*Nelec+kel;
          int kj = kel*Nelec+jel;
          U[jk] += u;
          U[kj] += u;
          dU[jk] += du_j;
          dU[kj] += du_k;
          d2U[jk] += d2u_j;
          d2U[kj] += d2u_k;
          // G[jel] +=  gr_ee - gradF[1]*r_Ij_inv * eI_table->dr(nn0+jel);
          // G[kel] -=  gr_ee + gradF[2]*r_Ik_inv * eI_table->dr(nn0+kel);
          // L[jel] -= (hessF(0,0) + 2.0*r_jk_inv*gradF[0] -
          // 	       2.0*hessF(0,1)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+jel))*r_jk_inv*r_Ij_inv
          // 	       + hessF(1,1) + 2.0*r_Ij_inv*gradF[1]);
          // L[kel] -= (hessF(0,0) + 2.0*r_jk_inv*gradF[0] +
          // 	       2.0*hessF(0,2)*dot(ee_table->dr(ee0+kel),eI_table->dr(nn0+kel))*r_jk_inv*r_Ik_inv
          // 	       + hessF(2,2) + 2.0*r_Ik_inv*gradF[2]);
        }
      }
    }
    int iat = 2;
    PosType G2 = 0.0;
    for (int jat=0; jat<Nelec; jat++)
      G2 -= dU[iat*Nelec+jat];
    // std::cerr << "G2   = " << G2   << std::endl;
    // std::cerr << "G[2] = " << G[2] << std::endl;
    // if (FirstTime) {
    // 	FirstTime = false;
    // 	ChiesaKEcorrection();
    // }
    // RealType dudr, d2udr2,u;
    // LogValue=0.0;
    // GradType gr;
    // for(int i=0; i<ee_table->size(SourceIndex); i++) {
    // 	for(int nn=ee_table->M[i]; nn<ee_table->M[i+1]; nn++) {
    // 	  int j = ee_table->J[nn];
    // 	  u = F[ee_table->PairID[nn]]->evaluate(ee_table->r(nn), dudr, d2udr2);
    // 	  LogValue -= u;
    // 	  dudr *= ee_table->rinv(nn);
    // 	  gr = dudr*ee_table->dr(nn);
    // 	  //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
    // 	  RealType lap = d2udr2+2.0*dudr;
    // 	  int ij = i*Nelec+j, ji=j*Nelec+i;
    // 	  U[ij]=u; U[ji]=u;
    // 	  //dU[ij] = gr; dU[ji] = -1.0*gr;
    // 	  dU[ij] = gr; dU[ji] = gr*-1.0;
    // 	  d2U[ij] = -lap; d2U[ji] = -lap;
    // 	  //add gradient and laplacian contribution
    // 	  dG[i] += gr;
    // 	  dG[j] -= gr;
    // 	  dL[i] -= lap;
    // 	  dL[j] -= lap;
    // 	}
    // }
  }

  inline void registerData(ParticleSet& P, WFBufferType& buf)
  {
    buf.add(U.begin(), U.end());
    buf.add(d2U.begin(), d2U.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
  }

  inline RealType updateBuffer(ParticleSet& P, WFBufferType& buf,
                               bool fromscratch=false)
  {
    evaluateLogAndStore(P,P.G,P.L);
    //RealType dudr, d2udr2,u;
    //LogValue=0.0;
    //GradType gr;
    //PairID.resize(ee_table->size(SourceIndex),ee_table->size(SourceIndex));
    //int nsp=P.groups();
    //for(int i=0; i<ee_table->size(SourceIndex); i++)
    //  for(int j=0; j<ee_table->size(SourceIndex); j++)
    //    PairID(i,j) = P.GroupID[i]*nsp+P.GroupID[j];
    //for(int i=0; i<ee_table->size(SourceIndex); i++) {
    //  for(int nn=ee_table->M[i]; nn<ee_table->M[i+1]; nn++) {
    //    int j = ee_table->J[nn];
    //    //ValueType sumu = F.evaluate(ee_table->r(nn));
    //    u = F[ee_table->PairID[nn]]->evaluate(ee_table->r(nn), dudr, d2udr2);
    //    LogValue -= u;
    //    dudr *= ee_table->rinv(nn);
    //    gr = dudr*ee_table->dr(nn);
    //    //(d^2 u \over dr^2) + (2.0\over r)(du\over\dr)\f$
    //    RealType lap = d2udr2+2.0*dudr;
    //    int ij = i*N+j, ji=j*N+i;
    //    U[ij]=u; U[ji]=u;
    //    //dU[ij] = gr; dU[ji] = -1.0*gr;
    //    dU[ij] = gr; dU[ji] = gr*-1.0;
    //    d2U[ij] = -lap; d2U[ji] = -lap;
    //    //add gradient and laplacian contribution
    //    P.G[i] += gr;
    //    P.G[j] -= gr;
    //    P.L[i] -= lap;
    //    P.L[j] -= lap;
    //  }
    //}
    U[NN]= LogValue;
    buf.put(U.begin(), U.end());
    buf.put(d2U.begin(), d2U.end());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    // for (int i=0; i<IonDataList.size(); i++) {
    // 	int n = IonDataList[i].elecs_inside.size();
    // 	buf.put ((RealType)n);
    // 	vector<RealType> elecs_inside(n);
    // 	for (int j=0; j<n; j++)
    // 	  elecs_inside[j] = IonDataList[i].elecs_inside[j];
    // 	buf.put(elecs_inside.begin(), elecs_inside.end());
    // }
    return LogValue;
  }

  inline void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    //      std::cerr << "Called copyFromBuffer.\n";
    buf.get(U.begin(), U.end());
    buf.get(d2U.begin(), d2U.end());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    // First, create lists of electrons within the sphere of each ion
    const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
    // for (int i=0; i<IonDataList.size(); i++) {
    // 	RealType nd;
    // 	buf.get(nd);
    // 	int n = (int)round(nd);
    // 	cerr << "n = " << n << std::endl;
    // 	vector<RealType> elecs_inside(n);
    // 	buf.get(elecs_inside.begin(), elecs_inside.end());
    // 	IonDataList[i].elecs_inside.resize(n);
    // 	for (int j=0; j<n; j++)
    // 	  IonDataList[i].elecs_inside[j] = (int)round(elecs_inside[j]);
    // }
    DiffValSum=0.0;
  }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    eeI_JastrowOrbital<FT>* eeIcopy=
      new eeI_JastrowOrbital<FT>(*IRef, tqp, false);
    std::map<const FT*,FT*> fcmap;
    for (int iG=0; iG<iGroups; iG++)
      for (int eG1=0; eG1<eGroups; eG1++)
        for (int eG2=0; eG2<eGroups; eG2++)
        {
          int ijk = iG*eGroups*eGroups + eG1*eGroups + eG2;
          if(F(iG,eG1,eG2)==0)
            continue;
          typename std::map<const FT*,FT*>::iterator fit=fcmap.find(F(iG,eG1,eG2));
          if(fit == fcmap.end())
          {
            FT* fc=new FT(*F(iG,eG1,eG2));
            eeIcopy->addFunc( iG, eG1, eG2, fc);
            //if (dPsi) (eeIcopy->dPsi)->addFunc(aname.str(),ig,jg,fc);
            fcmap[F(iG,eG1,eG2)]=fc;
          }
        }
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

  void copyFrom(const OrbitalBase& old)
  {
    //nothing to do
  }


  void
  finalizeOptimization()
  {
    ChiesaKEcorrection();
  }

  RealType KECorrection()
  {
    return KEcorr;
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
      const DistanceTableData* ee_table=P.DistTables[0];
      const DistanceTableData* eI_table=P.DistTables[myTableIndex];
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
            FT &func = *F.data()[TripletID(i, jel, kel)];
            int idx = J3UniqueIndex[F.data()[TripletID(i, jel, kel)]];
            func.evaluateDerivatives(r_jk, r_Ij, r_Ik, du_dalpha[idx],
                                     dgrad_dalpha[idx], dhess_dalpha[idx]);
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
  }
};




}
#endif

