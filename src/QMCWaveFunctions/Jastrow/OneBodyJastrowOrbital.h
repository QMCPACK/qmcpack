//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Bryan Clark, bclark@Princeton.edu, Princeton University
//                    Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_GENERIC_ONEBODYJASTROW_H
#define QMCPLUSPLUS_GENERIC_ONEBODYJASTROW_H
#include "Configuration.h"
#include "QMCWaveFunctions/OrbitalBase.h"
#include "QMCWaveFunctions/Jastrow/DiffOneBodyJastrowOrbital.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"

namespace qmcplusplus
{

/** @ingroup OrbitalComponent
 * @brief generic implementation of one-body Jastrow function.
 *
 *The One-Body Jastrow has the form
 \f$ \psi = \exp{\left(-J_1\right)} \f$
 \f[ J_{1}({\bf R}) = \sum_{I=1}^{N_I}
 \sum_{i=1}^{N_e} u(r_{iI}) \f]
 where \f[ r_{iI} = |{\bf r}_i - {\bf R}_I| \f]
 and the first summnation \f$ 1..N_I \f$ is over the centers,
 e.g., nulcei for molecular or solid systems,
 while the second summnation \f$ 1..N_e \f$ is over quantum particles.
 *
 *The template parameter FT is a functional \f$ u(r_{iI}), \f$ e.g.
 *class PadeJastrow<T> or NoCuspJastrow<T>
 *Requriement of the template function is
 *ValueType evaluate(ValueType r, ValueType& dudr, ValueType& d2udr2).
 *
 To calculate the Gradient use the identity
 \f[ {\bf \nabla}_i(r_{iI}) = \frac{{{\bf r_{iI}}}}{r_{iI}}. \f]
 \f[ {\bf \nabla}_i(J_{e}({\bf R})=
 \sum_{I=1}^{N_I} \frac{du}{dr_{iI}}{\bf \hat{r_{iI}}}.
 \f]
 To calculate The Laplacian use the identity
 \f[ \nabla^2_i(r_{iI})=\frac{2}{r_{iI}}, \f]
 and the vector product rule
 \f[
 \nabla^2 \cdot (f{\bf A}) =
 f({\bf \nabla} \cdot {\bf A})+{\bf A}\cdot({\bf \nabla} f)
 \f]
 \f[
 \nabla^2_i (J_{e}({\bf R})) = \sum_{I=1}^{N_I}
 \left(\frac{du}{dr_{iI}}\right) {\bf \nabla}_i
 \cdot {\bf \hat{r_{iI}}} - {\bf \hat{r_{iI}}} \cdot
 \left(\frac{d^2u}{dr_{iI}^2}\right){\bf \hat{r_{iI}}}
 \f]
 which can be simplified to
 \f[
 \nabla^2_i (J_{e}({\bf R})) = \sum_{I=1}^{N_I}
 \left(\frac{2}{r_{iI}}\frac{du}{dr_{iI}}
 + \frac{d^2u}{dr_{iI}^2}\right)
 \f]
 *
 *A generic OneBodyJastrowOrbital function uses a distance table.
 *The indices I(sources) and i(targets) are distinct. In general, the source
 *particles are fixed, e.g., the nuclei, while the target particles are updated
 *by MC methods.
 */
template<class FT>
class OneBodyJastrowOrbital: public OrbitalBase
{
protected:
  int myTableIndex;
  const ParticleSet& CenterRef;
  RealType curVal, curLap;
  PosType curGrad;
  ParticleAttrib<RealType> U,d2U;
  ParticleAttrib<PosType> dU;
  RealType *FirstAddressOfdU, *LastAddressOfdU;
  std::vector<FT*> Funique;

public:
  std::vector<FT*> Fs;

  typedef FT FuncType;


  ///constructor
  OneBodyJastrowOrbital(const ParticleSet& centers, ParticleSet& els)
    : CenterRef(centers), FirstAddressOfdU(0), LastAddressOfdU(0)
  {
    U.resize(els.getTotalNum());
    dU.resize(els.getTotalNum());
    d2U.resize(els.getTotalNum());

    myTableIndex=els.addTable(CenterRef,DT_AOS);
    //allocate vector of proper size  and set them to 0
    Funique.resize(CenterRef.getSpeciesSet().getTotalNum(),nullptr);
    Fs.resize(CenterRef.getTotalNum(),nullptr);
    OrbitalName = "OneBodyJastrow";
  }

  ~OneBodyJastrowOrbital() { }

  //evaluate the distance table with P
  void resetTargetParticleSet(ParticleSet& P)
  {
    if (dPsi)
      dPsi->resetTargetParticleSet(P);
  }

  void addFunc(int source_type, FT* afunc, int target_type=-1)
  {
    for (int i=0; i<Fs.size(); i++)
      if (CenterRef.GroupID[i] == source_type)
        Fs[i]=afunc;
    Funique[source_type]=afunc;
  }

  /** check in an optimizable parameter
   * @param o a super set of optimizable variables
   */
  void checkInVariables(opt_variables_type& active)
  {
    myVars.clear();
    for (int i=0; i<Funique.size(); ++i)
    {
      FT* fptr=Funique[i];
      if (fptr)
      {
        fptr->checkInVariables(active);
        fptr->checkInVariables(myVars);
      }
    }
  }

  /** check out optimizable variables
   */
  void checkOutVariables(const opt_variables_type& active)
  {
    myVars.getIndex(active);
    Optimizable=myVars.is_optimizable();
    for (int i=0; i<Funique.size(); ++i)
      if (Funique[i] != nullptr)
        Funique[i]->checkOutVariables(active);
    if (dPsi)
      dPsi->checkOutVariables(active);
  }

  ///reset the value of all the unique Two-Body Jastrow functions
  void resetParameters(const opt_variables_type& active)
  {
    if (!Optimizable)
      return;
    for (int i=0; i<Funique.size(); ++i)
      if (Funique[i])
        Funique[i]->resetParameters(active);
    for (int i=0; i<myVars.size(); ++i)
    {
      int ii=myVars.Index[i];
      if (ii>=0)
        myVars[i]= active[ii];
    }
    if (dPsi)
      dPsi->resetParameters(active);
  }

  /** print the state, e.g., optimizables */
  void reportStatus(std::ostream& os)
  {
    for (int i=0; i<Funique.size(); ++i)
    {
      if (Funique[i])
        Funique[i]->myVars.print(os);
    }
  }

  /**
   *@param P input configuration containing N particles
   *@param G a vector containing N gradients
   *@param L a vector containing N laplacians
   *@return The wavefunction value  \f$exp(-J({\bf R}))\f$
   *
   *Upon exit, the gradient \f$G[i]={\bf \nabla}_i J({\bf R})\f$
   *and the laplacian \f$L[i]=\nabla^2_i J({\bf R})\f$ are accumulated.
   *While evaluating the value of the Jastrow for a set of
   *particles add the gradient and laplacian contribution of the
   *Jastrow to G(radient) and L(aplacian) for local energy calculations
   *such that \f[ G[i]+={\bf \nabla}_i J({\bf R}) \f]
   *and \f[ L[i]+=\nabla^2_i J({\bf R}). \f]
   */
  RealType evaluateLog(ParticleSet& P,
                       ParticleSet::ParticleGradient_t& G,
                       ParticleSet::ParticleLaplacian_t& L)
  {
    LogValue=0.0;
    U=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    RealType dudr, d2udr2;
    for (int i=0; i<d_table->size(SourceIndex); i++)
    {
      FT* func=Fs[i];
      if (func == nullptr) continue;
      for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //ValueType uij= F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue -= uij;
        U[j] += uij;
        dudr *= d_table->rinv(nn);
        G[j] -= dudr*d_table->dr(nn);
        L[j] -= d2udr2+2.0*dudr;
      }
    }
    return LogValue;
  }

  ValueType evaluate(ParticleSet& P,
                     ParticleSet::ParticleGradient_t& G,
                     ParticleSet::ParticleLaplacian_t& L)
  {
    return std::exp(evaluateLog(P,G,L));
  }
  
  void evaluateHessian(ParticleSet& P, HessVector_t& grad_grad_psi)
  {
    LogValue=0.0;
    U=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    RealType dudr, d2udr2;

    Tensor<RealType,OHMMS_DIM> ident;
    grad_grad_psi=0.0;
    ident.diagonal(1.0);

    for (int i=0; i<d_table->size(SourceIndex); i++)
    {
      FT* func=Fs[i];
      if (func == nullptr) continue;
      for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        RealType rinv=d_table->rinv(nn);
        RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2);
        grad_grad_psi[j]-= rinv*rinv*outerProduct(d_table->dr(nn),d_table->dr(nn))*(d2udr2-dudr*rinv) + ident*dudr*rinv;
      }
    }
  } 
  
  /** evaluate the ratio \f$exp(U(iat)-U_0(iat))\f$
   * @param P active particle set
   * @param iat particle that has been moved.
   */
  inline ValueType ratio(ParticleSet& P, int iat)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    curVal=0.0;
    for (int i=0; i<d_table->size(SourceIndex); ++i)
      if (Fs[i] != nullptr)
        curVal += Fs[i]->evaluate(d_table->Temp[i].r1);
    return std::exp(U[iat]-curVal);
  }

  inline void evaluateRatios(VirtualParticleSet& VP, std::vector<ValueType>& ratios)
  {
    //vector<RealType> myr(ratios.size(),0.0);
    //const DistanceTableData* d_table=VP.DistTables[myTableIndex];
    //for (int i=0; i<d_table->size(SourceIndex); ++i)
    //  if (Fs[i])
    //    for (int nn=d_table->M[i],j=0; nn<d_table->M[i+1]; ++nn,++j)
    //      myr[j]+=Fs[i]->evaluate(d_table->r(nn));
    //RealType x=U[VP.activePtcl];
    //for(int k=0; k<ratios.size(); ++k)
    //  ratios[k]=std::exp(x-myr[k]);
    std::vector<RealType> myr(ratios.size(),U[VP.activePtcl]);
    const DistanceTableData* d_table=VP.DistTables[myTableIndex];
    for (int i=0; i<d_table->size(SourceIndex); ++i)
      if (Fs[i]!=nullptr)
        for (int nn=d_table->M[i],j=0; nn<d_table->M[i+1]; ++nn,++j)
          myr[j]-=Fs[i]->evaluate(d_table->r(nn));
    for(int k=0; k<ratios.size(); ++k)
      ratios[k]=std::exp(myr[k]);
  }

  inline void evaluateDerivRatios(VirtualParticleSet& VP, const opt_variables_type& optvars,
      std::vector<ValueType>& ratios, Matrix<ValueType>& dratios)
  {
    evaluateRatios(VP,ratios);
    if(dPsi)
      dPsi->evaluateDerivRatios(VP,optvars,dratios);
  }

  inline GradType evalGrad(ParticleSet& P, int iat)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int n=d_table->size(VisitorIndex);
    curGrad = 0.0;
    RealType ur,dudr, d2udr2;
    for (int i=0, nn=iat; i<d_table->size(SourceIndex); ++i,nn+= n)
    {
      if (Fs[i] != nullptr)
      {
        ur=Fs[i]->evaluate(d_table->r(nn),dudr,d2udr2);
        dudr *= d_table->rinv(nn);
        curGrad -= dudr*d_table->dr(nn);
      }
    }
    return curGrad;
  }

  inline GradType evalGradSource(ParticleSet& P,
                                 ParticleSet& source, int isrc)
  {
    if (&source != &CenterRef)
      return GradType();
    FT* func=Fs[isrc];
    if (func == 0)
      return GradType();
    GradType G(0.0);
    RealType dudr, d2udr2;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    for (int nn=d_table->M[isrc]; nn<d_table->M[isrc+1]; nn++)
    {
      RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2);
      dudr *= d_table->rinv(nn);
      G += dudr*d_table->dr(nn);
    }
    return G;
  }


  inline GradType
  evalGradSource(ParticleSet& P, ParticleSet& source, int isrc,
                 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
                 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
  {
    if (&source != &CenterRef)
      return GradType();
    FT* func=Fs[isrc];
    if (func == 0)
      return GradType();
    GradType G(0.0);
    RealType dudr, d2udr2, d3udr3=0.0;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    for (int nn=d_table->M[isrc],iel=0; nn<d_table->M[isrc+1]; nn++,iel++)
    {
      RealType rinv = d_table->rinv(nn);
      RealType uij= func->evaluate(d_table->r(nn), dudr, d2udr2, d3udr3);
      dudr *= rinv;
      d2udr2 *= rinv * rinv;
      G += dudr*d_table->dr(nn);
      for (int dim_ion=0; dim_ion < OHMMS_DIM; dim_ion++)
      {
        for (int dim_el=0; dim_el < OHMMS_DIM; dim_el++)
          grad_grad[dim_ion][iel][dim_el] += d2udr2 * d_table->dr(nn)[dim_ion] * d_table->dr(nn)[dim_el]
                                             - dudr * rinv * rinv * d_table->dr(nn)[dim_ion] * d_table->dr(nn)[dim_el];
        grad_grad[dim_ion][iel][dim_ion] += dudr;
        lapl_grad[dim_ion][iel] += (d3udr3*rinv + 2.0*d2udr2 - 2.0*rinv*rinv*dudr)*d_table->dr(nn)[dim_ion];
      }
    }
    return G;
  }

  inline ValueType ratioGrad(ParticleSet& P, int iat, GradType& grad_iat)
  {
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    int n=d_table->size(VisitorIndex);
    curVal=0.0;
    curGrad = 0.0;
    RealType dudr, d2udr2;
    for (int i=0; i<d_table->size(SourceIndex); ++i)
    {
      if (Fs[i] != nullptr)
      {
        curVal += Fs[i]->evaluate(d_table->Temp[i].r1,dudr,d2udr2);
        dudr *= d_table->Temp[i].rinv1;
        curGrad -= dudr*d_table->Temp[i].dr1;
      }
    }
    grad_iat += curGrad;
    return std::exp(U[iat]-curVal);
  }

  inline void restore(int iat) {}

  void acceptMove(ParticleSet& P, int iat)
  {
    U[iat] = curVal;
    dU[iat]=curGrad;
    d2U[iat]=curLap;
  }

  void evaluateLogAndStore(ParticleSet& P,
                           ParticleSet::ParticleGradient_t& dG,
                           ParticleSet::ParticleLaplacian_t& dL)
  {
    LogValue = 0.0;
    U=0.0;
    dU=0.0;
    d2U=0.0;
    RealType uij, dudr, d2udr2;
    const DistanceTableData* d_table=P.DistTables[myTableIndex];
    for (int i=0,ns=d_table->size(SourceIndex); i<ns; i++)
    {
      FT* func=Fs[i];
      if (func ==  nullptr) continue;
      for (int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++)
      {
        int j = d_table->J[nn];
        //uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
        uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
        LogValue-=uij;
        U[j]+=uij;
        dudr *= d_table->rinv(nn);
        dU[j] -= dudr*d_table->dr(nn);
        d2U[j] -= d2udr2+2.0*dudr;
        //add gradient and laplacian contribution
        dG[j] -= dudr*d_table->dr(nn);
        dL[j] -= d2udr2+2.0*dudr;
      }
    }
  }

  /** equivalent to evalaute with additional data management */
  void registerData(ParticleSet& P, WFBufferType& buf)
  {
    FirstAddressOfdU = &(dU[0][0]);
    LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    //add U, d2U and dU. Keep the order!!!
    DEBUG_PSIBUFFER(" OneBodyJastrow::registerData ",buf.current());
    buf.add(U.begin(), U.end());
    buf.add(d2U.begin(), d2U.end());
    buf.add(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodyJastrow::registerData ",buf.current());
  }

  void evaluateGL(ParticleSet& P)
  {
    evaluateLogAndStore(P,P.G,P.L);
  }

  RealType updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch=false)
  {
    evaluateLogAndStore(P,P.G,P.L);
    //LogValue = 0.0;
    //U=0.0; dU=0.0; d2U=0.0;
    //RealType uij, dudr, d2udr2;
    //for(int i=0; i<d_table->size(SourceIndex); i++) {
    //  FT* func=Fs[i];
    //  if(func == 0) continue;
    //  for(int nn=d_table->M[i]; nn<d_table->M[i+1]; nn++) {
    //    int j = d_table->J[nn];
    //    //uij = F[d_table->PairID[nn]]->evaluate(d_table->r(nn), dudr, d2udr2);
    //    uij = func->evaluate(d_table->r(nn), dudr, d2udr2);
    //    LogValue-=uij;
    //    U[j]+=uij;
    //    dudr *= d_table->rinv(nn);
    //    dU[j] -= dudr*d_table->dr(nn);
    //    d2U[j] -= d2udr2+2.0*dudr;
    //    //add gradient and laplacian contribution
    //    P.G[j] -= dudr*d_table->dr(nn);
    //    P.L[j] -= d2udr2+2.0*dudr;
    //  }
    //}
    //FirstAddressOfdU = &(dU[0][0]);
    //LastAddressOfdU = FirstAddressOfdU + dU.size()*DIM;
    DEBUG_PSIBUFFER(" OneBodyJastrow::updateBuffer ",buf.current());
    buf.put(U.first_address(), U.last_address());
    buf.put(d2U.first_address(), d2U.last_address());
    buf.put(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodyJastrow::updateBuffer ",buf.current());
    return LogValue;
  }

  /** copy the current data from a buffer
   *@param P the ParticleSet to operate on
   *@param buf PooledData which stores the data for each walker
   *
   *copyFromBuffer uses the data stored by registerData or evaluate(P,buf)
   */
  void copyFromBuffer(ParticleSet& P, WFBufferType& buf)
  {
    DEBUG_PSIBUFFER(" OneBodyJastrow::copyFromBuffer ",buf.current());
    buf.get(U.first_address(), U.last_address());
    buf.get(d2U.first_address(), d2U.last_address());
    buf.get(FirstAddressOfdU,LastAddressOfdU);
    DEBUG_PSIBUFFER(" OneBodyJastrow::copyFromBuffer ",buf.current());
  }

  OrbitalBasePtr makeClone(ParticleSet& tqp) const
  {
    OneBodyJastrowOrbital<FT>* j1copy=new OneBodyJastrowOrbital<FT>(CenterRef,tqp);
    j1copy->Optimizable=Optimizable;
    for (int i=0; i<Funique.size(); ++i)
    {
      if (Funique[i])
        j1copy->addFunc(i,new FT(*Funique[i]));
    }
    //j1copy->OrbitalName=OrbitalName+"_clone";
    if (dPsi)
    {
      j1copy->dPsi =  dPsi->makeClone(tqp);
    }
    return j1copy;
  }

  void copyFrom(const OrbitalBase& old)
  {
    //nothing to do
  }
};

}
#endif

