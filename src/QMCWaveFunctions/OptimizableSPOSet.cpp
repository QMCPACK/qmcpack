//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include "QMCWaveFunctions/OptimizableSPOSet.h"
#include "Numerics/OhmmsBlas.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

void
OptimizableSPOSet::addParameter ( std::string id, int iorb, int basis)
{
}

bool
OptimizableSPOSet::put (xmlNodePtr node, SPOPool_t &spo_pool)
{
  std::string gsName, basisName;
  bool same_k = false;
  bool mapped_k = false;
  bool same_orbital = false;
  bool lapl = false;
  derivScale=1.0;
  OhmmsAttributeSet attrib;
  attrib.add (gsName,    "gs_sposet");
  attrib.add (basisName, "basis_sposet");
  attrib.add (same_k,    "same_k");
  attrib.add (mapped_k,    "mapped_k");
  attrib.add (same_orbital,    "same_orbital");
  attrib.add (lapl,    "lapl");
  attrib.add (N,         "size");
  attrib.add (derivScale,"scale");
  attrib.add(thr, "thr");
  attrib.put (node);
  if (N == 0)
  {
    app_error() << "You must specify \"size\" attribute for linearopt sposet.\n";
    abort();
  }
  /////////////////////////////////////
  // First, find ground-state SPOSet //
  /////////////////////////////////////
  if (gsName == "")
  {
    app_error() << "You must supply \"gs_sposet\".  Aborting.\n";
    abort();
  }
  SPOPool_t::iterator iter = spo_pool.find(gsName);
  if (iter == spo_pool.end())
  {
    app_error() << "No sposet named \"" << gsName << "\" found.  Abort.\n";
    abort();
  }
  else
  {
    app_log() << "  Found ground-state SPOSet \"" << gsName << "\".\n";
    GSOrbitals = iter->second;
  }
  //////////////////////////////////////
  // Now, find basis SPOSet from pool //
  //////////////////////////////////////
  iter = spo_pool.find(basisName);
  if (iter != spo_pool.end())
  {
    BasisOrbitals = iter->second;
    app_log() << "  Found basis SPOSet \"" << basisName << "\".\n";
  }
  if (BasisOrbitals == GSOrbitals)
    BasisOrbitals = 0;
  if (BasisOrbitals)
    M = BasisOrbitals->getOrbitalSetSize();
  else
    M = GSOrbitals->getOrbitalSetSize() - N;
  resize(N,M);
  app_log() << "  linearopt sposet has " << N << " ground-state orbitals and "
            << M << " basis orbitals.\n";
  //Map K points to orbitals.
  ValueMatrix_t allowedOrbs;
  allowedOrbs.resize(N,M);
//    if (mapped_k)
//    {
//      std::cerr <<"Not Done"<< std::endl;
//    }
//    else
  if (same_k)
  {
    int off         = BasisOrbitals ? 0 : N;
    SPOSetBase* basOrbs = BasisOrbitals ? BasisOrbitals : GSOrbitals;
    for (int igs=0; igs<N; igs++)
    {
      PosType k_gs = GSOrbitals->get_k(igs);
      for (int ib=0; ib<M; ib++)
      {
        PosType k_b = basOrbs->get_k(ib+off);
        if (dot(k_gs-k_b, k_gs-k_b) < 1.0e-6)
        {
//             app_log()<<"Allowing ("<<igs<<","<<ib<<") :";
//             for(int ix=0;ix<3;ix++) app_log()<<k_gs[ix]<<"  ";
//             app_log()<<" : ";
//             for(int ix=0;ix<3;ix++) app_log()<<k_b[ix]<<"  ";
//             app_log()<< std::endl;
          allowedOrbs(igs,ib)=1;
        }
        else
        {
//             app_log()<<"Not allowing ("<<igs<<","<<ib<<") :";
//             for(int ix=0;ix<3;ix++) app_log()<<k_gs[ix]<<"  ";
//             app_log()<<" : ";
//             for(int ix=0;ix<3;ix++) app_log()<<k_b[ix]<<"  ";
//             app_log()<< std::endl;
          allowedOrbs(igs,ib)=0;
        }
      }
    }
  }
  else
    if (mapped_k)
    {
      xmlNodePtr xmlCoefs = node->xmlChildrenNode;
      while (xmlCoefs != NULL)
      {
        std::string cname((const char*)xmlCoefs->name);
        if (cname == "orbitalmap")
        {
          std::string type("0");
          OhmmsAttributeSet cAttrib;
          cAttrib.add(type, "type");
          cAttrib.put(xmlCoefs);
          if (type != "Array")
          {
            // app_error() << "Unknown correlation type " + type + " in OptimizableSPOSet." + "Resetting to \"Array\"\n";
            xmlNewProp(xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
          }
          std::vector<RealType> params;
          putContent(params, xmlCoefs);
          int indx(0);
          if(params.size()==N*M)
            for (int igs=0; igs<N; igs++)
            {
              for (int ib=0; ib<M; ib++)
              {
                allowedOrbs(igs,ib)=params[indx++];
              }
            }
          else
          {
            app_error()<<"Map size is incorrect. parameters given= "<<params.size()<<" parameters needed= "<<M*N<< std::endl;
          }
        }
        xmlCoefs = xmlCoefs->next;
      }
    }
    else
      if(same_orbital)
      {
        if(M<N)
        {
          app_error()<<"Map size is incorrect. N<M for same orbital mapping"<< std::endl;
        }
        allowedOrbs=0;
        for (int igs=0; igs<N; igs++)
          allowedOrbs(igs,igs)=1;
      }
      else
        if(lapl)
        {
          for (int igs=0; igs<N; igs++)
            for (int ib=0; ib<(M/N); ib++)
              allowedOrbs(igs,igs+N*ib)=1;
        }
        else
        {
          for (int igs=0; igs<N; igs++)
          {
            for (int ib=0; ib<M; ib++)
            {
              allowedOrbs(igs,ib)=1;
            }
          }
        }
  // Now, look for coefficients element
  xmlNodePtr xmlCoefs = node->xmlChildrenNode;
  while (xmlCoefs != NULL)
  {
    std::string cname((const char*)xmlCoefs->name);
    if (cname == "coefficients")
    {
      std::string type("0"), id("0");
      int state=-1;
      int asize(-1);
      OhmmsAttributeSet cAttrib;
      cAttrib.add(id, "id");
      cAttrib.add(type, "type");
      cAttrib.add(asize, "size");
      cAttrib.add(state, "state");
      cAttrib.put(xmlCoefs);
      if (state == -1)
      {
        app_error() << "You must specify the \"state\" attribute in <linearopt>'s <coefficient>.\n";
      }
      if (type != "Array")
      {
        // app_error() << "Unknown correlation type " + type + " in OptimizableSPOSet." + "Resetting to \"Array\"\n";
        xmlNewProp(xmlCoefs, (const xmlChar*) "type", (const xmlChar*) "Array");
      }
      std::vector<RealType> params;
      putContent(params, xmlCoefs);
      app_log() << "Coefficients for state" << state << ":\n";
      // std::cerr << "params.size() = " << params.size() << std::endl;
      //If params is missized resize and initialize to zero.
      if ((params.size()!=asize)&&(asize>0))
      {
        std::vector<RealType> t_params(params);
        params.resize(asize,0.0);
        int ipm(0);
        if (t_params.size()>0)
        {
          for (int ib=0; ib<M; ib++)
          {
#if defined(QMC_COMPLEX)
            if((allowedOrbs(state,ib).real()>0)&&(std::abs(t_params[ipm])>thr))
            {
              params[ib]=t_params[ipm++];
#else
            if((allowedOrbs(state,ib)>0)&&(std::abs(t_params[ipm])>thr))
            {
              params[ib]=t_params[ipm++];
#endif
            }
            else
            {
              params[ib]=0;
            }
          }
        }
      }
      else
      {
        if (params.size())
          asize=params.size();
      }
//    for (int i=0; i< params.size(); i++) {
      ValueMatrix_t& cref(*C);
      int indx=0;
      for (int i=0; i< M; i++)
      {
        std::stringstream sstr;
#ifndef QMC_COMPLEX
// 	  ParamPointers.push_back(&(C(state,i)));
// 	  ParamIndex.push_back(TinyVector<int,2>(state,i));
        sstr << id << "_" << indx;
        if (allowedOrbs(state,i))
        {
          ParamPointers.push_back(&(cref(state,i)));
          ParamIndex.push_back(TinyVector<int,2>(state,i));
          if ((indx<asize)&&(std::abs(params[indx])>thr))
            cref(state,i) = params[indx];
          else
            cref(state,i) = 0.0;
          myVars.insert(sstr.str(),cref(state,i),true,optimize::SPO_P);
          indx++;
        }
        else
        {
          cref(state,i) = 0.0;
//        myVars.insert(sstr.str(),C(state,i),false,optimize::LINEAR_P);
        }
#else
        //ParamPointers.push_back(&(C(state,i).real()));
        //ParamPointers.push_back(&(C(state,i).imag()));
        ParamPointers.push_back(reinterpret_cast<RealType*>(cref[state]+i));
        ParamPointers.push_back(reinterpret_cast<RealType*>(cref[state]+i)+1);
        ParamIndex.push_back(TinyVector<int,2>(state,i));
        ParamIndex.push_back(TinyVector<int,2>(state,i));
        sstr << id << "_" << 2*i+0;
        myVars.insert(sstr.str(),cref(state,i).real(),true,optimize::SPO_P);
        sstr << id << "_" << 2*i+1;
        myVars.insert(sstr.str(),cref(state,i).imag(),true,optimize::SPO_P);
#endif
      }
      for (int i=0; i<params.size(); i++)
      {
        char buf[100];
        snprintf (buf, 100, "  %12.5f\n", params[i]);
        app_log() << buf;
      }
    }
    xmlCoefs = xmlCoefs->next;
  }
  return SPOSetBase::put(node);
}


void
OptimizableSPOSet::resetTargetParticleSet(ParticleSet& P)
{
  GSOrbitals->resetTargetParticleSet(P);
  if (BasisOrbitals)
    BasisOrbitals->resetTargetParticleSet(P);
}

void
OptimizableSPOSet::setOrbitalSetSize(int norbs)
{
  OrbitalSetSize = norbs;
}


void
OptimizableSPOSet::checkInVariables(opt_variables_type& active)
{
  active.insertFrom(myVars);
}

void
OptimizableSPOSet::checkOutVariables(const opt_variables_type& active)
{
  myVars.getIndex(active);
}


void
OptimizableSPOSet::resetParameters(const opt_variables_type& active)
{
  for (int i=0; i<ParamPointers.size(); i++)
  {
    int loc=myVars.where(i);
    if (loc>=0)
      if(std::abs(active[loc])>thr)
        *(ParamPointers[i])=myVars[i]=active[loc];
  }
}

// Obsolete
void
OptimizableSPOSet::evaluateDerivatives
(ParticleSet& P, int iat, const opt_variables_type& active,
 ValueMatrix_t& d_phi, ValueMatrix_t& d_lapl_phi)
{
  for (int i=0; i<d_phi.size(); i++)
    for (int j=0; j<N; j++)
      d_phi[i][j] = d_lapl_phi[i][j] = ValueType();
  // Evaluate basis states
  if (BasisOrbitals)
  {
    BasisOrbitals->evaluate(P,iat,BasisVal);
    std::vector<TinyVector<int,2> >::iterator iter;
    std::vector<TinyVector<int,2> >& act = ActiveBasis[iat];
    for (iter=act.begin(); iter != act.end(); iter++)
    {
      int elem  = (*iter)[0];
      int param = (*iter)[1];
      int loc = myVars.where(param);
      if (loc >= 0)
        ;
    }
  }
  else
  {
  }
}

void
OptimizableSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
{
  GSOrbitals->evaluate(P,iat,GSVal);
  if (BasisOrbitals)
  {
    BasisOrbitals->evaluate(P,iat,BasisVal);
    BLAS::gemv_trans (N, M, C->data(), &(BasisVal[0]), &(psi[0]));
  }
  else
    BLAS::gemv_trans (N, M, C->data(), &(GSVal[N]), &(psi[0]));
//       for (int i=0; i<N; i++) {
// 	psi[i] = 0.0;
// 	for (int j=0; j<M; j++)
// 	  psi[i] += C(i,j)*GSVal[N+j];
//       }
  for (int i=0; i<N; i++)
    psi[i] += GSVal[i];
}

void
OptimizableSPOSet::evaluate(const ParticleSet& P, const PosType& r,
                            std::vector<RealType> &psi)
{
  app_error() << "OptimizableSPOSet::evaluate(const ParticleSet& P, const PosType& r, std::vector<RealType> &psi)\n  should not be called.  Abort.\n";
  abort();
}

void
OptimizableSPOSet::evaluate(const ParticleSet& P, int iat,
                            ValueVector_t& psi, GradVector_t& dpsi,
                            ValueVector_t& d2psi)
{
  GSOrbitals->evaluate(P,iat,GSVal,GSGrad,GSLapl);
  const ValueMatrix_t& cref(*C);
  if (BasisOrbitals)
  {
    BasisOrbitals->evaluate(P,iat,BasisVal,BasisGrad,BasisLapl);
    BLAS::gemv_trans (N, M, C->data(), &(BasisVal[0]),  &(psi[0]));
    BLAS::gemv_trans (N, M, C->data(), &(BasisLapl[0]), &(d2psi[0]));
    for (int iorb=0; iorb<N; iorb++)
    {
      psi  [iorb] += GSVal[iorb];
      dpsi [iorb] = GSGrad[iorb];
      d2psi[iorb] += GSLapl[iorb];
      for (int ibasis=0; ibasis<M; ibasis++)
        for (int dim=0; dim<OHMMS_DIM; dim++)
          dpsi[iorb][dim] += cref(iorb,ibasis) * BasisGrad[ibasis][dim];
    }
  }
  else
  {
    for (int iorb=0; iorb<N; iorb++)
    {
      psi  [iorb] = GSVal[iorb];
      dpsi [iorb] = GSGrad[iorb];
      d2psi[iorb] = GSLapl[iorb];
      for (int ibasis=0; ibasis<M; ibasis++)
      {
        psi[iorb] += cref(iorb,ibasis) * GSVal[N+ibasis];
        for (int dim=0; dim<OHMMS_DIM; dim++)
          dpsi[iorb][dim] += cref(iorb,ibasis) * GSGrad[N+ibasis][dim];
        d2psi[iorb] += cref(iorb,ibasis) * GSLapl[N+ibasis];
      }
    }
    // BLAS::gemv_trans (N, M, C.data(), &(GSVal[N]),  &(psi[0]));
    // BLAS::gemv_trans (N, M, C.data(), &(GSLapl[N]), &(d2psi[0]));
  }
  // for (int i=0; i<N; i++) {
  //   psi[i]  += GSVal[i];
  //   d2psi[i] += GSLapl[i];
  // }
}

void
OptimizableSPOSet::evaluate(const ParticleSet& P, int iat,
                            ValueVector_t& psi, GradVector_t& dpsi,
                            HessVector_t& gg_psi)
{
  app_error() << "Need specialization for OptimizableSPOSet::evaluate(HessVector_t)  Abort.\n";
  abort();
}

void
OptimizableSPOSet::evaluateBasis (const ParticleSet &P, int first, int last,
                                  ValueMatrix_t &basis_val, GradMatrix_t &basis_grad,
                                  ValueMatrix_t &basis_lapl)
{
  if (BasisOrbitals)
    BasisOrbitals->evaluate_notranspose(P, first, last, basis_val, basis_grad, basis_lapl);
  else
  {
    for (int iat=first; iat<last; iat++)
    {
      GSOrbitals->evaluate(P, iat, GSVal, GSGrad, GSLapl);
      for (int i=0; i<M; i++)
      {
        basis_val (iat-first,i) = GSVal[N+i];
        basis_grad(iat-first,i) = GSGrad[N+i];
        basis_lapl(iat-first,i) = GSLapl[N+i];
      }
    }
  }
}

void
OptimizableSPOSet::copyParamsFromMatrix (const opt_variables_type& active,
    const Matrix<RealType> &mat,
    std::vector<RealType> &destVec)
{
  for (int ip=0; ip<myVars.size(); ip++)
  {
    int loc = myVars.where(ip);
    if (loc >= 0)
    {
      TinyVector<int,2> idx = ParamIndex[ip];
      destVec[loc] += derivScale*mat(idx[0], idx[1]);
    }
  }
}

void
OptimizableSPOSet::copyParamsFromMatrix (const opt_variables_type& active,
    const Matrix<ComplexType> &mat,
    std::vector<RealType> &destVec)
{
  for (int ip=0; ip<myVars.size(); ip+=2)
  {
    int loc = myVars.where(ip);
    if (loc >= 0)
    {
      TinyVector<int,2> idx = ParamIndex[ip];
      assert (ParamIndex[ip+1][0] == idx[0] &&
              ParamIndex[ip+1][1] == idx[1]);
      destVec[loc] += derivScale*mat(idx[0], idx[1]).real();
      loc = myVars.where(ip+1);
      assert (loc >= 0);
      destVec[loc] += derivScale*mat(idx[0], idx[1]).imag();
    }
  }
}



void
OptimizableSPOSet::evaluate_notranspose
(const ParticleSet& P, int first, int last,
 ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
{
//     std::cerr << "GSValMatrix.size =(" << GSValMatrix.size(0) << ", " << GSValMatrix.size(1) << ")\n";
//     std::cerr << "GSGradMatrix.size =(" << GSGradMatrix.size(0) << ", " << GSGradMatrix.size(1) << ")\n";
//     std::cerr << "GSLaplMatrix.size =(" << GSLaplMatrix.size(0) << ", " << GSLaplMatrix.size(1) << ")\n";
//     std::cerr << "first=" << first << "  last=" << last << std::endl;
  GSOrbitals->evaluate_notranspose
  (P, first, last, GSValMatrix, GSGradMatrix, GSLaplMatrix);
  if (BasisOrbitals)
  {
    BasisOrbitals->evaluate_notranspose(P, first, last, BasisValMatrix, BasisGradMatrix, BasisLaplMatrix);
    //Note to Ken:
    //Use Numerics/MatrixOperators.h
    //for C=AB MatrixOperators::product(C,BasisValMatrix,logdet);
    //for C=AB^t MatrixOperators::ABt(C,BasisValMatrix,logdet);
    BLAS::gemm ('T', 'N', N, N, M, 1.0, C->data(),
                M, BasisValMatrix.data(), M, 0.0, logdet.data(), N);
    logdet += GSValMatrix;
    BLAS::gemm ('T', 'N', N, N, M, 1.0, C->data(),
                M, BasisLaplMatrix.data(), M, 0.0, d2logdet.data(), N);
    d2logdet += GSLaplMatrix;
    // Gradient part.
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)
          GradTmpSrc(i,j) = BasisGradMatrix(i,j)[dim];
      BLAS::gemm ('T', 'N', N, N, M, 1.0, C->data(), M,
                  GradTmpSrc.data(), M, 0.0, GradTmpDest.data(), N);
      for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
          dlogdet(i,j)[dim] = GradTmpDest(i,j) + GSGradMatrix(i,j)[dim];
    }
  }
  else
  {
    // HACK HACK HACK
//       BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
// 		  M, GSValMatrix.data()+N, M+N, 0.0, logdet.data(), N);
//       for (int i=0; i<N; i++)
// 	for (int j=0; j<N; j++)
// 	  logdet(i,j) += GSValMatrix(i,j);
    for (int iel=0; iel<N; iel++)
      for (int iorb=0; iorb<N; iorb++)
      {
        logdet(iel,iorb) = GSValMatrix(iel,iorb);
        const ValueType *restrict cptr=(*C)[iorb];
        for (int ibasis=0; ibasis<M; ibasis++)
          logdet(iel,iorb) += cptr[ibasis]*GSValMatrix(iel,N+ibasis);
      }
    //      logdet += GSValMatrix;
    // BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(),
    // 		  M, GSLaplMatrix.data()+N, M+N, 0.0, d2logdet.data(), N);
    for (int iel=0; iel<N; iel++)
      for (int iorb=0; iorb<N; iorb++)
      {
        const ValueType *restrict cptr=(*C)[iorb];
        d2logdet(iel,iorb) = GSLaplMatrix(iel,iorb);
        for (int ibasis=0; ibasis<M; ibasis++)
          d2logdet(iel,iorb) += cptr[ibasis]*GSLaplMatrix(iel,N+ibasis);
      }
    //d2logdet += GSLaplMatrix;
    // Gradient part.
    for (int i=0; i<N; i++)
      for (int iorb=0; iorb<N; iorb++)
      {
        const ValueType *restrict cptr=(*C)[iorb];
        dlogdet(i,iorb) = GSGradMatrix(i,iorb);
        for (int n=0; n<M; n++)
          dlogdet(i,iorb) += cptr[n]* GSGradMatrix(i,N+n);
      }
    // for (int dim=0; dim<OHMMS_DIM; dim++) {
    // 	for (int i=0; i<M; i++)
    // 	  for (int j=0; j<N; j++)
    // 	    GradTmpSrc(i,j) = GSGradMatrix(i,j+N)[dim];
    // 	BLAS::gemm ('T', 'N', N, N, M, 1.0, C.data(), M,
    // 		     GradTmpSrc.data(), M, 0.0, GradTmpDest.data(), N);
    // 	for (int i=0; i<N; i++)
    // 	  for (int j=0; j<N; j++)
    // 	    dlogdet(i,j)[dim] = GradTmpDest(i,j) + GSGradMatrix(i,j)[dim];
    //}
  }
}

SPOSetBase*
OptimizableSPOSet::makeClone() const
{
  SPOSetBase *gs, *basis(0);
  OptimizableSPOSet *clone;
  gs = GSOrbitals->makeClone();
  if (BasisOrbitals)
    basis = BasisOrbitals->makeClone();
  clone = new OptimizableSPOSet(N,gs,basis);
  clone->IsCloned=true;
  clone->C=C; //just pointer
  clone->thr=thr;
  clone->myVars=myVars;
  clone->derivScale=derivScale;
#if 0
  clone->ParamPointers.clear();
  clone->ParamIndex=ParamIndex;
  for(int i=0; i<ParamIndex.size() ; i++)
  {
#ifndef QMC_COMPLEX
    clone->ParamPointers.push_back(&(clone->C(ParamIndex[i][0],ParamIndex[i][1])));
#else
    int ci=ParamIndex[i][0];
    int cj=ParamIndex[i][1];
    //clone->ParamPointers.push_back(&(clone->C(ParamIndex[i][0],ParamIndex[i][1]).real()));
    //clone->ParamPointers.push_back(&(clone->C(ParamIndex[i][0],ParamIndex[i][1]).imag()));
    clone->ParamPointers.push_back(reinterpret_cast<RealType*>(clone->C[ci]+cj));
    clone->ParamPointers.push_back(reinterpret_cast<RealType*>(clone->C[ci]+cj)+1);
#endif
  }
#endif
//   for (int i=0; i< N; i++) {
//   for (int j=0; j< M; j++) {
//     std::stringstream sstr;
// #ifndef QMC_COMPLEX
//     clone->ParamPointers.push_back(&(clone->C(i,j)));
//     clone->ParamIndex.push_back(TinyVector<int,2>(i,j));
// #else
//     clone->ParamPointers.push_back(&(clone->C(i,j).real()));
//     clone->ParamPointers.push_back(&(clone->C(i,j).imag()));
//     clone->ParamIndex.push_back(TinyVector<int,2>(i,j));
//     clone->ParamIndex.push_back(TinyVector<int,2>(i,j));
// #endif
//   }
//   }
  return clone;
}


}
