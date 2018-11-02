//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"    
#include "QMCWaveFunctions/Fermion/MultiSlaterDeterminantFast.h"
#include "QMCWaveFunctions/Fermion/MultiDiracDeterminant.h"
#include "ParticleBase/ParticleAttribOps.h"

namespace qmcplusplus
{

MultiSlaterDeterminantFast::MultiSlaterDeterminantFast(ParticleSet& targetPtcl, MultiDiracDeterminant* up, MultiDiracDeterminant* dn):
  C2node_up(nullptr),C2node_dn(nullptr),C(nullptr),
  CSFcoeff(nullptr),DetsPerCSF(nullptr),CSFexpansion(nullptr),
  IsCloned(false),
  RatioTimer("MultiSlaterDeterminantFast::ratio"),
  RatioGradTimer("MultiSlaterDeterminantFast::ratioGrad"),
  RatioAllTimer("MultiSlaterDeterminantFast::ratio(all)"),
  Ratio1Timer("MultiSlaterDeterminantFast::detEval_ratio"),
  Ratio1GradTimer("MultiSlaterDeterminantFast::detEval_ratioGrad"),
  Ratio1AllTimer("MultiSlaterDeterminantFast::detEval_ratio(all)"),
  UpdateTimer("MultiSlaterDeterminantFast::updateBuffer"),
  EvaluateTimer("MultiSlaterDeterminantFast::evaluate"),
  AccRejTimer("MultiSlaterDeterminantFast::Accept_Reject")
{
  registerTimers();
  //Optimizable=true;
  Optimizable=true;
  ClassName="MultiSlaterDeterminantFast";
  usingCSF=false;
  NP = targetPtcl.getTotalNum();
  nels_up = targetPtcl.last(0)-targetPtcl.first(0);
  nels_dn = targetPtcl.last(1)-targetPtcl.first(1);
#if defined(ENABLE_SOA)
  m_nb_up = static_cast<qmcplusplus::LCAOrbitalSet*>(up->Phi)->BasisSetSize; 
  m_nb_dn = static_cast<qmcplusplus::LCAOrbitalSet*>(dn->Phi)->BasisSetSize; 
#else
  m_nb_up = (up->Phi)->BasisSetSize; 
  m_nb_dn = (dn->Phi)->BasisSetSize; 
#endif
  m_nmo_up= up->Phi->OrbitalSetSize;
  m_nmo_dn= dn->Phi->OrbitalSetSize;
  FirstIndex_up=targetPtcl.first(0);
  FirstIndex_dn=targetPtcl.first(1);
  Dets.resize(2);
  Dets[0]=up;
  Dets[1]=dn;
  myG.resize(NP);
  myL.resize(NP);
  myG_temp.resize(NP);
  myL_temp.resize(NP);
  myG_J.resize(NP);
  myL_J.resize(NP);


  usingBF=false;
  BFTrans=0;

#if defined(ENABLE_SOA)
  m_init_B_up = *(static_cast<qmcplusplus::LCAOrbitalSet*>(Dets[0]->Phi)->C); 
#else
  m_init_B_up = *((Dets[0]->Phi)->C); 
#endif
  Orbopt = false;
  CIopt  = false;

// NEED SOME SORT OF INTERNAL CHECK TO MAKE SURE THE REFERNCE MATRIX CORRESPONDS TO THE HARTREE FOCK MATRIX

// NEED TO INCLUDE AN INTERNAL CHECK TO MAKE SURE CSF AND ORBOPT AREN'T USED SIMULATANEOUSLY
//  if(usingCSF && Orbopt){
//    APP_ABORT("Orbital Optimization Does not work with CSFS only DETS")
//    app_log() << "Orbital Optimization Does not work with CSFS only DETS" << std::endl;
//  }

  if(m_nmo_up != m_nmo_dn || m_nb_up != m_nb_dn){
    APP_ABORT("OrbitalSetSize or BasisSetSize is not equal for spin up/dn Molecular orbitals. ")
    app_log() << "OrbitalSetSize or BasisSetSize is not equal for spin up/dn Molecular orbitals." << std::endl;
  }

}


void MultiSlaterDeterminantFast::buildOptVariables(std::vector<RealType>& input_params_up,
                                                   bool params_supplied_up,
                                                   std::vector<RealType>& input_params_dn,
                                                   bool params_supplied_dn)
{
  app_log() << "Multi-Slater class has CIopt set to "<< (CIopt ? "TRUE":"FALSE") <<" and Orbopt set to "<< (Orbopt ? "TRUE":"FALSE") << " which causes Optimizable to be set to "<< (Optimizable ? "TRUE":"FALSE") << std::endl;
  //a vector in which the element's index value correspond to Molecular Orbitals.
  //The element value at a index indicates how many times an electron is excited from or to that orbital in the Multi-Slater expansion i.e the indices with non-zero elements are active space orbitals
  std::vector<int> occupancy_vector_up (m_nmo_up,0);
  std::vector<int> occupancy_vector_dn (m_nmo_dn,0);

  // Function to fill occupancy_vectors, fill excitation_vectors and also return how many unique up/dn determinants are contained in the wavefunction 
  unique_ups =  this->build_occ_vec(Dets[0]->detData,
                                    nels_up,
                                    m_nmo_up,
                                    &occupancy_vector_up);

  unique_dns =  this->build_occ_vec(Dets[1]->detData,
                                    nels_dn,
                                    m_nmo_dn,
                                    &occupancy_vector_dn);

  // When calculating the parameter derivative of the Multi-Slater component of the wavefunction, each unique deterimant can contribute multiple times. 
  // The lookup_tbls are used so that a parameter derivative of a unique determinant is only done once and then scaled according to how many times it appears in the Multi-Slater expansion 
  lookup_tbls.resize(2);
  lookup_tbls[0].resize(unique_ups);
  lookup_tbls[1].resize(unique_dns);

  //construct lookup table for up spin
  for (int i(0); i < C2node_up->size(); i++)
  {
    lookup_tbls[0][(*C2node_up)[i]].push_back(i);
  } 
  //construct lookup table for dn spin
  for (int i(0); i < C2node_dn->size(); i++)
  {
    lookup_tbls[1][(*C2node_dn)[i]].push_back(i);
  } 

  // summing occupancy vectors to find all set of active rotations later on.
  for (int i=0; i<occupancy_vector_up.size(); i++){
    occupancy_vector_up[i] += occupancy_vector_dn[i];}

  if(CIopt)
  {
    if(usingCSF){
      m_first_var_pos = CSFcoeff->size() - 1;
    }
    else{
      // the -1 is because if there are n CI coefficients only n-1 of them are free variables
      m_first_var_pos = C->size() - 1;
    }  
  }
  else{
    m_first_var_pos = 0;
  }
  // create active rotations
  m_act_rot_inds_up.clear(); 

  int min_nel (nels_up < nels_dn ? nels_up : nels_dn);
  int max_nel (nels_up < nels_dn ? nels_dn : nels_up);

for(int i=0;i<m_nmo_up;i++)
  for(int j=i+1;j<m_nmo_up;j++)
   { 
    bool core_i(!occupancy_vector_up[i] and i <= min_nel-1); // true if orbital i is a 'core' orbital
    bool core_j(!occupancy_vector_up[j] and j <= min_nel-1); // true if orbital j is a 'core' orbital
    bool virt_i(!occupancy_vector_up[i] and i >  max_nel-1); // true if orbital i is a 'virtual' orbital
    bool virt_j(!occupancy_vector_up[j] and j >  max_nel-1); // true if orbital j is a 'virtual' orbital
      if( !( 
            ( core_i and core_j  ) 
                      or
            ( virt_i and virt_j ) 
           )    
        )
      {
        m_act_rot_inds_up.push_back(std::pair<int,int>(i,j)); // orbital rotation parameter accepted as long as rotation isn't core-core or virtual-virtual
      }
   }
 
    // This will add the orbital rotation parameters to myVars 
    // and will also read in initial parameter values supplied in input file     
    int p_up, q_up; 
    int nparams_active_up = m_act_rot_inds_up.size();

    for (int i=0; i< nparams_active_up; i++)
    {
      p_up = m_act_rot_inds_up[i].first;
      q_up = m_act_rot_inds_up[i].second;
      std::stringstream sstr_up;
      sstr_up << "SPO_UP"
              << "_orb_rot_"
              << ( p_up <   10 ? "0": "")
              << ( p_up <  100 ? "0": "")
              << ( p_up < 1000 ? "0": "")
              << p_up
              << "_"
              << ( q_up <   10 ? "0": "")
              << ( q_up <  100 ? "0": "")
              << ( q_up < 1000 ? "0": "")
              << q_up;

      // If the user input parameteres, use those. Otherwise, initialize the parameters to zero
      if (params_supplied_up) {
        myVars->insert(sstr_up.str(), input_params_up[i]);
      } else {
        myVars->insert(sstr_up.str(), 0.0);
      }
    }

    //Printing the parameters
    if(false){
      app_log() << std::string(16,' ') << "Parameter name" << std::string(15,' ') << "Value\n";
      myVars->print(app_log());
    }
 
    //the below code is very similar to  "Reset parameters function"
    // reading out the parameters that define the rotation into an antisymmetric matrix
    std::vector<RealType> rot_mat_up(m_nmo_up*m_nmo_up, 0.0);
    for (int i = 0; i < m_act_rot_inds_up.size(); i++) 
    {
      const int p = m_act_rot_inds_up[i].first;
      const int q = m_act_rot_inds_up[i].second;
      // m_first_var_pos is the index to the first parameter of the spin up electrons...
      rot_mat_up[p+q*m_nmo_up] =  (*myVars)[i + m_first_var_pos]; 
      rot_mat_up[q+p*m_nmo_up] = -(*myVars)[i + m_first_var_pos]; 

    }
    //exponentiate matrices and do BLAS command to perform rotation on m_b

    // exponentiate antisymmetric matrix to get the unitary rotation
    this->exponentiate_matrix(m_nmo_up, &rot_mat_up[0]);

#if defined(ENABLE_SOA)
    BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), static_cast<qmcplusplus::LCAOrbitalSet*>(Dets[0]->Phi)->C->data() , m_nb_up);
    BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), static_cast<qmcplusplus::LCAOrbitalSet*>(Dets[1]->Phi)->C->data() , m_nb_up);
#else
    BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), (Dets[0]->Phi)->C->data() , m_nb_up);
    BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), (Dets[1]->Phi)->C->data() , m_nb_up);
#endif    
    if(Orbopt)
    {
      T_up.resize(nels_up,m_nmo_up);
      T_dn.resize(nels_dn,m_nmo_dn);
    }
    if (!Orbopt)
    {
      //THIS ALLOWS FOR ORBITAL PARAMETERS TO BE READ IN EVEN WHEN THOSE PARAMETERS ARE NOT BEING OPTIMIZED
      //this assumes there are only CI coefficients ahead of the M_orb_coefficients 
      myVars->Index.erase(myVars->Index.begin()+m_first_var_pos,myVars->Index.end());
      myVars->NameAndValue.erase(myVars->NameAndValue.begin()+m_first_var_pos,myVars->NameAndValue.end());
      myVars->ParameterType.erase(myVars->ParameterType.begin()+m_first_var_pos,myVars->ParameterType.end());
      myVars->Recompute.erase(myVars->Recompute.begin()+m_first_var_pos,myVars->Recompute.end());
    }


}

void MultiSlaterDeterminantFast::initialize()
{
  if(C==nullptr)
  {
    C2node_up=new std::vector<size_t>;
    C2node_dn=new std::vector<size_t>;
    C=new std::vector<RealType>;
    CSFcoeff=new std::vector<RealType>;
    DetsPerCSF=new std::vector<size_t>;
    CSFexpansion=new std::vector<RealType>;
    myVars=new opt_variables_type;
    IsCloned=false;
  }
}

WaveFunctionComponentPtr MultiSlaterDeterminantFast::makeClone(ParticleSet& tqp) const
{
  MultiDiracDeterminant* up_clone = new MultiDiracDeterminant(*Dets[0]);
  MultiDiracDeterminant* dn_clone = new MultiDiracDeterminant(*Dets[1]);
  MultiSlaterDeterminantFast* clone = new MultiSlaterDeterminantFast(tqp,up_clone,dn_clone);
  if(usingBF)
  {
    BackflowTransformation *tr = BFTrans->makeClone(tqp);
    clone->setBF(tr);
  }
  clone->resetTargetParticleSet(tqp);

  //Set IsCloned so that only the main object handles the optimizable data
  clone->IsCloned=true;

  clone->C2node_up=C2node_up;
  clone->C2node_dn=C2node_dn;
  clone->C=C;
  clone->myVars=myVars;

  clone->Optimizable=Optimizable;
  clone->usingCSF=usingCSF;
  clone->usingBF=usingBF;

  if (usingCSF)
  {
    clone->CSFcoeff=CSFcoeff;
    clone->CSFexpansion=CSFexpansion;
    clone->DetsPerCSF=DetsPerCSF;
  }
  return clone;
}

MultiSlaterDeterminantFast::~MultiSlaterDeterminantFast() 
{ 
  if(!IsCloned)
  {
    delete myVars;
    delete CSFexpansion;
    delete DetsPerCSF;
    delete CSFcoeff;
    delete C;
    delete C2node_dn;
    delete C2node_up;
  }
  //clean up determinants too!
}

void MultiSlaterDeterminantFast::resetTargetParticleSet(ParticleSet& P)
{
  if(usingBF)
  {
    BFTrans->resetTargetParticleSet(P);
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(BFTrans->QP);
  }
  else
  {
    for(int i=0; i<Dets.size(); i++)
      Dets[i]->resetTargetParticleSet(P);
  }
}


void MultiSlaterDeterminantFast::testMSD(ParticleSet& P, int iat)
{
//     APP_ABORT("Testing disabled for safety");
  app_log() <<"Testing MSDFast. \n";
  int n = nels_up+nels_dn;
  ParticleSet::ParticleGradient_t G(n),G0(n);
  ParticleSet::ParticleLaplacian_t L(n),L0(n);
  ValueType log0;
  GradType G1;
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  /*
       app_log() <<"Testing evaluate(P,G,L). \n";
       std::cout << std::endl << std::endl;
       std::cout <<"Psi: " <<log <<"   " <<log0 <<"   " <<log/log0 << std::endl;

       for(int i=0; i<n; i++) {
         std::cout <<i  <<"\n"
             <<"  x: " <<G[i][0]-G0[i][0] <<"\n"
             <<"  y: " <<G[i][1]-G0[i][1] <<"\n"
             <<"  z: " <<G[i][2]-G0[i][2] <<"\n"
             <<"  d2: " <<L(i)-L0(i) <<"\n"
             << std::endl;
       }
       std::cout << std::endl << std::endl;
       APP_ABORT("end of test 1");
  */
  Walker_t::WFBuffer_t wbuffer;
  wbuffer.clear();
  registerData(P,wbuffer);
//     log = msd->evaluate(P,G,L);
  log0 = evaluate(P,G0,L0);
  PosType dr;
  dr[0] = 0.1;
  dr[1]=0.05;
  dr[2] = -0.01;
  PosType newpos(P.makeMove(iat,dr));
  app_log() <<"Testing ratio(P,dG,dL). \n";
  G=0;
  G0=0;
  L=0;
  log0 = ratioGrad(P,iat,G1);
  G0[iat]=G1;
  std::cout <<"Psi: " << log0 << std::endl;
  for(int i=0; i<n; i++)
  {
    std::cout <<i  <<"\n"
        <<"  x: " <<G[i][0]-G0[i][0] <<"  " <<G[i][0]   <<"\n"
        <<"  y: " <<G[i][1]-G0[i][1] <<"  " <<G[i][1] <<"\n"
        <<"  z: " <<G[i][2]-G0[i][2] <<"  " <<G[i][2] <<"\n"
        << std::endl;
  }
  std::cout << std::endl << std::endl;
  APP_ABORT("After MultiSlaterDeterminantFast::testMSD()");
}

/** Compute VGL of this MultiSlaterDeterminantFast
 *
 * THis is introduced to remove redundant code in 
 * - evaluate(P,G,L)
 * - evaluateLog(P,G,L,buf,fillbuffer)
 * Miguel's note: can this change over time??? I don't know yet
 */
WaveFunctionComponent::ValueType MultiSlaterDeterminantFast::evaluate_vgl_impl(ParticleSet& P
    , ParticleSet::ParticleGradient_t& g_tmp, ParticleSet::ParticleLaplacian_t& l_tmp)
{
  const ValueVector_t& detValues_up = Dets[0]->detValues;
  const ValueVector_t& detValues_dn = Dets[1]->detValues;
  const GradMatrix_t& grads_up = Dets[0]->grads;
  const GradMatrix_t& grads_dn = Dets[1]->grads;
  const ValueMatrix_t& lapls_up = Dets[0]->lapls;
  const ValueMatrix_t& lapls_dn = Dets[1]->lapls;
  const size_t N1 = Dets[0]->FirstIndex;
  const size_t N2 = Dets[1]->FirstIndex;
  const size_t NP1 = Dets[0]->NumPtcls;
  const size_t NP2 = Dets[1]->NumPtcls;
  CONSTEXPR ValueType czero(0);
  ValueType psi=czero;
  g_tmp=czero;
  l_tmp=czero;

  const RealType* restrict cptr=C->data();
  const size_t nc=C->size();
  const size_t* restrict upC=C2node_up->data();
  const size_t* restrict dnC=C2node_dn->data();
  for(size_t i=0; i<nc; ++i)
  {
    const RealType c=cptr[i];
    const size_t up=upC[i];
    const size_t down=dnC[i];
    psi += c*detValues_up[up]*detValues_dn[down];
    const ValueType c_up=c*detValues_dn[down];
    const ValueType c_dn=c*detValues_up[up];
    for(int k=0,n=N1; k<NP1; k++,n++)
    {
      g_tmp[n] += c_up*grads_up(up,k); //myG[n] += c*grads_up(up,k)*detValues_dn[down];
      l_tmp[n] += c_up*lapls_up(up,k); //myL[n] += c*lapls_up(up,k)*detValues_dn[down];
    }
    for(int k=0,n=N2; k<NP2; k++,n++)
    {
      g_tmp[n] += c_dn*grads_dn(down,k); //myG[n] += c*grads_dn(down,k)*detValues_up[up];
      l_tmp[n] += c_dn*lapls_dn(down,k); //myL[n] += c*lapls_dn(down,k)*detValues_up[up];
    }
  }
  ValueType psiinv = RealType(1)/psi;
  g_tmp *= psiinv;
  l_tmp *= psiinv;

  return psi;
}

WaveFunctionComponent::ValueType MultiSlaterDeterminantFast::evaluate(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  EvaluateTimer.start();
 
  Dets[0]->evaluateForWalkerMove(P);
  Dets[1]->evaluateForWalkerMove(P);

  psiCurrent=evaluate_vgl_impl(P,myG,myL);

  G += myG;
  for(size_t i=0; i<L.size(); i++)
    L[i] += myL[i] - dot(myG[i],myG[i]);

  EvaluateTimer.stop();
  return psiCurrent;
}

WaveFunctionComponent::RealType MultiSlaterDeterminantFast::evaluateLog(ParticleSet& P
    , ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  ValueType psi = evaluate(P,G,L);
  return LogValue = evaluateLogAndPhase(psi,PhaseValue);
}

WaveFunctionComponent::ValueType
MultiSlaterDeterminantFast::evalGrad_impl(ParticleSet& P, int iat, bool newpos, GradType& g_at)
{
  const bool upspin=(iat<FirstIndex_dn);
  const int spin0=(upspin)? 0: 1;
  const int spin1=(upspin)? 1: 0;

  if(newpos)
    Dets[spin0]->evaluateDetsAndGradsForPtclMove(P,iat);
  else
    Dets[spin0]->evaluateGrads(P,iat);

  const GradMatrix_t& grads = (newpos)? Dets[spin0]->new_grads:Dets[spin0]->grads;
  const ValueType *restrict detValues0 = (newpos)? Dets[spin0]->new_detValues.data(): Dets[spin0]->detValues.data();
  const ValueType *restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t *restrict det0=(upspin)? C2node_up->data():C2node_dn->data();
  const size_t *restrict det1=(upspin)? C2node_dn->data():C2node_up->data();
  const RealType *restrict cptr=C->data();
  const size_t nc=C->size();
  const size_t noffset=Dets[spin0]->FirstIndex;
  ValueType psi=ValueType(0);
  for(size_t i=0; i<nc; ++i)
  {
    const size_t d0=det0[i];
    //const size_t d1=det1[i];
    //psi +=  cptr[i]*detValues0[d0]        * detValues1[d1];
    //g_at += cptr[i]*grads(d0,iat-noffset) * detValues1[d1];
    const ValueType t=cptr[i]*detValues1[ det1[i] ];
    psi +=  t*detValues0[d0];
    g_at += t*grads(d0,iat-noffset);
  }
  return psi;
}

WaveFunctionComponent::GradType MultiSlaterDeterminantFast::evalGrad(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: evalGrad not implemented. \n");
  }
  CONSTEXPR RealType cone(1);
  GradType grad_iat;
  ValueType psi=evalGrad_impl(P,iat,false,grad_iat);;
  grad_iat*= (cone/psi);
  return grad_iat;
}

WaveFunctionComponent::ValueType MultiSlaterDeterminantFast::ratioGrad(ParticleSet& P
    , int iat, GradType& grad_iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratioGrad not implemented. \n");
  }
  UpdateMode=ORB_PBYP_PARTIAL;

  CONSTEXPR RealType cone(1);
  GradType dummy;
  ValueType psiNew=evalGrad_impl(P,iat,true,dummy);
  grad_iat+=(cone/psiNew)*dummy;
  curRatio=psiNew/psiCurrent;
  return curRatio;
}

WaveFunctionComponent::ValueType 
MultiSlaterDeterminantFast::ratio_impl(ParticleSet& P, int iat)
{
  const bool upspin=(iat<FirstIndex_dn);
  const int spin0=(upspin)? 0: 1;
  const int spin1=(upspin)? 1: 0;

  Dets[spin0]->evaluateDetsForPtclMove(P,iat);

  const ValueType *restrict detValues0 = Dets[spin0]->new_detValues.data(); //always new
  const ValueType *restrict detValues1 = Dets[spin1]->detValues.data();
  const size_t *restrict det0=(upspin)? C2node_up->data():C2node_dn->data();
  const size_t *restrict det1=(upspin)? C2node_dn->data():C2node_up->data();
  const RealType *restrict cptr=C->data();
  const size_t nc=C->size();

  ValueType psi=0;
  for(size_t i=0; i<nc; ++i)
    psi += cptr[i]*detValues0[ det0[i] ]*detValues1[ det1[i] ];
  return psi;
}

// use ci_node for this routine only
WaveFunctionComponent::ValueType MultiSlaterDeterminantFast::ratio(ParticleSet& P, int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: ratio not implemented. \n");
  }
  UpdateMode=ORB_PBYP_RATIO;
  ValueType psiNew=ratio_impl(P,iat);
  curRatio = psiNew/psiCurrent;
  return curRatio;
}

void MultiSlaterDeterminantFast::acceptMove(ParticleSet& P, int iat)
{
// this should depend on the type of update, ratio / ratioGrad
// for now is incorrect fot ratio(P,iat,dG,dL) updates
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: acceptMove not implemented. \n");
  }
// update psiCurrent,myG_temp,myL_temp
  AccRejTimer.start();
  psiCurrent *= curRatio;
  curRatio=1.0;

  Dets[iat>=nels_up]->acceptMove(P,iat);
  //Dets[DetID[iat]]->acceptMove(P,iat);

  AccRejTimer.stop();
}

void MultiSlaterDeterminantFast::restore(int iat)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }
  AccRejTimer.start();

  Dets[iat>=nels_up]->restore(iat);
  //Dets[DetID[iat]]->restore(iat);
  curRatio=1.0;
  AccRejTimer.stop();
}

void MultiSlaterDeterminantFast::registerData(ParticleSet& P, WFBufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: restore not implemented. \n");
  }

  Dets[0]->registerData(P,buf);
  Dets[1]->registerData(P,buf);

  buf.add(psiCurrent);
}

// FIX FIX FIX
WaveFunctionComponent::RealType MultiSlaterDeterminantFast::updateBuffer(ParticleSet& P, WFBufferType& buf, bool fromscratch)
{

  UpdateTimer.start();

  Dets[0]->updateBuffer(P,buf,fromscratch);
  Dets[1]->updateBuffer(P,buf,fromscratch);

  psiCurrent=evaluate_vgl_impl(P,myG,myL);
  
  P.G += myG;
  for(int i=0; i<P.L.size(); i++)
    P.L[i] += myL[i] - dot(myG[i],myG[i]);

  buf.put(psiCurrent);

  UpdateTimer.stop();
  return LogValue = evaluateLogAndPhase(psiCurrent,PhaseValue);
}

void MultiSlaterDeterminantFast::copyFromBuffer(ParticleSet& P, WFBufferType& buf)
{
  if(usingBF)
  {
    APP_ABORT("Fast MSD+BF: copyFromBuffer not implemented. \n");
  }
  Dets[0]->copyFromBuffer(P,buf);
  Dets[1]->copyFromBuffer(P,buf);

  buf.get(psiCurrent);
}


void MultiSlaterDeterminantFast::checkInVariables(opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
  {
    if(myVars->size())
      active.insertFrom(*myVars);
    else
      Optimizable=false;
  }
}

void MultiSlaterDeterminantFast::checkOutVariables(const opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
    myVars->getIndex(active);
}

/** resetParameters with optVariables
 *
 * USE_resetParameters
 */
void MultiSlaterDeterminantFast::resetParameters(const opt_variables_type& active)
{
  if(Optimizable && !IsCloned)
  {
    if(CIopt)
    {
      if(usingCSF)
      {
        RealType *restrict CSFcoeff_p=CSFcoeff->data();
        for(int i=0; i<CSFcoeff->size()-1; i++)
        {
          int loc=myVars->where(i);
          if(loc>=0)
          {
            CSFcoeff_p[i+1]= (*myVars)[i]=active[loc];
          }
        }
        int cnt=0;
        RealType *restrict C_p=C->data();
        const RealType *restrict CSFexpansion_p=CSFexpansion->data();
        for(int i=0; i<DetsPerCSF->size(); i++)
        {
          for(int k=0; k<(*DetsPerCSF)[i]; k++)
          {
            C_p[cnt] = CSFcoeff_p[i]*CSFexpansion_p[cnt];
            cnt++;
          }
        }
        //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
      }
      else
      {
        RealType *restrict C_p=C->data();
        for(int i=0; i<C->size()-1; i++)
        {
          int loc=myVars->where(i);
          if(loc>=0)
          {
            C_p[i+1]=(*myVars)[i]=active[loc];
          }
        }
        //for(int i=0; i<Dets.size(); i++) Dets[i]->resetParameters(active);
      }
    }
    if(Orbopt)
    {
      // read out the parameters for spin up electrons that define the rotation into an antisymmetric matrix
      std::vector<RealType> rot_mat_up(m_nmo_up*m_nmo_up, 0.0);
      for (int j=0, i = m_first_var_pos; i < m_first_var_pos + m_act_rot_inds_up.size(); j++, i++) 
      {
        int loc=myVars->where(i);

        const int p = m_act_rot_inds_up[j].first;
        const int q = m_act_rot_inds_up[j].second;
        // m_first_var_pos is the index to the first parameter of the spin up electrons...
        const RealType x = (*myVars)[i] = active[loc];
        //app_log() <<"active ["<< loc << "] = "<< x<< "\n";
        rot_mat_up[p+q*m_nmo_up] =  x;
        rot_mat_up[q+p*m_nmo_up] = -x;
      }

      this->exponentiate_matrix(m_nmo_up, &rot_mat_up[0]);

#if defined(ENABLE_SOA)
      BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), static_cast<qmcplusplus::LCAOrbitalSet*>(Dets[0]->Phi)->C->data() , m_nb_up);
      BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), static_cast<qmcplusplus::LCAOrbitalSet*>(Dets[1]->Phi)->C->data() , m_nb_up);
#else
      BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), (Dets[0]->Phi)->C->data() , m_nb_up);
      BLAS::gemm('N','T', m_nb_up, m_nmo_up, m_nmo_up, RealType(1.0), m_init_B_up.data(), m_nb_up, &rot_mat_up[0],m_nmo_up, RealType(0.0), (Dets[1]->Phi)->C->data() , m_nb_up);
#endif      
    }
  }
}
void MultiSlaterDeterminantFast::reportStatus(std::ostream& os)
{
}


void MultiSlaterDeterminantFast::evaluateDerivatives(ParticleSet& P,
    const opt_variables_type& optvars,
    std::vector<RealType>& dlogpsi,
    std::vector<RealType>& dhpsioverpsi)
{
  bool recalculate(false);
  for (int k=0; k<myVars->size(); ++k)
  {
    int kk=myVars->where(k);
    if (kk<0)
      continue;
    if (optvars.recompute(kk))
      recalculate=true;
  }
// need to modify for CSF later on, right now assume Slater Det basis
  if (recalculate)
  {
    if(Optimizable)
    {
      if(CIopt)
      {
        if(usingCSF)
        {
          if(laplSum_up.size() == 0)
            laplSum_up.resize(Dets[0]->detValues.size());
          if(laplSum_dn.size() == 0)
            laplSum_dn.resize(Dets[1]->detValues.size());
          // assume that evaluateLog has been called in opt routine before
          //   Dets[0]->evaluateForWalkerMove(P);
          //   Dets[1]->evaluateForWalkerMove(P);
          ValueVector_t& detValues_up = Dets[0]->detValues;
          ValueVector_t& detValues_dn = Dets[1]->detValues;
          GradMatrix_t& grads_up = Dets[0]->grads;
          GradMatrix_t& grads_dn = Dets[1]->grads;
          ValueMatrix_t& lapls_up = Dets[0]->lapls;
          ValueMatrix_t& lapls_dn = Dets[1]->lapls;
          const size_t N1 = Dets[0]->FirstIndex;
          const size_t N2 = Dets[1]->FirstIndex;
          const size_t NP1 = Dets[0]->NumPtcls;
          const size_t NP2 = Dets[1]->NumPtcls;
    // myG,myL should already be calculated
          const size_t n = P.getTotalNum();

          ValueType psiinv = (RealType)1.0/psiCurrent;
          ValueType lapl_sum=0.0;
          ValueType gg=0.0, ggP=0.0;
          myG_temp=0.0;
          int num=laplSum_up.size();
          ValueVector_t::iterator it(laplSum_up.begin());
          ValueVector_t::iterator last(laplSum_up.end());
          ValueType* ptr0 = lapls_up[0];
          while(it != last)
          {
            (*it)=0.0;
            for(int k=0; k<nels_up; k++,ptr0++)
              (*it) += *ptr0;
            it++;
          }
          it=laplSum_dn.begin();
          last=laplSum_dn.end();
          ptr0 = lapls_dn[0];
          while(it != last)
          {
            (*it)=0.0;
            for(int k=0; k<nels_dn; k++,ptr0++)
              (*it) += *ptr0;
            it++;
          }

          const RealType *restrict C_p=C->data();
          for(size_t i=0; i<C->size(); i++)
          {
            size_t upC = (*C2node_up)[i];
            size_t dnC = (*C2node_dn)[i];
            ValueType tmp1 = C_p[i]*detValues_dn[dnC]*psiinv;
            ValueType tmp2 = C_p[i]*detValues_up[upC]*psiinv;
            lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
            for(size_t k=0,j=N1; k<NP1; k++,j++)
              myG_temp[j] += tmp1*grads_up(upC,k);
            for(size_t k=0,j=N2; k<NP2; k++,j++)
              myG_temp[j] += tmp2*grads_dn(dnC,k);
          }
          gg=ggP=0.0;
          for(size_t i=0; i<n; i++)
          {
            gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
          }
    //       for(int i=0; i<C.size(); i++)
          num=CSFcoeff->size()-1;
          int cnt=0;
    //        this one is not optable
          cnt+=(*DetsPerCSF)[0];
          int ip(1);
          for(int i=0; i<num; i++,ip++)
          {
            int kk=myVars->where(i);
            if (kk<0)
            {
              cnt+=(*DetsPerCSF)[ip];
              continue;
            }
            ValueType cdet=0.0,q0=0.0,v1=0.0,v2=0.0;
            const RealType *restrict CSFexpansion_p=CSFexpansion->data();
            for(int k=0; k<(*DetsPerCSF)[ip]; k++)
            {
              size_t upC = (*C2node_up)[cnt];
              size_t dnC = (*C2node_dn)[cnt];
              ValueType tmp1=CSFexpansion_p[cnt]*detValues_dn[dnC]*psiinv;
              ValueType tmp2=CSFexpansion_p[cnt]*detValues_up[upC]*psiinv;
              cdet+=CSFexpansion_p[cnt]*detValues_up[upC]*detValues_dn[dnC]*psiinv;
              q0 += (tmp1*laplSum_up[upC] + tmp2*laplSum_dn[dnC]);
              for(size_t l=0,j=N1; l<NP1; l++,j++)
                v1 += tmp1*static_cast<ValueType>(dot(P.G[j],grads_up(upC,l))-dot(myG_temp[j],grads_up(upC,l)));
              for(size_t l=0,j=N2; l<NP2; l++,j++)
                v2 += tmp2*static_cast<ValueType>(dot(P.G[j],grads_dn(dnC,l))-dot(myG_temp[j],grads_dn(dnC,l)));
              cnt++;
            }
            convert(cdet,dlogpsi[kk]);
            ValueType dhpsi =  (RealType)-0.5*(q0-cdet*lapl_sum)
                               -cdet*gg-v1-v2;
            //ValueType dhpsi =  -0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
            //                         -cdet*lapl_sum)
            //                   -cdet*gg-(tmp1*v1+tmp2*v2);
            convert(dhpsi,dhpsioverpsi[kk]);
          }
        }
        else
          //usingCSF
        {
          if(laplSum_up.size() == 0)
            laplSum_up.resize(Dets[0]->detValues.size());
          if(laplSum_dn.size() == 0)
            laplSum_dn.resize(Dets[1]->detValues.size());
          // assume that evaluateLog has been called in opt routine before
          //   Dets[0]->evaluateForWalkerMove(P);
          //   Dets[1]->evaluateForWalkerMove(P);
          ValueVector_t& detValues_up = Dets[0]->detValues;
          ValueVector_t& detValues_dn = Dets[1]->detValues;
          GradMatrix_t& grads_up = Dets[0]->grads;
          GradMatrix_t& grads_dn = Dets[1]->grads;
          ValueMatrix_t& lapls_up = Dets[0]->lapls;
          ValueMatrix_t& lapls_dn = Dets[1]->lapls;
          int N1 = Dets[0]->FirstIndex;
          int N2 = Dets[1]->FirstIndex;
          int NP1 = Dets[0]->NumPtcls;
          int NP2 = Dets[1]->NumPtcls;
          int n = P.getTotalNum();
          ValueType psiinv = (RealType)1.0/psiCurrent;
          ValueType lapl_sum=0.0;
          ValueType gg=0.0, ggP=0.0;
          myG_temp=0.0;
          int num=laplSum_up.size();
          ValueVector_t::iterator it(laplSum_up.begin());
          ValueVector_t::iterator last(laplSum_up.end());
          ValueType* ptr0 = lapls_up[0];
          while(it != last)
          {
            (*it)=0.0;
            for(int k=0; k<nels_up; k++,ptr0++)
              (*it) += *ptr0;
            it++;
          }
          it=laplSum_dn.begin();
          last=laplSum_dn.end();
          ptr0 = lapls_dn[0];
          while(it != last)
          {
            (*it)=0.0;
            for(size_t k=0; k<nels_dn; k++,ptr0++)
              (*it) += *ptr0;
            it++;
          }
          const RealType *restrict C_p=C->data();
          for(size_t i=0; i<C->size(); i++)
          {
            size_t upC = (*C2node_up)[i];
            size_t dnC = (*C2node_dn)[i];
            ValueType tmp1 = C_p[i]*detValues_dn[dnC]*psiinv;
            ValueType tmp2 = C_p[i]*detValues_up[upC]*psiinv;
            lapl_sum += tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC];
            for(size_t k=0,j=N1; k<NP1; k++,j++)
              myG_temp[j] += tmp1*grads_up(upC,k);
            for(size_t k=0,j=N2; k<NP2; k++,j++)
              myG_temp[j] += tmp2*grads_dn(dnC,k);
          }
          gg=ggP=0.0;
          for(size_t i=0; i<n; i++)
          {
            gg += dot(myG_temp[i],myG_temp[i])-dot(P.G[i],myG_temp[i]);
          }
          for(size_t i=1; i<C->size(); i++)
          {
            int kk=myVars->where(i-1);
            if (kk<0)
              continue;
            const size_t upC = (*C2node_up)[i];
            const size_t dnC = (*C2node_dn)[i];
            ValueType cdet=detValues_up[upC]*detValues_dn[dnC]*psiinv;
            ValueType tmp1=detValues_dn[dnC]*psiinv;
            ValueType tmp2=detValues_up[upC]*psiinv;
            convert(cdet,dlogpsi[kk]);
            ValueType v1=0.0,v2=0.0;
            for(size_t k=0,j=N1; k<NP1; k++,j++)
              v1 += (dot(P.G[j],grads_up(upC,k))-dot(myG_temp[j],grads_up(upC,k)) );
            for(size_t k=0,j=N2; k<NP2; k++,j++)
              v2 += (dot(P.G[j],grads_dn(dnC,k))-dot(myG_temp[j],grads_dn(dnC,k)));
            ValueType dhpsi =  (RealType)-0.5*(tmp1*laplSum_up[upC]+tmp2*laplSum_dn[dnC]
                                     -cdet*lapl_sum)
                               -cdet*gg-(tmp1*v1+tmp2*v2);
            convert(dhpsi,dhpsioverpsi[kk]);
          }
        } // usingCSF
      }
      if(Orbopt)
      {
        const ValueVector_t& detValues_up = Dets[0]->detValues;
        const ValueVector_t& detValues_dn = Dets[1]->detValues;
        const GradMatrix_t& grads_up = Dets[0]->grads;
        const GradMatrix_t& grads_dn = Dets[1]->grads;
        const ValueMatrix_t& lapls_up = Dets[0]->lapls;
        const ValueMatrix_t& lapls_dn = Dets[1]->lapls;

        const size_t N1 = Dets[0]->FirstIndex;
        const size_t N2 = Dets[1]->FirstIndex;
        const size_t NP1 = Dets[0]->NumPtcls;
        const size_t NP2 = Dets[1]->NumPtcls;
  //      psiCurrent=0.0;
        myG_temp=0.0;
        myL_temp=0.0;
        myG_J=0.0;
        myL_J=0.0;

        const RealType *restrict C_p=C->data();
        for(int i=0; i<C->size(); i++)
        {
            const size_t upC = (*C2node_up)[i];
            const size_t dnC = (*C2node_dn)[i];
            const ValueType tmp1 = C_p[i]*detValues_dn[dnC];
            const ValueType tmp2 = C_p[i]*detValues_up[upC];
  //          psiCurrent+= (C[i]*detValues_up[upC]*detValues_dn[dnC]);
            for(size_t k=0,j=N1; k<NP1; k++,j++)
            {
              myG_temp[j] += tmp1*grads_up(upC,k);
              myL_temp[j] += tmp1*lapls_up(upC,k);
            }
            for(size_t k=0,j=N2; k<NP2; k++,j++)
            {
              myG_temp[j] += tmp2*grads_dn(dnC,k);
              myL_temp[j] += tmp2*lapls_dn(dnC,k);
            }
        }

        myG_temp *= (1/psiCurrent);
        myL_temp *= (1/psiCurrent);

        // calculation of myG_J which will be used to represent \frac{\nabla\psi_{J}}{\psi_{J}} 
        // calculation of myL_J will be used to represent \frac{\nabla^2\psi_{J}}{\psi_{J}}
        // IMPORTANT NOTE:  The value of P.L holds \nabla^2 ln[\psi] but we need  \frac{\nabla^2 \psi}{\psi} and this is what myL_J will hold 
        for(int iat=0; iat<(myL_temp.size()); iat++)
        {
          myG_J[iat] = (P.G[iat] - myG_temp[iat]);
          myL_J[iat] = (P.L[iat] + dot(P.G[iat],P.G[iat]) - myL_temp[iat]);
        }

//IMPORTANT NOTE: THE Dets[0]->psiMinv OBJECT DOES NOT HOLD THE INVERSE IF THE MULTIDIRACDETERMINANTBASE ONLY CONTAINES ONE ELECTRON. NEED A FIX FOR THIS CASE
        // The T matrix should be calculated and stored for use      
        // T = A^{-1} \widetilde A
        //REMINDER: that the ValueMatrix_t "matrix" stores data in a row major order and that BLAS commands assume column major
        BLAS::gemm('N','N',m_nmo_up, nels_up, nels_up, RealType(1.0), Dets[0]->psiM.data(), m_nmo_up, Dets[0]->psiMinv.data(), nels_up, RealType(0.0), T_up.data(), m_nmo_up);
        BLAS::gemm('N','N',m_nmo_dn, nels_dn, nels_dn, RealType(1.0), Dets[1]->psiM.data(), m_nmo_dn, Dets[1]->psiMinv.data(), nels_dn, RealType(0.0), T_dn.data(), m_nmo_dn);

        table_method_eval(myL_J,
                          myG_J,
                          nels_up,
                          m_nmo_up,
                          T_up.data(),
                          Dets[0]->psiM.data(),
                          Dets[0]->psiMinv.data(),
                          dlogpsi,
                          dhpsioverpsi,
                          m_first_var_pos,
                          m_act_rot_inds_up.size(),
                          &m_act_rot_inds_up,
                          0,
                          T_up,
                          Dets[0]->psiM);

        if (false)
        {
          for (int i=0; i<m_act_rot_inds_up.size(); i++){
           int kk= i+m_first_var_pos;
               const int p = m_act_rot_inds_up[i].first;
               const int q = m_act_rot_inds_up[i].second;
               std::vector<char> buff(1000, ' ');
               const int len = std::sprintf(&buff[0], " p = %4i   q = %4i     dlogpsi = %20.12f     dhpsioverpsi = %20.12f   kk = %4i", p, q,dlogpsi[kk], dhpsioverpsi[kk], kk);
               for (int k = 0; k < len; k++)
                 app_log() << buff[k];
               app_log() << std::endl;
          }
        }

        table_method_eval(myL_J,
                          myG_J,
                          nels_dn,
                          m_nmo_dn,
                          T_dn.data(),
                          Dets[1]->psiM.data(),
                          Dets[1]->psiMinv.data(),
                          dlogpsi,
                          dhpsioverpsi,
                          m_first_var_pos,
                          m_act_rot_inds_up.size(),
                          &m_act_rot_inds_up,
                          1,
                          T_dn,
                          Dets[1]->psiM);
        if (false)
        {
          for (int i=0; i<m_act_rot_inds_up.size(); i++){
           int kk= i+m_first_var_pos;
               const int p = m_act_rot_inds_up[i].first;
               const int q = m_act_rot_inds_up[i].second;
               std::vector<char> buff(1000, ' ');
               const int len = std::sprintf(&buff[0], " p = %4i   q = %4i     dlogpsi = %20.12f     dhpsioverpsi = %20.12f   kk = %4i", p, q,dlogpsi[kk], dhpsioverpsi[kk], kk);
               for (int k = 0; k < len; k++)
                 app_log() << buff[k];
               app_log() << std::endl;
          }
        }
      }
    }
  }
}

void MultiSlaterDeterminantFast::registerTimers()
{
  RatioTimer.reset();
  RatioGradTimer.reset();
  RatioAllTimer.reset();
  Ratio1Timer.reset();
  Ratio1GradTimer.reset();
  Ratio1AllTimer.reset();
  UpdateTimer.reset();
  EvaluateTimer.reset();
  AccRejTimer.reset();
  TimerManager.addTimer (&RatioTimer);
  TimerManager.addTimer (&RatioGradTimer);
  TimerManager.addTimer (&RatioAllTimer);
  TimerManager.addTimer (&Ratio1Timer);
  TimerManager.addTimer (&Ratio1GradTimer);
  TimerManager.addTimer (&Ratio1AllTimer);
  TimerManager.addTimer (&UpdateTimer);
  TimerManager.addTimer (&EvaluateTimer);
  TimerManager.addTimer (&AccRejTimer);
}


void MultiSlaterDeterminantFast::exponentiate_matrix(const int n, RealType * const mat)
{
  // save initial matrix and get some workspace
  std::vector<RealType> mat_a(mat, mat + n*n);
  std::vector<RealType> mat_b(mat, mat + n*n);
  std::vector<RealType> mat_c(mat, mat + n*n);
  RealType * ptr_a = &mat_a[0];
  RealType * ptr_b = &mat_b[0];
  RealType * ptr_c = &mat_c[0];

  // initialize output to identity matrix
  for (int j = 0; j < n; j++)
  for (int i = 0; i < n; i++)
    mat[i+n*j] = ( i == j ? 1.0 : 0.0 );

  // compute exponential of matrix
  for (int q = 1; q < 20; q++) {
    BLAS::axpy(n*n, RealType(1.0), ptr_b, mat);
    BLAS::gemm('N', 'N', n, n, n, RealType(1.0) / (q+1), ptr_a, n, ptr_b, n, RealType(0.0), ptr_c, n);
    std::swap(ptr_b, ptr_c);
    RealType max_elem = 0.0;
    for (int i = 0; i < n*n; i++)
      if ( std::abs(ptr_b[i]) > max_elem )
        max_elem = std::abs(ptr_b[i]);
    if ( max_elem < 1.0e-15 )
      break;
  }
}

int MultiSlaterDeterminantFast::build_occ_vec(std::vector<int> * data,
                                              const size_t& nel,
                                              const size_t& nmo,
                                              std::vector<int>* occ_vec)
{
  std::vector<int>::iterator it = (*data).begin();
  int count = 0; //number of determinants
  while(it != (*data).end())
  {
    int k = *it; // number of excitations with respect to the reference matrix  
    if(count == 0)
    {
      it += 3*k+1;
      count ++;
    }
    else
    {
      for (int i = 0; i<k; i++)
      {
      //for determining active orbitals
        (*occ_vec)[*(it+1+i)]++;
        (*occ_vec)[*(it+1+k+i)]++;
      }
      it += 3*k+1;
      count ++;
    }
  }
  return count;
}

void MultiSlaterDeterminantFast::table_method_eval(const ParticleSet::ParticleLaplacian_t& myL_J,
                                                   const ParticleSet::ParticleGradient_t& myG_J,
                                                   const size_t nel,
                                                   const size_t nmo,
                                                   double* T,
                                                   double* A,
                                                   double* Ainv,
                                                   std::vector<RealType>& dlogpsi , 
                                                   std::vector<RealType>& dhpsioverpsi, 
                                                   const int parameter_start_index,
                                                   const int parameters_size ,
                                                   const std::vector<std::pair<int,int>>* const m_act_rot_inds,
                                                   const int active_spin,
                                                   const SPOSetBase::ValueMatrix_t& Tr,
                                                   const SPOSetBase::ValueMatrix_t& Ar)
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
GUIDE TO THE MATICES BEING BUILT
----------------------------------------------
The idea here is that there is a loop over all unique determinants. For each determiant the table method is employed to calculate the contributions to the parameter derivatives (dhpsioverpsi/dlogpsi)

  loop through unquie determinants 
    loop through parameters
      evaluate contributaion to dlogpsi and dhpsioverpsi
\noindent 

  BLAS GUIDE  for matrix multiplication of  [  alpha * A.B + beta * C = C ]
  Matrix A is of dimensions a1,a2 and Matrix B is b1,b2   in which a2=b1
  The BLAS command is as follows...

 BLAS::gemm('N','N', b2, a1, a2 ,alpha, B, b2, A, a2, beta, C, b2);

Below is a human readable format for the matrix multiplications performed below...

This notation is inspired by http://dx.doi.org/10.1063/1.4948778
\newline
\hfill\break
$
    A_{i,j}=\phi_j(r_{i}) \\
    T = A^{-1} \widetilde{A} \\
    B_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i}) \\
    \hat{O_{I}} = \hat{O}D_{I} \\
    D_{I}=det(A_{I}) \newline 
    \psi_{MS} = \sum_{I=0} C_{I} D_{I\uparrow}D_{I\downarrow} \\
    \Psi_{total} = \psi_{J}\psi_{MS} \\
    \alpha_{I} = P^{T}_{I}TQ_{I} \\
    M_{I} = P^{T}_{I} \widetilde{M} Q_{I} = P^{T}_{I} (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )Q_{I} \\
$
\newline
There are three constants I use in the expressions for dhpsioverpsi and dlogpsi
\newline
\hfill\break
$
  const0 = C_{0}*det(A_{0\downarrow})+\sum_{I=1} C_{I}*det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{I=1} C_{I}*\hat{O}det(A_{I\downarrow})* det(\alpha_{I\uparrow}) \\
  const2 = \sum_{I=1} C_{I}*det(A_{I\downarrow})* Tr[\alpha_{I}^{-1}M_{I}]*det(\alpha_{I}) \\
$
\newline
Below is a translation of the shorthand I use to represent matrices independent of ``excitation matrix".
\newline
\hfill\break
$
    Y_{1} =  A^{-1}B   \\
    Y_{2} = A^{-1}BA^{-1}\widetilde{A} \\
    Y_{3} = A^{-1}\widetilde{B} \\
    Y_{4} = \widetilde{M} = (A^{-1}\widetilde{B} - A^{-1} B A^{-1}\widetilde{A} )\\
$
\newline
Below is a translation of the shorthand I use to represent matrices dependent on ``excitation" with respect to the reference Matrix and sums of matrices. Above this line I have represented these excitation matrices with a subscript ``I" but from this point on The subscript will be omitted and it is clear that whenever a matrix depends on $P^{T}_I$ and $Q_{I}$ that this is an excitation matrix. The reference matrix is always $A_{0}$ and is always the Hartree Fock Matrix.
\newline
\hfill\break
$
    Y_{5} = TQ \\
    Y_{6} = (P^{T}TQ)^{-1} = \alpha_{I}^{-1}\\
    Y_{7} = \alpha_{I}^{-1} P^{T} \\
    Y_{11} = \widetilde{M}Q \\
    Y_{23} = P^{T}\widetilde{M}Q \\
    Y_{24} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q \\
    Y_{25} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1} \\
    Y_{26} = \alpha_{I}^{-1}P^{T}\widetilde{M}Q\alpha_{I}^{-1}P^{T}\\
$
\newline
So far you will notice that I have not included up or down arrows to specify what spin the matrices are of. This is because we are calculating the derivative of all up or all down spin orbital rotation parameters at a time. If we are finding the up spin derivatives then any term that is down spin will be constant. The following assumes that we are taking up-spin MO rotation parameter derivatives. Of course the down spin expression can be retrieved by swapping the up and down arrows. I have dubbed any expression with lowercase p prefix as a "precursor" to an expression actually used...
\newline
\hfill\break
$
    \dot{C_{I}} = C_{I}*det(A_{I\downarrow})\\
    \ddot{C_{I}} = C_{I}*\hat{O}det(A_{I\downarrow}) \\
    pK1 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}) \\
    pK2 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK3 = \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK4 = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}) \\
    pK5 = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T}) \\
$
\newline
Now these p matrices will be used to make various expressions via BLAS commands.
\newline
\hfill\break
$
    K1T = const0^{-1}*pK1.T =const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (Q\alpha_{I}^{-1}P^{T}T) \\
    TK1T = T.K1T = const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) Tr[\alpha_{I}^{-1}M_{I}] (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    TK2AiB = T.K2AiB = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}A^{-1}\widetilde{B})\\
    K2XA =  const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}X\widetilde{A})\\
    TK2XA = T.K2XA = const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}X\widetilde{A})\\ \\
    K2T = \frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK2T = T.K2T =\frac{const1}{const0^{2}} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\
    MK2T = \frac{const0}{const1} Y_{4}.K2T= const0^{-1}  \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (\widetilde{M}Q\alpha_{I}^{-1}P^{T}T)\\ \\
    K3T = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK3T = T.K3T  = const0^{-1}  \sum_{I=1} \ddot{C_{I}} det(\alpha_{I}) (TQ\alpha_{I}^{-1}P^{T}T)\\ \\
    K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (Q\alpha_{I}^{-1}P^{T}T) \\
    TK4T = T.K4T = \sum_{I=1} \dot{C_{I}} det(A_{I}) (TQ\alpha_{I}^{-1}P^{T}T) \\ \\
    K5T =  const0^{-1} \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
    TK5T = T.K5T  = \sum_{I=1} \dot{C_{I}} det(\alpha_{I}) (T Q\alpha_{I}^{-1} M_{I} \alpha_{I}^{-1}P^{T} T)  \\
$
\newline
Now with all these matrices and constants the expressions of dhpsioverpsi and dlogpsi can be created.




In addition I will be using a special generalization of the kinetic operator which I will denote as O. Our Slater matrix with the special O operator applied to each element will be called B_bar

$
``Bbar"_{i,j} =\nabla^2 \phi_{j}(r_i) + \frac{\nabla_{i}J}{J} \cdot \nabla \phi_{j}(r_{i})  + \frac{\nabla^2_i J}{J} \phi_{j}(r_{i})
$
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
{
  const size_t num_unique_up_dets       (active_spin==0 ? unique_ups         : unique_dns);
  const size_t num_unique_dn_dets       (active_spin==0 ? unique_dns         : unique_ups);
  const ValueVector_t& detValues_up     (active_spin==0 ? Dets[0]->detValues : Dets[1]->detValues );
  const ValueVector_t& detValues_dn     (active_spin==0 ? Dets[1]->detValues : Dets[0]->detValues );
  const ValueMatrix_t& lapls_dn         (active_spin==0 ? Dets[1]->lapls     : Dets[0]->lapls );
  const GradMatrix_t& grads_dn          (active_spin==0 ? Dets[1]->grads     : Dets[0]->grads );
  const RealType* restrict cptr = C->data();
  const size_t nc = C->size();
  const size_t* restrict upC            (active_spin==0 ? C2node_up->data() : C2node_dn->data() );
  const size_t* restrict dnC            (active_spin==0 ? C2node_dn->data() : C2node_up->data() );
  //B_grad holds the gardient operator
  const SPOSetBase::GradMatrix_t& B_grad     (active_spin==0 ? Dets[0]->dpsiM : Dets[1]->dpsiM );
  //B_til holds the laplacian operator
  const SPOSetBase::ValueMatrix_t& B_til     (active_spin==0 ? Dets[0]->d2psiM : Dets[1]->d2psiM );
  //B_bar will hold our special O operator
  SPOSetBase::ValueMatrix_t Bbar;
  Bbar.resize(nel,nmo);

  const size_t N1  = Dets[0]->FirstIndex;
  const size_t N2  = Dets[1]->FirstIndex;  
  const size_t NP1 = Dets[0]->NumPtcls;
  const size_t NP2 = Dets[1]->NumPtcls;

  const int offset1 (active_spin==0 ? N1  : N2);
  const int offset2 (active_spin==0 ? N2  : N1);   
  const int NPother (active_spin==0 ? NP2 : NP1);

// possibly replace wit BLAS calls 
  for(int i=0; i<nel; i++)
    for(int j=0; j<nmo; j++)
      Bbar(i,j) = B_til(i,j) + 2*dot(myG_J[i+offset1], B_grad(i,j)) + myL_J[i+offset1]*Ar(i,j);

  const double* restrict B(Bbar.data());

  //Need to create the constants: (Oi, const0, const1, const2)to take advantage of minimal BLAS commands; 
  //Oi is the special operator applied to the slater matrix "A subscript i" from the total CI expansion
  //\hat{O_{i}} = \hat{O}D_{i} with D_{i}=det(A_{i}) and Multi-Slater component defined as \sum_{i=0} C_{i} D_{i\uparrow}D_{i\downarrow}
  std::vector<RealType> Oi(num_unique_dn_dets); 

  for (int index=0; index < num_unique_dn_dets; index++)
    for (int iat=0; iat < NPother; iat++)
      Oi[index]+= lapls_dn(index,iat) + 2*dot(grads_dn(index,iat),myG_J[offset2 + iat]) + myL_J[offset2 + iat]*detValues_dn[index];

  //const0 = C_{0}*det(A_{0\downarrow})+\sum_{i=1} C_{i}*det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const1 = C_{0}*\hat{O} det(A_{0\downarrow})+\sum_{i=1} C_{i}*\hat{O}det(A_{i\downarrow})* det(\alpha_{i\uparrow})
  //const2 = \sum_{i=1} C_{i}*det(A_{i\downarrow})* Tr[\alpha_{i}^{-1}M_{i}]*det(\alpha_{i})
  RealType const0(0.0),const1(0.0), const2(0.0);
  for(size_t i=0; i<nc; ++i)
  {
    const RealType c=cptr[i];
    const size_t up=upC[i];
    const size_t down=dnC[i];

    const0 += c * detValues_dn[down] * (detValues_up[up] / detValues_up[0]);
    const1 += c * Oi[down] * (detValues_up[up] / detValues_up[0]);
  }


  SPOSetBase::ValueMatrix_t Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y11,Y23,Y24,Y25,Y26;
  Y1.resize(nel,nel);
  Y2.resize(nel,nmo);
  Y3.resize(nel,nmo);
  Y4.resize(nel,nmo);

  SPOSetBase::ValueMatrix_t pK1,K1T,TK1T, pK2,K2AiB,TK2AiB,K2XA,TK2XA,K2T,TK2T,MK2T, pK3,K3T,TK3T, pK4,K4T,TK4T, pK5,K5T,TK5T;
  pK1.resize(nmo,nel);
  K1T.resize(nmo,nmo);
  TK1T.resize(nel,nmo);

  pK2.resize(nmo,nel);
  K2AiB.resize(nmo,nmo);
  TK2AiB.resize(nel,nmo);
  K2XA.resize(nmo,nmo);
  TK2XA.resize(nel,nmo);
  K2T.resize(nmo,nmo);
  TK2T.resize(nel,nmo);
  MK2T.resize(nel,nmo);

  pK3.resize(nmo,nel);
  K3T.resize(nmo,nmo);
  TK3T.resize(nel,nmo);

  pK4.resize(nmo,nel);
  K4T.resize(nmo,nmo);
  TK4T.resize(nel,nmo);

  pK5.resize(nmo,nel);
  K5T.resize(nmo,nmo);
  TK5T.resize(nel,nmo);

  BLAS::gemm('N','N', nel, nel, nel, RealType(1.0),    B, nmo,       Ainv, nel, RealType(0.0),  Y1.data(), nel);
  BLAS::gemm('N','N', nmo, nel, nel, RealType(1.0),    T, nmo,  Y1.data(), nel, RealType(0.0),  Y2.data(), nmo);
  BLAS::gemm('N','N', nmo, nel, nel, RealType(1.0),    B, nmo,       Ainv, nel, RealType(0.0),  Y3.data(), nmo);

  //possibly replace with BLAS call
  Y4 = Y3 - Y2;

  //Now we are going to loop through all unique determinants.
  //The few lines above are for the reference matrix contribution.
  //Although I start the loop below from index 0, the loop only performs actions when the index is => 1
  //the detData object contains all the information about the P^T and Q matrices (projection matrices) needed in the table method
  const int* restrict data_it = Dets[active_spin]->detData->data();
//  std::vector<int>::iterator data_it = Dets[active_spin]->detData->begin();
  for(int index=0, datum=0; index < num_unique_up_dets; index++)
  {
//    const int  k = *data_it;
    const int  k = data_it[datum];

    if (k==0)
    {
//      data_it += 3*k+1;
      datum += 3*k+1;
    }

    else
    {
      //Number of rows and cols of P^T
      const int prows=k;
      const int pcols=nel;
      //Number of rows and cols of Q
      const int qrows=nmo;
      const int qcols=k;
      
      Y5.resize(nel,k);
      Y6.resize(k,k);

      //Any matrix multiplication of P^T or Q is simply a projection
      //Explicit matrix multiplication can be avoided; instead column or row copying can be done
      //BlAS::copy(size of col/row being copied,
      //           Matrix pointer + place to begin copying,
      //           storage spacing (number of elements btw next row/col element),
      //           Pointer to resultant matrix + place to begin pasting,
      //           storage spacing of resultant matrix)
      //For example the next 4 lines is the matrix multiplication of T*Q = Y5      
      std::fill(Y5.begin(),Y5.end(),0.0);
      for ( int i=0; i<k; i++)
      {
//        BLAS::copy(nel, T + *(data_it+1+k+i), nmo, Y5.data()+i, k);   
        BLAS::copy(nel, T + data_it[datum+1+k+i], nmo, Y5.data()+i, k);   
      }

      std::fill(Y6.begin(),Y6.end(),0.0);
      for ( int i=0; i<k; i++)
      { 
//        BLAS::copy(k, Y5.data() + (*(data_it+1+i))*k, 1, (Y6.data() + i*k), 1);   
        BLAS::copy(k, Y5.data() + (data_it[datum+1+i])*k, 1, (Y6.data() + i*k), 1);   
      }

//      app_log() << "priting Y5 \n " << Y5 << std::endl;
//      app_log() << "priting Y6 \n " << Y6 << std::endl;

      Vector<ValueType> WS;
      Vector<IndexType> Piv;
      WS.resize(k);
      Piv.resize(k);
      // not sure if i need the following 2 lines
  //     std::fill(WS.begin(),WS.end(),0.0);
  //     std::fill(Piv.begin(),Piv.end(),0);
      RealType PhaseR=0.0;
      InvertWithLog(Y6.data(),k,k,WS.data(),Piv.data(),PhaseR);
//      app_log() << "priting inv Y6 \n " << Y6 << std::endl;

      Y11.resize(nel,  k);  
      Y23.resize(  k,  k);    
      Y24.resize(  k,  k);    
      Y25.resize(  k,  k);
      Y26.resize(  k,nel);

      std::fill(Y11.begin(),Y11.end(),0.0);
      for ( int i=0; i<k; i++)
      {
//        BLAS::copy(nel, Y4.data() + *(data_it+1+k+i), nmo, Y11.data()+i, k);   
        BLAS::copy(nel, Y4.data() + (data_it[datum+1+k+i]), nmo, Y11.data()+i, k);   
      }
//      app_log() << "priting Y11 \n " << Y11 << std::endl;

      std::fill(Y23.begin(),Y23.end(),0.0);
      for ( int i=0; i<k; i++)
      { 
//        BLAS::copy(k, Y11.data() + (*(data_it+1+i))*k, 1, (Y23.data() + i*k), 1);   
        BLAS::copy(k, Y11.data() + (data_it[datum+1+i])*k, 1, (Y23.data() + i*k), 1);   
      }
//      app_log() << "priting Y23 \n " << Y23 << std::endl;

      BLAS::gemm('N','N',   k,   k,   k, RealType(1.0), Y23.data(),   k,  Y6.data(),   k, RealType(0.0), Y24.data(),   k);
      BLAS::gemm('N','N',   k,   k,   k, RealType(1.0),  Y6.data(),   k, Y24.data(),   k, RealType(0.0), Y25.data(),   k);


      Y26.resize(  k,nel);

      std::fill(Y26.begin(),Y26.end(),0.0);
      for ( int i=0; i<k; i++)
      {
//        BLAS::copy(k, Y25.data() + i, k, Y26.data() + *(data_it+1+i), nel);   
        BLAS::copy(k, Y25.data() + i, k, Y26.data() + (data_it[datum+1+i]), nel);   
      }


      Y7.resize(  k,nel);

      std::fill(Y7.begin(),Y7.end(),0.0);
      for ( int i=0; i<k; i++)
      {
//        BLAS::copy(k, Y6.data() + i, k, Y7.data() + *(data_it+1+i), nel);   
        BLAS::copy(k, Y6.data() + i, k, Y7.data() + (data_it[datum+1+i]), nel);   
      }

//      app_log() << "priting Y7 \n " << Y7 << std::endl;
      
      // c_Tr_AlphaI_MI is a constant contributing to constant const2
      // c_Tr_AlphaI_MI = Tr[\alpha_{I}^{-1}(P^{T}\widetilde{M} Q)]
      RealType c_Tr_AlphaI_MI = 0.0;
      for(int i=0; i<k;i++){
        c_Tr_AlphaI_MI+=Y24(i,i); }  

      for(int p=0; p<lookup_tbls[active_spin][index].size(); p++)
      {
        //el_p is the element position that contains information about the CI coefficient, and det up/dn values associated with the current unique determinant
        const int el_p (lookup_tbls[active_spin][index][p]);
        const RealType c=cptr[el_p];
        const size_t up=upC[el_p];
        const size_t down=dnC[el_p];

        const RealType alpha_1  (  c * detValues_dn[ down ] * detValues_up[ up ]/detValues_up[0] * c_Tr_AlphaI_MI  ); 
        const RealType alpha_2  (  c * detValues_dn[ down ] * detValues_up[ up ]/detValues_up[0]                   );
        const RealType alpha_3  (  c * Oi[ down ]           * detValues_up[ up ]/detValues_up[0]                   );
        const RealType alpha_4  (  c * detValues_dn[ down ] * detValues_up[ up ]                 * (1/psiCurrent)  );

        const2 += alpha_1;

        for ( int i=0; i<k; i++)
        {
//          BLAS::axpy(nel, alpha_1, Y7.data() + i*nel,1, pK1.data() + (*(data_it+1+k+i))*nel,1);   
//          BLAS::axpy(nel, alpha_2, Y7.data() + i*nel,1, pK2.data() + (*(data_it+1+k+i))*nel,1);   
//          BLAS::axpy(nel, alpha_3, Y7.data() + i*nel,1, pK3.data() + (*(data_it+1+k+i))*nel,1);   
//          BLAS::axpy(nel, alpha_4, Y7.data() + i*nel,1, pK4.data() + (*(data_it+1+k+i))*nel,1);   
//          BLAS::axpy(nel, alpha_2,Y26.data() + i*nel,1, pK5.data() + (*(data_it+1+k+i))*nel,1);   
          BLAS::axpy(nel, alpha_1, Y7.data() + i*nel,1, pK1.data() + (data_it[datum+1+k+i])*nel,1);   
          BLAS::axpy(nel, alpha_2, Y7.data() + i*nel,1, pK2.data() + (data_it[datum+1+k+i])*nel,1);   
          BLAS::axpy(nel, alpha_3, Y7.data() + i*nel,1, pK3.data() + (data_it[datum+1+k+i])*nel,1);   
          BLAS::axpy(nel, alpha_4, Y7.data() + i*nel,1, pK4.data() + (data_it[datum+1+k+i])*nel,1);   
          BLAS::axpy(nel, alpha_2,Y26.data() + i*nel,1, pK5.data() + (data_it[datum+1+k+i])*nel,1);   
        }
      }
//      data_it += 3*k+1;
      datum += 3*k+1;
    }

  }


  BLAS::gemm('N','N', nmo, nmo, nel, 1.0/const0               , T           , nmo, pK1.data(), nel, RealType(0.0), K1T.data()   ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K1T.data()  , nmo, T         , nmo, RealType(0.0), TK1T.data()  ,  nmo);

  BLAS::gemm('N','N', nmo, nmo, nel, 1.0/const0               , Y3.data()   , nmo, pK2.data(), nel, RealType(0.0), K2AiB.data() ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K2AiB.data(), nmo, T         , nmo, RealType(0.0), TK2AiB.data(),  nmo);
  BLAS::gemm('N','N', nmo, nmo, nel, 1.0/const0               , Y2.data()   , nmo, pK2.data(), nel, RealType(0.0), K2XA.data()  ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K2XA.data() , nmo, T         , nmo, RealType(0.0), TK2XA.data() ,  nmo);

  BLAS::gemm('N','N', nmo, nmo, nel, const1/(const0*const0)   , T           , nmo, pK2.data(), nel, RealType(0.0), K2T.data()   ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K2T.data()  , nmo, T         , nmo, RealType(0.0), TK2T.data()  ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, const0/const1            , K2T.data()  , nmo, Y4.data() , nmo, RealType(0.0), MK2T.data()  ,  nmo);
 
  BLAS::gemm('N','N', nmo, nmo, nel, 1.0/const0               , T           , nmo, pK3.data(), nel, RealType(0.0), K3T.data()   ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K3T.data()  , nmo, T         , nmo, RealType(0.0), TK3T.data()  ,  nmo);

  BLAS::gemm('N','N', nmo, nmo, nel, RealType(1.0)            , T           , nmo, pK4.data(), nel, RealType(0.0), K4T.data()   ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K4T.data()  , nmo, T         , nmo, RealType(0.0), TK4T.data()  ,  nmo);

  BLAS::gemm('N','N', nmo, nmo, nel, 1.0/const0               , T           , nmo, pK5.data(), nel, RealType(0.0), K5T.data()   ,  nmo);
  BLAS::gemm('N','N', nmo, nel, nmo, RealType(1.0)            , K5T.data()  , nmo, T         , nmo, RealType(0.0), TK5T.data()  ,  nmo);

 

  for (int mu=0, kk=parameter_start_index; kk < (parameter_start_index + parameters_size) ; kk++, mu++)
  {
    const int i((*m_act_rot_inds)[mu].first), j((*m_act_rot_inds)[mu].second);
    if (i<=nel-1 && j>nel-1)
      { dlogpsi[kk] += detValues_up[0]*( Tr(i,j) )*const0*(1/psiCurrent)\
                       + ( K4T(i,j)-K4T(j,i)-TK4T(i,j)) ;

        dhpsioverpsi[kk] += -0.5 * Y4(i,j)\
                            -0.5*( -   K5T(i,j) +   K5T(j,i)    +   TK5T(i,j)  \
                                   + K2AiB(i,j) - K2AiB(j,i)    - TK2AiB(i,j)  \
                                   -  K2XA(i,j) +  K2XA(j,i)    +  TK2XA(i,j)  \
                                                                -   MK2T(i,j)  \
                                   +   K1T(i,j) -   K1T(j,i)    -   TK1T(i,j)  \
                                   -  const2/const1 * K2T(i,j) +  const2/const1 * K2T(j,i)    +   const2/const1 * TK2T(i,j)  \
                                   +   K3T(i,j) -   K3T(j,i)    -   TK3T(i,j)  \
                                   -   K2T(i,j) +     K2T(j,i)    +   TK2T(i,j) \
                                 );
      }
    else if (i<=nel-1 && j<=nel-1)
      { dlogpsi[kk] += detValues_up[0]* ( Tr(i,j) - Tr(j,i) )*const0*(1/psiCurrent)\
                       + (K4T(i,j)-TK4T(i,j)-K4T(j,i)+TK4T(j,i) );

        dhpsioverpsi[kk] += -0.5 *( Y4(i,j)-Y4(j,i) )\
                            -0.5*( -   K5T(i,j) +   K5T(j,i)    +   TK5T(i,j) -    TK5T(j,i) \
                                   + K2AiB(i,j) - K2AiB(j,i)    - TK2AiB(i,j) +  TK2AiB(j,i) \
                                   -  K2XA(i,j) +  K2XA(j,i)    +  TK2XA(i,j) -   TK2XA(j,i) \
                                                                -   MK2T(i,j) +    MK2T(j,i) \
                                   +   K1T(i,j) -   K1T(j,i)    -   TK1T(i,j) +    TK1T(j,i) \
                                   -   const2/const1 * K2T(i,j) +   const2/const1 * K2T(j,i)    +   const2/const1 * TK2T(i,j) -    const2/const1 * TK2T(j,i) \
                                   +   K3T(i,j) -   K3T(j,i)    -   TK3T(i,j) +    TK3T(j,i) \
                                   -   K2T(i,j) +   K2T(j,i)    +   TK2T(i,j) -    TK2T(j,i) \
                                 );


      }
    else  
      { dlogpsi[kk] +=  ( K4T(i,j)-K4T(j,i) );

        dhpsioverpsi[kk] += \
                            -0.5*( -   K5T(i,j) +   K5T(j,i)    \
                                   + K2AiB(i,j) - K2AiB(j,i)    \
                                   -  K2XA(i,j) +  K2XA(j,i)    \
                                                                \
                                   +   K1T(i,j) -   K1T(j,i)    \
                                   -   const2/const1 * K2T(i,j) +   const2/const1 * K2T(j,i)    \
                                   +   K3T(i,j) -   K3T(j,i)    \
                                   -   K2T(i,j) +   K2T(j,i)  \
                                 );

      }
  }

}

}
