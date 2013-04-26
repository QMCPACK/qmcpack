
/** constructor
 *@param spos the single-particle orbital set
 *@param first index of the first particle
 */
MultiDiracDeterminantBase::MultiDiracDeterminantBase(SPOSetBasePtr const &spos, int first):
  DiracDeterminantBase(spos,first)
{
  Optimizable=false;
  OrbitalName="MultiDiracDeterminantBase";
}

///default destructor
MultiDiracDeterminantBase::~MultiDiracDeterminantBase() {}

MultiDiracDeterminantBase& MultiDiracDeterminantBase::operator=(const MultiDiracDeterminantBase& s)
{
  NP=0;
  resize(s.NumPtcls, s.NumOrbitals);
  return *this;
}


MultiDiracDeterminantBase::RealType
MultiDiracDeterminantBase::registerData(ParticleSet& P, PooledData<RealType>& buf)
{
  if(NP == 0)
    //first time, allocate once
  {
    //int norb = cols();
    dpsiV.resize(NumOrbitals_total);
    d2psiV.resize(NumOrbitals_total);
    workV1.resize(NumOrbitals);
    workV2.resize(NumOrbitals);
    NP=P.getTotalNum();
    myG.resize(NP);
    myL.resize(NP);
    myG_temp.resize(NP);
    myL_temp.resize(NP);
    FirstAddressOfG = &myG[0][0];
    LastAddressOfG = FirstAddressOfG + NP*DIM;
    FirstAddressOfdV = &(dpsiM(0,0)[0]); //(*dpsiM.begin())[0]);
    LastAddressOfdV = FirstAddressOfdV + NumPtcls*NumOrbitals*DIM;
  }
  myG=0.0;
  myL=0.0;
  //ValueType x=evaluate(P,myG,myL);
  LogValue=evaluateLog(P,myG,myL);
  P.G += myG;
  P.L += myL;
  //add the data: determinant, inverse, gradient and laplacians
  buf.add(psiM.first_address(),psiM.last_address());
  buf.add(psiM_actual.first_address(),psiM_actual.last_address());
  buf.add(FirstAddressOfdV,LastAddressOfdV);
  buf.add(d2psiM.first_address(),d2psiM.last_address());
  buf.add(myL.first_address(), myL.last_address());
  buf.add(FirstAddressOfG,LastAddressOfG);
  buf.add(LogValue);
  buf.add(PhaseValue);
  return LogValue;
}


///reset the size: with the number of particles and number of orbtials
void MultiDiracDeterminantBase::resize(int nel, int morb)
{
  int norb=nel;
  if(norb <= 0)
    norb = nel; // for morb == -1 (default)
  psiM.resize(nel,norb);
  psiMinv.resize(nel,norb);
  WorkSpace.resize(nel);
  Pivot.resize(nel);
  LastIndex = FirstIndex + nel;
  NumPtcls=nel;
  NumOrbitals=norb;
  NumOrbitals_unoccupied=morb-norb;
  NumOrbitals_total=morb;
  psiV.resize(NumOrbitals_total);
  psiV_old.resize(NumOrbitals_total);
  dpsiM.resize(nel,NumOrbitals_total);
  d2psiM.resize(nel,NumOrbitals_total);
  psiM_temp.resize(nel,NumOrbitals);
  dpsiM_temp.resize(nel,NumOrbitals_total);
  d2psiM_temp.resize(nel,NumOrbitals_total);
  psiM_actual.resize(NumOrbitals_total,nel); //WARNING: Orbitals are first!
  Excitations.set_excitations(NumOrbitals,NumOrbitals_total-1,NumOrbitals-3,NumOrbitals-1,coefs);
  //    coef_vals.resize(coefs.size());
  //    coef_vals_temp.resize(coefs.size());
  // For forces
  grad_source_psiM.resize(nel,norb);
  grad_lapl_source_psiM.resize(nel,norb);
  grad_grad_source_psiM.resize(nel,norb);
  phi_alpha_Minv.resize(nel,norb);
  grad_phi_Minv.resize(nel,norb);
  lapl_phi_Minv.resize(nel,norb);
  grad_phi_alpha_Minv.resize(nel,norb);
}


/** dump the inverse to the buffer
*/
void MultiDiracDeterminantBase::dumpToBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  APP_ABORT("MultiDiracDeterminantBase::dumpToBuffer");
  buf.add(psiM.first_address(),psiM.last_address());
}

/** copy the inverse from the buffer
*/
void MultiDiracDeterminantBase::dumpFromBuffer(ParticleSet& P, PooledData<RealType>& buf)
{
  APP_ABORT("MultiDiracDeterminantBase::dumpFromBuffer");
  buf.get(psiM.first_address(),psiM.last_address());
}



MultiDiracDeterminantBase::GradType
MultiDiracDeterminantBase::evalGradSource(ParticleSet& P, ParticleSet& source,
    int iat)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex,
                           source, iat, grad_source_psiM);
//     Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
//     LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,WorkSpace.data(),Pivot.data(),PhaseValue);
  const ValueType* restrict yptr=psiM[0];
  const GradType* restrict dyptr=grad_source_psiM[0];
  GradType rv(0.0,0.0,0.0);
  for (int i=0; i<NumPtcls; i++)
    for(int j=0; j<NumOrbitals; j++)
      //rv += (*yptr++) *(*dyptr++);
      rv += grad_source_psiM(i,j) * psiM(i,j);
  // HACK HACK
  //return (grad_source_psiM(1,3));
  return rv;
}

MultiDiracDeterminantBase::GradType
MultiDiracDeterminantBase::evalGradSource
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat,
                           grad_source_psiM, grad_grad_source_psiM,
                           grad_lapl_source_psiM);
  Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
                         WorkSpace.data(),Pivot.data(),PhaseValue);
  GradMatrix_t &Phi_alpha(grad_source_psiM);
  GradMatrix_t &Grad_phi(dpsiM);
  ValueMatrix_t &Grad2_phi(d2psiM);
  HessMatrix_t &Grad_phi_alpha(grad_grad_source_psiM);
  GradMatrix_t &Grad2_phi_alpha(grad_lapl_source_psiM);
  //    vector<ValueType> grad_psi_over_psi_vector;
  //    vector<ValueType> grad_psi_alpha_over_psi_vector;
  GradType Psi_alpha_over_psi;
  Psi_alpha_over_psi = evalGradSource(P, source, iat);
  ofstream outfile;
  outfile.open("grad_psi_alpha_over_psi",ios::app);
  ValueMatrix_t toDet;
  ValueMatrix_t toDet_l;
  toDet.resize(2,2);
  toDet_l.resize(2,2);
  for (int ptcl=0; ptcl<NumPtcls; ptcl++)
  {
    ValueType Grad2_psi_over_psi(0.0);
    GradType Grad_psi_over_psi(0.0);
    HessType Grad_psi_alpha_over_psi(0.0);
    HessType one_row_change(0.0);
    HessType two_row_change(0.0);
    GradType one_row_change_l(0.0);
    GradType two_row_change_l(0.0);
    for (int el_dim=0; el_dim<3; el_dim++)
    {
      for (int orbital=0; orbital<NumOrbitals; orbital++)
      {
        Grad_psi_over_psi(el_dim)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl,orbital);
        if (el_dim==0)
          Grad2_psi_over_psi+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
      }
      for (int dim=0; dim<3; dim++)
      {
        one_row_change(dim,el_dim)=0.0;
        for (int orbital=0; orbital<NumOrbitals; orbital++)
        {
          one_row_change(dim,el_dim)+=Grad_phi_alpha(ptcl,orbital)(dim,el_dim)*psiM(ptcl,orbital);
          if (el_dim==0)
            one_row_change_l(dim)+=Grad2_phi_alpha(ptcl,orbital)(dim)*psiM(ptcl,orbital);
        }
        for (int ptcl2=0; ptcl2<NumPtcls; ptcl2++)
        {
          if (ptcl!=ptcl2)
          {
            toDet=0.0;
            toDet_l=0.0;
            for (int orbital=0; orbital<NumOrbitals; orbital++)
            {
              toDet(0,0)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl,orbital);
              toDet_l(0,0)+=Grad2_phi(ptcl,orbital)*psiM(ptcl,orbital);
              toDet(0,1)+=Grad_phi(ptcl,orbital)(el_dim)*psiM(ptcl2,orbital);
              toDet_l(0,1)+=Grad2_phi(ptcl,orbital)*psiM(ptcl2,orbital);
              toDet(1,0)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl,orbital);
              toDet_l(1,0)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl,orbital);
              toDet(1,1)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl2,orbital);
              toDet_l(1,1)+=Phi_alpha(ptcl2,orbital)(dim)*psiM(ptcl2,orbital);
            }
            two_row_change(dim,el_dim)+=toDet(0,0)*toDet(1,1)-toDet(1,0)*toDet(0,1);
            if (el_dim==0)
              two_row_change_l(dim)+=toDet_l(0,0)*toDet_l(1,1)-toDet_l(1,0)*toDet_l(0,1);
          }
        }
        Grad_psi_alpha_over_psi(dim,el_dim)=one_row_change(dim,el_dim)+two_row_change(dim,el_dim);
        outfile<<Grad_psi_alpha_over_psi(dim,el_dim)<<endl;
        grad_grad(dim)(ptcl)(el_dim)=one_row_change(dim,el_dim)+two_row_change(dim,el_dim)-
                                     Grad_psi_over_psi(el_dim)*Psi_alpha_over_psi(dim);
      }
    }
    for (int dim=0; dim<3; dim++)
    {
      lapl_grad(dim)(ptcl)=0.0;
      lapl_grad(dim)(ptcl)+=one_row_change_l(dim)+two_row_change_l(dim)- Psi_alpha_over_psi(dim)*Grad2_psi_over_psi;
      for (int el_dim=0; el_dim<3; el_dim++)
      {
        lapl_grad(dim)(ptcl)-= 2.0*Grad_psi_alpha_over_psi(dim,el_dim)*Grad_psi_over_psi(el_dim);
        lapl_grad(dim)(ptcl)+= 2.0*Psi_alpha_over_psi(dim)*(Grad_psi_over_psi(el_dim)*Grad_psi_over_psi(el_dim));
      }
    }
  }
  outfile.close();
  return Psi_alpha_over_psi;
}


MultiDiracDeterminantBase::GradType
MultiDiracDeterminantBase::evalGradSourcep
(ParticleSet& P, ParticleSet& source,int iat,
 TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> &grad_grad,
 TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> &lapl_grad)
{
  Phi->evaluateGradSource (P, FirstIndex, LastIndex, source, iat,
                           grad_source_psiM, grad_grad_source_psiM,
                           grad_lapl_source_psiM);
  // HACK HACK HACK
  // Phi->evaluate(P, FirstIndex, LastIndex, psiM, dpsiM, d2psiM);
  // LogValue=InvertWithLog(psiM.data(),NumPtcls,NumOrbitals,
  // 			   WorkSpace.data(),Pivot.data(),PhaseValue);
  // Compute matrices
  phi_alpha_Minv = 0.0;
  grad_phi_Minv = 0.0;
  lapl_phi_Minv = 0.0;
  grad_phi_alpha_Minv = 0.0;
  for (int i=0; i<NumPtcls; i++)
    for (int j=0; j<NumOrbitals; j++)
    {
      lapl_phi_Minv(i,j) = 0.0;
      for (int k=0; k<NumOrbitals; k++)
        lapl_phi_Minv(i,j) += d2psiM(i,k)*psiM(j,k);
    }
  for (int dim=0; dim<3; dim++)
  {
    for (int i=0; i<NumPtcls; i++)
      for (int j=0; j<NumOrbitals; j++)
      {
        for (int k=0; k<NumOrbitals; k++)
        {
          phi_alpha_Minv(i,j)[dim] += grad_source_psiM(i,k)[dim] * psiM(j,k);
          grad_phi_Minv(i,j)[dim] += dpsiM(i,k)[dim] * psiM(j,k);
          for (int dim_el=0; dim_el<3; dim_el++)
            grad_phi_alpha_Minv(i,j)(dim, dim_el) +=
              grad_grad_source_psiM(i,k)(dim,dim_el)*psiM(j,k);
        }
      }
  }
  GradType gradPsi;
  for(int i=0, iel=FirstIndex; i<NumPtcls; i++, iel++)
  {
    HessType dval (0.0);
    GradType d2val(0.0);
    for (int dim=0; dim<3; dim++)
      for (int dim_el=0; dim_el<3; dim_el++)
        dval(dim,dim_el) = grad_phi_alpha_Minv(i,i)(dim,dim_el);
    for(int j=0; j<NumOrbitals; j++)
    {
      gradPsi += grad_source_psiM(i,j) * psiM(i,j);
      for (int dim=0; dim<3; dim++)
        for (int k=0; k<3; k++)
          dval(dim,k) -= phi_alpha_Minv(j,i)[dim]*grad_phi_Minv(i,j)[k];
    }
    for (int dim=0; dim<OHMMS_DIM; dim++)
    {
      for (int k=0; k<3; k++)
        grad_grad[dim][iel][k] += dval(dim,k);
      for (int j=0; j<NumOrbitals; j++)
      {
        // First term, eq 9
        lapl_grad[dim][iel] += grad_lapl_source_psiM(i,j)[dim] *
                               psiM(i,j);
        // Second term, eq 9
        if (j == i)
          for (int dim_el=0; dim_el<3; dim_el++)
            lapl_grad[dim][iel] -= 2.0 * grad_phi_alpha_Minv(j,i)(dim,dim_el)
                                   * grad_phi_Minv(i,j)[dim_el];
        // Third term, eq 9
        // First term, eq 10
        lapl_grad[dim][iel] -= phi_alpha_Minv(j,i)[dim]*lapl_phi_Minv(i,j);
        // Second term, eq 11
        for (int dim_el=0; dim_el<3; dim_el++)
          lapl_grad[dim][iel] += 2.0*phi_alpha_Minv(j,i)[dim] *
                                 grad_phi_Minv(i,i)[dim_el]*grad_phi_Minv(i,j)[dim_el];
      }
    }
  }
  return gradPsi;
}



MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::logRatio(ParticleSet& P, int iat,
    ParticleSet::ParticleGradient_t& dG,
    ParticleSet::ParticleLaplacian_t& dL)
{
  APP_ABORT("  logRatio is not allowed");
  //THIS SHOULD NOT BE CALLED
  ValueType r=ratio(P,iat,dG,dL);
  return LogValue = evaluateLogAndPhase(r,PhaseValue);
}

MultiDiracDeterminantBase::RealType
MultiDiracDeterminantBase::evaluateLog(ParticleSet& P, PooledData<RealType>& buf)
{
  buf.put(psiM.first_address(),psiM.last_address());
  buf.put(psiM_actual.first_address(),psiM_actual.last_address());
  buf.put(FirstAddressOfdV,LastAddressOfdV);
  buf.put(d2psiM.first_address(),d2psiM.last_address());
  buf.put(myL.first_address(), myL.last_address());
  buf.put(FirstAddressOfG,LastAddressOfG);
  buf.put(LogValue);
  buf.put(PhaseValue);
  return LogValue;
}


MultiDiracDeterminantBase::MultiDiracDeterminantBase(const MultiDiracDeterminantBase& s):
  DiracDeterminantBase(s)
{
  this->resize(s.NumPtcls,s.NumOrbitals);
}

SPOSetBasePtr  MultiDiracDeterminantBase::clonePhi() const
{
  return Phi->makeClone();
}

/** return the ratio only for the  iat-th partcle move
 * @param P current configuration
 * @param iat the particle thas is being moved
 */
MultiDiracDeterminantBase::ValueType MultiDiracDeterminantBase::ratio_multi(ParticleSet& P, int iat)
{
  APP_ABORT("MultiDiracDeterminantBase::old_ratio");
  UpdateMode=ORB_PBYP_RATIO;
  WorkingIndex = iat-FirstIndex;
  Phi->evaluate(P, iat, psiV);
#ifdef DIRAC_USE_BLAS
  return curRatio = BLAS::dot(NumOrbitals,psiM[iat-FirstIndex],&psiV[0]);
#else
  return curRatio = DetRatio(psiM, psiV.begin(),iat-FirstIndex);
#endif
}


void MultiDiracDeterminantBase::update(ParticleSet& P,
                                       ParticleSet::ParticleGradient_t& dG,
                                       ParticleSet::ParticleLaplacian_t& dL,
                                       int iat)
{
  APP_ABORT("MultiDiracDeterminantBase::update..I didnt think this was ever called..if im wrong tell me and Ill fix it..Bryan");
  DetUpdate(psiM,psiV,workV1,workV2,WorkingIndex,curRatio);
  for(int j=0; j<NumOrbitals; j++)
  {
    dpsiM(WorkingIndex,j)=dpsiV[j];
    d2psiM(WorkingIndex,j)=d2psiV[j];
  }
  int kat=FirstIndex;
  for(int i=0; i<NumPtcls; i++,kat++)
  {
    GradType rv=dot(psiM[i],dpsiM[i],NumOrbitals);
    ValueType lap=dot(psiM[i],d2psiM[i],NumOrbitals);
    //GradType rv =psiM(i,0)*dpsiM(i,0);
    //ValueType lap=psiM(i,0)*d2psiM(i,0);
    //for(int j=1; j<NumOrbitals; j++) {
    //  rv += psiM(i,j)*dpsiM(i,j);
    //  lap += psiM(i,j)*d2psiM(i,j);
    //}
    lap -= dot(rv,rv);
    dG[kat] += rv - myG[kat];
    myG[kat]=rv;
    dL[kat] += lap -myL[kat];
    myL[kat]=lap;
  }
  PhaseValue += evaluatePhase(curRatio);
  LogValue +=std::log(std::abs(curRatio));
  curRatio=1.0;
}


/** Calculate the value of the Dirac determinant for particles
 *@param P input configuration containing N particles
 *@param G a vector containing N gradients
 *@param L a vector containing N laplacians
 *@return the value of the determinant
 *
 *\f$ (first,first+nel). \f$  Add the gradient and laplacian
 *contribution of the determinant to G(radient) and L(aplacian)
 *for local energy calculations.
 */
MultiDiracDeterminantBase::ValueType
MultiDiracDeterminantBase::evaluate(ParticleSet& P,
                                    ParticleSet::ParticleGradient_t& G,
                                    ParticleSet::ParticleLaplacian_t& L)
{
  APP_ABORT("  MultiDiracDeterminantBase::evaluate is distabled");
  Phi->evaluate(P, FirstIndex, LastIndex, psiM,dpsiM, d2psiM);
  ValueType CurrentDet;
  if(NumPtcls==1)
  {
    CurrentDet=psiM(0,0);
    ValueType y=1.0/CurrentDet;
    psiM(0,0)=y;
    GradType rv = y*dpsiM(0,0);
    G(FirstIndex) += rv;
    L(FirstIndex) += y*d2psiM(0,0) - dot(rv,rv);
  }
  else
  {
    CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals, WorkSpace.data(), Pivot.data());
    //CurrentDet = Invert(psiM.data(),NumPtcls,NumOrbitals);
    const ValueType* restrict yptr=psiM.data();
    const ValueType* restrict d2yptr=d2psiM.data();
    const GradType* restrict dyptr=dpsiM.data();
    for(int i=0, iat=FirstIndex; i<NumPtcls; i++, iat++)
    {
      GradType rv;
      ValueType lap=0.0;
      for(int j=0; j<NumOrbitals; j++,yptr++)
      {
        rv += *yptr * *dyptr++;
        lap += *yptr * *d2yptr++;
      }
      //Old index
      //    GradType rv = psiM(i,0)*dpsiM(i,0);
      //    ValueType lap=psiM(i,0)*d2psiM(i,0);
      //    for(int j=1; j<NumOrbitals; j++) {
      //      rv += psiM(i,j)*dpsiM(i,j);
      //      lap += psiM(i,j)*d2psiM(i,j);
      //    }
      G(iat) += rv;
      L(iat) += lap - dot(rv,rv);
    }
  }
  return CurrentDet;
}

