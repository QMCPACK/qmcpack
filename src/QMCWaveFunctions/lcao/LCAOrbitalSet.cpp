//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/lcao/LCAOrbitalSet.h"
#include <Numerics/MatrixOperators.h>
#include <boost/format.hpp>
#include <formic/utils/lapack_interface.h>

namespace qmcplusplus
{
  LCAOrbitalSet::LCAOrbitalSet(basis_type* bs,int rl):
    myBasisSet(nullptr), C(nullptr), ReportLevel(rl),
    BasisSetSize(0), Identity(true), IsCloned(false)
  {
    if(bs != nullptr) setBasisSet(bs);
  }

  LCAOrbitalSet::~LCAOrbitalSet()
  {
    if(!IsCloned && C!= nullptr) delete C;
  }

  void LCAOrbitalSet::setBasisSet(basis_type* bs)
  {
    myBasisSet=bs;
    BasisSetSize=myBasisSet->getBasisSetSize();
    Temp.resize(BasisSetSize);
  }

  bool LCAOrbitalSet::setIdentity(bool useIdentity)
  {
    Identity = useIdentity;
    if(Identity) return true;

    if ( C== nullptr && (OrbitalSetSize > 0) && (BasisSetSize > 0) )
    {
      C = new ValueMatrix_t(OrbitalSetSize, BasisSetSize);
    }
    else
    {
      app_error() << "either OrbitalSetSize or BasisSetSize has an invalid value !!\n";
      app_error() << "OrbitalSetSize = " << OrbitalSetSize << std::endl;
      app_error() << "BasisSetSize = " << BasisSetSize << std::endl;
      APP_ABORT("LCAOrbitalBuilder::setIdentiy ");
    }

    return true;
  }

  SPOSet* LCAOrbitalSet::makeClone() const
  {
    LCAOrbitalSet* myclone = new LCAOrbitalSet(*this);
    myclone->myBasisSet = myBasisSet->makeClone();
    myclone->IsCloned=true;
    return myclone;
  }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, 
      int iat, ValueVector_t& psi)
  {
    if(Identity)
    { //PAY ATTENTION TO COMPLEX
      myBasisSet->evaluateV(P,iat,psi.data());
    }
    else
    {
      Vector<ValueType> vTemp(Temp.data(0),BasisSetSize);
      myBasisSet->evaluateV(P,iat,vTemp.data());
      simd::gemv(*C,Temp.data(0),psi.data());
    }
  }

  /** Find a better place for other user classes, Matrix should be padded as well */
  template<typename T,unsigned D>
    inline void Product_ABt(const VectorSoaContainer<T,D>& A, const Matrix<T>& B, VectorSoaContainer<T,D>& C)
    {
      constexpr char transa = 't';
      constexpr char transb = 'n';
      constexpr T zone(1);
      constexpr T zero(0);
      BLAS::gemm(transa, transb, B.rows(), D, B.cols(),
          zone, B.data(), B.cols(), A.data(), A.capacity(),
          zero, C.data(), C.capacity());
    }

  inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp,
      ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi) const
  {
    simd::copy_n(temp.data(0),OrbitalSetSize,psi.data());
    const ValueType* restrict gx=temp.data(1);
    const ValueType* restrict gy=temp.data(2);
    const ValueType* restrict gz=temp.data(3);
    for(size_t j=0; j<OrbitalSetSize; j++)
    {
      dpsi[j][0]=gx[j];
      dpsi[j][1]=gy[j];
      dpsi[j][2]=gz[j];
    }
    simd::copy_n(temp.data(4),OrbitalSetSize,d2psi.data());
  }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, int iat,
      ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      //TAKE CARE OF IDENTITY
      myBasisSet->evaluateVGL(P,iat,Temp);
      if(Identity) 
        evaluate_vgl_impl(Temp,psi,dpsi,d2psi);
      else
      {
        Product_ABt(Temp,*C,Tempv);
        evaluate_vgl_impl(Tempv,psi,dpsi,d2psi);
      }
    }

  void LCAOrbitalSet::evaluateValues(const VirtualParticleSet& VP, ValueMatrix_t& psiM, ValueAlignedVector_t& SPOMem)
  {
    const int nVP = VP.getTotalNum();
    Matrix<ValueType> basisM(SPOMem.data(), nVP, BasisSetSize);
    for(size_t j=0; j<nVP; j++)
    {
      Vector<RealType> vTemp(basisM[j],BasisSetSize);
      myBasisSet->evaluateV(VP,j,vTemp.data());
    }
    MatrixOperators::product_ABt(basisM,*C,psiM);
  }

  size_t LCAOrbitalSet::estimateMemory(const int nP) { return BasisSetSize*nP; }

  void LCAOrbitalSet::evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi,
        HessVector_t& grad_grad_psi)
    {
#if 0
        myBasisSet->evaluateForPtclMoveWithHessian(P,iat);
        simd::gemv(C,myBasisSet->Phi.data(),psi.data());
        simd::gemv(C,myBasisSet->dPhi.data(),dpsi.data());
        simd::gemv(C,myBasisSet->grad_grad_Phi.data(),grad_grad_psi.data());
#endif
      }

  /* implement using gemm algorithm */
  inline void LCAOrbitalSet::evaluate_vgl_impl(const vgl_type& temp, int i,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet) const
  {
    simd::copy_n(temp.data(0),OrbitalSetSize,logdet[i]);
    const ValueType* restrict gx=temp.data(1);
    const ValueType* restrict gy=temp.data(2);
    const ValueType* restrict gz=temp.data(3);
    for(size_t j=0; j<OrbitalSetSize; j++)
    {
      dlogdet[i][j][0]=gx[j];
      dlogdet[i][j][1]=gy[j];
      dlogdet[i][j][2]=gz[j];
    }
    simd::copy_n(temp.data(4),OrbitalSetSize,d2logdet[i]);
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    if(Identity)
    {
      for(size_t i=0, iat=first; iat<last; i++,iat++)
      {
        myBasisSet->evaluateVGL(P,iat,Temp);
        evaluate_vgl_impl(Temp,i,logdet,dlogdet,d2logdet);
      }
    }
    else
    {
      for(size_t i=0, iat=first; iat<last; i++,iat++)
      {
        myBasisSet->evaluateVGL(P,iat,Temp);
        Product_ABt(Temp,*C,Tempv);
        evaluate_vgl_impl(Tempv,i,logdet,dlogdet,d2logdet);
      }
    }
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last,
      ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithHessian(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        for(int b=0; b<BasisSetSize; ++b,++jk)
        {
          dres +=  cptr[jk]*dptr[b];
          d2res +=  cptr[jk]*d2ptr[b];
        }
        dlogdet(ij)=dres;
        grad_grad_logdet(ij)=d2res;
        ++ij;
      }
    }
#endif
  }

  void LCAOrbitalSet::evaluate_notranspose(const ParticleSet& P, int first, int last
      , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, GGGMatrix_t& grad_grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateWithThirdDeriv(P,iat);
      MatrixOperators::product(C,myBasisSet->Phi,logdet[i]);
      const typename BS::GradType* restrict dptr=myBasisSet->dPhi.data();
      const typename BS::HessType* restrict d2ptr=myBasisSet->grad_grad_Phi.data();
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GradType dres;
        register HessType d2res;
        register GGGType gggres;
        for(int b=0; b<BasisSetSize; ++b)
        {
          dres +=  cptr[jk]*dptr[b];
          d2res +=  cptr[jk]*d2ptr[b];
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk++]*(gggptr[b])[2];
        }
        dlogdet(ij)=dres;
        grad_grad_logdet(ij)=d2res;
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
#endif
  }

  void LCAOrbitalSet::evaluateThirdDeriv(const ParticleSet& P, int first, int last
      , GGGMatrix_t& grad_grad_grad_logdet)
  {
#if 0
    const ValueType* restrict cptr=C.data();
#pragma ivdep
    for(int i=0,ij=0, iat=first; iat<last; i++,iat++)
    {
      myBasisSet->evaluateThirdDerivOnly(P,iat);
      const typename BS::GGGType* restrict gggptr=myBasisSet->grad_grad_grad_Phi.data();
      for(int j=0,jk=0; j<OrbitalSetSize; j++)
      {
        register GGGType gggres;
        for(int b=0; b<BasisSetSize; ++b)
        {
          gggres[0] +=  cptr[jk]*(gggptr[b])[0];
          gggres[1] +=  cptr[jk]*(gggptr[b])[1];
          gggres[2] +=  cptr[jk++]*(gggptr[b])[2];
        }
        grad_grad_grad_logdet(ij)=gggres;
        ++ij;
      }
    }
#endif
  }

  void LCAOrbitalSet::buildOptVariables(std::vector<RealType>& input_params, bool params_supplied, std::vector<int> * data, const size_t& nel, std::vector<size_t>& C2node, const int& spin)
  {
    const size_t& nmo = OrbitalSetSize;
    const size_t& nb = BasisSetSize;
    m_init_B = C;
    //a vector in which the element's index value correspond to Molecular Orbitals.
    //The element value at an index indicates how many times an electron is excited from or to that orbital in the Multi-Slater expansion i.e the indices with non-zero elements are active space orbitals
    std::vector<int> occupancy_vector (nmo,0);

    // Function to fill occupancy_vectors and also return number of unique determinants 
    const int unique_dets =  this->build_occ_vec(data, nel, nmo, &occupancy_vector);

    // When calculating the parameter derivative of the Multi-Slater component of the wavefunction, each unique deterimant can contribute multiple times. 
    // The lookup_tbls are used so that a parameter derivative of a unique determinant is only done once and then scaled according to how many times it appears in the Multi-Slater expansion 
    lookup_tbl.resize(unique_dets);
    //construct lookup table 
    for (int i(0); i < C2node.size(); i++)
    {
      lookup_tbl[C2node[i]].push_back(i);
    } 

    // create active rotations
    m_act_rot_inds.clear(); 


  for(int i=0;i<nmo;i++)
    for(int j=i+1;j<nmo;j++)
     { 
      bool core_i(!occupancy_vector[i] and i <= nel-1); // true if orbital i is a 'core' orbital
      bool core_j(!occupancy_vector[j] and j <= nel-1); // true if orbital j is a 'core' orbital
      bool virt_i(!occupancy_vector[i] and i >  nel-1); // true if orbital i is a 'virtual' orbital
      bool virt_j(!occupancy_vector[j] and j >  nel-1); // true if orbital j is a 'virtual' orbital
        if( !( 
              ( core_i and core_j  ) 
                        or
              ( virt_i and virt_j ) 
             )    
          )
        {
          m_act_rot_inds.push_back(std::pair<int,int>(i,j)); // orbital rotation parameter accepted as long as rotation isn't core-core or virtual-virtual
        }
     }
   
      // This will add the orbital rotation parameters to myVars 
      // and will also read in initial parameter values supplied in input file     
      int p, q; 
      int nparams_active = m_act_rot_inds.size();

      for (int i=0; i< nparams_active; i++)
      {
        p = m_act_rot_inds[i].first;
        q = m_act_rot_inds[i].second;
        std::stringstream sstr;
        sstr << "SPO_" << (spin==0 ? "UP":"DN") 
                << "_orb_rot_"
                << ( p <   10 ? "0": "")
                << ( p <  100 ? "0": "")
                << ( p < 1000 ? "0": "")
                << p
                << "_"
                << ( q <   10 ? "0": "")
                << ( q <  100 ? "0": "")
                << ( q < 1000 ? "0": "")
                << q;

        // If the user input parameteres, use those. Otherwise, initialize the parameters to zero
        if (params_supplied) {
          myVars.insert(sstr.str(), input_params[i]);
        } else {
          myVars.insert(sstr.str(), 0.0);
        }
      }

      //Printing the parameters
      if(true){
        app_log() << std::string(16,' ') << "Parameter name" << std::string(15,' ') << "Value\n";
        myVars.print(app_log());
      }
   
      //the below code is very similar to  "Reset parameters function"
      // reading out the parameters that define the rotation into an antisymmetric matrix
      std::vector<RealType> rot_mat(nmo*nmo, 0.0);
      for (int i = 0; i < m_act_rot_inds.size(); i++) 
      {
        const int p = m_act_rot_inds[i].first;
        const int q = m_act_rot_inds[i].second;

        rot_mat[p+q*nmo] =  myVars[i]; 
        rot_mat[q+p*nmo] = -myVars[i]; 

      }
      //exponentiate matrices and do BLAS command to perform rotation on m_b

      // exponentiate antisymmetric matrix to get the unitary rotation
      this->exponentiate_antisym_matrix(nmo, &rot_mat[0]);

      BLAS::gemm('N','T', nb, nmo, nmo, RealType(1.0), m_init_B->data(), nb, &rot_mat[0], nmo, RealType(0.0), C->data(), nb);
      if(Optimizable)
      {
        T.resize(nel,nmo);
        app_log()<< "VARIABLES HAVE NOT BE DELETED SDP\n";
      }
      else
      {
        //THIS ALLOWS FOR ORBITAL PARAMETERS TO BE READ IN EVEN WHEN THOSE PARAMETERS ARE NOT BEING OPTIMIZED
        //this assumes there are only CI coefficients ahead of the M_orb_coefficients 
        myVars.Index.erase(myVars.Index.begin(),myVars.Index.end());
        myVars.NameAndValue.erase(myVars.NameAndValue.begin(),myVars.NameAndValue.end());
        myVars.ParameterType.erase(myVars.ParameterType.begin(),myVars.ParameterType.end());
        myVars.Recompute.erase(myVars.Recompute.begin(),myVars.Recompute.end());
      }
  }

  int LCAOrbitalSet::build_occ_vec(std::vector<int> * data,
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


  // compute exponential of a real, antisymmetric matrix by diagonalizing and exponentiating eigenvalues
  void LCAOrbitalSet::exponentiate_antisym_matrix(const int n, RealType * const mat)
  {
    std::vector<std::complex<double> > mat_h(n*n, 0);
    std::vector<double> eval(n, 0);
    std::vector<std::complex<double> > work(2*n, 0);
    std::vector<double> rwork(3*n, 0);
    std::vector<std::complex<double> > mat_d(n*n, 0);
    std::vector<std::complex<double> > mat_t(n*n, 0);
    int info = 0;
    // exponentiating e^X = e^iY (Y hermitian)
    // i(-iX) = X, so -iX is hermitian
    // diagonalize -iX = UDU^T, exponentiate e^iD, and return U e^iD U^T
    // construct hermitian analogue of mat by multiplying by -i
    for(int i = 0; i < n; ++i)
    {
      for(int j = i; j < n; ++j)
      {
        mat_h[i+n*j] = std::complex<double>(0,-1.0*mat[i+n*j]);
        mat_h[j+n*i] = std::complex<double>(0,1.0*mat[i+n*j]);
      }
    }
    // diagonalize the matrix 
    formic::zheev('V', 'U', n, &mat_h.at(0), n, &eval.at(0), &work.at(0), 2*n, &rwork.at(0), info);
    if(info != 0)
    {
      std::ostringstream msg;
      msg << "zheev failed with info = " << info << " in MultiSlaterDeterminantFast::exponentiate_antisym_matrix";
      app_log() << msg.str() << std::endl;
      APP_ABORT(msg.str());
    }
    // iterate through diagonal matrix, exponentiate terms
    for(int i = 0; i < n; ++i)
    {
      for(int j = 0; j < n; ++j)
      {
        mat_d[i + j*n] = (i==j) ? std::exp( std::complex<double>(0.0,eval[i])) : std::complex<double>(0.0,0.0);
      }
    }
    // perform matrix multiplication 
    // assume row major
    formic::xgemm('N','C', n, n, n, std::complex<double>(1.0,0), &mat_d.at(0), n, &mat_h.at(0), n, std::complex<double>(0.0, 0.0), &mat_t.at(0), n);
    formic::xgemm('N','N', n, n, n, std::complex<double>(1.0,0), &mat_h.at(0), n, &mat_t.at(0), n, std::complex<double>(0.0, 0.0), &mat_d.at(0), n);
    for(int i = 0; i < n; ++i)
      for(int j = 0; j < n; ++j)
      {
        if(mat_d[i+n*j].imag() > 1e-12)
        {
          app_log() << "warning: large imaginary value in orbital rotation matrix: (i,j) = (" << i << "," << j << "), im = " << mat_d[i+n*j].imag() << std::endl;
        }
        mat[i+n*j] = mat_d[i+n*j].real();
      }
  }


}
