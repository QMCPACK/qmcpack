//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/CompositeSPOSet.h>
#include <Utilities/IteratorUtility.h>
#include <algorithm>

namespace qmcplusplus
{

  CompositeSPOSet::CompositeSPOSet()
  {
    className = "CompositeSPOSet";
    OrbitalSetSize = 0;
  }

  CompositeSPOSet::~CompositeSPOSet()
  {
    delete_iter(component_values.begin(),component_values.end());
    delete_iter(component_gradients.begin(),component_gradients.end());
    delete_iter(component_laplacians.begin(),component_laplacians.end());
  }


  void CompositeSPOSet::add(SPOSetBase* component)
  {
    int norbs = component->size();

    ValueVector_t* values     = new ValueVector_t;
    GradVector_t*  gradients  = new GradVector_t;
    ValueVector_t* laplacians = new ValueVector_t;

    values->resize(norbs);
    gradients->resize(norbs);
    laplacians->resize(norbs);

    components.push_back(component);
    component_values.push_back(values);
    component_gradients.push_back(gradients);
    component_laplacians.push_back(laplacians);

    OrbitalSetSize += norbs;
    BasisSetSize = OrbitalSetSize;
  }


  void CompositeSPOSet::report()
  {
    app_log()<<"CompositeSPOSet"<<endl;
    app_log()<<"  ncomponents = "<<components.size()<<endl;
    app_log()<<"  components"<<endl;
    for(int i=0;i<components.size();++i)
    {
      SPOSetBase& c = *components[i];
      app_log()<<"    "<<i<<endl;
      components[i]->basic_report("      ");
    }
  }

  void CompositeSPOSet::resetTargetParticleSet(ParticleSet& P)
  {
    for(int c=0;c<components.size();++c)
      components[c]->resetTargetParticleSet(P);
  }


  SPOSetBase* CompositeSPOSet::makeClone() const
  {
    CompositeSPOSet* clone = new CompositeSPOSet();
    for(int c=0;c<components.size();++c)
      clone->add(components[c]->makeClone());
    return clone;
  }

  
  void CompositeSPOSet::clone_from(CompositeSPOSet& other)
  {
    for(int c=0;c<other.components.size();++c)
      add(other.components[c]->makeClone());
  }


  void CompositeSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    int n=0;
    for(int c=0;c<components.size();++c)
    {
      SPOSetBase&    component = *components[c];
      ValueVector_t& values    = *component_values[c];
      component.evaluate(P,iat,values);
      copy(values.begin(),values.end(),psi.begin()+n);
      n += component.size();
    }
  }

  
  void CompositeSPOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, 
                GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    int n=0;
    for(int c=0;c<components.size();++c)
    {
      SPOSetBase&    component  = *components[c];
      ValueVector_t& values     = *component_values[c];
      GradVector_t&  gradients  = *component_gradients[c];
      ValueVector_t& laplacians = *component_laplacians[c];
      component.evaluate(P,iat,values,gradients,laplacians);
      copy(values.begin(),    values.end(),    psi.begin()+n  );
      copy(gradients.begin(), gradients.end(), dpsi.begin()+n );
      copy(laplacians.begin(),laplacians.end(),d2psi.begin()+n);
      n += component.size();
    }
  }


  //methods to be implemented later
  void CompositeSPOSet::resetParameters(const opt_variables_type& optVariables)
  {
    not_implemented("resetParameters");
  }
    
  void CompositeSPOSet::evaluate(
    const ParticleSet& P, PosType &r, ValueVector_t &psi)
  {
    not_implemented("evaluate(P,r,psi)");
  }

  void CompositeSPOSet::evaluate(
    const ParticleSet& P, int iat,ValueVector_t& psi, GradVector_t& dpsi, 
    HessVector_t& grad_grad_psi)
  {
    not_implemented("evaluate(P,iat,psi,dpsi,ddpsi)");
  }

  void CompositeSPOSet::evaluate(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    not_implemented("evaluate(P,first,last,logdet,dlogdet,d2logdet)");
  }

  void CompositeSPOSet::evaluate(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    not_implemented("evaluate(P,first,last,logdet,dlogdet,ddlogdet)");
  }

  void CompositeSPOSet::evaluate(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, 
    GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluate(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
  }

  void CompositeSPOSet::evaluateThirdDeriv(
    const ParticleSet& P,int first,int last,GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluateThirdDeriv(P,first,last,dddlogdet)");
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,d2logdet)");
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet)");
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, 
    GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
  }

  void CompositeSPOSet::evaluateGradSource(
    const ParticleSet &P, int first, int last, 
    const ParticleSet &source, int iat_src, GradMatrix_t &gradphi)
  {
    not_implemented("evaluateGradSource(P,first,last,source,iat,dphi)");
  }

  void CompositeSPOSet::evaluateGradSource(
    const ParticleSet &P, int first, int last, const ParticleSet &source, 
    int iat_src, GradMatrix_t &grad_phi, HessMatrix_t &grad_grad_phi, 
    GradMatrix_t &grad_lapl_phi)
  {
    not_implemented("evaluateGradSource(P,first,last,source,iat,dphi,ddphi,dd2phi)");
  }

  void CompositeSPOSet::evaluateBasis(
    const ParticleSet &P, int first, int last, ValueMatrix_t &basis_val,  
    GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl)
  {
    not_implemented("evaluateBasis");
  }

  void CompositeSPOSet::evaluateForDeriv(
    const ParticleSet &P, int first, int last, ValueMatrix_t &basis_val,  
    GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl)
  {
    not_implemented("evaluateForDeriv");
  }

  void CompositeSPOSet::copyParamsFromMatrix(
    const opt_variables_type& active, const ValueMatrix_t &mat, 
    vector<RealType> &destVec)
  {
    not_implemented("copyParamsFromMatrix");
  }


}
