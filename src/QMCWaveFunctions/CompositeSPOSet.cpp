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
    // base class and shallow copy
    CompositeSPOSet* clone = new CompositeSPOSet(*this);
    // remove component shallow copies then deep copy
    clone->clone_from(*this);
    return clone;
  }

  
  void CompositeSPOSet::clone_from(const CompositeSPOSet& master)
  {
    components.clear();
    component_values.clear();
    component_gradients.clear();
    component_laplacians.clear();
    OrbitalSetSize = 0;
    for(int c=0;c<master.components.size();++c)
      add(master.components[c]->makeClone());
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

  void CompositeSPOSet::evaluate(
      const ParticleSet& P, int iat,ValueVector_t& psi, GradVector_t& dpsi, 
      HessVector_t& grad_grad_psi)
  {
    not_implemented("evaluate(P,iat,psi,dpsi,ddpsi)");
  }

  void CompositeSPOSet::evaluate( const ParticleSet& P, PosType &r, ValueVector_t &psi)
  {
    not_implemented("evaluate(P,r,psi)");
  }




  //methods to be implemented later
  void CompositeSPOSet::resetParameters(const opt_variables_type& optVariables)
  {
    for(int c=0;c<components.size();++c)
      components[c]->resetParameters(optVariables);
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    int nat=last-first;
    int n=0;
    for(int c=0;c<components.size();++c)
    {
      int norb=components[c]->size();
      ValueMatrix_t v(nat,norb);
      GradMatrix_t g(nat,norb);
      ValueMatrix_t l(nat,norb);
      components[c]->evaluate_notranspose(P,first,last,v,g,l);
      for(int iat=0; iat<nat; ++iat)
        copy(v[iat],v[iat]+norb,logdet[iat]+n);
      for(int iat=0; iat<nat; ++iat)
        copy(g[iat],g[iat]+norb,dlogdet[iat]+n);
      for(int iat=0; iat<nat; ++iat)
        copy(l[iat],l[iat]+norb,d2logdet[iat]+n);
      n += norb;
    }
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    int nat=last-first;
    int n=0;
    for(int c=0;c<components.size();++c)
    {
      int norb=components[c]->size();
      ValueMatrix_t v(nat,norb);
      GradMatrix_t g(nat,norb);
      HessMatrix_t h(nat,norb);
      components[c]->evaluate_notranspose(P,first,last,v,g,h);
      for(int iat=0; iat<nat; ++iat)
        copy(v[iat],v[iat]+norb,logdet[iat]+n);
      for(int iat=0; iat<nat; ++iat)
        copy(g[iat],g[iat]+norb,dlogdet[iat]+n);
      for(int iat=0; iat<nat; ++iat)
        copy(h[iat],h[iat]+norb,grad_grad_logdet[iat]+n);
      n += norb;
    }
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, 
    GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
  }

}
