//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#include <QMCWaveFunctions/CompositeSPOSet.h>
#include <Utilities/IteratorUtility.h>
#include <algorithm>
#include <OhmmsData/AttributeSet.h>
#include <QMCWaveFunctions/BasisSetFactory.h>

namespace qmcplusplus
{

  namespace MatrixOperators
  {
    /** copy a small matrix (N, M1) to a big matrix (N, M2), M2>M1
     * @param small input matrix
     * @param big outout matrix
     * @param offset_c column offset
     *
     * @todo smater and more efficient matrix, move up for others
     * The columns [0,M1) are inserted into [offset_c,offset_c+M1).
     */
    template<typename MAT1, typename MAT2>
      inline void insert_columns(const MAT1& small, MAT2& big, int offset_c)
      {
        const int c=small.cols();
        for(int i=0; i<small.rows(); ++i) 
          copy(small[i],small[i]+c,big[i]+offset_c);
      }
  }

  CompositeSPOSet::CompositeSPOSet()
  {
    className = "CompositeSPOSet";
    OrbitalSetSize = 0;
    component_offsets.reserve(4);
  }

  CompositeSPOSet::~CompositeSPOSet()
  {
    delete_iter(component_values.begin(),component_values.end());
    delete_iter(component_gradients.begin(),component_gradients.end());
    delete_iter(component_laplacians.begin(),component_laplacians.end());
  }


  void CompositeSPOSet::add(SPOSetBase* component)
  {
    if(components.empty())
    {
      component_offsets.push_back(0); //add 0
    }

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

    component_offsets.push_back(OrbitalSetSize);
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
    component_offsets.clear(); //add 0

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
    const int nat=last-first;
    for(int c=0;c<components.size();++c)
    {
      int norb=components[c]->size();
      ValueMatrix_t v(nat,norb);
      GradMatrix_t g(nat,norb);
      ValueMatrix_t l(nat,norb);
      components[c]->evaluate_notranspose(P,first,last,v,g,l);
      int n=component_offsets[c];
      MatrixOperators::insert_columns(v,logdet,n);
      MatrixOperators::insert_columns(g,dlogdet,n);
      MatrixOperators::insert_columns(l,d2logdet,n);
    }
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    const int nat=last-first;
    for(int c=0;c<components.size();++c)
    {
      int norb=components[c]->size();
      ValueMatrix_t v(nat,norb);
      GradMatrix_t g(nat,norb);
      HessMatrix_t h(nat,norb);
      components[c]->evaluate_notranspose(P,first,last,v,g,h);
      int n=component_offsets[c];
      MatrixOperators::insert_columns(v,logdet,n);
      MatrixOperators::insert_columns(g,dlogdet,n);
      MatrixOperators::insert_columns(h,grad_grad_logdet,n);
    }
  }

  void CompositeSPOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, 
    GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
  }


  SPOSetBase* CompositeSPOSetBuilder::createSPOSetFromXML(xmlNodePtr cur)
  {
    vector<string> spolist;
    putContent(spolist,cur);
    if(spolist.empty())
    {
      return 0;
    }
    CompositeSPOSet* spo_now=new CompositeSPOSet;
    for(int i=0; i<spolist.size(); ++i)
    {
      SPOSetBase* spo=get_sposet(spolist[i]);
      if(spo) spo_now->add(spo);
    }
    return (spo_now->size())? spo_now:0;
  }

  SPOSetBase* CompositeSPOSetBuilder::createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input)
  {
    return createSPOSetFromXML(cur);
  }

  bool CompositeSPOSetBuilder::put(xmlNodePtr cur)
  {
    return true;
  }
}
