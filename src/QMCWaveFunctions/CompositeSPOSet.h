//////////////////////////////////////////////////////////////////
// (c) Copyright 2013-  by Jaron T. Krogel                      //
//////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COMPOSITE_SPOSET_H
#define QMCPLUSPLUS_COMPOSITE_SPOSET_H

#include <QMCWaveFunctions/SPOSetBase.h>
#include <QMCWaveFunctions/BasisSetBase.h>

namespace qmcplusplus
{

  class CompositeSPOSet : public SPOSetBase
  {
  public:
    ///component SPOSets
    vector<SPOSetBase*>    components;
    ///temporary storage for values
    vector<ValueVector_t*> component_values;
    ///temporary storage for gradients
    vector<GradVector_t*>  component_gradients;
    ///temporary storage for laplacians
    vector<ValueVector_t*> component_laplacians;
    ///store the precomputed offsets
    vector<int> component_offsets;

    CompositeSPOSet();
    ~CompositeSPOSet();

    ///add a sposet component to this composite sposet
    void add(SPOSetBase* component);

    ///print out component info
    void report();

    //SPOSetBase interface methods
    ///size is determined by component sposets and nothing else
    inline void setOrbitalSetSize(int norbs) { }

    void resetTargetParticleSet(ParticleSet& P);

    SPOSetBase* makeClone() const;

    /** add sposet clones from another Composite SPOSet
     *   should only be used in makeClone functions following shallow copy
     */
    void clone_from(const CompositeSPOSet& master);

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi);

    void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi, 
                  GradVector_t& dpsi, ValueVector_t& d2psi);

    ///unimplemented functions call this to abort
    inline void not_implemented(const string& method)
    {
      APP_ABORT("CompositeSPOSet::"+method+" has not been implemented");
    }


    //methods to be implemented in the future (possibly)
    void resetParameters(const opt_variables_type& optVariables);
    void evaluate(const ParticleSet& P, PosType &r, ValueVector_t &psi);
    void evaluate(const ParticleSet& P, int iat,ValueVector_t& psi, 
                  GradVector_t& dpsi, HessVector_t& ddpsi);
    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              ValueMatrix_t& d2logdet);
    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              HessMatrix_t& ddlogdet);
    void evaluate_notranspose(const ParticleSet& P, int first, int last, 
                              ValueMatrix_t& logdet, GradMatrix_t& dlogdet, 
                              HessMatrix_t& ddlogdet, GGGMatrix_t& dddlogdet);

  };

  struct CompositeSPOSetBuilder : public BasisSetBuilder
  {

    //BasisSetBuilder interface
    SPOSetBase* createSPOSetFromXML(xmlNodePtr cur);

    SPOSetBase* createSPOSet(xmlNodePtr cur,SPOSetInputInfo& input);
    
    bool put(xmlNodePtr cur);
  };
}

#endif
