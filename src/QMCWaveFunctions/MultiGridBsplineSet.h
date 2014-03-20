//////////////////////////////////////////////////////////////////
// (c) Copyright 2012-  by Jeongnim Kim and Ken Esler           //
//////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_EINSPLINE_MULTIGRID_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_MULTIGRID_ADOPTOR_H

#include <QMCWaveFunctions/EinsplineAdoptor.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus
{

  /** Class to test the multi-grid method without any truncation
  */
  template<typename SplineAdoptor>
    struct MultiGridBsplineSet: public SPOSetBase
  {
    ///type of BsplineSet
    typedef BsplineSet<SplineAdoptor> bspline_type;
    ///true with tiling
    bool is_complex;
    ///if true, either localized or extended
    bool IsExclusive;
    ///distance table index
    int myTableIndex;
    ///cutoff radius: use one box for all
    RealType Rcut;
    bspline_type* Extended;
    bspline_type* Localized;

    MultiGridBsplineSet()
      :myTableIndex(1), Rcut(0.5), Extended(0), Localized(0)
    {
      Extended=new bspline_type;
      Localized=new bspline_type;
      is_complex=Extended->is_complex;
    }

    SPOSetBase* makeClone() const
    {
      MultiGridBsplineSet* clone=new MultiGridBsplineSet<SplineAdoptor>(*this);
      clone->Extended=new bspline_type(*Extended);
      clone->Localized=new bspline_type(*Localized);
      return clone;
    }

    /** set spline table
     * @param spline_r real part
     * @param spline_i imaginary part
     * @param ispline index of the spline function
     * @param level 1 for extended/corase orbitals and 0 for localized/dense orbitals 
     */
    template<typename CT>
      void set_spline(CT* spline_r, CT* spline_i, int twist, int ispline, int level)
      {
        if(level)
          Extended->set_spline(spline_r,spline_i,twist,ispline,level);
        else
          Localized->set_spline(spline_r,spline_i,twist,ispline,level);
      }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
        Extended->evaluate_v(P.R[iat],psi);
      else
        Localized->evaluate_v(P.R[iat],psi);
    }

    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
      {
        //cout << iat << " Interstitial " << endl;
        Extended->evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
      }
      else
      {
        //cout << iat << " close to " << ic << " " << P.R[iat]-dt->origin().R[ic] << endl;
        Localized->evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
      }
    }

    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
      {
        //cout << iat << " Interstitial " << endl;
        Extended->evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
      }
      else
      {
        //cout << iat << " close to " << ic << " " << P.R[iat] << " " << dt->origin().R[ic] << endl;
        Localized->evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
      }
    }

    void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      typedef ValueMatrix_t::value_type value_type;
      typedef GradMatrix_t::value_type grad_type;
      for(int iat=first, i=0; iat<last; ++iat,++i)
      {
        VectorViewer<value_type> v(logdet[i],OrbitalSetSize);
        VectorViewer<grad_type> g(dlogdet[i],OrbitalSetSize);
        VectorViewer<value_type> l(d2logdet[i],OrbitalSetSize);
        int ic=dt->find_closest_source(iat,Rcut);
        if(ic<0)
        { 
          //cout << iat << " Interstitial " << endl;
          Extended->evaluate_vgl(P.R[iat],v,g,l);
        }
        else
        {
          //cout << iat << " close to " << ic << " " << P.R[iat]-dt->origin().R[ic] << endl;
          Localized->evaluate_vgl(P.R[iat],v,g,l);
        }
      }
    }

    virtual void evaluate_notranspose(const ParticleSet& P, int first, int last
        , ValueMatrix_t& logdet, GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      typedef ValueMatrix_t::value_type value_type;
      typedef GradMatrix_t::value_type grad_type;
      typedef HessMatrix_t::value_type hess_type;
      for(int iat=first, i=0; iat<last; ++iat,++i)
      {
        VectorViewer<value_type> v(logdet[i],OrbitalSetSize);
        VectorViewer<grad_type> g(dlogdet[i],OrbitalSetSize);
        VectorViewer<hess_type> h(grad_grad_logdet[i],OrbitalSetSize);
        int core=dt->find_closest_source(iat,Rcut);
        if(core<0)
          Extended->evaluate_vgh(P.R[iat],v,g,h);
        else
          Localized->evaluate_vgh(P.R[iat],v,g,h);
      }
    }

    /** implement virtual functions of SPOSetBase */
    void resetParameters(const opt_variables_type& active) { }

    void resetTargetParticleSet(ParticleSet& e) { }

    void setOrbitalSetSize(int norbs)
    {
      OrbitalSetSize = norbs;
      BasisSetSize=norbs;
    }


  };

}
#endif
