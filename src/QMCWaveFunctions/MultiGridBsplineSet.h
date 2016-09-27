//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef QMCPLUSPLUS_EINSPLINE_MULTIGRID_ADOPTOR_H
#define QMCPLUSPLUS_EINSPLINE_MULTIGRID_ADOPTOR_H

#include <QMCWaveFunctions/EinsplineAdoptor.h>
#include <Particle/DistanceTableData.h>

namespace qmcplusplus
{

  /** BsplineSet using a multigrid method
   *
   * Inspired by PRL 82, 5016 (1999) and generalization of the truncated orbitals for vacuum.
   */
  template<typename SplineAdoptor>
    struct MultiGridBsplineSet: public SPOSetBase
  {
    ///typedef of the spline data
    typedef typename SplineAdoptor::DataType DataType;
    typedef typename SplineAdoptor::SingleSplineType SingleSplineType;
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
    ///SPOSet on the global grid of the primitive cell
    bspline_type* Extended;
    ///SPOSet on the localized grid centered at the ions
    std::vector<bspline_type*> Localized;
    ///grid converter for each subdomain
    std::vector<GridConvert<DataType> > gTransform;
    ///copy of PCID of the source particleset
    std::vector<int> PCID;

    MultiGridBsplineSet()
      :myTableIndex(1), Rcut(0.5), Extended(0)
    {
      Extended=new bspline_type;
      is_complex=Extended->is_complex;
    }

    SPOSetBase* makeClone() const
    {
      MultiGridBsplineSet* clone=new MultiGridBsplineSet<SplineAdoptor>(*this);
      clone->Extended=new bspline_type(*Extended);
      for(int i=0; i<Localized.size(); ++i)
      {
        if(Localized[i])
          clone->Localized[i]=new bspline_type(*Localized[i]);
      }
      return clone;
    }

    /** create subdomains
     * @param ntot number of ions
     * @param n number of ions in a primitive cell
     */
    void resizeSubDomains(int ntot, int n)
    {
      Localized.resize(std::max(4,n),0);
      for(int i=0; i< n; ++i) Localized[i]=new bspline_type;
      gTransform.resize(n);
      PCID.resize(ntot);
    }

    /** allocate a subdomain and its GridConver
     * @param ic index of the subdomain
     * @param dense single spline which defines the origintal grid
     * @param lower Carteisan coordiate at the origin
     * @param upper Cartesian ccordiate at (1,1,1)  of the subdomain
     * @return size of the table
     */
    template <typename PT>
      size_t setSubDomain(int ic, typename SplineAdoptor::SingleSplineType* dense, PT& lower, PT& upper)
      {
        gTransform[ic].create(Localized[ic]->MultiSpline,dense,lower,upper
            ,Extended->MultiSpline->num_splines);

        char s[1024];
        sprintf(s," Offset = %d %d %d Ngrid = %d %d %d",
            gTransform[ic].Offset[0],gTransform[ic].Offset[1],gTransform[ic].Offset[2],
            gTransform[ic].N[0],     gTransform[ic].N[1],     gTransform[ic].N[2]);
        app_log() << "setSubDomain " << s << std::endl;
        return sizeof(DataType)*(Localized[ic]->MultiSpline->coefs_size);
      }

    size_t sizeOfExtended() const
    {
      return sizeof(DataType)*(Extended->MultiSpline->coefs_size);
    }

    /** set spline table
     * @param spline_r real part
     * @param spline_i imaginary part
     * @param ispline index of the spline function
     * @param center -1 for extended/corase orbitals and other positive index for localized/dense orbitals 
     */
    void set_spline(SingleSplineType* spline_r, SingleSplineType* spline_i, int twist, int ispline, int center)
    {
      const int dummy=0;
      if(center<0)
      {
        Extended->set_spline(spline_r,spline_i,twist,ispline,dummy);
      }
      else 
      {
        Localized[center]->set_spline_domain(spline_r,spline_i,twist,ispline, gTransform[center].Offset,gTransform[center].N);
      }
    }

    inline void evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
        Extended->evaluate_v(P.R[iat],psi);
      else
        Localized[PCID[ic]]->evaluate_v(P.R[iat],psi);
    }

    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, ValueVector_t& d2psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
        Extended->evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
      else
        Localized[PCID[ic]]->evaluate_vgl(P.R[iat],psi,dpsi,d2psi);
    }

    inline void evaluate(const ParticleSet& P, int iat,
        ValueVector_t& psi, GradVector_t& dpsi, HessVector_t& grad_grad_psi)
    {
      const DistanceTableData* dt=P.DistTables[myTableIndex];
      int ic=dt->find_closest_source(Rcut);
      if(ic<0)
        Extended->evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
      else
        Localized[PCID[ic]]->evaluate_vgh(P.R[iat],psi,dpsi,grad_grad_psi);
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
          Extended->evaluate_vgl(P.R[iat],v,g,l);
        else
          Localized[PCID[ic]]->evaluate_vgl(P.R[iat],v,g,l);
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
        int ic=dt->find_closest_source(iat,Rcut);
        if(ic<0)
          Extended->evaluate_vgh(P.R[iat],v,g,h);
        else
          Localized[PCID[ic]]->evaluate_vgh(P.R[iat],v,g,h);
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
