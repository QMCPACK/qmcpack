#ifndef QMCPLUSPLUS_HFDHE2MORONI1995_H
#define QMCPLUSPLUS_HFDHE2MORONI1995_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {
  /** @ingroup hamiltonian
   *@brief HFDHE2_Moroni1995_nonphysical for the indentical source and target particle sets. 
   */
  struct HFDHE2_Moroni1995_nonphysical: public QMCHamiltonianBase {

//    Return_t A,alpha,c1,c2,c3,D;
    Return_t rc,unsmoothen,tailcorr;
    // remember that the default units are Hartree and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;

    // epsilon = 3.42016039e-5, rm = 5.607384357
    // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
    HFDHE2_Moroni1995_nonphysical(ParticleSet& P): PtclRef(&P) {
      d_table = DistanceTable::add(P);
      Return_t rho = P.G.size()/P.Lattice.Volume, N0 = P.G.size();
      rc = P.Lattice.WignerSeitzRadius;
      tailcorr = 2.0*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
      unsmoothen = 0.0;

      cout<<"  HFDHE2_Moroni1995_nonphysical tail correction is  "<<tailcorr<<endl;
    }

    ~HFDHE2_Moroni1995_nonphysical() { }

    void resetTargetParticleSet(ParticleSet& P) {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
      Return_t rho = P.G.size()/P.Lattice.Volume, N0 = P.G.size(), rc = P.Lattice.WignerSeitzRadius;
      tailcorr = 2*M_PI*rho*N0*(-26.7433377905*std::pow(rc,-7.0) - 2.8440930339*std::pow(rc,-5.0)-0.486669351961 *std::pow(rc,-3.0)+ std::exp(-2.381392669*rc)*(2.75969257875+6.571911675726*rc+7.82515114293*rc*rc) );
      unsmoothen = 0.0;
    }

    inline void setValue(Return_t val) {
      unsmoothen = val;
      Value = unsmoothen + tailcorr;
    }

    inline Return_t evaluate(ParticleSet& P) {
//      Value = unsmoothen + tailcorr;
      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "HFDHE2_Moroni1995_nonphysical (T/S): " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return false;
    }
  };



  /** @ingroup hamiltonian
   *@brief HFDHE2_Moroni1995_physical for the indentical source and target particle sets. 
   */
  struct HFDHE2_Moroni1995_physical: public QMCHamiltonianBase {
    Return_t rc,A,alpha,c1,c2,c3,D,mirrorshift;
    // remember that the default units are Hartree and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;
    HFDHE2_Moroni1995_nonphysical* dep;

    // epsilon = 3.42016039e-5, rm = 5.607384357
    // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
    HFDHE2_Moroni1995_physical(ParticleSet& P): PtclRef(&P) {
      Dependants = 1;

      A = 18.63475757;
      alpha = -2.381392669;
      c1=1.460008056;
      c2=14.22016431;
      c3=187.2033646;
      D = 6.960524706;
      
      depName = "HFDHE2_np";
      
      d_table = DistanceTable::add(P);
      Return_t rho = P.G.size()/P.Lattice.Volume;
      Return_t N0 = P.G.size();
      rc = P.Lattice.WignerSeitzRadius;

      Return_t rcm2 = 1.0/(rc*rc), rcm6 = rcm2*rcm2*rcm2, rcm8 = rcm6*rcm2, rcm10 = rcm8*rcm2;
      mirrorshift = -2.0*(A*std::exp(alpha*rc) - (c1*rcm6+c2*rcm8+c3*rcm10)*dampF(rc));

    }

    ~HFDHE2_Moroni1995_physical() { }

    void resetTargetParticleSet(ParticleSet& P)  {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
      Return_t rc = P.Lattice.WignerSeitzRadius;

      Return_t rcm2 = 1.0/(rc*rc), rcm6 = rcm2*rcm2*rcm2, rcm8 = rcm6*rcm2, rcm10 = rcm8*rcm2;
      mirrorshift = -2.0*(A*std::exp(alpha*rc) - (c1*rcm6+c2*rcm8+c3*rcm10)*dampF(rc));
    }

    inline Return_t evaluate(ParticleSet& P) {
      Return_t smoothen = 0.0;
      Value = 0.0;
      
      for(int i=0; i<d_table->getTotNadj(); i++) {
	Return_t r1 = d_table->r(i);
        if (r1 < rc) {
          Return_t rm2 = 1.0/(r1*r1), rm6 = rm2*rm2*rm2, rm8 = rm6*rm2, rm10 = rm8*rm2;
          Return_t rd1 = 2.0*rc - r1, rdm2 = 1.0/(rd1*rd1), rdm6 = rdm2*rdm2*rdm2, rdm8 = rdm6*rdm2, rdm10 = rdm8*rdm2;
	  Value += (A*std::exp(alpha*r1) - (c1*rm6+c2*rm8+c3*rm10)*dampF(r1));
	  smoothen += A*std::exp(alpha*rd1) - (c1*rdm6+c2*rdm8+c3*rdm10)*dampF(rd1) + mirrorshift;
	}
      }
      Value += smoothen;
      dep->setValue(-smoothen);

      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    inline Return_t dampF(Return_t r) {
      if (r < D){
        Return_t t1=(D/r - 1.0);
	return std::exp(-t1*t1);
      }
      else
	return 1.0;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "HFDHE2_Moroni1995_physical (T/S): " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new HFDHE2_Moroni1995_physical(qp);
    }

    QMCHamiltonianBase* makeDependants(ParticleSet& qp) {
      dep = new HFDHE2_Moroni1995_nonphysical(qp);
      return dep;
    }
  };
}
#endif
