#ifndef QMCPLUSPLUS_HFDBHE_SMOOTHED_H
#define QMCPLUSPLUS_HFDBHE_SMOOTHED_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {
  /** @ingroup hamiltonian
   *@brief HFDBHE_smoothed_np for the indentical source and target particle sets. 
   */
  struct HFDBHE_smoothed_np: public QMCHamiltonianBase {

    Return_t tailcorr, unsmooth;
    // remember that the default units are Kelvins and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;

    // epsilon = 3.42016039e-5, rm = 5.607384357
    // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
    HFDBHE_smoothed_np(ParticleSet& P): PtclRef(&P) {
      d_table = DistanceTable::add(P);
//      tailcorr = -0.000274387876825;	// N = 32, RHO = 0.022 /Angstrom^3 : 11x
      // Note the 2 powers of N
      // Specific value for:
      // Equilibrium density 0.022 Angstroms^-3
      // N = 32
      // Will be different for other densities.
      // Generalize later.
//      tailcorr = -0.0002704899426;	// N = 32, RHO = 0.3648 sigma^2 : 01x
//      tailcorr = -0.0002650846591;	// N = 64, RHO = 0.3648 sigma^3 : 02x
      tailcorr = -0.0002618259317;	// N = 128, RHO = 0.3648 sigma^2 : 03x
//      tailcorr = -0.0002598311735;	// N = 256, RHO = 0.3648 sigma^2 : 04x
      unsmooth = 0.0;

      app_log() << "  HFDBHE_smoothed_np tail correction is " << tailcorr << endl;
    }

    ~HFDBHE_smoothed_np() { }

    void resetTargetParticleSet(ParticleSet& P)  {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
//      tailcorr = -0.000274387876825;	// N = 32, RHO = 0.022 /Angstroms^3 : 11x
//      tailcorr = -0.0002704899426;	// N = 32, RHO = 0.3648 sigma^2 : 01x
//      tailcorr = -0.0002650846591;	// N = 64, RHO = 0.3648 sigma^3 : 02x
      tailcorr = -0.0002618259317;	// N = 128, RHO = 0.3648 sigma^2 : 03x
//      tailcorr = -0.0002598311735;	// N = 256, RHO = 0.3648 sigma^2 : 04x
      unsmooth = 0.0;
    }


    inline void setValue(Return_t val) {
      unsmooth = val;
      Value = unsmooth + tailcorr;
    }

    inline Return_t evaluate(ParticleSet& P) {
//      Value = unsmooth + tailcorr;
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
      os << "HFD-B(HE)_smoothed_np: " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return false;
    }
  };



  /** @ingroup hamiltonian
   *@brief HFDBHE_smoothed for the indentical source and target particle sets. 
   */
  struct HFDBHE_smoothed: public QMCHamiltonianBase {
    Return_t rc, mirrorshift;
    // remember that the default units are Kelvins and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;
    HFDBHE_smoothed_np* dep;

    // epsilon = 3.42016039e-5, rm = 5.607384357
    // C6 = 1.460008056, C8 = 14.22016431, C10 = 187.2033646;
    HFDBHE_smoothed(ParticleSet& P): PtclRef(&P) {
      Dependants = 1;
      depName = "HFDBHE_np";

      rc = P.Lattice.WignerSeitzRadius;

      d_table = DistanceTable::add(P);

      mirrorshift = -2.0*(6.394277071*std::exp(-1.863335173*rc-0.072712207*rc*rc) - (1.461*std::pow(rc,-6.0)+14.11*std::pow(rc,-8.0)+183.5*std::pow(rc,-10.0))*dampF(rc));
    }

    ~HFDBHE_smoothed() { }

    void resetTargetParticleSet(ParticleSet& P)  {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
      rc = P.Lattice.WignerSeitzRadius;

      mirrorshift = -2.0*(6.394277071*std::exp(-1.863335173*rc-0.072712207*rc*rc) - (1.461*std::pow(rc,-6.0)+14.11*std::pow(rc,-8.0)+183.5*std::pow(rc,-10.0))*dampF(rc));
    }

    inline Return_t evaluate(ParticleSet& P) {
      Return_t smooth = 0.0;
      Value = 0.0;

      for(int i=0; i<d_table->getTotNadj(); i++) {
	Return_t r1 = d_table->r(i);
	if (r1 < rc) {
          Return_t r2i = 1.0/(r1*r1), r6i = r2i*r2i*r2i, r8i = r6i*r2i, r10i = r8i*r2i;
          Return_t rd1 = 2.0*rc - r1, rd2i = 1.0/(rd1*rd1), rd6i = rd2i*rd2i*rd2i, rd8i = rd6i*rd2i, rd10i = rd8i*rd2i;

	  Value += (6.394277071*std::exp(-1.863335173*r1-0.072712207*r1*r1) - (1.461*r6i+14.11*r8i+183.5*r10i)*dampF(r1));
	  smooth += (6.394277071*std::exp(-1.863335173*rd1-0.072712207*rd1*rd1) - (1.461*rd6i+14.11*rd8i+183.5*rd10i)*dampF(rd1) + mirrorshift);
	}
      }
      Value += smooth;
      dep->setValue(-smooth);

      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    inline Return_t dampF(Return_t r) {
      const Return_t D = 8.301460704;
      if (r < D)
	return std::exp(-(D/r - 1.0)*(D/r - 1.0));
      else
	return 1.0;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "HFDB(HE)Potential (T/S): " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new HFDBHE_smoothed(qp);
    }

    QMCHamiltonianBase* makeDependants(ParticleSet& qp)
    {
      dep = new HFDBHE_smoothed_np(qp);
      return dep;
    }
  };
}
#endif
