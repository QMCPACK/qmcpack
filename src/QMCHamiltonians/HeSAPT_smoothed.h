#ifndef QMCPLUSPLUS_HESAPT_SMOOTHED_H
#define QMCPLUSPLUS_HESAPT_SMOOTHED_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"

namespace qmcplusplus {
  /** @ingroup hamiltonian
   *@brief HeSAPT_smoothed_np for the indentical source and target particle sets. 
   *
   * \f[ H = \sum_i \frac{q^2}{r} \f] 
   * where \f$ q \f$ is the charge of the set of quantum 
   * particles.  For instance, \f$ q = -1 \f$ for electrons 
   * and \f$ q = 1 \f$ for positrons.
   */
  struct HeSAPT_smoothed_np: public QMCHamiltonianBase {
    Return_t tailcorr, unsmooth;
    // remember that the default units are Hartrees and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;

    HeSAPT_smoothed_np(ParticleSet& P): d_table(NULL), PtclRef(&P) {
      d_table = DistanceTable::add(P);
//      tailcorr = -0.000274497151179;	// N = 32, RHO = 0.022 /Angstrom^3 : 11x
      // N = 32, rho = 0.022 Angstroms^-3 = 0.003260063604 a0^-3
      // Will be different for other densities.
      // Generalize later.
//      tailcorr = -0.0002705961159;	// N = 32, RHO = 0.3648 sigma^2 : 01x
      tailcorr = -0.0002651066072;	// N = 64, RHO = 0.3648 sigma^2 : 02x
//      tailcorr = -0.0002618281770;	// N = 128, RHO = 0.3648 sigma^2 : 03x
//      tailcorr = -0.0002598285275;	// N = 256, RHO = 0.3648 sigma^2 : 04x
//      tailcorr = -0.0002178620081;	// N = 32, RHO = 0.3280 sigma^2 : 21x
//      tailcorr = -0.0002137701388;	// N = 64, RHO = 0.3280 sigma^2 : 22x
//      tailcorr = -0.0002113168481;	// N = 128, RHO = 0.3280 sigma^2 : 23x
//      tailcorr = -0.0002098169307;	// N = 256, RHO = 0.3280 sigma^2 : 24x
//      tailcorr = -0.0002443783556;	// N = 32, RHO = 0.3470 sigma^2 : 31x
//      tailcorr = -0.0002395972888;	// N = 64, RHO = 0.3470 sigma^2 : 32x
//      tailcorr = -0.0002367366424;	// N = 128, RHO = 0.3470 sigma^2 : 33x
//      tailcorr = -0.0002349898273;	// N = 256, RHO = 0.3470 sigma^2 : 34x
//      tailcorr = -0.0003280694504;	// N = 32, RHO = 0.4009 sigma^2 : 41x
//      tailcorr = -0.0003209410359;	// N = 64, RHO = 0.4009 sigma^2 : 42x
//      tailcorr = -0.0003166990570;	// N = 128, RHO = 0.4009 sigma^2 : 43x
//      tailcorr = -0.0003141175814;	// N = 256, RHO = 0.4009 sigma^2 : 44x
//      tailcorr = -0.0003926875054;	// N = 32, RHO = 0.4378 sigma^2 : 51x
//      tailcorr = -0.0003835898243;	// N = 64, RHO = 0.4378 sigma^2 : 52x
//      tailcorr = -0.0003781940744;	// N = 128, RHO = 0.4378 sigma^2 : 53x
//      tailcorr = -0.0003749179381;	// N = 256, RHO = 0.4378 sigma^2 : 54x
//      tailcorr = -0.0004941981444;	// N = 32, RHO = 0.4899 sigma^2 : 61x
//      tailcorr = -0.0004817718362;	// N = 64, RHO = 0.4899 sigma^2 : 62x
//      tailcorr = -0.0004744316477;	// N = 128, RHO = 0.4899 sigma^2 : 63x
//      tailcorr = -0.0004699888991;	// N = 256, RHO = 0.4899 sigma^2 : 64x
      unsmooth = 0.0;

      app_log() << "  HeSAPT_smoothed_np tail correction is " << tailcorr << endl;
    }

    ~HeSAPT_smoothed_np() { }

    void resetTargetParticleSet(ParticleSet& P)  {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
//      tailcorr = -0.000274497151179;	// N = 32, RHO = 0.022 /Angstrom^3 : 11x
//      tailcorr = -0.0002705961159;	// N = 32, RHO = 0.3648 sigma^2 : 01x
      tailcorr = -0.0002651066072;	// N = 64, RHO = 0.3648 sigma^2 : 02x
//      tailcorr = -0.0002618281770;	// N = 128, RHO = 0.3648 sigma^2 : 03x
//      tailcorr = -0.0002598285275;	// N = 256, RHO = 0.3648 sigma^2 : 04x
//      tailcorr = -0.0002178620081;	// N = 32, RHO = 0.3280 sigma^2 : 21x
//      tailcorr = -0.0002137701388;	// N = 64, RHO = 0.3280 sigma^2 : 22x
//      tailcorr = -0.0002113168481;	// N = 128, RHO = 0.3280 sigma^2 : 23x
//      tailcorr = -0.0002098169307;	// N = 256, RHO = 0.3280 sigma^2 : 24x
//      tailcorr = -0.0002443783556;	// N = 32, RHO = 0.3470 sigma^2 : 31x
//      tailcorr = -0.0002395972888;	// N = 64, RHO = 0.3470 sigma^2 : 32x
//      tailcorr = -0.0002367366424;	// N = 128, RHO = 0.3470 sigma^2 : 33x
//      tailcorr = -0.0002349898273;	// N = 256, RHO = 0.3470 sigma^2 : 34x
//      tailcorr = -0.0003280694504;	// N = 32, RHO = 0.4009 sigma^2 : 41x
//      tailcorr = -0.0003209410359;	// N = 64, RHO = 0.4009 sigma^2 : 42x
//      tailcorr = -0.0003166990570;	// N = 128, RHO = 0.4009 sigma^2 : 43x
//      tailcorr = -0.0003141175814;	// N = 256, RHO = 0.4009 sigma^2 : 44x
//      tailcorr = -0.0003926875054;	// N = 32, RHO = 0.4378 sigma^2 : 51x
//      tailcorr = -0.0003835898243;	// N = 64, RHO = 0.4378 sigma^2 : 52x
//      tailcorr = -0.0003781940744;	// N = 128, RHO = 0.4378 sigma^2 : 53x
//      tailcorr = -0.0003749179381;	// N = 256, RHO = 0.4378 sigma^2 : 54x
//      tailcorr = -0.0004941981444;	// N = 32, RHO = 0.4899 sigma^2 : 61x
//      tailcorr = -0.0004817718362;	// N = 64, RHO = 0.4899 sigma^2 : 62x
//      tailcorr = -0.0004744316477;	// N = 128, RHO = 0.4899 sigma^2 : 63x
//      tailcorr = -0.0004699888991;	// N = 256, RHO = 0.4899 sigma^2 : 64x
      unsmooth = 0.0;
    }


    inline void setValue (Return_t val) {
      unsmooth = val;
      Value = unsmooth + tailcorr;
    }

    inline Return_t evaluate(ParticleSet& P) {
      // Value = unsmooth + tailcorr;
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
      os << "HeSAPT_smoothed_np: " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return false;
    }
  };


  /** @ingroup hamiltonian
   *@brief HeSAPT_smoothed for the indentical source and target particle sets. 
   *
   * \f[ H = \sum_i \frac{q^2}{r} \f] 
   * where \f$ q \f$ is the charge of the set of quantum 
   * particles.  For instance, \f$ q = -1 \f$ for electrons 
   * and \f$ q = 1 \f$ for positrons.
   */
  struct HeSAPT_smoothed: public QMCHamiltonianBase {
    Return_t rc, mirrorshift;
    // remember that the default units are Hartrees and Bohrs
    DistanceTableData* d_table;
    ParticleSet* PtclRef;
    HeSAPT_smoothed_np* dep;

    HeSAPT_smoothed(ParticleSet& P): d_table(NULL), PtclRef(&P) {
      Dependants = 1;
      depName = "HeSAPT_np";
      d_table = DistanceTable::add(P);
      rc = P.Lattice.WignerSeitzRadius;

      mirrorshift = -2.0*((-4.27467067939-32.3451957988*rc-37.5775397337*rc*rc+4.0/rc)*std::exp(-5.72036885319*rc) + (18.2962387439-6.16632555293*rc+6.91482730781*rc*rc)*std::exp(-2.80857770752*rc) - damp(7,rc)*1.460977837725/std::pow(rc,6.) - damp(9,rc)*14.11785737/std::pow(rc,8.) - damp(11,rc)*183.691075/std::pow(rc,10.) + damp(12,rc)*76.70/std::pow(rc,11.) - damp(13,rc)*3372./std::pow(rc,12.) + damp(14,rc)*3806./std::pow(rc,13.) - damp(15,rc)*85340./std::pow(rc,14.) + damp(16,rc)*171000./std::pow(rc,15.) - damp(17,rc)*2860000./std::pow(rc,16.));
    }

    ~HeSAPT_smoothed() { }

    void resetTargetParticleSet(ParticleSet& P)  {
      d_table = DistanceTable::add(P);
      PtclRef=&P;
      rc = P.Lattice.WignerSeitzRadius;

      mirrorshift = -2.0*((-4.27467067939-32.3451957988*rc-37.5775397337*rc*rc+4.0/rc)*std::exp(-5.72036885319*rc) + (18.2962387439-6.16632555293*rc+6.91482730781*rc*rc)*std::exp(-2.80857770752*rc) - damp(7,rc)*1.460977837725/std::pow(rc,6.) - damp(9,rc)*14.11785737/std::pow(rc,8.) - damp(11,rc)*183.691075/std::pow(rc,10.) + damp(12,rc)*76.70/std::pow(rc,11.) - damp(13,rc)*3372./std::pow(rc,12.) + damp(14,rc)*3806./std::pow(rc,13.) - damp(15,rc)*85340./std::pow(rc,14.) + damp(16,rc)*171000./std::pow(rc,15.) - damp(17,rc)*2860000./std::pow(rc,16.));
    }

    inline Return_t evaluate(ParticleSet& P) {
      Return_t smooth = 0.0;
      Value = 0.0;

      for(int i=0; i < d_table->getTotNadj(); i++) {
	Return_t r1 = d_table->r(i);
	if (r1 < rc) {
          Return_t r1i = 1.0/r1, r2i = r1i*r1i, r6i = r2i*r2i*r2i, r8i = r6i*r2i, r10i = r8i*r2i, r11i = r10i*r1i, r12i = r10i*r2i, r13i = r11i*r2i, r14i = r12i*r2i, r15i = r13i*r2i, r16i = r14i*r2i;
          Return_t rd1 = 2.0*rc - r1, rd1i = 1.0/rd1, rd2i = rd1i*rd1i, rd6i = rd2i*rd2i*rd2i, rd8i = rd6i*rd2i, rd10i = rd8i*rd2i, rd11i = rd10i*rd1i, rd12i = rd10i*rd2i, rd13i = rd11i*rd2i, rd14i = rd12i*rd2i, rd15i = rd13i*rd2i, rd16i = rd14i*rd2i;

	  Value += ((-4.27467067939-32.3451957988*r1-37.5775397337*r1*r1+4.0*r1i)*std::exp(-5.72036885319*r1) + (18.2962387439-6.16632555293*r1+6.91482730781*r1*r1)*std::exp(-2.80857770752*r1) - damp(7,r1)*1.460977837725*r6i - damp(9,r1)*14.11785737*r8i - damp(11,r1)*183.691075*r10i + damp(12,r1)*76.70*r11i - damp(13,r1)*3372.*r12i + damp(14,r1)*3806.*r13i - damp(15,r1)*85340.*r14i + damp(16,r1)*171000.*r15i - damp(17,r1)*2860000.*r16i);
	  smooth += ((-4.27467067939-32.3451957988*rd1-37.5775397337*rd1*rd1+4.0*rd1i)*std::exp(-5.72036885319*rd1) + (18.2962387439-6.16632555293*rd1+6.91482730781*rd1*rd1)*std::exp(-2.80857770752*rd1) - damp(7,rd1)*1.460977837725*rd6i - damp(9,rd1)*14.11785737*rd8i - damp(11,rd1)*183.691075*rd10i + damp(12,rd1)*76.70*rd11i - damp(13,rd1)*3372.*rd12i + damp(14,rd1)*3806.*rd13i - damp(15,rd1)*85340.*rd14i + damp(16,rd1)*171000.*rd15i - damp(17,rd1)*2860000.*rd16i + mirrorshift);
	}
      }
      Value += smooth;
      dep->setValue(-smooth);
      // temporary tail corr.  See Mathematica notebook.
      return Value;
    }

    inline Return_t evaluate(ParticleSet& P, vector<NonLocalData>& Txy) {
      return evaluate(P);
    }

    Return_t damp(int n, Return_t r) {
      Return_t wert, br, sum, term;
      int i;

      br = 2.41324077320*r;

      if (r >= 1.0) {
	sum = term = 1.0;
	for (i=1;i<=n;i++) {
	  term *= (br/i);
	  sum += term;
	}
	wert = 1.0 - std::exp(-br)*sum;
      } else {
	sum = 0.0;
	term = 1.0;
	for (i=1;i<=n;i++)
	  term *= (br/i);
	for (i=n+1;i<=30;i++) {
	  term *= (br/i);
	  sum += term;
	}
	wert = std::exp(-br)*sum;
      }

      return wert;
    }

    /** Do nothing */
    bool put(xmlNodePtr cur) {
      return true;
    }

    bool get(std::ostream& os) const {
      os << "HeSAPT_smoothed: " << PtclRef->getName();
      return true;
    }

    QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
    {
      return new HeSAPT_smoothed(qp);
    }
    
    QMCHamiltonianBase* makeDependants(ParticleSet& qp)
    {
      dep = new HeSAPT_smoothed_np(qp);
      return dep;
    }
  };
}
#endif
