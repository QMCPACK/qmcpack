//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#ifndef OHMMS_EWALDSUM_H
#define OHMMS_EWALDSUM_H

#include <map>
#include "ParticleBase/ParticleBase.h"
#include "ParticleBase/PtclPairListBase.h"

template<class T, bool EVUNIT, bool AAUNIT>
struct EwaldSumTraits {};

template<class T>
struct EwaldSumTraits<T,true,true>
{

  T ConvE;
  T ConvF;

  EwaldSumTraits()
  {
    const T bohr = 0.529177e0; // in AA
    const T rydberg = 13.6058e0; // in eV
    ConvE = rydberg*bohr;     // default energy conversion
    ConvF = rydberg*bohr*bohr;// default force conversion
  }
};

template<class T>
struct EwaldSumTraits<T,false,false>
{

  T ConvE;
  T ConvF;

  EwaldSumTraits()
  {
    ConvE = 1.0;
    ConvF = 1.0;
  }
};

/** class to perform EwalkSum
 *
 * The template parameters
 * - PT is the trait class for a aprticle
 *  Usage: EwaldSum<PtclOnLatticeTraits> in a QMCHamiltonianBase
 * - EVUNIT false, the energy is in hartree
 * - AAUNIT false, the length unit is AA
 */
template<class PT>
struct EwaldSum: public EwaldSumTraits<typename PT::Scalar_t,false,false>
{

  typedef ParticleBase<PT> Particle_t;
  typedef typename Particle_t::ParticleLayout_t    ParticleLayout_t;
  typedef typename Particle_t::Scalar_t            Scalar_t;
  typedef typename Particle_t::SingleParticleIndex_t SingleParticleIndex_t;
  typedef typename Particle_t::SingleParticlePos_t SingleParticlePos_t;
  typedef typename Particle_t::Tensor_t            Tensor_t;

  typedef typename Particle_t::ParticlePos_t       ParticlePos_t;
  typedef typename Particle_t::ParticleScalar_t    ParticleScalar_t;

  typedef PtclPairListBase<double,3>  PtclPairList_t;

  EwaldSum()
  {
    Eps = -1;
    Tr_over_Tg = 10.0/1.0; // PIV with icc, need to be modified on different machines
  }

  ~EwaldSum() { }

  void init(ParticleLayout_t& lat, const ParticleScalar_t& q);

  void makecells(const ParticleLayout_t& lattice);

  Scalar_t energyG(const ParticleLayout_t& lattice,
                   const ParticlePos_t& r,
                   const ParticleScalar_t& q);

  Scalar_t energyR(const ParticleLayout_t& lattice,
                   const ParticlePos_t& r,
                   const ParticleScalar_t& q);

  Scalar_t energyR(const PtclPairList_t& nnlist,
                   int first, int last,
                   const ParticleScalar_t& q);

  Scalar_t forcestress(const ParticleLayout_t& lattice,
                       const ParticlePos_t& r,
                       const ParticleScalar_t& q,
                       ParticlePos_t& f,
                       Tensor_t& stress);

  void print(std::ostream& os)
  {
    os.setf(std::ios::scientific, std::ios::floatfield);
    os << "Eps                 = " << Eps << std::endl;
    os << "Reciprocal     sum  = " << gamrec << std::endl;
    os << "Direct         sum  = " << gamdir << std::endl;
    os << "Point self-senergy  = " << s3 << std::endl;
    os << "Charged system term = " << s4 << std::endl;
    os << "Maximum K components= " << maxG << std::endl;
    os << "Maximum R components= " << maxR << std::endl;
    os << "Rmax                = " << Rmax << std::endl;
  }

  ///Eps = alpha^2
  Scalar_t Eps;

  /// ratio to determine an optimal eps
  Scalar_t Tr_over_Tg;

  SingleParticleIndex_t maxG;
  SingleParticleIndex_t maxR;

  Scalar_t Rmax;

  //!< tolerance
  Scalar_t acclog;
  //!< |G_max|^2 for a given Eps
  Scalar_t g2max;
  //!< |R_n|^2   for a given Eps
  Scalar_t x2max;
  //!< alpha = sqrt(eps)
  Scalar_t sqeps;
  //!< sum over reciprocal terms
  Scalar_t gamrec;
  //!< sum over direct terms
  Scalar_t gamdir;

  //!< point self-energy; charged system term
  Scalar_t s3, s4;

  Scalar_t Volume;
  Scalar_t E0;
  Scalar_t Cg;
  Scalar_t Cr;
  std::vector<SingleParticleIndex_t> RShell;
};


template<class PT>
void
EwaldSum<PT>::init(ParticleLayout_t& lat, const ParticleScalar_t& q)
{
  int nat = q.size();
  const Scalar_t accur =1.0e-10;
  acclog = log(accur);
  Eps =
    M_PI*(pow(static_cast<Scalar_t>(nat)*Tr_over_Tg/(lat.Volume*lat.Volume),1.0/3.0));
  g2max = 4.0*Eps*std::abs(acclog);
  x2max = std::abs(acclog)/Eps;
  Rmax = sqrt(x2max);
  //Using fixed cutoff
  //Rmax = 6.1;
  //x2max = Rmax*Rmax;
  //Eps = std::abs(acclog)/x2max;
  //g2max = 4.0*Eps*std::abs(acclog);
  sqeps = sqrt(Eps);
  maxG[0] = static_cast<int>(sqrt( g2max/dot(lat.b(0),lat.b(0)))) + 1;
  maxG[1] = static_cast<int>(sqrt( g2max/dot(lat.b(1),lat.b(1)))) + 1;
  maxG[2] = static_cast<int>(sqrt( g2max/dot(lat.b(2),lat.b(2)))) + 1;
  maxR[0] = static_cast<int>(sqrt( x2max/dot(lat.a(0),lat.a(0)))) + 0;
  maxR[1] = static_cast<int>(sqrt( x2max/dot(lat.a(1),lat.a(1)))) + 0;
  maxR[2] = static_cast<int>(sqrt( x2max/dot(lat.a(2),lat.a(2)))) + 0;
//   std::map<int,std::vector<SingleParticleIndex_t>*> gs;
//   std::map<int,std::vector<SingleParticleIndex_t>*> rs;
//   int ic=0;
//   for(int ix1=-maxR[0]; ix1<=maxR[0]; ix1++) {
//     for(int ix2=-maxR[1]; ix2<=maxR[1]; ix2++) {
//       for(int ix3=-maxR[2]; ix3<=maxR[2]; ix3++) {
// 	ic++;
// 	int ih=ix1*ix1+ix2*ix2+ix3*ix3;
// 	std::map<int,std::vector<SingleParticleIndex_t>*>::iterator it = rs.find(ih);
// 	if(it == rs.end()) {
// 	  std::vector<SingleParticleIndex_t>* ns = new std::vector<SingleParticleIndex_t>;
// 	  ns->push_back(SingleParticleIndex_t(ix1,ix2,ix3));
//           rs[ih] = ns;
// 	} else {
// 	  (*it).second->push_back(SingleParticleIndex_t(ix1,ix2,ix3));
// 	}
//       }
//     }
//   }
//   RShell.resize(ic);
//   int ir=0;
//   std::map<int,std::vector<SingleParticleIndex_t>*>::const_iterator cit = rs.begin();
//   while(cit != rs.end()) {
//     const std::vector<SingleParticleIndex_t>& iv = *((*cit).second);
//     for(int i=0; i<iv.size(); i++) RShell[ir++] = iv[i];
//     cit++;
//   }
//   std::map<int,std::vector<SingleParticleIndex_t>*>::iterator it = rs.begin();
//   while(it != rs.end()) {
//     delete (*it).second;it++;
//   }
  //the third and fourth sums
  s3 = 0.0;
  s4 = 0.0;
  for(int iat=0; iat<nat; iat++)
  {
    s3 += q[iat]*q[iat];
    s4 += q[iat];
  }
  Volume = lat.Volume;
  Cg = M_PI/(Volume*Eps)*ConvE;
  Cr = sqeps*ConvE;
  const Scalar_t spi = sqrt(M_PI);
  E0 =  - ConvE*2.0*sqeps/spi*s3 - Cg*s4*s4;
  //R.setUnit(PosUnit::CartesianUnit);
}

/** evaluate the Ewald-Sum
 *@param lat the lattice
 *@param Rcart the input Cartesian Position
 *@param q the input point charge
 */
template<class PT>
typename EwaldSum<PT>::Scalar_t
EwaldSum<PT>::energyG(const ParticleLayout_t& lat,
                      const ParticlePos_t& Rcart,
                      const ParticleScalar_t& q)
{
  //  int nat = rin.size();
  int nat = Rcart.size();
  //R.resize(nat);//make sure the size is consistent
  //convert the current positions to Cartesian coordinates in a super cell
  //wrapAroundBox(lat, rin, R);
  Scalar_t enorm = 1.0/(4.0*Eps);
  gamrec = 0.0;
  for(int ig1=-maxG[0]; ig1<= maxG[0]; ig1++)
  {
    for(int ig2=-maxG[1]; ig2<=maxG[1]; ig2++)
    {
      for(int ig3=-maxG[2]; ig3<=maxG[2]; ig3++)
      {
        if(ig1 == 0 && ig2 == 0 && ig3 == 0)
          continue;	// exclude G = 0;
        // this G
        SingleParticlePos_t tau =
          static_cast<Scalar_t>(ig1)*lat.b(0) +
          static_cast<Scalar_t>(ig2)*lat.b(1) +
          static_cast<Scalar_t>(ig3)*lat.b(2);
        tau *= TWOPI; // multiply 2*pi
        Scalar_t tau2 = dot(tau,tau);
        Scalar_t t2e  = tau2*enorm;// |2\pi G|^2/(4*Eps)
        if ( -t2e < acclog)
          continue;	//  contribution neglegible
        Scalar_t expfac = exp( - t2e)/t2e; // 4\eps correction later
        //gamma-ewald: sum over kappa
        Scalar_t sumr = 0.0;
        Scalar_t sumi = 0.0;
        for(int iat=0; iat<nat; iat++)
        {
          Scalar_t kdotr = dot(Rcart[iat], tau);
          sumr += q[iat]*cos(kdotr);
          sumi += q[iat]*sin(kdotr);
        }
        Scalar_t term = expfac*(sumr*sumr + sumi*sumi);
        gamrec += term; // unit = 1/(1/A)^2 = A^2
      }//loop-ig3
    }//loop-ig2
  }//loop-ig1
  return Cg*gamrec;
}

template<class PT>
typename EwaldSum<PT>::Scalar_t
EwaldSum<PT>::energyR(const ParticleLayout_t& lat,
                      const ParticlePos_t& Rcart,
                      const ParticleScalar_t& q)
{
  //summation over the direct lattice
  int nat = q.size();
  gamdir = 0.0;
  for(int i=0; i<RShell.size(); i++)
  {
    SingleParticlePos_t xlp = lat.toCart(RShell[i]);
    for(int iat=0; iat<nat; iat++)
    {
      for(int jat=0; jat<nat; jat++)
      {
        if(i==0 && iat==jat)
          continue;
        SingleParticlePos_t xxx = Rcart[iat]-Rcart[jat]-xlp;
        Scalar_t xxx2 = dot(xxx,xxx);
        // (l-prim,kappa) = (0,kappa0) is excluded
        //if (xxx2 < 1.0e-8) continue;
        Scalar_t arg = sqrt(xxx2) * sqeps;
        // neglegible contribution
        if ( -arg < acclog)
          continue;
        // h(x) = erfc(x)/x, sqrt(eps) is corrected later
        Scalar_t hx = erfc(arg)/arg;
        //.....gamma-ewald:
        gamdir += q[iat]*q[jat]*hx;
      }//jat
    }//iat
  }
  return Cr*gamdir;
}

template<class PT>
typename EwaldSum<PT>::Scalar_t
EwaldSum<PT>::energyR(const PtclPairList_t& nnlist,
                      int first, int last,
                      const ParticleScalar_t& q)
{
  int nnloc = 0;
  gamdir = 0.0;
  for(int iat=first, iL = 0; iat<last; iat++, iL++)
  {
    for(int nni = 0; nni< nnlist.nadj(iL); nni++, nnloc++)
    {
      Scalar_t arg = nnlist.R[nnloc] * sqeps;
      // neglegible contribution
      if ( -arg < acclog)
        continue;
      // h(x) = erfc(x)/x, sqrt(eps) is corrected later
      Scalar_t hx = erfc(arg)/arg;
      gamdir += q[iat]*q[nnlist.iadj(iL,nni)]*hx;
    }
  }
  return Cr*gamdir;
}


template<class PT>
typename EwaldSum<PT>::Scalar_t
EwaldSum<PT>::forcestress(const ParticleLayout_t& lat,
                          const ParticlePos_t& Rcart,
                          const ParticleScalar_t& q,
                          ParticlePos_t& f0,
                          Tensor_t& stress)
{
  const Scalar_t pi = M_PI;
  const Scalar_t twopi = 2*M_PI;
  const Scalar_t spi = sqrt(M_PI);
  int nat = Rcart.size();
//   if(Eps < 0) {
//     // Eps = \pi *(t_R/T_G N/V^2)^{1/3}
//     // t_R/t_G ~ 1/9 optimal on PIV with icc
//     Eps = pi*(pow(static_cast<Scalar_t>(nat)*Tr_over_Tg/(lat.Volume*lat.Volume),1.0/3.0));
//     init();
//   }
//   maxg1 = static_cast<int>(sqrt( g2max/dot(lat.b(0),lat.b(0)))) + 1;
//   maxg2 = static_cast<int>(sqrt( g2max/dot(lat.b(1),lat.b(1)))) + 1;
//   maxg3 = static_cast<int>(sqrt( g2max/dot(lat.b(2),lat.b(2)))) + 1;
//   maxx1 = static_cast<int>(sqrt( x2max/dot(lat.a(0),lat.a(0)))) + 1;
//   maxx2 = static_cast<int>(sqrt( x2max/dot(lat.a(1),lat.a(1)))) + 1;
//   maxx3 = static_cast<int>(sqrt( x2max/dot(lat.a(2),lat.a(2)))) + 1;
  //R.resize(nat);//make sure the size is consistent
  //convert the current positions to Cartesian coordinates in a super cell
  //wrapAroundBox(lat, rin, R);
  Scalar_t vol = lat.Volume;
  //!< pi/(V*eps) factor for rec force
  Scalar_t rec_fac = 2.0*ConvF*pi/(vol*Eps);
  Tensor_t stress_rec, stress_dir;//stress for rec and dir terms
  Scalar_t enorm = 1.0/(4.0*Eps);
  gamrec = 0.0;
  for(int ig1=-maxG[0]; ig1<= maxG[0]; ig1++)
  {
    for(int ig2=-maxG[1]; ig2<=maxG[1]; ig2++)
    {
      for(int ig3=-maxG[2]; ig3<=maxG[2]; ig3++)
      {
        if(ig1 == 0 && ig2 == 0 && ig3 == 0)
          continue;
        // this G
        SingleParticlePos_t tau =
          static_cast<Scalar_t>(ig1)*lat.b(0) +
          static_cast<Scalar_t>(ig2)*lat.b(1) +
          static_cast<Scalar_t>(ig3)*lat.b(2);
        tau *= twopi; // multiply 2*pi
        Scalar_t tau2 = dot(tau,tau);   // G^2
        Scalar_t t2e  = tau2*enorm; //G^2/4/Eps
        if ( -t2e < acclog)
          continue;	//contribution neglegible
        Scalar_t expfac = exp( - t2e)/t2e;
        //gamma-ewald: sum over kappa
        Scalar_t sumr = 0.0;
        Scalar_t sumi = 0.0;
        for(int iat=0; iat<nat; iat++)
        {
          Scalar_t kdotr = dot(Rcart[iat], tau);
          sumr += q[iat]*cos(kdotr);
          sumi += q[iat]*sin(kdotr);
        }
        Scalar_t term   = expfac*(sumr*sumr + sumi*sumi);
        gamrec += term;
        //stress:
        Scalar_t factor = term * 2.0*(t2e + 1.0)/tau2;
        for(int idir=0; idir<3; idir++)
        {
          for(int jdir=0; jdir<3; jdir++)
          {
            Scalar_t sum = tau[idir]*tau[jdir]*factor;
            if(idir == jdir)
              sum -= term;
            stress_rec(idir,jdir) += sum;
          }
        }
        //forces: loop over the atoms kappa-0
        for(int iat=0; iat<nat; iat++)
        {
          Scalar_t sum = 0.0;
          // summation over kappa
          for(int jat=0; jat<nat; jat++)
          {
            SingleParticlePos_t dr =Rcart[iat]-Rcart[jat];
            sum += q[jat]*sin(dot(tau,dr));
          }
          sum *= expfac*rec_fac;
          f0[iat] += sum*tau;
        }//loop-iat
      }//loop-ig3
    }//loop-ig2
  }//loop-ig1
  // Eps factor for direct force x2
  Scalar_t dir_fac = 2.0*ConvF*Eps;
  gamdir = 0.0;
  //summation over the direct lattice
  for(int i=0; i<RShell.size(); i++)
  {
    SingleParticlePos_t xlp = lat.toCart(RShell[i]);
    for(int iat=0; iat<nat; iat++)
    {
      for(int jat=0; jat<nat; jat++)
      {
        if(i==0 && iat==jat)
          continue;
        SingleParticlePos_t xxx = Rcart[iat]-Rcart[jat]-xlp;
        Scalar_t xxx2 = dot(xxx,xxx);
        // (l-prim,kappa) = (0,kappa0) is excluded
        if (xxx2 < 1.0e-8)
          continue;
        Scalar_t arg = sqrt(xxx2) * sqeps;
        // neglegible contribution
        if ( -arg < acclog)
          continue;
        // function h(x) = erfc(x)/x
        Scalar_t hx = erfc(arg)/arg;
        //.....gamma-ewald:
        gamdir += q[iat]*q[jat]*hx;
        //
        //stress:
        // function h'(x)*x
        Scalar_t  arg2 = arg*arg, hprim;
        if(-arg2>acclog)
          hprim = - hx - 2.0/spi*exp( - arg2);
        else
          hprim = - hx;
        Scalar_t factor = q[iat]*q[jat]*hprim/xxx2;
        stress_dir += factor*outerProduct(xxx,xxx);
        //forces
        hprim = hprim/arg;//h'(x)
        factor = hprim*q[jat]/sqrt(xxx2)*dir_fac;
        // Important: f0 contains conversion factor + sign
        f0[iat] -= factor*xxx;
      }//jat
    }//iat
  }
  //the third and fourth sums
  s3 = 0.0;
  s4 = 0.0;
  for(int iat=0; iat<nat; iat++)
  {
    s3 += q[iat]*q[iat];
    s4 += q[iat];
  }
  rec_fac = ConvE*pi/(vol*Eps);//!< recip term including 2
  dir_fac = ConvE*sqeps;       //!< direct term including 2
  for(int idir=0; idir<OHMMS_DIM; idir++)
  {
    for(int jdir=0; jdir<OHMMS_DIM; jdir++)
    {
      stress(idir,jdir) +=
        rec_fac*stress_rec(idir,jdir) + dir_fac*stress_dir(idir,jdir);
      if(idir==jdir)
        stress(idir,jdir) += rec_fac*s4*s4;
    }
  }
  /** reciprocal term:
    \f$ \frac{2\pi}{V} \frac{1}{4\epsilon} V_G= \frac{\pi}{2 V \epsilon}V_G \f$
    direct term:
    \f$\frac{\sqrt{\epsilon}}{2} V_R\f$
    unit = (au)^2/A
  */
  return
    rec_fac*gamrec+dir_fac*gamdir-2.0*ConvE*sqeps/spi*s3 - rec_fac*s4*s4;
}

#endif


