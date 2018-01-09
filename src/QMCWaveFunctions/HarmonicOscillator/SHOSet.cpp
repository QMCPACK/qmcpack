//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    

#include <QMCWaveFunctions/HarmonicOscillator/SHOSet.h>
#include <Utilities/string_utils.h>


namespace qmcplusplus
{

  SHOSet::SHOSet(RealType l, PosType c, const std::vector<SHOState*>& sho_states)
    : length(l), center(c)
  {
    state_info.resize(sho_states.size());
    for(int s=0;s<sho_states.size();++s)
      state_info[s] = *sho_states[s];
    initialize();
  }


  void SHOSet::initialize()
  {
    using std::sqrt;

    className = "SHOSet";

    OrbitalSetSize = state_info.size();
    BasisSetSize = OrbitalSetSize;

    qn_max = -1;
    for(int s=0;s<state_info.size();++s)
      for(int d=0;d<DIM;++d)
        qn_max[d] = std::max(qn_max[d],state_info[s].quantum_number[d]);
    qn_max+=1;

    nmax = -1;
    for(int d=0;d<DIM;++d)
      nmax = std::max(nmax,qn_max[d]);

    prefactors.resize(nmax);
    hermite.resize(DIM,nmax);
    bvalues.resize(DIM,nmax);

    //d0_values.resize(DIM,nmax);
    //d1_values.resize(DIM,nmax);
    //d2_values.resize(DIM,nmax);

    if(nmax>0)
    {
      prefactors[0] = 1.0/(sqrt(sqrt(M_PI)*length));
      for(int n=1;n<nmax;++n)
        prefactors[n] = prefactors[n-1]/sqrt(2.*n);
    }

  }


  SHOSet::~SHOSet()
  {
  }


  SPOSetBase* SHOSet::makeClone() const
  {
    return new SHOSet(*this);
  }


  void SHOSet::report(const std::string& pad)
  {
    app_log()<<pad<<"SHOSet report"<< std::endl;
    app_log()<<pad<<"  length    = "<< length    << std::endl;
    app_log()<<pad<<"  center    = "<< center    << std::endl;
    app_log()<<pad<<"  nmax      = "<< nmax      << std::endl;
    app_log()<<pad<<"  qn_max    = "<< qn_max    << std::endl;
    app_log()<<pad<<"  # states  = "<< state_info.size() << std::endl;
    app_log()<<pad<<"  states"<< std::endl;
    for(int s=0;s<state_info.size();++s)
      state_info[s].sho_report(pad+"    "+int2string(s)+" ");
    app_log()<<pad<<"end SHOSet report"<< std::endl;
    app_log().flush();
  }


  void SHOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi)
  {
    const PosType &r(P.activeR(iat));
    ValueViewer_t p(&psi[0],size());
    evaluate_v(r,p);
  }


  void SHOSet::evaluate(const ParticleSet& P, int iat, ValueVector_t& psi,
                        GradVector_t& dpsi, ValueVector_t& d2psi)
  {
    const PosType &r(P.activeR(iat));
    ValueViewer_t   p(  &psi[0],size());
    GradViewer_t   dp( &dpsi[0],size());
    ValueViewer_t d2p(&d2psi[0],size());
    evaluate_vgl(r,p,dp,d2p);
  }
  

  void SHOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, ValueMatrix_t& d2logdet)
  {
    for(int iat=first, i=0; iat<last; ++iat,++i)
    {
      ValueViewer_t   p(  logdet[i],size());
      GradViewer_t   dp( dlogdet[i],size());
      ValueViewer_t d2p(d2logdet[i],size());
      evaluate_vgl(P.R[iat],p,dp,d2p);
    } 
  }


  void SHOSet::evaluate_v(PosType r, ValueViewer_t& psi)
  {
    PosType x = (r-center)/length;
    evaluate_hermite(x);
    evaluate_d0(x,psi);
  }
  

  void SHOSet::evaluate_vgl(PosType r, ValueViewer_t& psi, 
                            GradViewer_t& dpsi, ValueViewer_t& d2psi)
  {
    PosType x = (r-center)/length;
    evaluate_hermite(x);
    evaluate_d0(x,psi);
    evaluate_d1(x,psi,dpsi);
    evaluate_d2(x,psi,d2psi);
  }


  void SHOSet::evaluate_hermite(const PosType& xpos)
  {
    int shift=0;
    for(int d=0;d<DIM;++d,shift+=nmax)
    {
      int nh = qn_max[d];
      if(nh>0)
      {
        RealType x = xpos[d];
        hermite(0+shift) = 1.0;
        RealType Hnm2 = 0.0;
        RealType Hnm1 = 1.0;
        for(int n=1;n<nh;++n)
        {
          RealType Hn = 2*(x*Hnm1-(n-1)*Hnm2);
          hermite(n+shift) = Hn;
          Hnm2 = Hnm1;
          Hnm1 = Hn;
        }
      }
    }
  }


  void SHOSet::evaluate_d0(const PosType& xpos, ValueViewer_t& psi)
  {
    using std::exp;
    int shift=0;
    for(int d=0;d<DIM;++d,shift+=nmax)
    {
      RealType x = xpos[d];
      RealType g = exp(-.5*x*x);
      for(int n=0;n<qn_max[d];++n)
      {
        int ns = n+shift;
        bvalues(ns) = prefactors[n]*g*hermite(ns);
        //d0_values(ns) = bvalues(ns);
      }
    }
    for(int s=0;s<state_info.size();++s)
    {
      const SHOState& state = state_info[s];
      RealType phi = 1.0;
      int shift = 0;
      for(int d=0;d<DIM;++d,shift+=nmax)
        phi *= bvalues(shift+state.quantum_number[d]); 
      psi[s] = phi;
    }
  }


  void SHOSet::evaluate_d1(const PosType& xpos, ValueViewer_t& psi, GradViewer_t& dpsi)
  {
    RealType ol = 1.0/length;
    int shift=0;
    for(int d=0;d<DIM;++d,shift+=nmax)
    {
      RealType x = xpos[d];
      RealType Hnm1 = 0.0;
      for(int n=0;n<qn_max[d];++n)
      {
        int ns = n+shift;
        RealType Hn = hermite(ns);
        bvalues(ns) = (-x + 2*n*Hnm1/Hn)*ol;
        //d1_values(ns) = bvalues(ns);
        Hnm1 = Hn;
      }
    }
    for(int s=0;s<state_info.size();++s)
    {
      const SHOState& state = state_info[s];
      TinyVector<ValueType,DIM> dphi;
      int shift = 0;
      for(int d=0;d<DIM;++d,shift+=nmax)
        dphi[d] = bvalues(shift+state.quantum_number[d]); 
      dphi *= psi[s];
      dpsi[s] = dphi;
    }
  }


  void SHOSet::evaluate_d2(const PosType& xpos, ValueViewer_t& psi, ValueViewer_t& d2psi)
  {
    RealType ol2 = 1.0/(length*length);
    int shift = 0;
    for(int d=0;d<DIM;++d,shift+=nmax)
    {
      RealType x = xpos[d];
      RealType x2  = x*x;
      RealType Hnm1 = 0.0;
      RealType Hnm2 = 0.0;
      for(int n=0;n<qn_max[d];++n)
      {
        int ns = n+shift;
        RealType Hn = hermite(ns);
        bvalues(ns) = (-1.0 + x2 - 2*n)*ol2;
        //d2_values(ns) = bvalues(ns);
        Hnm2 = Hnm1;
        Hnm1 = Hn;
      }
    }
    for(int s=0;s<state_info.size();++s)
    {
      const SHOState& state = state_info[s];
      ValueType d2phi=0.0;
      int shift = 0;
      for(int d=0;d<DIM;++d,shift+=nmax)
        d2phi += bvalues(shift+state.quantum_number[d]); 
      d2phi *= psi[s];
      d2psi[s] = d2phi;
    }
  }


  void SHOSet::evaluate_check(PosType r, ValueViewer_t& psi, 
                              GradViewer_t& dpsi, ValueViewer_t& d2psi)
  {
    using std::exp;
    using std::sqrt;

    evaluate_vgl(r,psi,dpsi,d2psi);

    const int N = 6;
    RealType H[N],dH[N],d2H[N],pre[N];
    RealType p[N],dp[N],d2p[N];

    pre[0] = 1.0/(sqrt(sqrt(M_PI)*length));
    for(int n=1;n<N;++n)
      pre[n] = pre[n-1]/sqrt(2.*n);

    for(int d=0;d<DIM;++d)
    {
      RealType x = (r[d]-center[d])/length;
      RealType x2 = x*x,x3=x*x*x,x4=x*x*x*x,x5=x*x*x*x*x;
      H[0] = 1;                  dH[0] = 0;                 d2H[0] = 0;
      H[1] = 2*x;                dH[1] = 2;                 d2H[1] = 0;
      H[2] = 4*x2-2;             dH[2] = 8*x;               d2H[2] = 8;
      H[3] = 8*x3-12*x;          dH[3] = 24*x2-12;          d2H[3] = 48*x;
      H[4] = 16*x4-48*x2+12;     dH[4] = 64*x3-96*x;        d2H[4] = 192*x2-96;
      H[5] = 32*x5-160*x3+120*x; dH[5] = 160*x4-480*x2+120; d2H[5] = 640*x3-960*x;
      RealType g = exp(-x2/2);
      for(int n=0;n<N;++n)
      {
        p[n]   = pre[n]*g*H[n];
        dp[n]  = pre[n]*g*(-x*H[n]+dH[n]);
        d2p[n] = pre[n]*g*((x2-1)*H[n]-2*x*dH[n]+d2H[n]);
      }
      app_log()<<"eval check dim = "<<d<<"  x = "<<x<< std::endl;
      app_log()<<"  hermite check"<< std::endl;
      for(int n=0;n<qn_max[d];++n)
      {
        app_log()<<"    "<<n<<" "<<H[n]<< std::endl;
        app_log()<<"    "<<n<<" "<<hermite(d,n)<< std::endl;
      }
      app_log()<<"  phi d0 check"<< std::endl;
      for(int n=0;n<qn_max[d];++n)
      {
        app_log()<<"    "<<n<<" "<<p[n]<< std::endl;
        app_log()<<"    "<<n<<" "<<d0_values(d,n)<< std::endl;
      }
      app_log()<<"  phi d1 check"<< std::endl;
      for(int n=0;n<qn_max[d];++n)
      {
        app_log()<<"    "<<n<<" "<<dp[n]/p[n]<< std::endl;
        app_log()<<"    "<<n<<" "<<d1_values(d,n)<< std::endl;
      }
      app_log()<<"  phi d2 check"<< std::endl;
      for(int n=0;n<qn_max[d];++n)
      {
        app_log()<<"    "<<n<<" "<<d2p[n]/p[n]<< std::endl;
        app_log()<<"    "<<n<<" "<<d2_values(d,n)<< std::endl;
      }
        
    }
    
  }


  void SHOSet::test_derivatives()
  {
    int     n  = 3;
    PosType c  = 5.123;
    PosType L  = 1.0;
    PosType drg = L/n;
    PosType dr  = L/1000;
    int nphi   = state_info.size();

    PosType o2dr,odr2;

    ValueVector_t   vpsi,vpsitmp; 
    GradVector_t   vdpsi,vdpsin;
    ValueVector_t vd2psi,vd2psin;


    vpsi.resize(nphi);
    vdpsi.resize(nphi);
    vd2psi.resize(nphi);

    vpsitmp.resize(nphi);
    vdpsin.resize(nphi);
    vd2psin.resize(nphi);

    
    ValueViewer_t   psi(  &vpsi[0],size());
    GradViewer_t   dpsi( &vdpsi[0],size());
    ValueViewer_t d2psi(&vd2psi[0],size());
    
    ValueViewer_t   psitmp(&vpsitmp[0],size());
    GradViewer_t   dpsin(   &vdpsin[0],  size());
    ValueViewer_t d2psin(  &vd2psin[0],  size());


    app_log()<<" loading dr"<< std::endl;

    RealType odr2sum = 0.0;
    for(int d=0;d<DIM;++d)
    {
      RealType odr = 1.0/dr[d];
      o2dr[d] = .5*odr;
      odr2[d] = odr*odr;
      odr2sum += odr2[d];
    }

    app_log()<<"SHOSet::test_derivatives"<< std::endl;

    ParticleSet Ps;

    int p=0;
    PosType r,rtmp;
    for(int i=0;i<n;++i)
    {
      r[0] = c[0]+i*drg[0];
      for(int j=0;j<n;++j)
      {
        r[1] = c[1]+j*drg[1];
        for(int k=0;k<n;++k)
        {
          r[2] = c[2]+k*drg[2];

          //evaluate_check(r,psi,dpsi,d2psi);
          //APP_ABORT("SHOSet eval check");

          evaluate_vgl(r,psi,dpsi,d2psi);

          for(int m=0;m<nphi;++m)
            d2psin[m] = -2*odr2sum*psi[m];
          for(int d=0;d<DIM;++d)
          {
            rtmp=r;
            rtmp[d]+=dr[d];
            evaluate_v(rtmp,psitmp);
            for(int m=0;m<nphi;++m)
            {
              ValueType phi = psitmp[m];
              dpsin[m][d] = phi*o2dr[d];
              d2psin[m]  += phi*odr2[d];
            }
            rtmp=r;
            rtmp[d]-=dr[d];
            evaluate_v(rtmp,psitmp);
            for(int m=0;m<nphi;++m)
            {
              ValueType phi = psitmp[m];
              dpsin[m][d]-= phi*o2dr[d];
              d2psin[m]  += phi*odr2[d];
            }
          }

          RealType dphi_diff = 0.0;
          RealType d2phi_diff = 0.0;
          for(int m=0;m<nphi;++m)
            for(int d=0;d<DIM;++d)
              dphi_diff = std::max<RealType>(dphi_diff,std::abs(dpsi[m][d]-dpsin[m][d])/std::abs(dpsin[m][d]));
          for(int m=0;m<nphi;++m)
            d2phi_diff = std::max<RealType>(d2phi_diff,std::abs(d2psi[m]-d2psin[m])/std::abs(d2psin[m]));
          app_log()<<"  "<<p<<" "<<dphi_diff<<" "<<d2phi_diff<< std::endl;
          app_log()<<"    derivatives"<< std::endl;
          for(int m=0;m<nphi;++m)
          {
            std::string qn="";
            for(int d=0;d<DIM;++d)
              qn+=int2string(state_info[m].quantum_number[d])+" ";
            app_log()<<"    "<<qn;
            for(int d=0;d<DIM;++d)
              app_log()<<real(dpsi[m][d])<<" ";
            app_log()<< std::endl;
            app_log()<<"    "<<qn;
            for(int d=0;d<DIM;++d)
              app_log()<<real(dpsin[m][d])<<" ";
            app_log()<< std::endl;
          }
          app_log()<<"    laplacians"<< std::endl;
          PosType x=r/length;
          RealType x2 = dot(x,x);
          for(int m=0;m<nphi;++m)
          {
            std::string qn="";
            for(int d=0;d<DIM;++d)
              qn+=int2string(state_info[m].quantum_number[d])+" ";
            app_log()<<"    "<<qn<<real(d2psi[m]/psi[m])<< std::endl;
            app_log()<<"    "<<qn<<real(d2psin[m]/psi[m])<< std::endl;
          }
          p++;
        }
      }
    }

    app_log()<<"end SHOSet::test_derivatives"<< std::endl;
  }


  void SHOSet::test_overlap()
  {
    app_log()<<"SHOSet::test_overlap"<< std::endl;


    //linear
    int d=0;

    app_log()<<"  length = "<<length<< std::endl;
    app_log()<<"  prefactors"<< std::endl;
    for(int n=0;n<qn_max[d];++n)
      app_log()<<"    "<<n<<" "<<prefactors[n]<< std::endl;

    app_log()<<"  1d overlap"<< std::endl;

    ValueVector_t vpsi;
    vpsi.resize(size());
    ValueViewer_t psi(  &vpsi[0],size());

    double xmax = 4.0;
    double dx = .1;
    double dr = length*dx;
    
    int nphi = qn_max[d];
    Array<double,2> omat;
    omat.resize(nphi,nphi);
    for(int i=0;i<nphi;++i)
      for(int j=0;j<nphi;++j)
        omat(i,j) = 0.0;

    PosType xp = 0.0;
    for(double x=-xmax;x<xmax;x+=dx)
    {
      xp[d] = x;
      evaluate_hermite(xp);
      evaluate_d0(xp,psi);
      
      for(int i=0;i<nphi;++i)
        for(int j=0;j<nphi;++j)
          omat(i,j) += bvalues(d,i)*bvalues(d,j)*dr;
    }

    for(int i=0;i<nphi;++i){
      app_log()<< std::endl;
      for(int j=0;j<nphi;++j)
        app_log()<<omat(i,j)<<" ";
    }
    app_log()<< std::endl;



    //volumetric
    app_log()<<"  3d overlap"<< std::endl;
    double dV = dr*dr*dr;
    nphi = size();
    omat.resize(nphi,nphi);
    for(int i=0;i<nphi;++i)
      for(int j=0;j<nphi;++j)
        omat(i,j) = 0.0;
    for(double x=-xmax;x<xmax;x+=dx)
      for(double y=-xmax;y<xmax;y+=dx)
        for(double z=-xmax;z<xmax;z+=dx)
        {
          xp[0]=x; xp[1]=y; xp[2]=z;
          evaluate_hermite(xp);
          evaluate_d0(xp,psi);
          
          for(int i=0;i<nphi;++i)
            for(int j=0;j<nphi;++j)
              omat(i,j) += std::abs(psi[i]*psi[j])*dV;
        }
    for(int i=0;i<nphi;++i){
      app_log()<< std::endl;
      for(int j=0;j<nphi;++j)
        app_log()<<omat(i,j)<<" ";
    }
    app_log()<< std::endl;


    app_log()<<"end SHOSet::test_overlap"<< std::endl;
  }




  //methods to be implemented later
  void SHOSet::resetParameters(const opt_variables_type& optVariables)
  {
    not_implemented("resetParameters");
  }
    
  void SHOSet::evaluate(
    const ParticleSet& P, PosType &r, ValueVector_t &psi)
  {
    not_implemented("evaluate(P,r,psi)");
  }

  void SHOSet::evaluate(
    const ParticleSet& P, int iat,ValueVector_t& psi, 
    GradVector_t& dpsi, HessVector_t& ddpsi)
  {
    not_implemented("evaluate(P,iat,psi,dpsi,ddpsi)");
  }

  void SHOSet::evaluateThirdDeriv(
    const ParticleSet& P,int first,int last,GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluateThirdDeriv(P,first,last,dddlogdet)");
  }

  void SHOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet)");
  }

  void SHOSet::evaluate_notranspose(
    const ParticleSet& P, int first, int last, ValueMatrix_t& logdet, 
    GradMatrix_t& dlogdet, HessMatrix_t& grad_grad_logdet, 
    GGGMatrix_t& grad_grad_grad_logdet)
  {
    not_implemented("evaluate_notranspose(P,first,last,logdet,dlogdet,ddlogdet,dddlogdet)");
  }

  void SHOSet::evaluateGradSource(
    const ParticleSet &P, int first, int last, 
    const ParticleSet &source, int iat_src, GradMatrix_t &gradphi)
  {
    not_implemented("evaluateGradSource(P,first,last,source,iat,dphi)");
  }

  void SHOSet::evaluateGradSource(
    const ParticleSet &P, int first, int last, const ParticleSet &source, 
    int iat_src, GradMatrix_t &grad_phi, HessMatrix_t &grad_grad_phi, 
    GradMatrix_t &grad_lapl_phi)
  {
    not_implemented("evaluateGradSource(P,first,last,source,iat,dphi,ddphi,dd2phi)");
  }

  void SHOSet::evaluateBasis(
    const ParticleSet &P, int first, int last, ValueMatrix_t &basis_val,  
    GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl)
  {
    not_implemented("evaluateBasis");
  }

  void SHOSet::evaluateForDeriv(
    const ParticleSet &P, int first, int last, ValueMatrix_t &basis_val,  
    GradMatrix_t  &basis_grad, ValueMatrix_t &basis_lapl)
  {
    not_implemented("evaluateForDeriv");
  }

  void SHOSet::copyParamsFromMatrix(
    const opt_variables_type& active, const ValueMatrix_t &mat, 
    std::vector<RealType> &destVec)
  {
    not_implemented("copyParamsFromMatrix");
  }


}
