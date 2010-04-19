#include<iostream>

  template<class BS>
  void CuspCorr<BS>::execute(int curOrb_, int curCenter_, double Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<double>& xgrid, Vector<double>& rad_orb, string thisFile)
  {
    Psi1=Phi;  
    Psi2=Eta; 
    int norb=Psi1->OrbitalSetSize;
    ValueVector_t X(5);
    curOrb=curOrb_;
    curCenter=curCenter_;

    C=0.0;
    Z=Zion;

    ELideal.resize(nElms);
    ELorig.resize(nElms);
    ELcurr.resize(nElms);
    pos.resize(nElms);
    val1.resize(norb);
    val2.resize(norb);
    grad1.resize(norb);
    grad2.resize(norb);
    lapl1.resize(norb);
    lapl2.resize(norb);

    TinyVector<RealType,3> ddr=0;
    Rc=Rc_init;
    RealType chi2,chi2Min,phi0,phiMin;
    RealType tmpo;
    evaluate(Psi2,ddr,val2,grad2,lapl2);  // eta(0)
    eta0=val2[curOrb];

//    Rc=Rc_init;
//    RealType tmpo=getELorig();
//    getELideal(tmpo);
//    valAtZero = phi(0.0);
//    RealType chi2,chi2Min,phi0=valAtZero,phiMin;
///*
    double xa, xb, xc, xd;
    double fa,fb,fc, fd,delx;
    int cnt=0;
    // FIX FIX FIX
    // bracket the minimum and use golden search to 
    // find optimum Rc 
      cnt=0;
      tmpo=getELorig();
      getELideal(tmpo);
    // try to bracket the minimum
      xa=phi(0.0);
      fa=loop(xa,X);
      xb=xa+0.005;
      fb=loop(xb,X);
      if( fb > fa ) {
        xc=xa;xa=xb;xb=xc;
        fc=fa;fa=fb;fb=fc;
      }

      delx = 1.61*(xb-xa);
      xd = xb + delx;
      fd=loop(xd,X);
      while( fb > fd ) {
        xd += delx; // dumb, fix later 
        fd=loop(xd,X); 
        if(++cnt == 50) { 
          APP_ABORT("Problems bracketing minimum. \n");  
        }
      } 
      xc=xb + 0.4*(xd-xb);
      fc=loop(xc,X);
         
    // find minimum
      while( std::fabs(xa-xd) > 1e-5*(std::fabs(xb)+std::fabs(xc)) ) // from Num Rec. 
      {
      //app_log()<<"xa,fa: " <<xa <<"  " <<loop(xa,X) <<endl
      //   <<"xb,fb: " <<xb <<"  " <<fb <<endl
      //   <<"xc,fc: " <<xc <<"  " <<fc <<endl
      //   <<"xd,fd: " <<xd <<"  " <<loop(xd,X) <<endl <<endl;
     
        if(fb > fc) {
          xa=xb;
          //xb=xc; fb=fc;
          xb=xa+0.4*(xd-xa); 
          fb=loop(xb,X);
          xc=xa+0.6*(xd-xa); 
          fc=loop(xc,X);
        } else {
          xd=xc;
          //xc=xb; fc=fb;
          xb=xa+0.4*(xd-xa); 
          fb=loop(xb,X);
          xc=xa+0.6*(xd-xa); 
          fc=loop(xc,X);
        }  
      }  
      app_log()<<"xa,fa: " <<xa <<"  " <<loop(xa,X) <<endl
         <<"xb,fb: " <<xb <<"  " <<fb <<endl
         <<"xc,fc: " <<xc <<"  " <<fc <<endl
         <<"xd,fd: " <<xd <<"  " <<loop(xd,X) <<endl;
      if(fb < fc) {
        chi2=loop(xb,X);
      } else {
        chi2=loop(xc,X);
      }

//*/
/*
      for(int i=0; i<70; i++) {
        phi0 += sg*i*0.001; 
        chi2=loop(phi0,X);
        if(i==0) { 
          chi2Min=chi2;
          phiMin=phi0;
        } else if(chi2 < chi2Min) {
          chi2Min=chi2;
          phiMin=phi0;
        } 
        app_log()<<"i , chi2: " <<i <<"  " <<chi2 <<endl;
      }
      chi2=loop(phiMin,X);
*/
      app_log()<<setprecision(12) <<"alpha: " <<alpha[0] <<"  "
             <<alpha[1] <<"  "
             <<alpha[2] <<"  "
             <<alpha[3] <<"  "
             <<alpha[4] <<"  " <<endl;
      app_log()<<"C, sg: " <<C <<"  " <<sg <<endl;
      app_log()<<"chi2: " <<chi2 <<endl; 

    for(int i=0; i<xgrid.size(); i++) {
      rad_orb[i] = phiBar(xgrid[i])*std::sqrt(4.0*3.14159265358979);
      //rad_orb[i] = phi(xgrid[i])*std::sqrt(4.0*3.14159265358979);
    }
// write report file if desired 
    std::ofstream out(thisFile.c_str());
    for(int i=0; i<nElms; i++) {
      out<<pos[i] <<"  " 
         <<ELideal[i] <<"  "
         <<ELorig[i] <<"  "
         <<ELcurr[i] <<"  "
         <<phi(pos[i]) <<"  "
         <<phiBar(pos[i]) <<"  "
         <<pr(pos[i]) <<"  "
         <<endl;
    }
    out.close();

  }

