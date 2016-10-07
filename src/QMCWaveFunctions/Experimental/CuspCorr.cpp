//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include<iostream>

void fillRadFunWithPhiBar(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, RealType* data)
{
  Psi1=Phi;
  Psi2=Eta;
  int norb=Psi1->OrbitalSetSize;
  curOrb=curOrb_;
  curCenter=curCenter_;
  C=data[1];
  sg=data[2];
  Rc=data[3];
  alpha[0]=data[4];
  alpha[1]=data[5];
  alpha[2]=data[6];
  alpha[3]=data[7];
  alpha[4]=data[8];
  Z=Zion;
  val1.resize(norb);
  grad1.resize(norb);
  lapl1.resize(norb);
  for(int i=0; i<xgrid.size(); i++)
  {
    rad_orb[i] = phiBar(xgrid[i])*std::sqrt(4.0*3.14159265358979);
  }
}

void fillRadFunWithPhi(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb)
{
  Psi1=Phi;
  Psi2=Eta;
  int norb=Psi1->OrbitalSetSize;
  curOrb=curOrb_;
  curCenter=curCenter_;
  Z=Zion;
  val1.resize(norb);
  grad1.resize(norb);
  lapl1.resize(norb);
  for(int i=0; i<xgrid.size(); i++)
  {
    rad_orb[i] = phi(xgrid[i])*std::sqrt(4.0*3.14159265358979);
  }
}

void executeWithRCLoop(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, std::string file, RealType cutoff, RealType* data)
{
  RealType bestRc = cutoff, smallX2;
  RealType xa, xb, xc, xd;
  RealType fa,fb,fc, fd,delx;
  int cnt=0;
  // FIX FIX FIX
  // bracket the minimum and use golden search to
  // find optimum Rc
  cnt=0;
  xa=bestRc;
  fa=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xa,data);
  xb=xa+0.005;
  fb=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xb,data);
  if( fb > fa )
  {
    xc=xa;
    xa=xb;
    xb=xc;
    fc=fa;
    fa=fb;
    fb=fc;
  }
  bestRc = xb;
  smallX2=fb;
  delx = 1.61*(xb-xa);
  xd = xb + delx;
  fd=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xd,data);
  while( fb > fd )
  {
    int xx = xb, ff = fb;
    xb = xd;
    fb = fd;
    xa=xx;
    fa=ff;
    xd += delx; // dumb, fix later
    if(xd < 0.0 || xd > Rc_max)
    {
      smallX2==execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,file,bestRc,data);
      return;
    }
    fd=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xd,data);
    if(++cnt == 500)
    {
//          APP_ABORT("Problems bracketing minimum. \n");
      smallX2==execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,file,bestRc,data);
      app_log() <<"Problems bracketing minimum, using best rc so far: " <<bestRc <<"  " <<smallX2 << std::endl;
      return;
    }
  }
  xc=xb + 0.4*(xd-xb);
  fc=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xc,data);
  // find minimum
  //while( std::abs(xa-xd) > 1e-4*(std::abs(xb)+std::abs(xc)) ) // from Num Rec.
  while( std::abs(xa-xd) > 0.01 ) // from Num Rec.
  {
    if(fb > fc)
    {
      xa=xb;
      //xb=xc; fb=fc;
      xb=xa+0.4*(xd-xa);
      fb=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xb,data);
      xc=xa+0.6*(xd-xa);
      fc=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xc,data);
    }
    else
    {
      xd=xc;
      //xc=xb; fc=fb;
      xb=xa+0.4*(xd-xa);
      fb=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xb,data);
      xc=xa+0.6*(xd-xa);
      fc=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,"NULL",xc,data);
    }
  }
  if(fb < fc)
  {
    fb=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,file,xb,data);
  }
  else
  {
    fc=execute(curOrb_,curCenter_,Zion,Phi,Eta,xgrid,rad_orb,file,xc,data);
  }
}

RealType execute(int curOrb_, int curCenter_, RealType Zion, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta, Vector<RealType>& xgrid, Vector<RealType>& rad_orb, std::string thisFile, RealType cutoff, RealType* data)
{
  bool writeout=(thisFile!="NULL");
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
  //Rc=Rc_init;
  Rc=cutoff;
  RealType chi2,chi2Min,phi0,phiMin;
  RealType tmpo;
  evaluate(Psi2,ddr,val2,grad2,lapl2);  // eta(0)
  eta0=val2[curOrb];
  ddr=0;
  ddr[0]=Rc;
  evaluate(Psi1,ddr,val1,grad1,lapl1);
  valAtRc = val1[curOrb];
//    Rc=Rc_init;
//    RealType tmpo=getELorig();
//    getELideal(tmpo);
//    valAtZero = phi(0.0);
//    RealType chi2,chi2Min,phi0=valAtZero,phiMin;
///*
  RealType xa, xb, xc, xd;
  RealType fa,fb,fc, fd,delx;
  int cnt=0;
  // FIX FIX FIX
  // bracket the minimum and use golden search to
  // find optimum Rc
  cnt=0;
  ELorigAtRc=getELorig();
  getELideal(ELorigAtRc);
  // try to bracket the minimum
  xa=phi(0.0);
  fa=loop(xa,X);
  xb=xa+0.005;
  fb=loop(xb,X);
  if( fb > fa )
  {
    xc=xa;
    xa=xb;
    xb=xc;
    fc=fa;
    fa=fb;
    fb=fc;
  }
  delx = 1.61*(xb-xa);
  xd = xb + delx;
  fd=loop(xd,X);
  while( fb > fd )
  {
    int xx = xb, ff = fb;
    xb = xd;
    fb = fd;
    xa=xx;
    fa=ff;
    xd += delx; // dumb, fix later
    fd=loop(xd,X);
//std::cout <<"x,f: " <<xd <<"  " <<fd << std::endl;
//std::cout.flush();
    if(cnt == 50)
      delx*=5;
    if(cnt == 100)
      delx*=5;
    if(++cnt == 1000)
    {
      APP_ABORT("Problems bracketing minimum. \n");
    }
  }
  xc=xb + 0.4*(xd-xb);
  fc=loop(xc,X);
  // find minimum
  while( std::abs(xa-xd) > 1e-5*(std::abs(xb)+std::abs(xc)) ) // from Num Rec.
  {
    //app_log()<<"xa,fa: " <<xa <<"  " <<loop(xa,X) << std::endl
    //   <<"xb,fb: " <<xb <<"  " <<fb << std::endl
    //   <<"xc,fc: " <<xc <<"  " <<fc << std::endl
    //   <<"xd,fd: " <<xd <<"  " <<loop(xd,X) << std::endl << std::endl;
    if(fb > fc)
    {
      xa=xb;
      //xb=xc; fb=fc;
      xb=xa+0.4*(xd-xa);
      fb=loop(xb,X);
      xc=xa+0.6*(xd-xa);
      fc=loop(xc,X);
    }
    else
    {
      xd=xc;
      //xc=xb; fc=fb;
      xb=xa+0.4*(xd-xa);
      fb=loop(xb,X);
      xc=xa+0.6*(xd-xa);
      fc=loop(xc,X);
    }
  }
  //app_log()<<"xa,fa: " <<xa <<"  " <<loop(xa,X) << std::endl
  //   <<"xb,fb: " <<xb <<"  " <<fb << std::endl
  //   <<"xc,fc: " <<xc <<"  " <<fc << std::endl
  //   <<"xd,fd: " <<xd <<"  " <<loop(xd,X) << std::endl;
  if(fb < fc)
  {
    chi2=loop(xb,X);
  }
  else
  {
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
          app_log()<<"i , chi2: " <<i <<"  " <<chi2 << std::endl;
        }
        chi2=loop(phiMin,X);
  */
  if(writeout)
  {
    app_log()<< std::setprecision(12) <<"alpha: " <<alpha[0] <<"  "
             <<alpha[1] <<"  "
             <<alpha[2] <<"  "
             <<alpha[3] <<"  "
             <<alpha[4] <<"  " << std::endl;
    app_log()<<"mo,ion,C, sg, rc, chi2: " <<curOrb_ <<"  " <<curCenter_ <<"  " <<C <<"  " <<sg <<"  " <<Rc <<" " <<chi2 << std::endl;
  }
  for(int i=0; i<xgrid.size(); i++)
  {
    rad_orb[i] = phiBar(xgrid[i])*std::sqrt(4.0*3.14159265358979);
    //rad_orb[i] = phi(xgrid[i])*std::sqrt(4.0*3.14159265358979);
  }
// write report file if desired
  if(writeout && printOrbs)
  {
    std::ofstream out(thisFile.c_str());
    for(int i=0; i<nElms; i++)
    {
      out<<pos[i] <<"  "
         <<ELideal[i] <<"  "
         <<ELorig[i] <<"  "
         <<ELcurr[i] <<"  "
         <<phi(pos[i]) <<"  "
         <<phiBar(pos[i]) <<"  "
         <<pr(pos[i]) <<"  "
         << std::endl;
    }
    out.close();
    std::string full;
    full = thisFile + ".full.dat";
    out.open(full.c_str());
    int nx = 500;
    for(int i=0; i<nx; i++)
    {
      RealType xp = i*5.0/(nx-1.0);
      out<<xp <<"  "
         <<phi(xp) << std::endl;
    }
    out.close();
  }
  data[1]=C;
  data[2]=sg;
  data[3]=Rc;
  data[4]=alpha[0];
  data[5]=alpha[1];
  data[6]=alpha[2];
  data[7]=alpha[3];
  data[8]=alpha[4];
  return chi2;
}

