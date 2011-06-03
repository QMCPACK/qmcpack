//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//   Tel:    217-244-6319 (NCSA) 217-333-3324 (MCC)
//
// Supported by
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//   Department of Physics, Ohio State University
//   Ohio Supercomputer Center
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
using namespace std;
#include "Particle/MCWalkerConfiguration.h"
#include "Particle/DistanceTable.h"
#include "ParticleBase/ParticleUtility.h"
#include "ParticleBase/RandomSeqGenerator.h"
#include "Message/Communicate.h"
#include "QMCDrivers/WaveFunctionTester.h"
#include "QMCDrivers/DriftOperators.h"
#include "Utilities/OhmmsInform.h"
#include "LongRange/StructFact.h"
#include "OhmmsData/AttributeSet.h"
#include "OhmmsData/ParameterSet.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "QMCWaveFunctions/Fermion/SlaterDet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "Numerics/DeterminantOperators.h"
#include "Numerics/SymmetryOperations.h"
#include "Numerics/Blasf.h"

namespace qmcplusplus
  {


  WaveFunctionTester::WaveFunctionTester(MCWalkerConfiguration& w,
                                         TrialWaveFunction& psi,
                                         QMCHamiltonian& h,
                                         ParticleSetPool &ptclPool, WaveFunctionPool& ppool):
      QMCDriver(w,psi,h,ppool),checkRatio("no"),checkClone("no"), checkHamPbyP("no"),
      PtclPool(ptclPool), wftricks("no"),checkEloc("no")
  {
    m_param.add(checkRatio,"ratio","string");
    m_param.add(checkClone,"clone","string");
    m_param.add(checkHamPbyP,"hamiltonianpbyp","string");
    m_param.add(sourceName,"source","string");
    m_param.add(wftricks,"orbitalutility","string");
    m_param.add(checkEloc,"printEloc","string");

    char fname[16];
    sprintf(fname,"wftest.%03d",OHMMS::Controller->rank());

    fout=new ofstream(fname);
    fout->precision(15);

  }

  WaveFunctionTester::~WaveFunctionTester()
  {
  }

  /*!
   * \brief Test the evaluation of the wavefunction, gradient and laplacian
   by comparing to the numerical evaluation.
   *
   Use the finite difference formulas formulas
   \f[
   \nabla_i f({\bf R}) = \frac{f({\bf R+\Delta r_i}) - f({\bf R})}{2\Delta r_i}
   \f]
   and
   \f[
   \nabla_i^2 f({\bf R}) = \sum_{x,y,z} \frac{f({\bf R}+\Delta x_i)
   - 2 f({\bf R}) + f({\bf R}-\Delta x_i)}{2\Delta x_i^2},
   \f]
   where \f$ f = \ln \Psi \f$ and \f$ \Delta r_i \f$ is a
   small displacement for the ith particle.
  */

  bool
  WaveFunctionTester::run()
  {

    app_log() << "Starting a Wavefunction tester" << endl;

    //DistanceTable::create(1);

    put(qmcNode);

    if (checkRatio == "yes")
      {
        runRatioTest();
        runRatioTest2();
      }
    else if (checkClone == "yes")
      runCloneTest();
    else if(checkEloc != "no")   
      printEloc();
    else if (sourceName.size() != 0)
      {
        runGradSourceTest();
        runZeroVarianceTest();
      }
      else if (checkRatio =="deriv")
        runDerivTest();
      else if (checkRatio =="derivclone")
        runDerivCloneTest();
    else if (wftricks !="no")
      runwftricks();
    else
      runBasicTest();

    RealType ene = H.evaluate(W);

    *fout << " Energy " << ene << endl;

    return true;
  }

  void WaveFunctionTester::runCloneTest()
  {

    for (int iter=0; iter<4; ++iter)
      {
        ParticleSet* w_clone = new MCWalkerConfiguration(W);
        TrialWaveFunction *psi_clone = Psi.makeClone(*w_clone);
        QMCHamiltonian *h_clone = H.makeClone(*w_clone,*psi_clone);
        h_clone->setPrimary(false);

        IndexType nskipped = 0;
        RealType sig2Enloc=0, sig2Drift=0;
        RealType delta = 0.00001;
        RealType delta2 = 2*delta;
        ValueType c1 = 1.0/delta/2.0;
        ValueType c2 = 1.0/delta/delta;

        int nat = W.getTotalNum();

        ParticleSet::ParticlePos_t deltaR(nat);
        MCWalkerConfiguration::PropertyContainer_t Properties;
        //pick the first walker
        MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

        //copy the properties of the working walker
        Properties = awalker->Properties;

        W.R = awalker->R;
        W.update();
        ValueType logpsi1 = Psi.evaluateLog(W);
        RealType eloc1  = H.evaluate(W);

        w_clone->R=awalker->R;
        w_clone->update();
        ValueType logpsi2 = psi_clone->evaluateLog(*w_clone);
        RealType eloc2  = h_clone->evaluate(*w_clone);

        app_log() << "Testing walker-by-walker functions " << endl;
        app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << endl;
        app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << endl;


        app_log() << "Testing pbyp functions " << endl;
        Walker_t::Buffer_t &wbuffer(awalker->DataSet);
        wbuffer.clear();
        app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;
        //W.registerData(wbuffer);
        logpsi1=Psi.registerData(W,wbuffer);
        eloc1= H.evaluate(W);
        app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;

        wbuffer.clear();
        app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;
        //w_clone->registerData(wbuffer);
        logpsi2=psi_clone->registerData(W,wbuffer);
        eloc2= H.evaluate(*w_clone);
        app_log() << "  Walker Buffer State current=" << wbuffer.current() << " size=" << wbuffer.size() << endl;

        app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << endl;
        app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << endl;

        delete h_clone;
        delete psi_clone;
        delete w_clone;
      }
  }

  void WaveFunctionTester::printEloc()
  {
    ParticleSetPool::PoolType::iterator p;
    for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
      app_log() << "ParticelSet = " << p->first << endl;

    // Find source ParticleSet
    ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
    ParticleSet& source = *((*pit).second);
    app_log() << "Source = " <<sourceName <<"  " <<(*pit).first << endl;

    int nel = W.getTotalNum();
    int ncenter = source.getTotalNum();
//    int ncenter = 3;
//    cout<<"number of centers: " <<source.getTotalNum() <<endl;
//    cout<<"0: " <<source.R[0] <<endl;
//    cout<<"1: " <<source.R[1] <<endl;
//    cout<<"2: " <<source.R[2] <<endl;

    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

    //copy the properties of the working walker
    Properties = awalker->Properties;

    W.R = awalker->R;

    W.update();
    //ValueType psi = Psi.evaluate(W);
    ValueType logpsi = Psi.evaluateLog(W);
    RealType eloc=H.evaluate(W);

    app_log() << "  Logpsi: " <<logpsi  << endl;
    app_log() << "  HamTest " << "  Total " <<  eloc << endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

    //RealType psi = Psi.evaluateLog(W);

    //int iat=0;
    double maxR = 1000000.0;
    vector<int> closestElectron(ncenter);

    for(int iat=0; iat<ncenter; iat++) {
      maxR=10000000;
      for(int k=0; k<nel; k++)
      {
        double dx = std::sqrt( (W.R[k][0]-source.R[iat][0])*(W.R[k][0]-source.R[iat][0])  
                              +(W.R[k][1]-source.R[iat][1])*(W.R[k][1]-source.R[iat][1]) 
                              +(W.R[k][2]-source.R[iat][2])*(W.R[k][2]-source.R[iat][2])); 
        if(dx < maxR) { maxR = dx; closestElectron[iat]=k; }
      }
    }

//    closestElectron[iat]=1;   
    ofstream out("eloc.dat");
    double x,dx=1.0/499.0;
    for (int k=0; k<500; k++)
    {
      x=-0.5+k*dx; 
      out<<x <<"  ";
      for(int iat=0; iat<ncenter; iat++) {
        PosType tempR = W.R[closestElectron[iat]];
        W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;

        W.R[closestElectron[iat]][0] += x;
        W.update();
        ValueType logpsi_p = Psi.evaluateLog(W); 
        ValueType ene = H.evaluate(W);
        out<<ene <<"  ";
        W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;

        W.R[closestElectron[iat]][1] += x;
        W.update();
        logpsi_p = Psi.evaluateLog(W);
        ene = H.evaluate(W);
        out<<ene <<"  ";
        W.R[closestElectron[iat]]=source.R[iat];
//        W.R[closestElectron[iat]]=0.0;

        W.R[closestElectron[iat]][2] += x;
        W.update();
        logpsi_p = Psi.evaluateLog(W);
        ene = H.evaluate(W);
        out<<ene <<"  ";

        W.R[closestElectron[iat]] = tempR;
      } 
      out<<endl;

    }
    out.close();
  }

  void WaveFunctionTester::runBasicTest()
  {
    IndexType nskipped = 0;
    RealType sig2Enloc=0, sig2Drift=0;
    RealType delta = 1e-5;
    RealType delta2 = 2*delta;
    ValueType c1 = 1.0/delta/2.0;
    ValueType c2 = 1.0/delta/delta;

    int nat = W.getTotalNum();

    ParticleSet::ParticlePos_t deltaR(nat);
    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

    //copy the properties of the working walker
    Properties = awalker->Properties;

    //sample a new walker configuration and copy to ParticleSet::R
    //makeGaussRandom(deltaR);

    W.R = awalker->R;

    //W.R += deltaR;

    W.update();
    //ValueType psi = Psi.evaluate(W);
    ValueType logpsi = Psi.evaluateLog(W);
    RealType eloc=H.evaluate(W);

    app_log() << "  Logpsi: " <<logpsi  << endl;
    app_log() << "  HamTest " << "  Total " <<  eloc << endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

    //RealType psi = Psi.evaluateLog(W);
    ParticleSet::ParticleGradient_t G(nat), G1(nat);
    ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
    G = W.G;
    L = W.L;

//#if defined(QMC_COMPLEX)
//  ValueType logpsi(std::log(psi));
//#else
//  ValueType logpsi(std::log(std::abs(psi)));
//#endif

/*
{
   int iat=0;
   ofstream out("eloc.dat");
   double dx=7.0/999.0;
   W.R[iat] = 0.0;
   W.R[iat][0] = -3.5;
   for (int k=0; k<1000; k++)
   {
       W.R[iat][0] += dx;
       W.update();
       ValueType logpsi_p = Psi.evaluateLog(W); 
       ValueType ene = H.evaluate(W);

       out<<W.R[iat][0] <<"  "
           <<std::exp(logpsi_p) <<"  "
           <<logpsi_p <<"  "
           <<ene <<endl;
   }
   out.close();
}
*/

    for (int iat=0; iat<nat; iat++)
      {
        PosType r0 = W.R[iat];
        GradType g0;
        ValueType lap = 0.0;
        for (int idim=0; idim<OHMMS_DIM; idim++)
          {

            W.R[iat][idim] = r0[idim]+delta;
            W.update();
            ValueType logpsi_p = Psi.evaluateLog(W);

            W.R[iat][idim] = r0[idim]-delta;
            W.update();
            ValueType logpsi_m = Psi.evaluateLog(W);
            lap += logpsi_m+logpsi_p;
            g0[idim] = logpsi_p-logpsi_m;
//#if defined(QMC_COMPLEX)
//      lap += std::log(psi_m) + std::log(psi_p);
//      g0[idim] = std::log(psi_p)-std::log(psi_m);
//#else
//      lap += std::log(std::abs(psi_m)) + std::log(abs(psi_p));
//      g0[idim] = std::log(std::abs(psi_p/psi_m));
//#endif
            W.R[iat] = r0;
          }
        G1[iat] = c1*g0;
        L1[iat] = c2*(lap-2.0*OHMMS_DIM*logpsi);
        *fout << "G1 = " << G1[iat] << endl;
        *fout << "L1 = " << L1[iat] << endl;
      }

    for (int iat=0; iat<nat; iat++)
      {
        *fout << "For particle #" << iat << " at " << W.R[iat] << endl;
        *fout << "Gradient      = " << setw(12) << G[iat] << endl
             << "  Finite diff = " << setw(12) << G1[iat] << endl
             << "  Error       = " << setw(12) << G[iat]-G1[iat] << endl << endl;
        *fout << "Laplacian     = " << setw(12) << L[iat] << endl
             << "  Finite diff = " << setw(12) << L1[iat] << endl
             << "  Error       = " << setw(12) << L[iat]-L1[iat] << endl << endl;
      }


    makeGaussRandom(deltaR);
    //testing ratio alone
    for (int iat=0; iat<nat; iat++)
      {
        W.update();
        //ValueType psi_p = log(fabs(Psi.evaluate(W)));
        RealType psi_p = Psi.evaluateLog(W);
        RealType phase_p=Psi.getPhase();

        W.makeMove(iat,deltaR[iat]);
        //W.update();
        RealType aratio = Psi.ratio(W,iat);
        W.rejectMove(iat);
        Psi.rejectMove(iat);

        W.R[iat] += deltaR[iat];
        W.update();
        //ValueType psi_m = log(fabs(Psi.evaluate(W)));
        RealType psi_m = Psi.evaluateLog(W);
        RealType phase_m=Psi.getPhase();

        RealType ratDiff=std::exp(psi_m-psi_p)*std::cos(phase_m-phase_p) ;
        *fout << iat << " ratio " << aratio/ratDiff << " " << ratDiff << endl;
      }
  }

  void WaveFunctionTester::runRatioTest()
  {

    int nat = W.getTotalNum();
    ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
    ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);

    bool checkHam=(checkHamPbyP == "yes");
    Tau=0.025;
    MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
    while (it != it_end)
      {
        makeGaussRandom(deltaR);
        Walker_t::Buffer_t tbuffer;
        W.R = (**it).R+Tau*deltaR;
        (**it).R=W.R;
        W.update();
        RealType logpsi=Psi.registerData(W,tbuffer);
        RealType ene;
        if (checkHam)
          ene = H.registerData(W,tbuffer);
        else
          ene = H.evaluate(W);
        (*it)->DataSet=tbuffer;

        //RealType ene = H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
        H.saveProperty((*it)->getPropertyBase());
        ++it;

        app_log() << "  HamTest " << "  Total " <<  ene << endl;
        for (int i=0; i<H.sizeOfObservables(); i++)
          app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;
      }

    *fout << "  Update using drift " << endl;
    bool pbyp_mode=true;
    for (int iter=0; iter<4;++iter)
      {
        int iw=0;
        it=W.begin();
        while (it != it_end)
          {

            *fout << "\nStart Walker " << iw++ << endl;
            Walker_t& thisWalker(**it);
            W.loadWalker(thisWalker,pbyp_mode);
            Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);

            Psi.copyFromBuffer(W,w_buffer);
            H.copyFromBuffer(W,w_buffer);
//             Psi.evaluateLog(W);

            RealType eold(thisWalker.Properties(LOCALENERGY));
            RealType logpsi(thisWalker.Properties(LOGPSI));
            RealType emixed(eold), enew(eold);

            makeGaussRandom(deltaR);

            //mave a move
            RealType ratio_accum(1.0);
            for (int iat=0; iat<nat; iat++)
              {
                PosType dr(Tau*deltaR[iat]);

                PosType newpos(W.makeMove(iat,dr));

                //RealType ratio=Psi.ratio(W,iat,dGp,dLp);
                W.dG=0;
                W.dL=0;
                RealType ratio=Psi.ratio(W,iat,W.dG,W.dL);
                Gp = W.G + W.dG;
                //RealType enew = H.evaluatePbyP(W,iat);
                if (checkHam) enew = H.evaluatePbyP(W,iat);

                if (ratio > Random())
                  {
                    *fout << " Accepting a move for " << iat << endl;
                    *fout << " Energy after a move " << enew << endl;
                    W.G += W.dG;
                    W.L += W.dL;
                    W.acceptMove(iat);
                    Psi.acceptMove(W,iat);
                    if (checkHam) H.acceptMove(iat);
                    ratio_accum *= ratio;
                  }
                else
                  {
                    *fout << " Rejecting a move for " << iat << endl;
                    W.rejectMove(iat);
                    Psi.rejectMove(iat);
                    //H.rejectMove(iat);
                  }
              }

            *fout << " Energy after pbyp = " << H.getLocalEnergy() << endl;
            RealType newlogpsi_up = Psi.evaluateLog(W,w_buffer);
            W.saveWalker(thisWalker);
            RealType ene_up;
            if (checkHam)
              ene_up= H.evaluate(W,w_buffer);
            else
              ene_up = H.evaluate(W);

            Gp=W.G;
            Lp=W.L;
            W.R=thisWalker.R;
            W.update();
            RealType newlogpsi=Psi.updateBuffer(W,w_buffer,false);
            RealType ene = H.evaluate(W);
            thisWalker.resetProperty(newlogpsi,Psi.getPhase(),ene);
            //thisWalker.resetProperty(std::log(psi),Psi.getPhase(),ene);

            *fout << iter << "  Energy by update = "<< ene_up << " " << ene << " "  << ene_up-ene << endl;
            *fout << iter << " Ratio " << ratio_accum*ratio_accum
                 << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " "
                 << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << endl
                 << " new log(psi) updated " << newlogpsi_up
                 << " new log(psi) calculated " << newlogpsi
                 << " old log(psi) " << logpsi << endl;

            *fout << " Gradients " << endl;
            for (int iat=0; iat<nat; iat++)
              *fout << W.G[iat]-Gp[iat] << W.G[iat] << endl; //W.G[iat] << G[iat] << endl;
            *fout << " Laplacians " << endl;
            for (int iat=0; iat<nat; iat++)
              *fout << W.L[iat]-Lp[iat] << " " << W.L[iat] << endl;

            ++it;
          }
      }

    *fout << "  Update without drift : for VMC useDrift=\"no\"" << endl;
    for (int iter=0; iter<4;++iter)
      {
        it=W.begin();
        int iw=0;
        while (it != it_end)
          {

            *fout << "\nStart Walker " << iw++ << endl;
            Walker_t& thisWalker(**it);
            W.loadWalker(thisWalker,pbyp_mode);
            Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
            //Psi.updateBuffer(W,w_buffer,true);
            Psi.copyFromBuffer(W,w_buffer);
            RealType eold(thisWalker.Properties(LOCALENERGY));
            RealType logpsi(thisWalker.Properties(LOGPSI));
            RealType emixed(eold), enew(eold);

            //mave a move
            RealType ratio_accum(1.0);

            for (int substep=0; substep<3; ++substep)
              {
                makeGaussRandom(deltaR);

                for (int iat=0; iat<nat; iat++)
                  {
                    PosType dr(Tau*deltaR[iat]);

                    PosType newpos(W.makeMove(iat,dr));

                    RealType ratio=Psi.ratio(W,iat);
                    RealType prob = ratio*ratio;
                    if (prob > Random())
                      {
                        *fout << " Accepting a move for " << iat << endl;
                        W.acceptMove(iat);
                        Psi.acceptMove(W,iat);
                        ratio_accum *= ratio;
                      }
                    else
                      {
                        *fout << " Rejecting a move for " << iat << endl;
                        W.rejectMove(iat);
                        Psi.rejectMove(iat);
                      }
                  }

                RealType logpsi_up = Psi.updateBuffer(W,w_buffer,false);
                W.saveWalker(thisWalker);

                RealType ene = H.evaluate(W);
                thisWalker.resetProperty(logpsi_up,Psi.getPhase(),ene);
              }

            Gp=W.G;
            Lp=W.L;

            W.update();
            RealType newlogpsi=Psi.evaluateLog(W);
            *fout << iter << " Ratio " << ratio_accum*ratio_accum
                 << " | " << std::exp(2.0*(newlogpsi-logpsi)) << " "
                 << ratio_accum*ratio_accum/std::exp(2.0*(newlogpsi-logpsi)) << endl
                 << " new log(psi) " << newlogpsi
                 << " old log(psi) " << logpsi << endl;

            *fout << " Gradients " << endl;
            for (int iat=0; iat<nat; iat++)
              {
                *fout << W.G[iat]-Gp[iat] << W.G[iat] << endl; //W.G[iat] << G[iat] << endl;
              }
            *fout << " Laplacians " << endl;
            for (int iat=0; iat<nat; iat++)
              {
                *fout << W.L[iat]-Lp[iat] << " " << W.L[iat] << endl;
              }
            ++it;
          }
      }

    //for(it=W.begin();it != it_end; ++it)
    //{
    //  Walker_t& thisWalker(**it);
    //  Walker_t::Buffer_t& w_buffer((*it)->DataSet);
    //  w_buffer.rewind();
    //  W.updateBuffer(**it,w_buffer);
    //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
    //}


  }

  void WaveFunctionTester::runRatioTest2()
  {

    int nat = W.getTotalNum();
    ParticleSet::ParticleGradient_t Gp(nat), dGp(nat);
    ParticleSet::ParticleLaplacian_t Lp(nat), dLp(nat);

    Tau=0.025;
    MCWalkerConfiguration::iterator it(W.begin()), it_end(W.end());
    for (; it != it_end; ++it)
      {
        makeGaussRandom(deltaR);
        Walker_t::Buffer_t tbuffer;
        (**it).R  +=  Tau*deltaR;
        W.loadWalker(**it,true);
        RealType logpsi=Psi.registerData(W,tbuffer);
        RealType ene = H.evaluate(W);
        (*it)->DataSet=tbuffer;

        //RealType ene = H.evaluate(W);
        (*it)->resetProperty(logpsi,Psi.getPhase(),ene,0.0,0.0,1.0);
        H.saveProperty((*it)->getPropertyBase());

        app_log() << "  HamTest " << "  Total " <<  ene << endl;
        for (int i=0; i<H.sizeOfObservables(); i++)
          app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;
      }

    for (int iter=0; iter<20;++iter)
      {
        int iw=0;
        it=W.begin();
        //while(it != it_end)
        for (; it != it_end; ++it)
          {
            *fout << "\nStart Walker " << iw++ << endl;
            Walker_t& thisWalker(**it);
            W.loadWalker(thisWalker,true);
            Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
            Psi.copyFromBuffer(W,w_buffer);

            RealType eold(thisWalker.Properties(LOCALENERGY));
            RealType logpsi(thisWalker.Properties(LOGPSI));
            RealType emixed(eold), enew(eold);
            Psi.evaluateLog(W);
            ParticleSet::ParticleGradient_t realGrad(W.G);

            makeGaussRandom(deltaR);

            //mave a move
            RealType ratio_accum(1.0);
            for (int iat=0; iat<nat; iat++)
              {
                GradType grad_now=Psi.evalGrad(W,iat), grad_new;
                for(int sds=0;sds<3;sds++) *fout<< realGrad[iat][sds]-grad_now[sds]<<" ";
                
                PosType dr(Tau*deltaR[iat]);
                PosType newpos(W.makeMove(iat,dr));

                RealType ratio2 = Psi.ratioGrad(W,iat,grad_new);
                W.rejectMove(iat);
                Psi.rejectMove(iat);

                newpos=W.makeMove(iat,dr);
                RealType ratio1 = Psi.ratio(W,iat);
                //Psi.rejectMove(iat);
                W.rejectMove(iat);
                *fout << "  ratio1 = " << ratio1 << " ration2 = " << ratio2 << endl;
              }
              
          }
      }

    //for(it=W.begin();it != it_end; ++it)
    //{
    //  Walker_t& thisWalker(**it);
    //  Walker_t::Buffer_t& w_buffer((*it)->DataSet);
    //  w_buffer.rewind();
    //  W.updateBuffer(**it,w_buffer);
    //  RealType logpsi=Psi.updateBuffer(W,w_buffer,true);
    //}


  }


  void WaveFunctionTester::runGradSourceTest()
  {
    ParticleSetPool::PoolType::iterator p;
    for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
      app_log() << "ParticelSet = " << p->first << endl;

    // Find source ParticleSet
    ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
    app_log() << pit->first << endl;
    // if(pit == PtclPool.getPool().end())
    //   APP_ABORT("Unknown source \"" + sourceName + "\" WaveFunctionTester.");

    ParticleSet& source = *((*pit).second);

    IndexType nskipped = 0;
    RealType sig2Enloc=0, sig2Drift=0;
    RealType delta = 0.00001;
    RealType delta2 = 2*delta;
    ValueType c1 = 1.0/delta/2.0;
    ValueType c2 = 1.0/delta/delta;

    int nat = W.getTotalNum();

    ParticleSet::ParticlePos_t deltaR(nat);
    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

    //copy the properties of the working walker
    Properties = awalker->Properties;

    //sample a new walker configuration and copy to ParticleSet::R
    //makeGaussRandom(deltaR);

    W.R = awalker->R;

    //W.R += deltaR;

    W.update();
    //ValueType psi = Psi.evaluate(W);
    ValueType logpsi = Psi.evaluateLog(W);
    RealType eloc=H.evaluate(W);

    app_log() << "  HamTest " << "  Total " <<  eloc << endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

    //RealType psi = Psi.evaluateLog(W);
    ParticleSet::ParticleGradient_t G(nat), G1(nat);
    ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
    G = W.G;
    L = W.L;

    for (int isrc=0; isrc < 1/*source.getTotalNum()*/; isrc++)
      {
        TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad;
        TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad;
        TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad_FD;
        TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_FD;
        for (int dim=0; dim<OHMMS_DIM; dim++)
          {
            grad_grad[dim].resize(nat);
            lapl_grad[dim].resize(nat);
            grad_grad_FD[dim].resize(nat);
            lapl_grad_FD[dim].resize(nat);
          }
        Psi.evaluateLog(W);
        GradType grad_log = Psi.evalGradSource(W, source, isrc, grad_grad, lapl_grad);
        ValueType log = Psi.evaluateLog(W);
        //grad_log = Psi.evalGradSource (W, source, isrc);

        for (int iat=0; iat<nat; iat++)
          {
            PosType r0 = W.R[iat];
            GradType gFD[OHMMS_DIM];
            GradType lapFD = ValueType();
            for (int eldim=0; eldim<3; eldim++)
              {
                W.R[iat][eldim] = r0[eldim]+delta;
                W.update();
                ValueType log_p = Psi.evaluateLog(W);
                GradType gradlogpsi_p =  Psi.evalGradSource(W, source, isrc);
                W.R[iat][eldim] = r0[eldim]-delta;
                W.update();
                ValueType log_m = Psi.evaluateLog(W);
                GradType gradlogpsi_m = Psi.evalGradSource(W, source, isrc);
                lapFD    += gradlogpsi_m + gradlogpsi_p;
                gFD[eldim] = gradlogpsi_p - gradlogpsi_m;
                W.R[iat] = r0;
                W.update();
                //Psi.evaluateLog(W);
              }
            for (int iondim=0; iondim<OHMMS_DIM; iondim++)
              {
                for (int eldim=0; eldim<OHMMS_DIM; eldim++)
                  grad_grad_FD[iondim][iat][eldim] = c1*gFD[eldim][iondim];
                lapl_grad_FD[iondim][iat] = c2*(lapFD[iondim]-6.0*grad_log[iondim]);
              }
          }
        for (int dimsrc=0; dimsrc<OHMMS_DIM; dimsrc++)
          {
            for (int iat=0; iat<nat; iat++)
              {
                *fout << "For particle #" << iat << " at " << W.R[iat] << endl;
                *fout << "Gradient      = " << setw(12) << grad_grad[dimsrc][iat] << endl
                     << "  Finite diff = " << setw(12) << grad_grad_FD[dimsrc][iat] << endl
                     << "  Error       = " << setw(12)
                     <<  grad_grad_FD[dimsrc][iat] - grad_grad[dimsrc][iat] << endl << endl;
                *fout << "Laplacian     = " << setw(12) << lapl_grad[dimsrc][iat] << endl
                     << "  Finite diff = " << setw(12) << lapl_grad_FD[dimsrc][iat] << endl
                     << "  Error       = " << setw(12)
                     << lapl_grad_FD[dimsrc][iat] - lapl_grad[dimsrc][iat] << endl << endl;
              }
          }
      }
  }


  void WaveFunctionTester::runZeroVarianceTest()
  {
    ParticleSetPool::PoolType::iterator p;
    for (p=PtclPool.getPool().begin(); p != PtclPool.getPool().end(); p++)
      app_log() << "ParticelSet = " << p->first << endl;

    // Find source ParticleSet
    ParticleSetPool::PoolType::iterator pit(PtclPool.getPool().find(sourceName));
    app_log() << pit->first << endl;
    // if(pit == PtclPool.getPool().end())
    //   APP_ABORT("Unknown source \"" + sourceName + "\" WaveFunctionTester.");

    ParticleSet& source = *((*pit).second);

    int nat = W.getTotalNum();

    ParticleSet::ParticlePos_t deltaR(nat);
    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

    //copy the properties of the working walker
    Properties = awalker->Properties;

    //sample a new walker configuration and copy to ParticleSet::R
    //makeGaussRandom(deltaR);

    W.R = awalker->R;

    //W.R += deltaR;

    W.update();
    //ValueType psi = Psi.evaluate(W);
    ValueType logpsi = Psi.evaluateLog(W);
    RealType eloc=H.evaluate(W);

    //RealType psi = Psi.evaluateLog(W);
    ParticleSet::ParticleGradient_t G(nat), G1(nat);
    ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
    G = W.G;
    L = W.L;

    PosType r1(5.0, 2.62, 2.55);
    W.R[1] = PosType(4.313, 5.989, 4.699);
    W.R[2] = PosType(5.813, 4.321, 4.893);
    W.R[3] = PosType(4.002, 5.502, 5.381);
    W.R[4] = PosType(5.901, 5.121, 5.311);
    W.R[5] = PosType(5.808, 4.083, 5.021);
    W.R[6] = PosType(4.750, 5.810, 4.732);
    W.R[7] = PosType(4.690, 5.901, 4.989);
    for (int i=1; i<8; i++)
      W.R[i] -= PosType(2.5, 2.5, 2.5);

    char fname[32];
    sprintf(fname,"ZVtest.%03d.dat",OHMMS::Controller->rank());
    FILE *fzout = fopen(fname, "w");

    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad;
    TinyVector<ParticleSet::ParticleGradient_t, OHMMS_DIM> grad_grad_FD;
    TinyVector<ParticleSet::ParticleLaplacian_t,OHMMS_DIM> lapl_grad_FD;
    for (int dim=0; dim<OHMMS_DIM; dim++)
      {
        grad_grad[dim].resize(nat);
        lapl_grad[dim].resize(nat);
        grad_grad_FD[dim].resize(nat);
        lapl_grad_FD[dim].resize(nat);
      }

    for (r1[0]=0.0; r1[0]<5.0; r1[0]+=1.0e-4)
      {
        W.R[0] = r1;
        fprintf(fzout, "%1.8e %1.8e %1.8e ", r1[0], r1[1], r1[2]);
        ValueType log = Psi.evaluateLog(W);
        ValueType psi = std::cos(Psi.getPhase())*std::exp(log);//*W.PropertyList[SIGN];
        double E = H.evaluate(W);
        //double KE = E - W.PropertyList[LOCALPOTENTIAL];
        double KE = -0.5*(Sum(W.L) + Dot(W.G,W.G));
#if defined(QMC_COMPLEX)
        fprintf(fzout, "%16.12e %16.12e %16.12e ", psi.real(), psi.imag(),KE);
#else
        fprintf(fzout, "%16.12e %16.12e ", psi, KE);
#endif
        for (int isrc=0; isrc < source.getTotalNum(); isrc++)
          {
            GradType grad_log = Psi.evalGradSource(W, source, isrc, grad_grad, lapl_grad);
            for (int dim=0; dim<OHMMS_DIM; dim++)
              {
                double ZV = 0.5*Sum(lapl_grad[dim]) + Dot(grad_grad[dim], W.G);
#if defined(QMC_COMPLEX)
                fprintf(fzout, "%16.12e %16.12e %16.12e ", ZV, grad_log[dim].real(), grad_log[dim].imag());
#else
                fprintf(fzout, "%16.12e %16.12e ", ZV, grad_log[dim]);
#endif
              }
          }
        fprintf(fzout, "\n");
      }
    fclose(fzout);
  }



  bool
  WaveFunctionTester::put(xmlNodePtr q)
  {
    myNode = q;
    return putQMCInfo(q);
  }

  void WaveFunctionTester::runDerivTest()
  {
    app_log()<<" Testing derivatives"<<endl;
    IndexType nskipped = 0;
    RealType sig2Enloc=0, sig2Drift=0;
    RealType delta = 0.00001;
    RealType delta2 = 2*delta;
    ValueType c1 = 1.0/delta/2.0;
    ValueType c2 = 1.0/delta/delta;

    int nat = W.getTotalNum();

    ParticleSet::ParticlePos_t deltaR(nat);
    MCWalkerConfiguration::PropertyContainer_t Properties;
    //pick the first walker
    MCWalkerConfiguration::Walker_t* awalker = *(W.begin());

    //copy the properties of the working walker
    Properties = awalker->Properties;

    //sample a new walker configuration and copy to ParticleSet::R
    //makeGaussRandom(deltaR);

    W.R = awalker->R;

    //W.R += deltaR;

    W.update();
    //ValueType psi = Psi.evaluate(W);
    ValueType logpsi = Psi.evaluateLog(W);
    RealType eloc=H.evaluate(W);

    app_log() << "  HamTest " << "  Total " <<  eloc << endl;
    for (int i=0; i<H.sizeOfObservables(); i++)
      app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;

    //RealType psi = Psi.evaluateLog(W);
    ParticleSet::ParticleGradient_t G(nat), G1(nat);
    ParticleSet::ParticleLaplacian_t L(nat), L1(nat);
    G = W.G;
    L = W.L;
    *fout<<"Gradients"<<endl;
    for (int iat=0;iat<W.R.size();iat++)
      {
        for (int i=0; i<3 ; i++) *fout<<W.G[iat][i]<<"  ";
        *fout<<endl;
      }

    opt_variables_type wfVars,wfvar_prime;
//build optimizables from the wavefunction
    wfVars.clear();
    Psi.checkInVariables(wfVars);

    wfVars.resetIndex();
    Psi.checkOutVariables(wfVars);
    
    wfvar_prime= wfVars;
    wfVars.print(*fout);
    int Nvars= wfVars.size();
    vector<RealType> Dsaved(Nvars);
    vector<RealType> HDsaved(Nvars);
    vector<RealType> PGradient(Nvars);
    vector<RealType> HGradient(Nvars);
    Psi.resetParameters(wfVars);
    logpsi = Psi.evaluateLog(W);
    eloc=H.evaluate(W);
    Psi.evaluateDerivatives(W, wfVars, Dsaved, HDsaved);
    RealType FiniteDiff = 1e-5;

    QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
    for (int i=0; i<Nvars ; i++)
      {
        for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
        wfvar_prime[i] = wfVars[i]+ FiniteDiff;
//     Psi.checkOutVariables(wfvar_prime);
        Psi.resetParameters(wfvar_prime);
        Psi.reset();
        W.update();
        W.G=0;
        W.L=0;
        RealType logpsiPlus = Psi.evaluateLog(W);
        H.evaluate(W);
        RealType elocPlus=H.getLocalEnergy()-H.getLocalPotential();
        wfvar_prime[i] = wfVars[i]- FiniteDiff;
//     Psi.checkOutVariables(wfvar_prime);
        Psi.resetParameters(wfvar_prime);
        Psi.reset();
        W.update();
        W.G=0;
        W.L=0;
        RealType logpsiMinus = Psi.evaluateLog(W);
        H.evaluate(W);
        RealType elocMinus = H.getLocalEnergy()-H.getLocalPotential();
        
        PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
        HGradient[i]= (elocPlus-elocMinus)*dh;
      }
    *fout<<endl<<"Deriv  Numeric Analytic"<<endl;
    for (int i=0; i<Nvars ; i++) *fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i]) <<endl;
    *fout<<endl<<"Hderiv  Numeric Analytic"<<endl;
    for (int i=0; i<Nvars ; i++) *fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i]) <<endl;

  }



void WaveFunctionTester::runDerivCloneTest()
{
  app_log()<<" Testing derivatives clone"<<endl;
  RandomGenerator_t* Rng1= new RandomGenerator_t();
  RandomGenerator_t* Rng2= new RandomGenerator_t();
  (*Rng1) = (*Rng2);
  MCWalkerConfiguration* w_clone = new MCWalkerConfiguration(W);
  TrialWaveFunction *psi_clone = Psi.makeClone(*w_clone);
  QMCHamiltonian *h_clone = H.makeClone(*w_clone,*psi_clone);
  h_clone->setRandomGenerator(Rng2);
  H.setRandomGenerator(Rng1);
  h_clone->setPrimary(true);
  
  
  int nat = W.getTotalNum();
  
  ParticleSet::ParticlePos_t deltaR(nat);
  //pick the first walker
  MCWalkerConfiguration::Walker_t* awalker = *(W.begin());
//   MCWalkerConfiguration::Walker_t* bwalker = *(w_clone->begin());
//   bwalker->R = awalker->R;
  
  
  W.R = awalker->R;
  W.update();

  w_clone->R=awalker->R;
  w_clone->update(); 

  opt_variables_type wfVars;
  //build optimizables from the wavefunction
//   wfVars.clear();
  Psi.checkInVariables(wfVars);
  wfVars.resetIndex();
  Psi.checkOutVariables(wfVars);
  wfVars.print(*fout);
  int Nvars= wfVars.size();
  
  
  opt_variables_type wfvar_prime;
//   wfvar_prime.insertFrom(wfVars);
//   wfvar_prime.clear();
  psi_clone->checkInVariables(wfvar_prime);
  wfvar_prime.resetIndex();
  for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
  psi_clone->checkOutVariables(wfvar_prime);
  wfvar_prime.print(*fout);
  
  psi_clone->resetParameters(wfvar_prime);
  Psi.resetParameters(wfVars);
  
  
  vector<RealType> Dsaved(Nvars,0), og_Dsaved(Nvars,0);
  vector<RealType> HDsaved(Nvars,0), og_HDsaved(Nvars,0);
  vector<RealType> PGradient(Nvars,0), og_PGradient(Nvars,0);
  vector<RealType> HGradient(Nvars,0), og_HGradient(Nvars,0);
  
  
  ValueType logpsi2 = psi_clone->evaluateLog(*w_clone);
  RealType eloc2  = h_clone->evaluate(*w_clone); 
  psi_clone->evaluateDerivatives(*w_clone, wfvar_prime, Dsaved, HDsaved);
  
  ValueType logpsi1 = Psi.evaluateLog(W);
  RealType eloc1  = H.evaluate(W);
  Psi.evaluateDerivatives(W, wfVars, og_Dsaved, og_HDsaved);
  
  app_log() << "log (original) = " << logpsi1 << " energy = " << eloc1 << endl;
  for (int i=0; i<H.sizeOfObservables(); i++)
    app_log() << "  HamTest " << H.getObservableName(i) << " " << H.getObservable(i) << endl;
  
  app_log() << "log (clone)    = " << logpsi2 << " energy = " << eloc2 << endl;
  for (int i=0; i<h_clone->sizeOfObservables(); i++)
    app_log() << "  HamTest " << h_clone->getObservableName(i) << " " << h_clone->getObservable(i) << endl;
  
  
//   app_log()<<" Saved quantities:"<<endl;
//   for(int i=0;i<Nvars;i++) app_log()<<Dsaved[i]<<" "<<og_Dsaved[i]<<"         "<<HDsaved[i]<<" "<<og_HDsaved[i]<<endl;
  
  RealType FiniteDiff = 1e-6;
  QMCTraits::RealType dh=1.0/(2.0*FiniteDiff);
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
    
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
    psi_clone->resetParameters(wfvar_prime);
    psi_clone->reset();
    w_clone->update();
    w_clone->G=0;
    w_clone->L=0;
    RealType logpsiPlus = psi_clone->evaluateLog(*w_clone);
    h_clone->evaluate(*w_clone);
    RealType elocPlus=h_clone->getLocalEnergy()-h_clone->getLocalPotential();
    
    wfvar_prime[i] = wfVars[i]- FiniteDiff;
    psi_clone->resetParameters(wfvar_prime);
    psi_clone->reset();
    w_clone->update();
    w_clone->G=0;
    w_clone->L=0;
    RealType logpsiMinus = psi_clone->evaluateLog(*w_clone);
    h_clone->evaluate(*w_clone);
    RealType elocMinus = h_clone->getLocalEnergy()-h_clone->getLocalPotential();
    
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }
  *fout<<"CLONE"<<endl;
  *fout<<endl<<"   Deriv  Numeric Analytic"<<endl;
  for (int i=0; i<Nvars ; i++) *fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i])/PGradient[i] <<endl;
  *fout<<endl<<"   Hderiv  Numeric Analytic"<<endl;
  for (int i=0; i<Nvars ; i++) *fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i])/HGradient[i] <<endl;
  
  
  for (int i=0; i<Nvars ; i++)
  {
    for (int j=0; j<Nvars; j++) wfvar_prime[j]=wfVars[j];
    
    wfvar_prime[i] = wfVars[i]+ FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiPlus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocPlus=H.getLocalEnergy()-H.getLocalPotential();
    
    wfvar_prime[i] = wfVars[i]- FiniteDiff;
    Psi.resetParameters(wfvar_prime);
    Psi.reset();
    W.update();
    W.G=0;
    W.L=0;
    RealType logpsiMinus = Psi.evaluateLog(W);
    H.evaluate(W);
    RealType elocMinus = H.getLocalEnergy()-H.getLocalPotential();
    
    PGradient[i]= (logpsiPlus-logpsiMinus)*dh;
    HGradient[i]= (elocPlus-elocMinus)*dh;
  }
  *fout<<"ORIGINAL"<<endl;
  *fout<<endl<<"   Deriv  Numeric Analytic"<<endl;
  for (int i=0; i<Nvars ; i++) *fout<<i<<"  "<<PGradient[i]<<"  "<<Dsaved[i] <<"  " <<(PGradient[i]-Dsaved[i])/PGradient[i] <<endl;
  *fout<<endl<<"   Hderiv  Numeric Analytic"<<endl;
  for (int i=0; i<Nvars ; i++) *fout<<i <<"  "<<HGradient[i]<<"  "<<HDsaved[i] <<"  " <<(HGradient[i]-HDsaved[i])/HGradient[i] <<endl;
  
}
    void WaveFunctionTester::runwftricks()
    {
      vector<OrbitalBase*>& Orbitals=Psi.getOrbitals();
      app_log()<<" Total of "<<Orbitals.size()<<" orbitals."<<endl;
      int SDindex(0);
      for (int i=0;i<Orbitals.size();i++)
        if ("SlaterDet"==Orbitals[i]->OrbitalName) SDindex=i;
      
      SPOSetBasePtr Phi= dynamic_cast<SlaterDet *>(Orbitals[SDindex])->getPhi();
      int NumOrbitals=Phi->getBasisSetSize();
      app_log()<<"Basis set size: "<<NumOrbitals<<endl;
      
      vector<int> SPONumbers(0,0);
      vector<int> irrepRotations(0,0);
      vector<int> Grid(0,0);
      
      xmlNodePtr kids=myNode->children;
      
      string doProj("yes");
      string doRotate("yes");
      string sClass("C2V");
      
      ParameterSet aAttrib;
      aAttrib.add(doProj,"projection","string");
      aAttrib.add(doRotate,"rotate","string");
      aAttrib.put(myNode);
     
      while(kids != NULL) 
      {
        string cname((const char*)(kids->name));
        if(cname == "orbitals") 
        {
          putContent(SPONumbers,kids);
        }
        else if(cname == "representations") 
        {
          putContent(irrepRotations,kids);
        }
        else if(cname=="grid")
          putContent(Grid,kids);
        
      kids=kids->next;
      }
      
      ParticleSet::ParticlePos_t R_cart(1);
      R_cart.setUnit(PosUnit::CartesianUnit);
      ParticleSet::ParticlePos_t R_unit(1);
      R_unit.setUnit(PosUnit::LatticeUnit);
      
//       app_log()<<" My crystals basis set is:"<<endl;
//       vector<vector<RealType> > BasisMatrix(3, vector<RealType>(3,0.0));
//       
//       for (int i=0;i<3;i++)
//       {
//         R_unit[0][0]=0;
//         R_unit[0][1]=0;
//         R_unit[0][2]=0;
//         R_unit[0][i]=1;
//         W.convert2Cart(R_unit,R_cart);
//         app_log()<<"basis_"<<i<<":  ("<<R_cart[0][0]<<", "<<R_cart[0][1]<<", "<<R_cart[0][2]<<")"<<endl;
//         for (int j=0;j<3;j++) BasisMatrix[j][i]=R_cart[0][j];
//       }

      int Nrotated(SPONumbers.size());
      app_log()<<" Projected orbitals: ";
      for(int i=0;i<Nrotated;i++) app_log()<< SPONumbers[i] <<" ";
      app_log()<<endl;
      //indexing trick
//       for(int i=0;i<Nrotated;i++) SPONumbers[i]-=1;
      
      SymmetryBuilder SO;
      SO.put(myNode);
      SymmetryGroup symOp(*SO.getSymmetryGroup());
      
//       SO.changeBasis(InverseBasisMatrix);
      
      OrbitalSetTraits<ValueType>::ValueVector_t values;
      values.resize(NumOrbitals);
      
      RealType overG0(1.0/Grid[0]);
      RealType overG1(1.0/Grid[1]);
      RealType overG2(1.0/Grid[2]);
      RealType overNpoints=  overG0*overG1*overG2;

      vector<RealType> NormPhi(Nrotated, 0.0);
      
      int totsymops = symOp.getSymmetriesSize();
      Matrix<RealType> SymmetryOrbitalValues;
      SymmetryOrbitalValues.resize(Nrotated,totsymops);
       
      int ctabledim = symOp.getClassesSize();
      Matrix<double> projs(Nrotated,ctabledim);
      Matrix<double> orthoProjs(Nrotated,Nrotated);
      
      vector<RealType> brokenSymmetryCharacter(totsymops);
      for(int k=0;k<Nrotated;k++) for(int l=0;l<totsymops;l++) 
        brokenSymmetryCharacter[l] += irrepRotations[k]*symOp.getsymmetryCharacter(l,irrepRotations[k]-1);
//       app_log()<<"bsc: ";
//       for(int l=0;l<totsymops;l++) app_log()<<brokenSymmetryCharacter[l]<<" ";
//       app_log()<<endl;
//       for(int l=0;l<totsymops;l++) brokenSymmetryCharacter[l]+=0.5;
      
      if ((doProj=="yes")||(doRotate=="yes"))
      {
        OrbitalSetTraits<ValueType>::ValueVector_t identityValues(values.size());
        //Loop over grid
      for(int i=0;i<Grid[0];i++)
        for(int j=0;j<Grid[1];j++)
          for(int k=0;k<Grid[2];k++)
          {
            //Loop over symmetry classes and small group operators
            for(int l=0;l<totsymops;l++)
            {
                R_unit[0][0]=overG0*RealType(i);// R_cart[0][0]=0;
                R_unit[0][1]=overG1*RealType(j);// R_cart[0][1]=0;
                R_unit[0][2]=overG2*RealType(k);// R_cart[0][2]=0;
                
//                 for(int a=0; a<3; a++) for(int b=0;b<3;b++) R_cart[0][a]+=BasisMatrix[a][b]*R_unit[0][b];
                W.convert2Cart(R_unit,R_cart);

                symOp.TransformSinglePosition(R_cart,l);
                W.R[0]=R_cart[0];
                values=0.0;
                //evaluate orbitals
                
//                 Phi->evaluate(W,0,values);
                Psi.evaluateLog(W);
                for(int n=0;n<NumOrbitals;n++) values[n] = Phi->t_logpsi(0,n);
                if (l==0){
                  identityValues=values;
                  #if defined(QMC_COMPLEX)
                  for(int n=0;n<Nrotated;n++) NormPhi[n] += totsymops*real(values[SPONumbers[n]]*values[SPONumbers[n]]);
                  #else
                  for(int n=0;n<Nrotated;n++) NormPhi[n] += totsymops*(values[SPONumbers[n]]*values[SPONumbers[n]]);
                  #endif
                }
                
                //now we have phi evaluated at the rotated/inverted/whichever coordinates
                for(int n=0;n<Nrotated;n++) 
                {
                  int N=SPONumbers[n];
                  #if defined(QMC_COMPLEX)
                  RealType phi2 = real(values[N]*identityValues[N]);
                  #else
                  RealType phi2 = (values[N]*identityValues[N]);
                  #endif
                  SymmetryOrbitalValues(n,l) += phi2;
                }
                
                for(int n=0;n<Nrotated;n++) for(int p=0;p<Nrotated;p++)
                {
                  int N=SPONumbers[n];
                  int P=SPONumbers[p];
                  #if defined(QMC_COMPLEX)
                  orthoProjs(n,p) += 0.5*real(identityValues[N]*values[P]+identityValues[P]*values[N])*brokenSymmetryCharacter[l];
                  #else
                  orthoProjs(n,p) +=0.5*(identityValues[N]*values[P]+identityValues[P]*values[N])*brokenSymmetryCharacter[l];
                  #endif
                }
          }
      }
      for(int n=0;n<Nrotated;n++) for(int l=0;l<totsymops;l++) 
        SymmetryOrbitalValues(n,l)/= NormPhi[n];
      
      for(int n=0;n<Nrotated;n++) for(int l=0;l<Nrotated;l++) 
        orthoProjs(n,l) /= std::sqrt(NormPhi[n]*NormPhi[l]);
      
//       if (true){
//         app_log()<<endl;
//         for(int n=0;n<Nrotated;n++) {
//           for(int l=0;l<totsymops;l++) app_log()<<SymmetryOrbitalValues(n,l)<<" ";
//           app_log()<<endl;
//         }
//       app_log()<<endl;
//       }
      

      
      for(int n=0;n<Nrotated;n++)
      {
        if (false) app_log()<<" orbital #"<<SPONumbers[n]<<endl;
        for(int i=0;i<ctabledim;i++)
        {
          double proj(0);
          for(int j=0;j<totsymops;j++) proj+=symOp.getsymmetryCharacter(j,i)*SymmetryOrbitalValues(n,j);
          if (false) app_log()<<"  Rep "<<i<< ": "<<proj;
          projs(n,i)=proj<1e-4?0:proj;
        }
        if (false) app_log()<<endl;
      }
      
      if (true){
        app_log()<<"Printing Projection Matrix"<<endl;
        for(int n=0;n<Nrotated;n++) {
          for(int l=0;l<ctabledim;l++) app_log()<<projs(n,l)<<" ";
          app_log()<<endl;
        }
      app_log()<<endl;
      }
      if (true){
        app_log()<<"Printing Coefficient Matrix"<<endl;
        for(int n=0;n<Nrotated;n++) {
          for(int l=0;l<ctabledim;l++) app_log()<<std::sqrt(projs(n,l))<<" ";
          app_log()<<endl;
        }
      app_log()<<endl;
      }
      
      if (doRotate=="yes")
      {

//         app_log()<<"Printing Broken Symmetry Projection Matrix"<<endl;
//           for(int n=0;n<Nrotated;n++) {
//             for(int l=0;l<Nrotated;l++) app_log()<<orthoProjs(n,l)<<" ";
//             app_log()<<endl;
//           } 
        
        char JOBU('A');
        char JOBVT('A');
        int vdim=Nrotated;
        Vector<double> Sigma(vdim);
        Matrix<double> U(vdim,vdim);
        Matrix<double> VT(vdim,vdim);
        int lwork=8*Nrotated;
        vector<double> work(lwork,0);
        int info(0);
        
        dgesvd(&JOBU, &JOBVT, &vdim, &vdim,
            orthoProjs.data(), &vdim, Sigma.data(), U.data(),
            &vdim, VT.data(), &vdim, &(work[0]),
            &lwork, &info);
        
        app_log()<<"Printing Rotation Matrix"<<endl;
          for(int n=0;n<vdim;n++) {
            for(int l=0;l<vdim;l++) app_log()<<VT(l,n)<<" ";
            app_log()<<endl;
          }  
          app_log()<<endl<<"Printing Eigenvalues"<<endl;
          for(int n=0;n<vdim;n++) app_log()<<Sigma[n]<<" ";
          app_log()<<endl;  
      }      
      }

    }
    
}

/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$
 ***************************************************************************/
