#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "QMCWaveFunctions/LCOrbitalSet.h"

namespace qmcplusplus {

  template<class BS>
  BS* LCOrbitalSetWithCorrection<BS,false>::extractHighYLM( std::vector<bool> &rmv)
  {
     int nUniqCenters = myBasisSet->LOBasisSet.size(), cnt=0;
     int cntEta;
     vector<int> nSorbs;
     vector<vector<bool> > sOrbs, tags;

     rmv.resize(myBasisSet->BasisSetSize);
     nSorbs.resize(nUniqCenters);
     sOrbs.resize(nUniqCenters);
     tags.resize(nUniqCenters);
     for(int i=0; i<nUniqCenters; i++)
     {
        cnt=0;
        int bss = myBasisSet->LOBasisSet[i]->BasisSetSize;
        int nrf = myBasisSet->LOBasisSet[i]->Rnl.size();
        tags[i].resize(bss);
        sOrbs[i].resize(nrf);

        for(int k=0; k<nrf; k++)
        {
           sOrbs[i][k] = (myBasisSet->LOBasisSet[i]->RnlID[k])[1]==0;
           if( sOrbs[i][k] ) cnt++;
        }
        for(int k=0; k<bss; k++)
        {
           tags[i][k] = (myBasisSet->LOBasisSet[i]->RnlID[myBasisSet->LOBasisSet[i]->NL[k]])[1]==0;
        }
        nSorbs[i]=cnt;
     }

     cnt=0;
     for(int i=0; i<myBasisSet->CenterSys.getTotalNum(); i++)
     {
        int bss = myBasisSet->LOBasis[i]->BasisSetSize;
        for(int k=0; k<bss; k++)
           rmv[cnt++] = (myBasisSet->LOBasis[i]->RnlID[myBasisSet->LOBasis[i]->NL[k]])[1]==0;
     }

     BS* Eta = new BS(*myBasisSet);
     for(int i=0; i<nUniqCenters; i++)
     {
       Eta->LOBasisSet[i]=myBasisSet->LOBasisSet[i]->makeClone();
       for(int j=0; j<myBasisSet->CenterSys.getTotalNum(); ++j)
       {
         if(myBasisSet->CenterSys.GroupID[j]==i)
         {
           Eta->LOBasis[j]=Eta->LOBasisSet[i];
         }
       }
     }

// modify objects
     vector<int> LM_eta, NL_eta;
     vector<typename BS::ThisCOT_t::RadialOrbital_t*> Rnl_eta;
//     vector<COT::RadialOrbital_t*> Rnl_eta;
     vector<QuantumNumberType> RnlID_eta;
     for(int i=0; i<nUniqCenters; i++)
     {
        int bss = myBasisSet->LOBasisSet[i]->BasisSetSize;
        int nrf = myBasisSet->LOBasisSet[i]->Rnl.size();
        LM_eta.resize(bss-nSorbs[i]);
        NL_eta.resize(bss-nSorbs[i]);
        Rnl_eta.resize(nrf-nSorbs[i]);
        RnlID_eta.resize(nrf-nSorbs[i]);

        map<int,int> old2eta;

        cntEta=0;
        for(int k=0; k<nrf; k++)
        {
           if(sOrbs[i][k])
           {
             delete Eta->LOBasisSet[i]->Rnl[k]; // not used any more
           }
           else
           {
             old2eta[k] = cntEta;
             Rnl_eta[cntEta] = Eta->LOBasisSet[i]->Rnl[k];
             RnlID_eta[cntEta] = Eta->LOBasisSet[i]->RnlID[k];
             (RnlID_eta[cntEta])[0] = cntEta; // same as Rnl 
             (RnlID_eta[cntEta])[3] = 0; 
             if(i==0 && cntEta==0)
               Rnl_eta[cntEta]->setGridManager(true);
             else
               Rnl_eta[cntEta]->setGridManager(false);
             cntEta++;
           }
        }

        cntEta=0;
        for(int k=0; k<bss; k++)
        {
          if(!tags[i][k])
          {
            LM_eta[cntEta] = myBasisSet->LOBasisSet[i]->LM[k];
            NL_eta[cntEta] = old2eta[myBasisSet->LOBasisSet[i]->NL[k]];
            cntEta++;
          }
        }

        Eta->LOBasisSet[i]->LM = LM_eta;
        Eta->LOBasisSet[i]->NL = NL_eta;
        Eta->LOBasisSet[i]->Rnl = Rnl_eta;
        Eta->LOBasisSet[i]->RnlID = RnlID_eta;

        Eta->LOBasisSet[i]->setBasisSetSize(-1);
     }

// now set basis set sizes
     Eta->setBasisSetSize(-10);
     Eta->resize(Eta->NumTargets);

     return Eta;
  }

  // phi --> S component of the basis set, eta--> the rest
  template<class BS>
  void LCOrbitalSetWithCorrection<BS,false>::createLCOSets(int centr, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta) 
  {
     int numCtr = myBasisSet->NumCenters; 
     int numUniq = myBasisSet->LOBasisSet.size();
     vector<bool> rmv;

     int cnt=0;
     rmv.resize(myBasisSet->BasisSetSize);
     for(int i=0; i<myBasisSet->CenterSys.getTotalNum(); i++)
     {
       int bss = myBasisSet->LOBasis[i]->BasisSetSize;
       if(centr == i)
       {
// WARNING, assuming that COT has to be NGOrbital, otherwise problems
         COT* cur = (COT*) myBasisSet->LOBasis[i];
         for(int k=0; k<bss; k++)
           rmv[cnt++] = (cur->RnlID[cur->NL[k]])[1]==0;
       }
       else
       {
         cnt+=bss;
       }
     }
     int nOrbs = OrbitalSetSize, bss=BasisSetSize;
     for(int i=0; i<bss; i++)
     {
       if(rmv[i])
         for(int k=0; k<nOrbs; k++)
           Eta->C(k,i) = 0.0;
       else
         for(int k=0; k<nOrbs; k++)
           Phi->C(k,i) = 0.0;
     }
  }

  template<class BS>
  LCOrbitalSet<BS,false>* LCOrbitalSetWithCorrection<BS,false>::clone2LCOrbitalSet()
  {
 
     BS* newBS = (BS*) myBasisSet->makeClone();
     LCOrbitalSet<BS,false>* newSPO = new LCOrbitalSet<BS,false>(newBS,ReportLevel);

     newSPO->setOrbitalSetSize(OrbitalSetSize);
     newSPO->TotalOrbitalSize=TotalOrbitalSize;
     newSPO->setIdentity(Identity);  
     newSPO->C = C;
     newSPO->Occ.resize(Occ.size());
     newSPO->Occ = Occ;
     newSPO->className = "LCOrbitalSet";

     return newSPO;
  }

  template<class BS>  
  bool LCOrbitalSetWithCorrection<BS,false>::transformSPOSet()
  {
     app_log()<<" Transforming Single Particle Orbital Set with cusp correcting algorithm. \n";

      SpeciesSet& tspecies(sourcePtcl->getSpeciesSet());
      int iz = tspecies.addAttribute("charge");
      Z.resize(sourcePtcl->getTotalNum());
      for(int iat=0; iat<Z.size();iat++) {
        Z[iat] = tspecies(iz,sourcePtcl->GroupID[iat]);
      }

     if(Rcut < 0) Rcut=0.1;
     //int indexRc=-1;

     originalSPOSet = clone2LCOrbitalSet();
     corrBasisSet = new CorrectingOrbitalBasisSet<COT>(myBasisSet->myTable,*sourcePtcl,myBasisSet->NumTargets);

     TinyVector<double,3> dr=0, storeR = targetPtcl->R[0];
     int numCentr = myBasisSet->NumCenters;
     LCOrbitalSet<BS,false> *dummyLO1, *dummyLO2;
//       GridType_ *mygrid = myBasisSet->LOBasis[0]->Grids[0]; 
     ValueVector_t psi,d2y;
     GradVector_t dy;
     psi.resize(OrbitalSetSize);
     d2y.resize(OrbitalSetSize);
     dy.resize(OrbitalSetSize);

     GridType_ *mygrid = new LogGrid<RealType>;
// this should depend on smallest exponent of Gaussians in the original basis
// FIX FIX FIX: control from input file
     mygrid->set(0.000001,100.0,1001);

     dummyLO1 = new LCOrbitalSet<BS,false>(myBasisSet,ReportLevel);
     dummyLO1->OrbitalSetSize = OrbitalSetSize;
     dummyLO1->TotalOrbitalSize=TotalOrbitalSize;
     dummyLO1->BasisSetSize = BasisSetSize;
     dummyLO1->setIdentity(Identity);
     dummyLO1->C = C;
     dummyLO1->Occ.resize(Occ.size());
     dummyLO1->Occ = Occ;
     dummyLO2 = (LCOrbitalSet<BS,false>*) dummyLO1->makeClone();

     targetPtcl->R[0]=0;
     CuspCorr<BS> myCorr(Rcut,500,targetPtcl,sourcePtcl); 
     char buff1[10],buff2[10];
     Vector<double> rad_orb, xgrid;
     std::vector<bool> rmv;
     rmv.resize(myBasisSet->BasisSetSize);
     xgrid.resize(mygrid->size());
     rad_orb.resize(mygrid->size());
     for(int ig=0; ig<mygrid->size(); ++ig) 
     {
       xgrid[ig] = mygrid->r(ig);
       //if(xgrid[ig]>2*Rcut && indexRc < 0) indexRc=ig;
     }

     for(int i=0; i<numCentr; i++ )
     {
       app_log()<<"Transforming orbitals of center " <<i <<endl;
       dummyLO1->C = C;
       dummyLO2->C = C;
       createLCOSets(i,dummyLO1,dummyLO2);

       COT *myCOT = new COT(0,true);  // assuming gaussian package
       myCOT->Grids.resize(1);
//       myCOT->Grids[0] = myBasisSet->LOBasis[0]->Grids[0]; 
       myCOT->LM.resize(OrbitalSetSize);
       myCOT->NL.resize(OrbitalSetSize);
       myCOT->Rnl.resize(OrbitalSetSize);
       myCOT->RnlID.resize(OrbitalSetSize);
       std::sprintf(buff1,"%d",i);

       for(int k=0; k<OrbitalSetSize; k++ )
       {
          bool corrO = false;
          for(int ip=0; ip<dummyLO1->C.cols(); ip++) { 
            if(std::fabs(dummyLO1->C(k,ip)) > 1e-8) {
              corrO = true;
              break; 
            }
          }   
          if(corrO) {
// this should be Zion = charge of ion i
            std::sprintf(buff2,"%d",k);
            myCorr.execute(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb,"newOrbs.C"+string(buff1)+".MO"+string(buff2)); 
          }
          else {
            for(int ip=0; ip<rad_orb.size(); ip++)
              rad_orb[ip]=0.0; 
          } 
/*
          for(int ig=0; ig<xgrid.size(); ++ig)
          {
            (targetPtcl->R[0])[2] = xgrid[ig];
            TinyVector<RealType,3> ddr2=targetPtcl->makeMove(0,dr);
            dummyLO1->evaluate(*targetPtcl,0,psi);
// FIX FIX FIX: nasty now, do this correctly later
            rad_orb[ig] = psi[k]*std::sqrt(4.0*3.14159265358979);
          }
*/
          myCOT->LM[k] = 0;  // they are all s-type
          myCOT->NL[k] = k;
          myCOT->Rnl[k] = new NGOrbital(mygrid,rad_orb);
          RealType yprime_i = (rad_orb[1]-rad_orb[0])/((mygrid->r(1)-mygrid->r(0)));
          myCOT->Rnl[k]->spline(0,yprime_i,rad_orb.size()-1,0.0);
          if(k==0)
            myCOT->Rnl[k]->setGridManager(true);
          else
            myCOT->Rnl[k]->setGridManager(false);
          (myCOT->RnlID[k])[0] = k;
          (myCOT->RnlID[k])[1] = 0;
          (myCOT->RnlID[k])[2] = 0;
          (myCOT->RnlID[k])[3] = 0;
       }
       myCOT->setBasisSetSize(-1);
       corrBasisSet->add(i,myCOT);
     }
     corrBasisSet->setBasisSetSize(-1);

     BS *dum3 = extractHighYLM(rmv);
     dum3->setBasisSetSize(-1);
     int norb=0, cnt=0;
     for(int i=0; i<rmv.size(); i++) if(!rmv[i]) norb++;
     app_log()<<"Found " <<rmv.size()-norb <<" spherically symmetric basis functions and " <<norb <<" non symmetric ones. \n";
/*
     C.resize(OrbitalSetSize,norb);
     for(int k=0; k<OrbitalSetSize; k++)
     {
       cnt=0;
       for(int i=0; i<rmv.size(); i++)
          if(!rmv[i]) C(k,cnt++) = originalSPOSet->C(k,i);
     }
// FIX FIX FIX
// you can't delete myBasisSet because the other spin
// determinant might use a copy of it
//       delete myBasisSet;
     this->setBasisSet(dum3);
*/
// for now
     for(int i=0; i<rmv.size(); i++)
       if(rmv[i]) { 
 app_log()<<"removing basis element i: " <<i <<endl;
         for(int k=0; k<OrbitalSetSize; k++)
            C(k,i) = 0.0;
       }
     this->checkObject();

/*
     ofstream out("test.txt");
     ofstream out9("test9.txt");
     ValueVector_t psi_2,d2y_2;
     GradVector_t dy_2;
     psi_2.resize(OrbitalSetSize);
     d2y_2.resize(OrbitalSetSize);
     dy_2.resize(OrbitalSetSize);
     vector<double> v1(OrbitalSetSize,0.0),v2(OrbitalSetSize,0.0); 
     vector<double> v3(OrbitalSetSize,0.0),v4(OrbitalSetSize,0.0); 
     vector<double> v5(OrbitalSetSize,0.0),v6(OrbitalSetSize,0.0); 
     vector<double> v7(OrbitalSetSize,0.0),v8(OrbitalSetSize,0.0); 
     vector<double> v9(OrbitalSetSize,0.0),v10(OrbitalSetSize,0.0); 
     vector<double> v11(OrbitalSetSize,0.0),v12(OrbitalSetSize,0.0); 
     (targetPtcl->R[0]) = 0.0;
     for(int ig=0; ig<xgrid.size(); ++ig)
     {
       (targetPtcl->R[0])[0] = xgrid[ig];
       TinyVector<RealType,3> ddr2=targetPtcl->makeMove(0,dr);
       this->evaluate(*targetPtcl,0,psi);
       originalSPOSet->evaluate(*targetPtcl,0,psi_2);
       for(int i=0; i<OrbitalSetSize; i++)
          v1[i] += std::fabs(psi[i]-psi_2[i]);

       this->evaluate(*targetPtcl,0,psi,dy,d2y);
       originalSPOSet->evaluate(*targetPtcl,0,psi_2,dy_2,d2y_2);
       out<<xgrid[ig] <<"  ";
       out9<<xgrid[ig] <<"  ";
       for(int i=0; i<OrbitalSetSize; i++) {
          out9<<psi[i] <<"  " <<psi_2[i] <<"  ";
          out<<dy[i][0] <<"  " <<dy_2[i][0]
             <<d2y[i] <<"  " <<d2y_2[i] <<"  ";
          v2[i] += std::fabs(psi[i]-psi_2[i]);
          v3[i] += std::fabs(dy[i][0]-dy_2[i][0]);
          v4[i] += std::fabs(dy[i][1]-dy_2[i][1]);
          v5[i] += std::fabs(dy[i][2]-dy_2[i][2]);
          v6[i] += std::fabs(d2y[i]-d2y_2[i]);
       } 
       out<<endl;
       out9<<endl;
     }
     app_log()<<"Integrated differences after transformation: \n";
     for(int i=0; i<OrbitalSetSize; i++) {
       app_log()<<i <<"  "
                <<v1[i] <<"  "
                <<v2[i] <<"  "
                <<v3[i] <<"  "
                <<v4[i] <<"  "
                <<v5[i] <<"  "
                <<v6[i] <<"\n";
     }
     for(int i=0; i<OrbitalSetSize; i++) { 
       if(v1[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"1: evaluate ";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
       if(v2[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"2: evaluate \n";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
       if(v3[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"3: evaluate \n";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
       if(v4[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"4: evaluate \n";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
       if(v5[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"5: evaluate \n";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
       if(v6[i] > 1.0e-8) {
            app_log()<<"transformSPOSet: Orbitals do not agree. \n";
            app_log()<<"6: evaluate \n";
            APP_ABORT("transformSPOSet: Orbita ls do not agree");
       }
     }
     out.close();
     out9.close();

     (targetPtcl->R[0]) = 0;
     (targetPtcl->R[0])[0] = xgrid[400];
     TinyVector<RealType,3> ddr2=targetPtcl->makeMove(0,dr);
     this->evaluate(*targetPtcl,0,psi,dy,d2y);
     originalSPOSet->evaluate(*targetPtcl,0,psi_2,dy_2,d2y_2);

     ofstream out1("test2.txt");
     out1<<xgrid[400] <<endl <<endl;
     out1<<corrBasisSet->d2Phi[0] <<endl <<endl;
     double tmp1=0.0,tmp2=0.0;
     for(int i=0; i<myBasisSet->Phi.size(); i++) {
        out1<<i <<"  " <<C(0,i) <<"  " <<myBasisSet->d2Phi[i] <<"  " <<originalSPOSet->myBasisSet->d2Phi[i] <<"  " <<myBasisSet->d2Phi[i]-originalSPOSet->myBasisSet->d2Phi[i] <<endl;
        tmp1+=C(0,i)*myBasisSet->d2Phi[i];
        tmp2+=originalSPOSet->C(0,i)*originalSPOSet->myBasisSet->d2Phi[i];
     }
     out1<<endl <<endl;
     out1<<tmp1 <<"   "  <<tmp2 <<endl;
     out1.close();
*/

     targetPtcl->R[0]=storeR;
     delete dummyLO1;
     delete dummyLO2;
     app_log()<<"Done transforming SPOSet. \n";

//APP_ABORT("ABORTING");

     return true;
  }


}
