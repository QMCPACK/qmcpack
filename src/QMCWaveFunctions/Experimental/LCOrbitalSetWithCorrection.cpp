//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/LCOrbitalSet.h"
#include "ParticleIO/XMLParticleIO.h"
#include "OhmmsData/Libxml2Doc.h"


namespace qmcplusplus
{

template<class BS>
BS* LCOrbitalSetWithCorrection<BS,false>::extractHighYLM( std::vector<bool> &rmv)
{
  int nUniqCenters = myBasisSet->LOBasisSet.size(), cnt=0;
  int cntEta;
  std::vector<int> nSorbs;
  std::vector<std::vector<bool> > sOrbs, tags;
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
      if( sOrbs[i][k] )
        cnt++;
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
      rmv[cnt++] = (((myBasisSet->LOBasis[i]->RnlID[myBasisSet->LOBasis[i]->NL[k]])[1]==0)&&(corrCenter[i]));
  }
  BS* Eta=0;
  /* not used any more, restore if you want to test and propagate changes due to corrCenter
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
       std::vector<int> LM_eta, NL_eta;
       std::vector<typename BS::ThisCOT_t::RadialOrbital_t*> Rnl_eta;
  //     std::vector<COT::RadialOrbital_t*> Rnl_eta;
       std::vector<QuantumNumberType> RnlID_eta;
       for(int i=0; i<nUniqCenters; i++)
       {
          int bss = myBasisSet->LOBasisSet[i]->BasisSetSize;
          int nrf = myBasisSet->LOBasisSet[i]->Rnl.size();
          LM_eta.resize(bss-nSorbs[i]);
          NL_eta.resize(bss-nSorbs[i]);
          Rnl_eta.resize(nrf-nSorbs[i]);
          RnlID_eta.resize(nrf-nSorbs[i]);

          std::map<int,int> old2eta;

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
  */
  return Eta;
}

// phi --> S component of the basis set, eta--> the rest
template<class BS>
void LCOrbitalSetWithCorrection<BS,false>::createLCOSets(int centr, LCOrbitalSet<BS,false>* Phi, LCOrbitalSet<BS,false>* Eta)
{
  int numCtr = myBasisSet->NumCenters;
  int numUniq = myBasisSet->LOBasisSet.size();
  std::vector<bool> rmv;
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
        rmv[cnt++] = (((cur->RnlID[cur->NL[k]])[1]==0)&&(corrCenter[i]));;
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
    {
      auto &cref(*(Eta->C));
      for(int k=0; k<nOrbs; k++)
        cref(k,i)=0.0; //Eta->C(k,i) = 0.0;
    }
    else
    {
      auto &cref(*(Phi->C));
      for(int k=0; k<nOrbs; k++)
        cref(k,i)=0.0; //Phi->C(k,i) = 0.0;
    }
  }
}

template<class BS>
LCOrbitalSet<BS,false>* LCOrbitalSetWithCorrection<BS,false>::clone2LCOrbitalSet()
{
  BS* newBS = (BS*) myBasisSet->makeClone();
  LCOrbitalSet<BS,false>* newSPO = new LCOrbitalSet<BS,false>(newBS,ReportLevel);
  newSPO->IsCloned=true;
  newSPO->setOrbitalSetSize(OrbitalSetSize);
  newSPO->setIdentity(Identity);
  newSPO->C = C;
  newSPO->Occ.resize(Occ.size());
  newSPO->Occ = Occ;
  newSPO->className = "LCOrbitalSet";
  return newSPO;
}

template<class BS>
bool LCOrbitalSetWithCorrection<BS,false>::readCuspInfo(Matrix<TinyVector<RealType,9> > &info)
{
  bool success=true;
  std::string cname;
  int ncenter = info.rows();
  int nOrbs = info.cols();
  app_log() <<"Reading cusp info from : " <<cuspInfoFile << std::endl;
  Libxml2Document adoc;
  if(!adoc.parse(cuspInfoFile))
  {
    app_log()<<"Could not find precomputed cusp data for spo set: " <<objectName << std::endl;
    app_log() <<"Recalculating data.\n";
    return false;
  }
  xmlNodePtr head=adoc.getRoot();
  head=head->children;
  xmlNodePtr cur=NULL,ctr;
  while (head != NULL)
  {
    getNodeName(cname,head);
    if(cname == "sposet")
    {
      std::string name;
      OhmmsAttributeSet spoAttrib;
      spoAttrib.add (name, "name");
      spoAttrib.put(head);
      if(name == objectName)
      {
        cur=head;
        break;
      }
    }
    head=head->next;
  }
  if(cur==NULL)
  {
    app_log()<<"Could not find precomputed cusp data for spo set: " <<objectName << std::endl;
    app_log() <<"Recalculating data.\n";
    return false;
  }
  else
  {
    app_log() <<"Found precomputed cusp data for spo set: " <<objectName << std::endl;
  }
  cur=cur->children;
  while (cur != NULL)
  {
    getNodeName(cname,cur);
    if(cname == "center")
    {
      int num=-1;
      OhmmsAttributeSet Attrib;
      Attrib.add (num, "num");
      Attrib.put(cur);
      if(num < 0 || num >= ncenter )
      {
        APP_ABORT("Error with cusp info xml block. incorrect center number. \n");
      }
      ctr=cur->children;
      while (ctr != NULL)
      {
        getNodeName(cname,ctr);
        if(cname == "orbital")
        {
          int orb=-1;
          OhmmsAttributeSet orbAttrib;
          RealType a1,a2,a3,a4,a5,a6,a7,a8,a9;
          orbAttrib.add (orb,"num");
          orbAttrib.add (a1, "redo");
          orbAttrib.add (a2, "C");
          orbAttrib.add (a3, "sg");
          orbAttrib.add (a4, "rc");
          orbAttrib.add (a5, "a1");
          orbAttrib.add (a6, "a2");
          orbAttrib.add (a7, "a3");
          orbAttrib.add (a8, "a4");
          orbAttrib.add (a9, "a5");
          orbAttrib.put(ctr);
          if(orb < OrbitalSetSize)
          {
            info(num,orb)[0] = a1;
            info(num,orb)[1] = a2;
            info(num,orb)[2] = a3;
            info(num,orb)[3] = a4;
            info(num,orb)[4] = a5;
            info(num,orb)[5] = a6;
            info(num,orb)[6] = a7;
            info(num,orb)[7] = a8;
            info(num,orb)[8] = a9;
          }
          /*
          std::cout <<" Found: num,orb:" <<num <<"  " <<orb << std::endl
              <<info(num,orb)[0] <<"\n"
              <<info(num,orb)[1] <<"\n"
              <<info(num,orb)[2] <<"\n"
              <<info(num,orb)[3] <<"\n"
              <<info(num,orb)[4] <<"\n"
              <<info(num,orb)[5] <<"\n"
              <<info(num,orb)[6] <<"\n"
              <<info(num,orb)[7] <<"\n"
              <<info(num,orb)[8] <<"\n"
               << std::endl;cout.flush();
          */
        }
        ctr=ctr->next;
      }
    }
    cur=cur->next;
  }
  return success;
}

template<class BS>
bool LCOrbitalSetWithCorrection<BS,false>::transformSPOSet()
{
  app_log()<<" Transforming Single Particle Orbital Set with cusp correcting algorithm. \n";
  app_log() <<"cuspFile: " <<cuspInfoFile << std::endl;
  SpeciesSet& tspecies(sourcePtcl->getSpeciesSet());
  int iz = tspecies.addAttribute("charge");
  Z.resize(sourcePtcl->getTotalNum());
  for(int iat=0; iat<Z.size(); iat++)
  {
    Z[iat] = tspecies(iz,sourcePtcl->GroupID[iat]);
  }
  int numAttrib = tspecies.numAttributes();
  int cctag = tspecies.addAttribute("cuspCorr");
  corrCenter.resize(Z.size(),"true");
  if(cctag < numAttrib)
    // parameter is not new
  {
    for(int iat=0; iat<Z.size(); iat++)
    {
      if( tspecies(cctag,sourcePtcl->GroupID[iat]) != 0)
      {
        corrCenter[iat] = false;
        app_log() <<"Not using cusp correction algorithm in atoms of type: " <<tspecies.speciesName[sourcePtcl->GroupID[iat]] << std::endl;
      }
    }
  }
  if(Rcut < 0)
    Rcut=0.1;
  //int indexRc=-1;
  originalSPOSet = clone2LCOrbitalSet();
  corrBasisSet = new CorrectingOrbitalBasisSet<COT>(myBasisSet->myTable,*sourcePtcl,myBasisSet->NumTargets);
  TinyVector<RealType,3> dr=0, storeR = targetPtcl->R[0];
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
  dummyLO1->setOrbitalSetSize(OrbitalSetSize);
  dummyLO1->BasisSetSize = BasisSetSize;
  dummyLO1->setIdentity(Identity);
  (*dummyLO1->C) = *C;
  dummyLO1->Occ.resize(Occ.size());
  dummyLO1->Occ = Occ;
  dummyLO2 = new LCOrbitalSet<BS,false>(myBasisSet,ReportLevel);
  dummyLO2->setOrbitalSetSize(OrbitalSetSize);
  dummyLO2->BasisSetSize = BasisSetSize;
  dummyLO2->setIdentity(Identity);
  (*dummyLO2->C) = *C;
  dummyLO2->Occ.resize(Occ.size());
  dummyLO2->Occ = Occ;
  Matrix<TinyVector<RealType,9> > info;
  info.resize(numCentr,OrbitalSetSize);
  info=0;
  bool readCuspCoeff=false;
  if(cuspInfoFile != "")
    readCuspCoeff = readCuspInfo(info);
  targetPtcl->R[0]=0;
  CuspCorr<BS> myCorr(0.2,500,targetPtcl,sourcePtcl,true);
  Vector<RealType> rad_orb, xgrid;
  std::vector<bool> rmv;
  rmv.resize(myBasisSet->BasisSetSize);
  xgrid.resize(mygrid->size());
  rad_orb.resize(mygrid->size());
  for(int ig=0; ig<mygrid->size(); ++ig)
  {
    xgrid[ig] = mygrid->r(ig);
    //if(xgrid[ig]>2*Rcut && indexRc < 0) indexRc=ig;
  }
  xmlNodePtr spo = xmlNewNode(NULL,(const xmlChar*)"sposet");
  xmlNewProp(spo,(const xmlChar*)"name",(const xmlChar*)objectName.c_str());
  for(int i=0; i<numCentr; i++ )
  {
    app_log()<<"Transforming orbitals of center " << i << std::endl;
    (*dummyLO1->C) = *C;
    (*dummyLO2->C) = *C;
    createLCOSets(i,dummyLO1,dummyLO2);
    COT *myCOT = new COT(0,true);  // assuming gaussian package
    myCOT->Grids.resize(1);
//       myCOT->Grids[0] = myBasisSet->LOBasis[0]->Grids[0];
    myCOT->LM.resize(OrbitalSetSize);
    myCOT->NL.resize(OrbitalSetSize);
    myCOT->Rnl.resize(OrbitalSetSize);
    myCOT->RnlID.resize(OrbitalSetSize);
    xmlNodePtr ctr = xmlNewNode(NULL,(const xmlChar*)"center");
    std::ostringstream num;
    num <<i;
    xmlNewProp(ctr,(const xmlChar*)"num",(const xmlChar*)num.str().c_str());
    for(int k=0; k<OrbitalSetSize; k++ )
    {
      xmlNodePtr orb = xmlNewNode(NULL,(const xmlChar*)"orbital");
      const std::string fileprefix("newOrbs."+objectName+".C"+std::to_string(i)+".MO"+std::to_string(k));
      std::ostringstream num0,redo,C,sg,rc,a1,a2,a3,a4,a5;
      num0<<k;
      xmlNewProp(orb,(const xmlChar*)"num",(const xmlChar*)num0.str().c_str());
      bool corrO = false;
      if(corrCenter[i])
      {
        auto& cref(*(dummyLO1->C));
        for(int ip=0; ip<cref.cols(); ip++)
        {
          if(std::abs(cref(k,ip)) > 1e-8)
          {
            corrO = true;
            break;
          }
        }
      }
      if(corrO)
      {
        if(readCuspCoeff)
        {
          RealType redo = info(i,k)[0];
          if(redo < -1)
            // no correction to this orbital
          {
            myCorr.fillRadFunWithPhi(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb);
          }
          else
            if(redo > 10)
              // recompute with rc loop
            {
              myCorr.executeWithRCLoop(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb,fileprefix,Rcut,info(i,k).data());
            }
            else
              if(redo > 1)
                // no rc loop, read rc from file
              {
                RealType rc = info(i,k)[3];
                myCorr.execute(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb,fileprefix,rc,info(i,k).data());
              }
              else
                // read from file
              {
                myCorr.fillRadFunWithPhiBar(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb,info(i,k).data());
              }
        }
        else
        {
          myCorr.executeWithRCLoop(k,i,Z[i],dummyLO1,dummyLO2,xgrid,rad_orb,fileprefix,Rcut,info(i,k).data());
          info(i,k)[0]=0;
        }
      }
      else
      {
        for(int ip=0; ip<rad_orb.size(); ip++)
          rad_orb[ip]=0.0;
      }
      C.setf(std::ios::scientific, std::ios::floatfield);
      C.precision(14);
      C<<info(i,k)[1];
      sg.setf(std::ios::scientific, std::ios::floatfield);
      sg.precision(14);
      sg<<info(i,k)[2];
      rc.setf(std::ios::scientific, std::ios::floatfield);
      rc.precision(14);
      rc<<info(i,k)[3];
      a1.setf(std::ios::scientific, std::ios::floatfield);
      a1.precision(14);
      a1<<info(i,k)[4];
      a2.setf(std::ios::scientific, std::ios::floatfield);
      a2.precision(14);
      a2<<info(i,k)[5];
      a3.setf(std::ios::scientific, std::ios::floatfield);
      a3.precision(14);
      a3<<info(i,k)[6];
      a4.setf(std::ios::scientific, std::ios::floatfield);
      a4.precision(14);
      a4<<info(i,k)[7];
      a5.setf(std::ios::scientific, std::ios::floatfield);
      a5.precision(14);
      a5<<info(i,k)[8];
      xmlNewProp(orb,(const xmlChar*)"C",(const xmlChar*)C.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"sg",(const xmlChar*)sg.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"rc",(const xmlChar*)rc.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"a1",(const xmlChar*)a1.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"a2",(const xmlChar*)a2.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"a3",(const xmlChar*)a3.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"a4",(const xmlChar*)a4.str().c_str());
      xmlNewProp(orb,(const xmlChar*)"a5",(const xmlChar*)a5.str().c_str());
      xmlAddChild(ctr,orb);
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
    xmlAddChild(spo,ctr);
  }
  corrBasisSet->setBasisSetSize(-1);
  BS *dum3 = extractHighYLM(rmv);
  //dum3->setBasisSetSize(-1);
  int norb=0, cnt=0;
  for(int i=0; i<rmv.size(); i++)
    if(!rmv[i])
      norb++;
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
    if(rmv[i])
    {
      app_log()<<"removing basis element i: " <<i << std::endl;
      for(int k=0; k<OrbitalSetSize; k++)
        (*C)(k,i) = 0.0;
    }
  this->checkObject();
  if(!readCuspCoeff)
  {
    xmlDocPtr doc = xmlNewDoc((const xmlChar*)"1.0");
    xmlNodePtr cuspRoot = xmlNewNode(NULL, BAD_CAST "qmcsystem");
    xmlAddChild(cuspRoot,spo);
    xmlDocSetRootElement(doc, cuspRoot);
    std::string fname = objectName+".cuspInfo.xml";
    app_log() <<"Saving resulting cusp Info xml block to: " <<fname << std::endl;
    xmlSaveFormatFile(fname.c_str(),doc,1);
    xmlFreeDoc(doc);
  }
  /*
       std::ofstream out("test.txt");
       std::ofstream out9("test9.txt");
       ValueVector_t psi_2,d2y_2;
       GradVector_t dy_2;
       psi_2.resize(OrbitalSetSize);
       d2y_2.resize(OrbitalSetSize);
       dy_2.resize(OrbitalSetSize);
       std::vector<RealType> v1(OrbitalSetSize,0.0),v2(OrbitalSetSize,0.0);
       std::vector<RealType> v3(OrbitalSetSize,0.0),v4(OrbitalSetSize,0.0);
       std::vector<RealType> v5(OrbitalSetSize,0.0),v6(OrbitalSetSize,0.0);
       std::vector<RealType> v7(OrbitalSetSize,0.0),v8(OrbitalSetSize,0.0);
       std::vector<RealType> v9(OrbitalSetSize,0.0),v10(OrbitalSetSize,0.0);
       std::vector<RealType> v11(OrbitalSetSize,0.0),v12(OrbitalSetSize,0.0);
       (targetPtcl->R[0]) = 0.0;
       for(int ig=0; ig<xgrid.size(); ++ig)
       {
         (targetPtcl->R[0])[0] = xgrid[ig];
         TinyVector<RealType,3> ddr2=targetPtcl->makeMove(0,dr);
         this->evaluate(*targetPtcl,0,psi);
         originalSPOSet->evaluate(*targetPtcl,0,psi_2);
         for(int i=0; i<OrbitalSetSize; i++)
            v1[i] += std::abs(psi[i]-psi_2[i]);

         this->evaluate(*targetPtcl,0,psi,dy,d2y);
         originalSPOSet->evaluate(*targetPtcl,0,psi_2,dy_2,d2y_2);
         out<<xgrid[ig] <<"  ";
         out9<<xgrid[ig] <<"  ";
         for(int i=0; i<OrbitalSetSize; i++) {
            out9<<psi[i] <<"  " <<psi_2[i] <<"  ";
            out<<dy[i][0] <<"  " <<dy_2[i][0]
               <<d2y[i] <<"  " <<d2y_2[i] <<"  ";
            v2[i] += std::abs(psi[i]-psi_2[i]);
            v3[i] += std::abs(dy[i][0]-dy_2[i][0]);
            v4[i] += std::abs(dy[i][1]-dy_2[i][1]);
            v5[i] += std::abs(dy[i][2]-dy_2[i][2]);
            v6[i] += std::abs(d2y[i]-d2y_2[i]);
         }
         out<< std::endl;
         out9<< std::endl;
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

       std::ofstream out1("test2.txt");
       out1<<xgrid[400] << std::endl << std::endl;
       out1<<corrBasisSet->d2Phi[0] << std::endl << std::endl;
       RealType tmp1=0.0,tmp2=0.0;
       for(int i=0; i<myBasisSet->Phi.size(); i++) {
          out1<<i <<"  " <<C(0,i) <<"  " <<myBasisSet->d2Phi[i] <<"  " <<originalSPOSet->myBasisSet->d2Phi[i] <<"  " <<myBasisSet->d2Phi[i]-originalSPOSet->myBasisSet->d2Phi[i] << std::endl;
          tmp1+=C(0,i)*myBasisSet->d2Phi[i];
          tmp2+=originalSPOSet->C(0,i)*originalSPOSet->myBasisSet->d2Phi[i];
       }
       out1<< std::endl << std::endl;
       out1<<tmp1 <<"   "  <<tmp2 << std::endl;
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
