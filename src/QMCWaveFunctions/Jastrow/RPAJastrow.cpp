// -*- C++ -*-
/** @file RPAJastrow.cpp
 * @brief Definitions of RPAJastrow
 */

#include "QMCWaveFunctions/Jastrow/RPAJastrow.h"
#include "QMCWaveFunctions/OrbitalBuilderBase.h"
#include "QMCWaveFunctions/Jastrow/TwoBodyJastrowOrbital.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/LRTwoBodyJastrow.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRRPAHandlerTemp.h"
#include "Utilities/IteratorUtility.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/ParameterSet.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

RPAJastrow::RPAJastrow(ParticleSet& target, bool is_manager)
  :IsManager(is_manager), targetPtcl(target)
{
  Optimizable=true;
  OrbitalName="RPAJastrow";
}

RPAJastrow::~RPAJastrow()
{
  delete_iter(Psi.begin(), Psi.end());
  delete myHandler;
}

bool RPAJastrow::put(xmlNodePtr cur)
{
  ReportEngine PRE("RPAJastrow","put");
  xmlNodePtr myNode=xmlCopyNode(cur,1);
  //capture attribute jastrow/@name
  MyName="RPA_Jee";
  string useL="yes";
  string useS="yes";
  rpafunc="rpa";
  OhmmsAttributeSet a;
  a.add(MyName,"name");
  a.add(useL,"longrange");
  a.add(useS,"shortrange");
  a.add(rpafunc,"function");
  a.put(cur);
  Rs=-1.0;
  Kc=-1.0;
  string ID_Rs="RPA_rs";
  ParameterSet params;
  params.add(Rs,"rs","double");
  params.add(Kc,"kc","double");
  params.put(cur);
  buildOrbital(MyName, useL, useS, rpafunc, Rs, Kc);
//     app_log() <<endl<<"   LongRangeForm is "<<rpafunc<<endl;
//
//     DropLongRange = (useL == "no");
//     DropShortRange = (useS=="no");
//
//     app_log() << "    Rs can be optimized using ID=" << ID_Rs << endl;
//     RealType tlen = std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/ static_cast<RealType>(targetPtcl.getTotalNum()) ,1.0/3.0);
//
//     if(Rs<0) {
//       if(targetPtcl.Lattice.SuperCellEnum) {
//         Rs=tlen;
//       } else {
//         cout<<"  Error finding rs. Is this an open system?!"<<endl;
//         Rs=100.0;
//       }
//     }
//
//     //Add Rs to optimizable list
//     myVars.insert(ID_Rs,Rs,true);
//
//     int indx = targetPtcl.SK->KLists.ksq.size()-1;
//     double Kc_max=std::pow(targetPtcl.SK->KLists.ksq[indx],0.5);
//
//     if(Kc<0){
//       Kc = 2.0*  std::pow(2.25*M_PI,1.0/3.0)/tlen ;
//     }
//
//     if(Kc>Kc_max){
//       Kc=Kc_max;
//       app_log() << "    Kc set too high. Resetting to the maximum value"<<endl;
//     }
//
//     app_log() << "    RPAJastrowBuilder::addTwoBodyPart Rs = " << Rs <<  "  Kc= " << Kc << endl;
//
//     if (rpafunc=="Yukawa" || rpafunc=="breakup"){
//       myHandler= new LRHandlerTemp<YukawaBreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
//     } else if (rpafunc=="RPA"){
//       myHandler= new LRRPAHandlerTemp<RPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
//     } else if (rpafunc=="dYukawa"){
//       myHandler= new LRHandlerTemp<DerivYukawaBreakup<RealType>,LPQHIBasis >(targetPtcl,Kc);
//     } else if (rpafunc=="dRPA"){
//       myHandler= new LRRPAHandlerTemp<DerivRPABreakup<RealType>,LPQHIBasis >(targetPtcl,Kc);
//     }
//
//
//     myHandler->Breakup(targetPtcl,Rs);
//
//     app_log() << "  Maximum K shell " << myHandler->MaxKshell << endl;
//     app_log() << "  Number of k vectors " << myHandler->Fk.size() << endl;
//
//     if(!DropLongRange) makeLongRange();
//     if(!DropShortRange) makeShortRange();
  return true;
}

void RPAJastrow::buildOrbital(const string& name, const string& UL
                              , const string& US, const string& RF, RealType R, RealType K)
{
  string ID_Rs="RPA_rs";
  MyName = name;
  string useL=UL;
  string useS=US;
  rpafunc=RF;
  Rs=R;
  Kc=K;
  app_log() <<endl<<"   LongRangeForm is "<<rpafunc<<endl;
  DropLongRange = (useL == "no");
  DropShortRange = (useS=="no");
  app_log() << "    Rs can be optimized using ID=" << ID_Rs << endl;
  RealType tlen = std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/ static_cast<RealType>(targetPtcl.getTotalNum()) ,1.0/3.0);
  if(Rs<0)
  {
    if(targetPtcl.Lattice.SuperCellEnum)
    {
      Rs=tlen;
    }
    else
    {
      cout<<"  Error finding rs. Is this an open system?!"<<endl;
      Rs=100.0;
    }
  }
  //Add Rs to optimizable list
//     myVars.insert(ID_Rs,Rs,true);
  int indx = targetPtcl.SK->KLists.ksq.size()-1;
  double Kc_max=std::pow(targetPtcl.SK->KLists.ksq[indx],0.5);
  if(Kc<0)
  {
    Kc = 2.0*  std::pow(2.25*M_PI,1.0/3.0)/tlen ;
  }
  if(Kc>Kc_max)
  {
    Kc=Kc_max;
    app_log() << "    Kc set too high. Resetting to the maximum value"<<endl;
  }
  app_log() << "    RPAJastrowBuilder::addTwoBodyPart Rs = " << Rs <<  "  Kc= " << Kc << endl;
  if (rpafunc=="yukawa" || rpafunc=="breakup")
  {
    myHandler= new LRHandlerTemp<YukawaBreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
  }
  else
    if (rpafunc=="rpa")
    {
      myHandler= new LRRPAHandlerTemp<RPABreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
    }
    else
      if (rpafunc=="dyukawa")
      {
        myHandler= new LRHandlerTemp<DerivYukawaBreakup<RealType>,LPQHIBasis >(targetPtcl,Kc);
      }
      else
        if (rpafunc=="drpa")
        {
          myHandler= new LRRPAHandlerTemp<DerivRPABreakup<RealType>,LPQHIBasis >(targetPtcl,Kc);
        }
  myHandler->Breakup(targetPtcl,Rs);
  app_log() << "  Maximum K shell " << myHandler->MaxKshell << endl;
  app_log() << "  Number of k vectors " << myHandler->Fk.size() << endl;
  if(!DropLongRange)
  {
    makeLongRange();
    app_log()<<"  Using LongRange part"<<endl;
  }
  if(!DropShortRange)
  {
    makeShortRange();
    app_log()<<"  Using ShortRange part"<<endl;
  }
}

void RPAJastrow::makeLongRange()
{
  LongRangeRPA = new LRTwoBodyJastrow(targetPtcl);
  LongRangeRPA->resetByHandler(myHandler);
  Psi.push_back(LongRangeRPA);
}

void RPAJastrow::makeShortRange()
{
//     app_log()<< "  Adding Short Range part of RPA function"<<endl;
  //short-range uses realHandler
  Rcut = myHandler->get_rc()-0.1;
  myGrid = new GridType;
  int npts=static_cast<int>(Rcut/0.01)+1;
  myGrid->set(0,Rcut,npts);
  //create the numerical functor
  nfunc = new FuncType;
  SRA = new ShortRangePartAdapter<RealType>(myHandler);
  SRA->setRmax(Rcut);
  nfunc->initialize(SRA, myGrid);
  //Do not write the table
  //static  int counter=0;
  //if(IsManager && counter==0)
  //{
  //  char fname[32];
  //  sprintf(fname,"%s.%d.dat",MyName.c_str(),counter++);
  //  ofstream fout(fname);
  //  for (int i = 0; i < myGrid->size(); i++) {
  //    RealType r=(*myGrid)(i);
  //    fout << r << "   " << nfunc->evaluate(r) << "   "
  //      << myHandler->evaluate(r,1.0/r) << " "
  //      << myHandler->evaluateLR(r) << endl;
  //  }
  //}
  TwoBodyJastrowOrbital<FuncType> *j2 = new TwoBodyJastrowOrbital<FuncType>(targetPtcl,IsManager);
  j2->addFunc(0,0,nfunc);
  ShortRangeRPA=j2;
  Psi.push_back(ShortRangeRPA);
}

void RPAJastrow::resetParameters(const opt_variables_type& active)
{
  /*
   int loc=myVars.where(0);
   if(loc>=0) {
     Rs=myVars[0]=active[loc];
     ///Insert breakup etc.
     myHandler->Breakup(targetPtcl,Rs);

     if(!DropLongRange){
       delete LongRangeRPA;
       makeLongRange();
     }
     if(!DropShortRange){
       delete ShortRangeRPA;
       makeShortRange();
     }
   };*/
}

void RPAJastrow::checkOutVariables(const opt_variables_type& active)
{
//     myVars.getIndex(active);
}

void RPAJastrow::checkInVariables(opt_variables_type& active)
{
//     active.insertFrom(myVars);
}

void RPAJastrow::reportStatus(ostream& os)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->reportStatus(os);
}

void RPAJastrow::resetTargetParticleSet(ParticleSet& P)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->resetTargetParticleSet(P);
}

RPAJastrow::RealType
RPAJastrow::evaluateLog(ParticleSet& P,
                        ParticleSet::ParticleGradient_t& G, ParticleSet::ParticleLaplacian_t& L)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->evaluateLog(P,G,L);
  return LogValue;
}

RPAJastrow::ValueType
RPAJastrow::ratio(ParticleSet& P, int iat,
                  ParticleSet::ParticleGradient_t& dG,
                  ParticleSet::ParticleLaplacian_t& dL)
{
  ValueType r(1.0);
  for(int i=0; i<Psi.size(); i++)
    r *= Psi[i]->ratio(P,iat,dG,dL);
  return r;
}

RPAJastrow::ValueType
RPAJastrow::ratio(ParticleSet& P, int iat)
{
  ValueType r(1.0);
  for(int i=0; i<Psi.size(); i++)
    r *= Psi[i]->ratio(P,iat);
  return r;
}

//RPAJastrow::ValueType
//  RPAJastrow::logRatio(ParticleSet& P, int iat,
//      ParticleSet::ParticleGradient_t& dG,
//      ParticleSet::ParticleLaplacian_t& dL) {
//    ValueType r(0.0);
//    for(int i=0; i<Psi.size(); i++)
//      r += Psi[i]->logRatio(P,iat,dG,dL);
//    return r;
//  }

void RPAJastrow::acceptMove(ParticleSet& P, int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->acceptMove(P,iat);
}

void RPAJastrow::restore(int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->restore(iat);
}

void RPAJastrow::update(ParticleSet& P,
                        ParticleSet::ParticleGradient_t& dG,
                        ParticleSet::ParticleLaplacian_t& dL,
                        int iat)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->update(P,dG,dL,iat);
}

RPAJastrow::RealType
RPAJastrow::registerData(ParticleSet& P, BufferType& buf)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->registerData(P,buf);
  return LogValue;
}

RPAJastrow::RealType
RPAJastrow::updateBuffer(ParticleSet& P, BufferType& buf, bool fromscratch)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->updateBuffer(P,buf,fromscratch);
  return LogValue;
}

void
RPAJastrow::copyFromBuffer(ParticleSet& P, BufferType& buf)
{
  for(int i=0; i<Psi.size(); i++)
    Psi[i]->copyFromBuffer(P,buf);
}

RPAJastrow::RealType
RPAJastrow::evaluateLog(ParticleSet& P,BufferType& buf)
{
  LogValue=0.0;
  for(int i=0; i<Psi.size(); i++)
    LogValue += Psi[i]->evaluateLog(P,buf);
  return LogValue;
}

OrbitalBase* RPAJastrow::makeClone(ParticleSet& tpq) const
{
  HandlerType* tempHandler;
  if (rpafunc=="yukawa" || rpafunc=="breakup")
  {
    tempHandler= new LRHandlerTemp<YukawaBreakup<RealType>,LPQHIBasis>(dynamic_cast<const LRHandlerTemp<YukawaBreakup<RealType>,LPQHIBasis>& > (*myHandler),tpq);
  }
  else
    if (rpafunc=="rpa")
    {
      tempHandler= new LRRPAHandlerTemp<RPABreakup<RealType>,LPQHIBasis>(dynamic_cast<const LRRPAHandlerTemp<RPABreakup<RealType>,LPQHIBasis>& > (*myHandler),tpq);
    }
    else
      if (rpafunc=="dyukawa")
      {
        tempHandler= new LRHandlerTemp<DerivYukawaBreakup<RealType>,LPQHIBasis >(dynamic_cast<const LRHandlerTemp<DerivYukawaBreakup<RealType>,LPQHIBasis>& > (*myHandler),tpq);
      }
      else
        if (rpafunc=="drpa")
        {
          tempHandler= new LRRPAHandlerTemp<DerivRPABreakup<RealType>,LPQHIBasis >(dynamic_cast<const LRRPAHandlerTemp<DerivRPABreakup<RealType>,LPQHIBasis>& > (*myHandler),tpq);
        }
  RPAJastrow* myClone = new RPAJastrow(tpq,IsManager);
  myClone->setHandler(tempHandler);
  if(!DropLongRange)
    myClone->makeLongRange();
  if(!DropShortRange)
    myClone->makeShortRange();
  return myClone;
}
};
