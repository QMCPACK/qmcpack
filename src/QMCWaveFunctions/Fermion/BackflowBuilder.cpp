//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "QMCWaveFunctions/Fermion/BackflowBuilder.h"
#include "Utilities/ProgressReportEngine.h"
#include "OhmmsData/AttributeSet.h"
#include "QMCWaveFunctions/TrialWaveFunction.h"
#include "QMCWaveFunctions/Fermion/BackflowTransformation.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee.h"
#include "QMCWaveFunctions/Fermion/Backflow_ee_kSpace.h"
#include "QMCWaveFunctions/Fermion/Backflow_eI.h"
#include "QMCWaveFunctions/Fermion/Backflow_eI_spin.h"
#include "QMCWaveFunctions/Fermion/GaussianFunctor.h"
#include "QMCWaveFunctions/Jastrow/BsplineFunctor.h"
#include "QMCWaveFunctions/Jastrow/LRBreakupUtilities.h"
#include "QMCWaveFunctions/Jastrow/SplineFunctors.h"
#include "LongRange/LRHandlerTemp.h"
#include "LongRange/LRRPABFeeHandlerTemp.h"
#include "Particle/ParticleSet.h"
#include "Configuration.h"
#include <map>
#include <cmath>
#include "OhmmsPETE/OhmmsArray.h"
#include "OhmmsData/ParameterSet.h"
#include "Numerics/LinearFit.h"

namespace qmcplusplus
{

BackflowBuilder::BackflowBuilder(ParticleSet& els, PtclPoolType& pool, TrialWaveFunction& psi)
  : OrbitalBuilderBase(els,psi), ptclPool(pool), BFTrans(0), cutOff(-1.0)
{
  ClassName="BackflowBuilder";
}

BackflowBuilder::~BackflowBuilder()
{
}

bool BackflowBuilder::put(xmlNodePtr cur)
{
  bool first=true;
  bool success=true;
  xmlNodePtr curRoot=cur;
  std::string cname;
  BFTrans = new BackflowTransformation(targetPtcl);
  cur = curRoot->children;
  while (cur != NULL)
  {
    getNodeName(cname,cur);
    if (cname == "transf" || cname == "transformation")
    {
      OhmmsAttributeSet spoAttrib;
      std::string source("none");
      std::string name("bf0");
      std::string type("none");
      spoAttrib.add (name, "name");
      spoAttrib.add (type, "type");
      spoAttrib.add (source, "source");
      spoAttrib.put(cur);
      BFTrans->sources[source] = BFTrans->names.size();
      BFTrans->names.push_back(name);
      if(type == "e-e")
      {
        addTwoBody(cur);
      }
      else
        if(type == "e-e-I")
        {
          APP_ABORT("e-e-I backflow is not implemented yet. \n");
        }
        else
          if(type == "e-I")
          {
            addOneBody(cur);
          }
          else
            if(type == "rpa")
            {
              app_log() <<"Adding RPA backflow functions. \n";
              addRPA(cur);
            }
            else
            {
              APP_ABORT("Unknown backflow type. \n");
            }
    }
    cur = cur->next;
  }
  BFTrans->cutOff = cutOff;
  return success;
}

void BackflowBuilder::addOneBody(xmlNodePtr cur)
{
  OhmmsAttributeSet spoAttrib;
  std::string source("none");
  std::string name("bf0");
  std::string type("none");
  std::string funct("Gaussian");
  std::string unique("no");
  std::string spin("no"); //add spin attribute, with default spin="no"
  spoAttrib.add (name, "name");
  spoAttrib.add (type, "type");
  spoAttrib.add (source, "source");
  spoAttrib.add (funct, "function");
  spoAttrib.add (unique, "unique");
  spoAttrib.add (spin, "spin");
  spoAttrib.put(cur);
  ParticleSet* ions=0;
  PtclPoolType::iterator pit(ptclPool.find(source));
  if(pit == ptclPool.end())
  {
    APP_ABORT("Missing backflow/@source.");
  }
  else
  {
    ions=(*pit).second;
  }
  app_log() <<"Adding electron-Ion backflow for source:"
            <<source <<" \n";
  BackflowFunctionBase *tbf=0;
  int nIons = ions->getTotalNum();
  SpeciesSet &sSet = ions->getSpeciesSet();
  SpeciesSet &tSet = targetPtcl.getSpeciesSet();
  int numSpecies = sSet.getTotalNum();
  if(spin=="yes")
  {
    if(funct!="Bspline")
      APP_ABORT("DON'T KNOW WHAT TO DO YET"); //will template this
    Backflow_eI_spin<BsplineFunctor<RealType> >*  tbf1
    =new Backflow_eI_spin<BsplineFunctor<RealType> >(*ions,targetPtcl);
    tbf1->numParams=0;
    xmlNodePtr cur1=cur->children;
    std::string cname;
    while (cur1 != NULL)
    {
      getNodeName(cname,cur1);
      if(cname =="correlation")
      {
        RealType my_cusp=0.0;
        std::string speciesA("0");
        std::string speciesB("e"); //assume electrons
        OhmmsAttributeSet anAttrib;
        anAttrib.add(speciesA, "elementType");
        anAttrib.add (speciesA, "speciesA");
        anAttrib.add(speciesB, "speciesB");
        anAttrib.add(my_cusp, "cusp");
        anAttrib.put(cur1);
        BsplineFunctor<RealType> *afunc = new BsplineFunctor<RealType>(my_cusp);
        afunc->elementType = speciesA;
        int ig = sSet.findSpecies (speciesA);
        afunc->cutoff_radius = ions->Lattice.WignerSeitzRadius;
        int jg=-1;
        if(speciesB.size())
          jg=tSet.findSpecies(speciesB);
        if (ig < numSpecies)
        {
          //ignore
          afunc->put (cur1);
          tbf1->addFunc (ig,afunc,jg);
          //WHAT IS THIS
          afunc->myVars.setParameterType(optimize::BACKFLOW_P);
          tbf1->myVars.insertFrom(afunc->myVars);
          tbf1->myVars.resetIndex();
          afunc->myVars.getIndex(tbf1->myVars);
          tbf1->numParams=tbf1->myVars.size();
          int offset_a = tbf1->myVars.getIndex(afunc->myVars.name(0));
          if(jg<0)
          {
            for(int jjg=0; jjg<targetPtcl.groups(); ++jjg)
              tbf1->offsetPrms(ig,jjg)=offset_a;
          }
          else
            tbf1->offsetPrms(ig,jg)=offset_a;
        }
      }
      cur1 = cur1->next;
    }
    //synch parameters
    tbf1->resetParameters(tbf1->myVars);
    tbf1->derivs.resize(tbf1->numParams);
    tbf=tbf1;
  }
  else //keep original implementation
  {
    std::vector<xmlNodePtr> funs;
    std::vector<int> ion2functor(nIons,-1);
    std::vector<RealType> cusps;
    xmlNodePtr curRoot=cur;
    std::string cname;
    cur = curRoot->children;
    while (cur != NULL)
    {
      getNodeName(cname,cur);
      if (cname == "correlation")
      {
        RealType my_cusp=0.0;
        std::string elementType("none");
        OhmmsAttributeSet anAttrib;
        anAttrib.add (elementType, "elementType");
        anAttrib.add (my_cusp, "cusp");
        anAttrib.put(cur);
        funs.push_back(cur);
        cusps.push_back(my_cusp);
        if(unique == "yes") // look for <index> block, and map based on that
        {
          xmlNodePtr kids=cur;
          std::string aname;
          kids = cur->children;
          while (kids != NULL)
          {
            getNodeName(aname,kids);
            if (aname == "index")
            {
              std::vector<int> pos;
              putContent(pos, kids);
              for(int i=0; i<pos.size(); i++)
              {
                app_log() << "Adding backflow transformation of type " << funs.size()-1 << " for atom " << pos[i] << ".\n";
                ion2functor[pos[i]]=funs.size()-1;
              }
            }
            kids = kids->next;
          }
        }
        else
          // map based on elementType
        {
          int ig = sSet.findSpecies (elementType);
          if (ig < numSpecies)
          {
            for (int i=0; i<ion2functor.size(); i++)
              if (ions->GroupID[i] == ig)
              {
                ion2functor[i]=funs.size()-1;
                app_log() << "Adding backflow transformation of element type " << elementType << " for atom " << i << ".\n";
              }
          }
        }
      }
      cur = cur->next;
    }
    if(funct == "Bspline")
    {
      app_log() <<"Using BsplineFunctor type. \n";
      tbf = (BackflowFunctionBase *) new Backflow_eI<BsplineFunctor<RealType> >(*ions,targetPtcl);
      Backflow_eI<BsplineFunctor<RealType> > *dum = (Backflow_eI<BsplineFunctor<RealType> >*) tbf;
      tbf->numParams=0;
      std::vector<int> offsets;
      for(int i=0; i<funs.size(); i++)
      {
        //           BsplineFunctor<RealType> *bsp = new BsplineFunctor<RealType>(cusps[i]);
        BsplineFunctor<RealType> *bsp = new BsplineFunctor<RealType>();
        bsp->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
        bsp->put(funs[i]);
        if(bsp->cutoff_radius > cutOff)
          cutOff = bsp->cutoff_radius;
        bsp->myVars.setParameterType(optimize::BACKFLOW_P);
        bsp->print();
        dum->uniqueRadFun.push_back(bsp);
        offsets.push_back(tbf->numParams);
        tbf->numParams += bsp->NumParams;
      }
      tbf->derivs.resize(tbf->numParams);
      dum->offsetPrms.resize(nIons);
      dum->RadFun.resize(nIons);
      for(int i=0; i<ion2functor.size(); i++)
      {
        if(ion2functor[i] < 0 || ion2functor[i] >= funs.size())
        {
          APP_ABORT("backflowTransformation::put() ion not mapped to radial function.\n");
        }
        dum->RadFun[i] = dum->uniqueRadFun[ion2functor[i]];
        dum->offsetPrms[i] = offsets[ion2functor[i]];
      }
    }
    else
    {
      APP_ABORT("Unknown function type in e-I BF Transformation.\n");
    }
  }//spin="no"
  BFTrans->bfFuns.push_back(tbf);
//      tbf->reportStatus(cerr);
}

void BackflowBuilder::addTwoBody(xmlNodePtr cur)
{
  app_log() <<"Adding electron-electron backflow. \n";
  OhmmsAttributeSet trAttrib;
  std::string source("none");
  std::string name("bf0");
  std::string type("none");
  std::string funct("Bspline");
  trAttrib.add (name, "name");
  trAttrib.add (funct, "function");
  trAttrib.put(cur);
  xmlNodePtr curRoot=cur;
  //BackflowFunctionBase *tbf = (BackflowFunctionBase *) new Backflow_ee<BsplineFunctor<RealType> >(targetPtcl,targetPtcl);
  Backflow_ee<BsplineFunctor<RealType> > *tbf = new Backflow_ee<BsplineFunctor<RealType> >(targetPtcl,targetPtcl);
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  std::vector<int> offsets;
  if(funct == "Gaussian")
  {
    APP_ABORT("Disabled GaussianFunctor for now, \n");
  }
  else
    if(funct == "Bspline")
    {
      app_log() <<"Using BsplineFunctor type. \n";
//         BsplineFunctor<RealType> *bsp = new BsplineFunctor<RealType>(cusp);
      std::string cname;
      cur = curRoot->children;
      while (cur != NULL)
      {
        getNodeName(cname,cur);
        if (cname == "correlation")
        {
          RealType cusp=0;
          std::string spA(species.speciesName[0]);
          std::string spB(species.speciesName[0]);
          OhmmsAttributeSet anAttrib;
          anAttrib.add (cusp, "cusp");
          anAttrib.add(spA,"speciesA");
          anAttrib.add(spB,"speciesB");
          anAttrib.add(cusp,"cusp");
          anAttrib.put(cur);
          int ia = species.findSpecies(spA);
          int ib = species.findSpecies(spB);
          if(ia==species.size() || ib == species.size())
          {
            APP_ABORT("Failed. Species are incorrect in e-e backflow.");
          }
          app_log() <<"Adding radial component for species: " <<spA <<" " <<spB <<" " <<ia <<"  " <<ib << std::endl;
          BsplineFunctor<RealType> *bsp = new BsplineFunctor<RealType>();
          bsp->cutoff_radius = targetPtcl.Lattice.WignerSeitzRadius;
          bsp->put(cur);
          if(bsp->cutoff_radius > cutOff)
            cutOff = bsp->cutoff_radius;
          bsp->myVars.setParameterType(optimize::BACKFLOW_P);
          tbf->addFunc(ia,ib,bsp);
          offsets.push_back(tbf->numParams);
          tbf->numParams += bsp->NumParams;
//            if(OHMMS::Controller->rank()==0)
//            {
//              char fname[64];
//              sprintf(fname,"BFe-e.%s.dat",(spA+spB).c_str());
//              std::ofstream fout(fname);
//              fout.setf(std::ios::scientific, std::ios::floatfield);
//              fout << "# Backflow radial function \n";
//              bsp->print(fout);
//              fout.close();
//            }
        }
        cur = cur->next;
      }
      tbf->derivs.resize(tbf->numParams);
      // setup offsets
      // could keep a std::map<std::pair<>,int>
      for(int i=0; i<tbf->RadFun.size(); i++)
      {
        bool done = false;
        for(int k=0; k<tbf->uniqueRadFun.size(); k++)
        {
          if(tbf->RadFun[i] == tbf->uniqueRadFun[k])
          {
            done=true;
            tbf->offsetPrms[i] = offsets[k];
            break;
          }
        }
        if(!done)
        {
          APP_ABORT("Error creating Backflow_ee object. \n");
        }
      }
      BFTrans->bfFuns.push_back((BackflowFunctionBase *)tbf);
    }
    else
    {
      APP_ABORT("Unknown function type in e-e BF Transformation.\n");
    }
}

void BackflowBuilder::addRPA(xmlNodePtr cur)
{
  ReportEngine PRE("BackflowBuilder","addRPA");
  /*
      std::string useL="yes";
      std::string useS="yes";

      OhmmsAttributeSet a;
      a.add(useL,"longrange");
      a.add(useS,"shortrange");
      a.put(cur);
  */
  Rs=-1.0;
  Kc=-1.0;
  RealType my_cusp=0.0;
  OhmmsAttributeSet anAttrib0;
  anAttrib0.add (Rs,"rs");
  anAttrib0.add (Kc,"kc");
  anAttrib0.put(cur);
  //ParameterSet params;
  //params.add(Rs,"rs","RealType");
  //params.add(Kc,"kc","RealType");
  //params.add(Kc,"Kc","RealType");
  //params.put(cur);
  RealType tlen = std::pow(3.0/4.0/M_PI*targetPtcl.Lattice.Volume/ static_cast<RealType>(targetPtcl.getTotalNum()) ,1.0/3.0);
  if(Rs<0)
  {
    if(targetPtcl.Lattice.SuperCellEnum)
    {
      Rs=tlen;
    }
    else
    {
      std::cout <<"  Error finding rs. Is this an open system?!"<< std::endl;
      Rs=100.0;
    }
  }
  int indx = targetPtcl.SK->KLists.ksq.size()-1;
  RealType Kc_max=std::pow(targetPtcl.SK->KLists.ksq[indx],0.5);
  if(Kc<0)
  {
    Kc = 2.0*  std::pow(2.25*M_PI,1.0/3.0)/tlen ;
  }
  if(Kc>Kc_max)
  {
    Kc=Kc_max;
    app_log() << " BackflowBuilderBuilder   Kc set too high. Resetting to the maximum value"<< std::endl;
  }
  app_log() << "    BackflowBuilderBuilder   Rs = " << Rs <<  "  Kc= " << Kc << std::endl;
  // mmorales: LRRPABFeeHandlerTemp is a copy of LRRPAHandlerTemp for now,
  // in case I need to specialize it later on
  myHandler= new LRRPABFeeHandlerTemp<RPABFeeBreakup<RealType>,LPQHIBasis>(targetPtcl,Kc);
  myHandler->Breakup(targetPtcl,Rs);
  app_log() << "  Maximum K shell " << myHandler->MaxKshell << std::endl;
  app_log() << "  Number of k vectors " << myHandler->Fk.size() << std::endl;
  std::vector<int> offsetsSR;
  std::vector<int> offsetsLR;
  Backflow_ee<BsplineFunctor<RealType> > *tbf = 0;
  Backflow_ee_kSpace *tbfks = 0;
  // now look for components
  std::string cname;
  xmlNodePtr curRoot = cur;
  cur = cur->children;
  while (cur != NULL)
  {
    getNodeName(cname,cur);
    if (cname == "correlation")
    {
      std::string type = "none";
      OhmmsAttributeSet anAttrib;
      anAttrib.add(type,"type");
      anAttrib.put(cur);
      if(type == "shortrange")
      {
        if(tbf==0)
        {
          tbf = new Backflow_ee<BsplineFunctor<RealType> >(targetPtcl,targetPtcl);
        }
        makeShortRange_twoBody(cur,tbf,offsetsSR);
      }
      else
        if(type == "longrange")
        {
          if(tbfks==0)
          {
            tbfks = new Backflow_ee_kSpace(targetPtcl,targetPtcl);
          }
          else
          {
            APP_ABORT("Only a single LongRange RPAbackflow allowed for now. ");
          }
          makeLongRange_twoBody(cur,tbfks,offsetsLR);
        }
        else
        {
          APP_ABORT("Unknown rpa backflow type in <correlation/>.");
        }
    }
    cur = cur->next;
  }
  if(tbf != 0)
  {
    tbf->derivs.resize(tbf->numParams);
    // setup offsets
    // could keep a std::map<std::pair<>,int>
    for(int i=0; i<tbf->RadFun.size(); i++)
    {
      bool done = false;
      for(int k=0; k<tbf->uniqueRadFun.size(); k++)
      {
        if(tbf->RadFun[i] == tbf->uniqueRadFun[k])
        {
          done=true;
// how do I calculate the offset now???
          tbf->offsetPrms[i] = offsetsSR[k];
          break;
        }
      }
      if(!done)
      {
        APP_ABORT("Error creating Backflow_ee object in addRPA. \n");
      }
    }
    BFTrans->bfFuns.push_back((BackflowFunctionBase *)tbf);
  }
}

void BackflowBuilder::makeLongRange_oneBody() {}
void BackflowBuilder::makeLongRange_twoBody(xmlNodePtr cur, Backflow_ee_kSpace *tbfks, std::vector<int>& offsets)
{
  int size=-1;
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  std::string spA(species.speciesName[0]);
  std::string spB(species.speciesName[0]);
  std::string init = "yes";
  OhmmsAttributeSet anAttrib;
  anAttrib.add(size,"size");
  anAttrib.add(init,"init");
  anAttrib.add(init,"initialize");
  anAttrib.add(spA,"speciesA");
  anAttrib.add(spB,"speciesB");
  anAttrib.put(cur);
  int ia = species.findSpecies(spA);
  int ib = species.findSpecies(spB);
  if(ia==species.size() || ib == species.size())
  {
    APP_ABORT("Failed. Species are incorrect in longrange RPA backflow.");
  }
  app_log() <<"Adding RPABackflow longrange component for species: " <<spA <<" " <<spB <<" " <<ia <<"  " <<ib << std::endl;
  // Now read coefficents
  xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
  while (xmlCoefs != NULL)
  {
    std::string cname((const char*)xmlCoefs->name);
    if (cname == "coefficients")
    {
      std::string type("0"), id("0");
      std::string optimize("no");
      OhmmsAttributeSet cAttrib;
      cAttrib.add(id, "id");
      cAttrib.add(type, "type");
      cAttrib.add(optimize, "optimize");
      cAttrib.put(xmlCoefs);
      if (type != "Array")
      {
        APP_ABORT("Unknown correlation type " + type + " in Backflow.");
      }
      if(optimize == "true" || optimize == "yes")
        tbfks->Optimize = true;
      std::vector<RealType> yk;
      if(init == "true" || init == "yes")
      {
        app_log() <<"Initializing k-space backflow function with RPA form.";
        yk.resize(myHandler->Fk_symm.size());
        for(int i=0; i<yk.size(); i++)
          yk[i] = myHandler->Fk_symm[i];
      }
      else
      {
        app_log() <<"Reading k-space backflow function from xml.";
        putContent(yk, xmlCoefs);
      }
      tbfks->initialize(targetPtcl,yk);
      tbfks->myVars.setParameterType(optimize::BACKFLOW_P);
      tbfks->addFunc(ia,ib);
      offsets.push_back(tbfks->numParams);
      if(OHMMS::Controller->rank()==0)
      {
        char fname[16];
        sprintf(fname,"RPABFee-LR.%s.dat",(spA+spB).c_str());
        std::ofstream fout(fname);
        fout.setf(std::ios::scientific, std::ios::floatfield);
        fout << "# Backflow longrange  \n";
        for(int i=0; i<tbfks->NumKShells; i++)
        {
          fout<<std::pow(targetPtcl.SK->KLists.ksq[targetPtcl.SK->KLists.kshell[i]],0.5) <<" " <<yk[i] << std::endl;
        }
        fout.close();
      }
    }
    xmlCoefs=xmlCoefs->next;
  }
}

void BackflowBuilder::makeShortRange_oneBody() {}
void BackflowBuilder::makeShortRange_twoBody(xmlNodePtr cur, Backflow_ee<BsplineFunctor<RealType> > *tbf, std::vector<int>& offsets)
{
  int size=-1;
  SpeciesSet& species(targetPtcl.getSpeciesSet());
  std::string spA(species.speciesName[0]);
  std::string spB(species.speciesName[0]);
  RealType cusp=0.0;
  std::string init = "yes";
  OhmmsAttributeSet anAttrib;
  anAttrib.add(cusp,"cusp");
  anAttrib.add(size,"size");
  anAttrib.add(init,"init");
  anAttrib.add(init,"initialize");
  anAttrib.add(spA,"speciesA");
  anAttrib.add(spB,"speciesB");
  anAttrib.put(cur);
  int ia = species.findSpecies(spA);
  int ib = species.findSpecies(spB);
  if(ia==species.size() || ib == species.size())
  {
    APP_ABORT("Failed. Species are incorrect in e-e backflow.");
  }
  app_log() <<"Adding radial component for species: " <<spA <<" " <<spB <<" " <<ia <<"  " <<ib << std::endl;
  // Now read coefficents
  xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
  while (xmlCoefs != NULL)
  {
    std::string cname((const char*)xmlCoefs->name);
    if (cname == "coefficients")
    {
      std::string type("0"), id("0");
      std::string optimize("no");
      OhmmsAttributeSet cAttrib;
      cAttrib.add(id, "id");
      cAttrib.add(type, "type");
      cAttrib.add(optimize, "optimize");
      cAttrib.put(xmlCoefs);
      if (type != "Array")
      {
        APP_ABORT("Unknown correlation type " + type + " in Backflow.");
      }
      BsplineFunctor<RealType> *bsp = new BsplineFunctor<RealType>();
      if(init == "true" || init == "yes")
      {
        app_log() <<"Initializing backflow radial functions with RPA.";
        Rcut = myHandler->get_rc()-0.01;
        GridType* myGrid = new GridType;
        int npts=static_cast<int>(Rcut/0.01)+3;
        myGrid->set(0,Rcut-0.01,npts);
        //create the numerical functor
        std::vector<RealType> x(myGrid->size()),y(myGrid->size());
//              x[0]=(*myGrid)(0);
//              RealType x0=x[0];
//              if(x0 < 1.0e-6) x0 = 1.0e-6;
////              y[0]=-myHandler->srDf(x0,1.0/x0)/x0;
//              y[0]=myHandler->evaluate(x0,1.0/x0);
        for  (int i = 0; i < myGrid->size(); i++)
        {
          x[i]=(*myGrid)(i);
//                y[i]=-myHandler->srDf(x[i],1.0/x[i])/x[i];
          y[i]=myHandler->evaluate(x[i],0.0);
        }
        app_log() <<"Rcut,npts:" <<Rcut <<"  " <<npts <<"  " <<x[myGrid->size()-1] << std::endl;
//fit potential to gaussians
        int nfitgaussians(3);
        Matrix<RealType> basis(myGrid->size(),nfitgaussians);
        for (int i=0; i<npts; i++)
        {
          RealType r = x[i];
          for (int j=0; j<nfitgaussians; j++)
            basis(i,j) = std::exp(-r*r/((j+1)*Rcut*Rcut))-std::exp(-1/(j+1));
        }
        std::vector<RealType> gb(nfitgaussians);
        LinearFit(y, basis, gb);
//spline deriv of gaussian fit
        for  (int i = 0; i < myGrid->size(); i++)
        {
          RealType r = x[i];
          y[i]=0;
          for (int j=0; j<nfitgaussians; j++)
            y[i]+=gb[j]*(2.0/((j+1)*Rcut*Rcut)*std::exp(-r*r/((j+1)*Rcut*Rcut)));
        }
//              RealType y1_c(0),y2_c(0);
//              for (int j=0; j<nfitgaussians; j++) y1_c+=gb[j]*(std::exp(-x[1]*x[1]/((j+1)*Rcut*Rcut)));
//              for (int j=0; j<nfitgaussians; j++) y2_c+=gb[j]*(std::exp(-x[2]*x[2]/((j+1)*Rcut*Rcut)));
//make a temp functor to ensure right BC's (Necessary?)
        BsplineFunctor<RealType> *tmp_bsp = new BsplineFunctor<RealType>();
        tmp_bsp->initialize(12,x,y,cusp,Rcut,id,optimize);
//              tmp_bsp->print(app_log());
        for  (int i = 0; i < myGrid->size(); i++)
        {
          y[i]=tmp_bsp->f(x[i]);
        }
        delete tmp_bsp;
//make functor for backflow
        bsp->initialize(size,x,y,cusp,Rcut,id,optimize);
      }
      else
      {
        bsp->put(cur);
      }
      if(bsp->cutoff_radius > cutOff)
        cutOff = bsp->cutoff_radius;
      bsp->myVars.setParameterType(optimize::BACKFLOW_P);
      tbf->addFunc(ia,ib,bsp);
      offsets.push_back(tbf->numParams);
      tbf->numParams += bsp->NumParams;
//            if(OHMMS::Controller->rank()==0)
//            {
//              char fname[64];
//              sprintf(fname,"RPABFee-SR.%s.dat",(spA+spB).c_str());
//              std::ofstream fout(fname);
//              fout.setf(std::ios::scientific, std::ios::floatfield);
//              fout << "# Backflow radial function \n";
//              bsp->print(fout);
//              fout.close();
//            }
    }
    xmlCoefs=xmlCoefs->next;
  }
}




}

