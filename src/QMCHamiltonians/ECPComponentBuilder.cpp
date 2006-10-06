//////////////////////////////////////////////////////////////////
// (c) Copyright 2003-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
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
//////////////////////////////////////////////////////////////////
#include "QMCHamiltonians/ECPComponentBuilder.h"
#include "Numerics/GaussianTimesRN.h"
#include "Numerics/Transform2GridFunctor.h"
#include "QMCFactory/OneDimGridFactory.h"

namespace qmcplusplus {

  ECPComponentBuilder::ECPComponentBuilder(const string& aname):
    NumNonLocal(0), Lmax(0), Zeff(0), Species(aname),
    pp_loc(0), pp_nonloc(0) {
    angMon["s"]=0; angMon["p"]=1; angMon["d"]=2; angMon["f"]=3; angMon["g"]=4;
    angMon["0"]=0; angMon["1"]=1; angMon["2"]=2; angMon["3"]=3; angMon["4"]=4;
  }

  bool ECPComponentBuilder::parse(const string& fname) {

    // build an XML tree from a the file;
    xmlDocPtr m_doc = xmlParseFile(fname.c_str());
    if (m_doc == NULL) {
      ERRORMSG("File " << fname << " is invalid")
        xmlFreeDoc(m_doc);
    }    
    // Check the document is of the right kind
    xmlNodePtr cur = xmlDocGetRootElement(m_doc);
    if (cur == NULL) {
      ERRORMSG("Empty document");
      xmlFreeDoc(m_doc);
    }

    bool success=put(cur);
    xmlFreeDoc(m_doc);

    return success;
  }

  bool ECPComponentBuilder::put(xmlNodePtr cur) {
    int nk=0;
    vector<RealType> kpts;
    vector<xmlNodePtr> semiPtr;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "header") {
        Zeff = atoi((const char*)xmlGetProp(cur,(const xmlChar*)"zval"));
      } else if(cname == "semilocal") {
        semiPtr.push_back(cur);//save the pointer
      } else if(cname == "local") {
        buildLocal(cur);
      } else if(cname == "sphericalGrid") {
        nk=atoi((const char*)xmlGetProp(cur,(const xmlChar*)"size"));
        kpts.resize(nk*4);
        putContent(kpts,cur);
      }
      cur=cur->next;
    }

    if(semiPtr.size()) {
      if(pp_nonloc==0) pp_nonloc=new NonLocalECPComponent; 
      if(pp_loc) {
        for(int i=0; i<semiPtr.size(); i++) addSemiLocal(semiPtr[i]);
      } else {
        buildSemiLocalAndLocal(semiPtr);
      }
      app_log() << "    Maximum angular momentum of NonLocalECP " << pp_nonloc->lmax << endl;
      app_log() << "    Maximum cutoff of NonLocalECP " << pp_nonloc->Rmax << endl;
    } 

    if(pp_nonloc) {
      pp_nonloc->lmax=Lmax;
      if(nk == 0) { //assign default
        RealType w=1.0/12.0;
        pp_nonloc->addknot(PosType( 0.000000, 0.525731, 0.850651),w);
        pp_nonloc->addknot(PosType( 0.000000,-0.525731,-0.850651),w);
        pp_nonloc->addknot(PosType( 0.000000,-0.525731, 0.850651),w);
        pp_nonloc->addknot(PosType( 0.000000, 0.525731,-0.850651),w);
        pp_nonloc->addknot(PosType(-0.850651, 0.000000, 0.525731),w);
        pp_nonloc->addknot(PosType( 0.850651, 0.000000,-0.525731),w);
        pp_nonloc->addknot(PosType( 0.850651, 0.000000, 0.525731),w);
        pp_nonloc->addknot(PosType(-0.850651, 0.000000,-0.525731),w);
        pp_nonloc->addknot(PosType(-0.525731, 0.850651, 0.000000),w);
        pp_nonloc->addknot(PosType( 0.525731,-0.850651, 0.000000),w);
        pp_nonloc->addknot(PosType( 0.525731, 0.850651, 0.000000),w);
        pp_nonloc->addknot(PosType(-0.525731,-0.850651, 0.000000),w);
        pp_nonloc->resize_warrays(12,NumNonLocal,Lmax);
      } else {
        for(int ik=0, kk=0; ik<nk; ik++, kk+=4) {
          pp_nonloc->addknot(PosType(kpts[kk],kpts[kk+1],kpts[kk+2]),kpts[kk+3]);
        }
        pp_nonloc->resize_warrays(nk,NumNonLocal,Lmax);
      }
    }

    return true;
  }

  void ECPComponentBuilder::addSemiLocal(xmlNodePtr cur) {

    GridType* grid_semilocal=0;
    RealType rmax= pp_nonloc->Rmax;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_semilocal=OneDimGridFactory::createGrid(cur);
        rmax=grid_semilocal->rmax();
      } else if(cname == "vps") {
        //should be able to overwrite rmax
        int l=angMon[(const char*)xmlGetProp(cur,(const xmlChar*)"l")];
        Lmax = std::max(l,Lmax);
        xmlNodePtr cur1=cur->children;
        while(cur1 != NULL) {
          if(xmlStrEqual(cur1->name,(const xmlChar*)"basisGroup")) {
            pp_nonloc->add(l,createVr(cur1,grid_semilocal));
          }
          cur1=cur1->next;
        }
        NumNonLocal++;
      }
      cur=cur->next;
    }

    pp_nonloc->lmax=Lmax;
    pp_nonloc->Rmax=rmax;

    //ofstream fout("C.semilocal.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //int ig=0;
    //while(ig<grid_semilocal->size()) {
    //  double r=(*grid_semilocal)[ig++];
    //  fout << setw(20) << r;
    //  for(int i=0; i<pp_nonloc->nlpp_m.size(); i++) {
    //    fout << setw(20) << pp_nonloc->nlpp_m[i]->operator()(ig);
    //  }
    //  fout << endl;
    //}
  }

  ECPComponentBuilder::RadialPotentialType*
  ECPComponentBuilder::createVr(xmlNodePtr cur, GridType* agrid) {
    //todo rcut should be reset if necessary
    typedef GaussianTimesRN<RealType> InFuncType;
    InFuncType a;
    a.putBasisGroup(cur);
    RadialPotentialType *app=new RadialPotentialType(agrid);
    Transform2GridFunctor<InFuncType,RadialPotentialType> transform(a,*app);
    transform.generate(agrid->rmin(),agrid->rmax(), agrid->size());
    //int ig=0;
    //while(ig<agrid->size()) {
    //  (*app)(ig)=a.f((*agrid)[ig]);
    //  ++ig;
    //}
    //app->spline(0,0.0,app->size()-1,0.0);
    return app;
  }

  void ECPComponentBuilder::buildLocal(xmlNodePtr cur) {

    if(pp_loc) return; //something is wrong

    typedef GaussianTimesRN<RealType> InFuncType;
    GridType* grid_local=0;
    InFuncType vr;

    bool bareCoulomb=true;
    cur=cur->children;
    while(cur != NULL) {
      string cname((const char*)cur->name);
      if(cname == "grid") {
        grid_local=OneDimGridFactory::createGrid(cur);
      } else if(cname == "basisGroup") {
        vr.putBasisGroup(cur);
        bareCoulomb=false;
      }
      cur=cur->next;
    }

    int ng=grid_local->size();
    vector<RealType> v(ng);
    if(bareCoulomb) {
      for(int ig=0; ig<ng; ig++) {
        v[ig]=1.0/(*grid_local)[ig];
      }
    } else {
      RealType zinv=1.0/Zeff;
      for(int ig=0; ig<ng; ig++) {
        double r=(*grid_local)[ig];
        v[ig]=1.0/r-zinv*vr.f(r);
      }
    }
    //if(bareCoulomb) {
    //  for(int ig=0; ig<ng; ig++) {
    //    v[ig]=-Zeff/(*grid_local)[ig];
    //  }
    //} else {
    //  for(int ig=0; ig<ng; ig++) {
    //    double r=(*grid_local)[ig];
    //    v[ig]=vr.f(r)-Zeff/r;
    //  }
    //}

    pp_loc=new RadialPotentialType(grid_local,v);
    int imin = 0;
    RealType yprime_i = ((*pp_loc)(imin+1)-(*pp_loc)(imin))/pp_loc->dr(imin);
    pp_loc->spline(imin,yprime_i,ng-1,0.0);

    //cout << " Effective Z = " << Zeff << endl;
    //ofstream fout("C.local.dat");
    //fout.setf(std::ios::scientific, std::ios::floatfield);
    //fout.precision(12);
    //int ig=0;
    //while(ig<grid_semilocal->size()) {
    //  double r=(*grid_semilocal)[ig++];
    //  fout << setw(20) << r << setw(20) << localpp->f(r)-Zeff/r<< endl;
    //}
  }

  void ECPComponentBuilder::buildSemiLocalAndLocal(vector<xmlNodePtr>& semiPtr) {
    //this is abinit/siesta format

  }
}
/***************************************************************************
 * $RCSfile$   $Author$
 * $Revision$   $Date$
 * $Id$ 
 ***************************************************************************/
