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
    
    
#ifndef QMCPLUSPLUS_NUMBERFLUC_H
#define QMCPLUSPLUS_NUMBERFLUC_H
#include "Particle/ParticleSet.h"
#include "Particle/WalkerSetRef.h"
#include "Particle/DistanceTableData.h"
#include "Particle/DistanceTable.h"
#include "QMCHamiltonians/QMCHamiltonianBase.h"
#include <numeric>
#include <OhmmsData/AttributeSet.h>

namespace qmcplusplus
{

struct NumberFluctuations: public QMCHamiltonianBase
{

  std::vector<Return_t> regions;
  std::string shape;
  ParticleSet::ParticlePos_t *L, *C;
  bool PBCType;
  int nreg,grps,pwrs,ngrid;
  std::vector<int> grp_v;
  std::vector<Return_t>  Values, voffsets;
  DistanceTableData* d_table;


  NumberFluctuations(ParticleSet& P): ngrid(2), pwrs(2), shape("sphere")
  {
//       C = new ParticleSet::ParticlePos_t(P.R.size());
//       C->setUnit(PosUnit::LatticeUnit);
    L = new ParticleSet::ParticlePos_t(P.R.size());
    L->setUnit(PosUnit::CartesianUnit);
    grps=P.groups();
    grp_v.resize(P.R.size(),0);
    int i(P.R.size()),g(grps-1);
    while((i>0)&&(g>0))
    {
      if(i<P.first(g))
        g--;
      grp_v[i]=g;
      i--;
    }
    PBCType=(P.Lattice.SuperCellEnum>1);
    d_table = DistanceTable::add(P,DT_AOS);
  }

  ~NumberFluctuations() { }

  void resetTargetParticleSet(ParticleSet& P)
  {
    grps=P.groups();
    grp_v.resize(P.R.size(),0);
    int i(P.R.size()),g(grps-1);
    while((i>0)&&(g>0))
    {
      if(i<P.first(g))
        g--;
      grp_v[i]=g;
      i--;
    }
    Values.resize( (grps+1)*nreg*pwrs,0 );
  }

//     inline Return_t
//     evaluateSphere(ParticleSet& P)
//     {
//      for (int g=0; g<Values.size(); ++g) Values[g]=0;
// #if OHMMS_DIM==3
//      int npts=std::pow(ngrid,3);
// #endif
// #if  OHMMS_DIM==2
//      int npts=std::pow(ngrid,2);
// #endif
//      Return_t grdnrm=1.0/npts;
//      for(int j(0);j<npts;j++)
//      {
//       C->setUnit(PosUnit::LatticeUnit);
//       L->setUnit(PosUnit::CartesianUnit);
//       int l(j);
//       for(int k(DIM-1);k>0;k--)
//       {
//         int dvi=(std::pow(ngrid,k));
//         C[k]=(l/dvi)*grdnrm;
//         l-=(l/dvi)*std::pow(ngrid,k);
//       }
// //       for(int k(0);k<DIM;k++) app_log()<<C[k]<<" ";
// //       app_log()<< std::endl;
//       P.convert2Cart((*C));
// //       for(int k(0);k<DIM;k++) app_log()<<C[k]<<" ";
// //       app_log()<< std::endl;
//       (*L)=P.R;
//       for(int k(0);k<(*L).size();k++) (*L)[k]+= -1.0*(*C)[0];
// //       for(int k(0);k<DIM;k++) app_log()<<L[0][k]<<" ";
// //       app_log()<< std::endl;
//       //       P.applyBC(*L);
//       P.convert2CartInBox(*L, *C);
//       Matrix<int> np(grps+1,nreg);
//       np=0;
//       for (int g=0; g<grps; ++g)
//         for (int iat=P.first(g); iat<P.last(g); ++iat)
//         {
//           Return_t z=dot((*C)[iat],(*C)[iat]);
//           int a(nreg-1);
//           while((z<regions[a])&&(a>=0))
//           {
//             np(g,a)+=1;
//             np(grps,a)+=1;
//             a--;
//           }
//         }
//
//       int indx(0);
//       for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++)
//       {
//         Return_t vx=np(g,h)-voffsets[g*nreg+h];
//         Return_t vy(1);
//         for (int i=0; i<pwrs; i++,++indx)
//           Values[indx]+=(vy*=vx);
//       }
//      }
//       int indx(0);
//       for (int g=0; g<grps+1; ++g)
//         for (int h=0; h<nreg; h++)
//           for (int i=0; i<pwrs; i++,++indx)
//             Values[indx]*=grdnrm;
//        return 0;
//     }

  inline Return_t
  evaluateSphere(ParticleSet& P)
  {
    for (int g=0; g<Values.size(); ++g)
      Values[g]=0;
    Return_t nnrm(1.0/P.R.size());
    Return_t pnrm(1.0/P.Lattice.Volume);
    std::vector<Matrix<int> > np(P.R.size(), Matrix<int>(grps+1,nreg));
    for(int i(0); i<P.R.size(); i++)
      np[i]=0;
    //Center
    for(int i(0); i<P.R.size(); i++)
      for (int h=0; h<nreg; h++)
      {
        np[i](grp_v[i],h)+=1;
        np[i](grps,h)+=1;
      }
    d_table=P.DistTables[0];
    for (int g=0; g<grps; ++g)
      for (int iat=P.first(g); iat<P.last(g); ++iat)
      {
        for(int nn=0; nn<d_table->nadj(iat); ++nn)
        {
          int jat=d_table->iadj(iat,nn);
          Return_t z = d_table->r(d_table->loc(iat,nn));
//             app_log()<<iat<<" "<<jat<<" "<<z<< std::endl;
          z*=z;
          int a(nreg-1);
          while((z<regions[a])&&(a>=0))
          {
            np[iat](grp_v[jat],a)+=1;
            np[iat](grps,a)+=1;
            np[jat](grp_v[iat],a)+=1;
            np[jat](grps,a)+=1;
            a--;
          }
        }
      }
    for(int i(0); i<P.R.size(); i++)
    {
      int indx(0);
      for (int g=0; g<grps+1; ++g)
        for (int h=0; h<nreg; h++)
        {
          Return_t vx=np[i](g,h)-voffsets[g*nreg+h];
          Return_t vy(1.0);
          for (int k=0; k<pwrs; k++,++indx)
            Values[indx]+=(vy*=vx);
        }
    }
//       int indx(0);
//       for (int g=0; g<grps+1; ++g) for (int h=0; h<nreg; h++) for (int k=0; k<pwrs; k++,++indx) Values[indx]*=pnrm;
    for (int g=0; g<Values.size(); ++g)
      Values[g]*=nnrm;
    return 0;
  }



  inline Return_t
  evaluateBox(ParticleSet& P)
  {
    P.convert2CartInBox(P.R,*L);
    Matrix<int> np(grps+1,nreg);
    np=0;
    for (int g=0; g<grps; ++g)
      for (int iat=P.first(g); iat<P.last(g); ++iat)
      {
        Return_t z=(*L)[iat][0];
        for (int x=1; x<DIM; x++)
          z=std::max(z,static_cast<Return_t>((*L)[iat][x]));
        int a(nreg-1);
        while((z<regions[a])&&(a>=0))
        {
          np(g,a)+=1;
          np(grps,a)+=1;
          a--;
        }
      }
    for (int g=0; g<Values.size(); ++g)
      Values[g]=0;
    int indx(0);
    for (int g=0; g<grps+1; ++g)
      for (int h=0; h<nreg; h++)
        for (int i=0; i<pwrs; i++,++indx)
          Values[indx]+=std::pow(np(g,h)-voffsets[g*nreg+h],i+1);
    return 0;
  }

  inline Return_t
  evaluateHalfSpace(ParticleSet& P)
  {
    if(PBCType)
      P.convert2CartInBox(P.R,*L);
    Matrix<int> np(grps+1,nreg);
    np=0;
    for (int g=0; g<grps; ++g)
      for (int iat=P.first(g); iat<P.last(g); ++iat)
      {
        Return_t z=(*L)[iat][DIM-1];
        int a(nreg-1);
        while((z<regions[a])&&(a>=0))
        {
          np(g,a)+=1;
          np(grps,a)+=1;
          a--;
        }
      }
    for (int g=0; g<Values.size(); ++g)
      Values[g]=0;
    int indx(0);
    for (int g=0; g<grps+1; ++g)
      for (int h=0; h<nreg; h++)
        for (int i=0; i<pwrs; i++,++indx)
          Values[indx]+=std::pow(np(g,h)-voffsets[g*nreg+h],i+1);
    return 0;
  }

  void addObservables(PropertySetType& plist, BufferType& collectables)
  {
    myIndex=plist.size();
    for (int g=0; g<grps+1; ++g)
      for (int h=0; h<nreg; h++)
        for (int i=0; i<pwrs; i++)
        {
          std::stringstream sstr;
          sstr << "N"<<g<<"_R"<<h<<"_P"<<i+1;
          int id=plist.add(sstr.str());
        }
  }

  void setObservables(PropertySetType& plist)
  {
    copy(Values.begin(),Values.end(),plist.begin()+myIndex);
  }

  void setParticlePropertyList(PropertySetType& plist, int offset)
  {
    copy(Values.begin(),Values.end(),plist.begin()+myIndex+offset);
  }

  inline Return_t evaluate(ParticleSet& P)
  {
    if (shape=="sphere")
      return evaluateSphere(P);
    else
      if (shape=="box")
        return evaluateBox(P);
      else
        if (shape=="half")
          return evaluateHalfSpace(P);
    return 0.0;
  }

  inline Return_t evaluate(ParticleSet& P, std::vector<NonLocalData>& Txy)
  {
    return evaluate(P);
  }

  /** Do nothing */
  bool put(xmlNodePtr cur)
  {
    OhmmsAttributeSet aAttrib;
    aAttrib.add(pwrs, "max");
    aAttrib.add(ngrid,"grid"),
                aAttrib.put(cur);
    xmlNodePtr xmlCoefs = cur->xmlChildrenNode;
    while (xmlCoefs != NULL)
    {
      std::string cname((const char*)xmlCoefs->name);
      if (cname == "regions")
      {
        OhmmsAttributeSet cAttrib;
        cAttrib.add(shape, "shape");
        cAttrib.put(xmlCoefs);
        putContent(regions, xmlCoefs);
//               std::sort(regions.begin(),regions.end());
      }
      if (cname == "rho")
        putContent(voffsets, xmlCoefs);
      xmlCoefs = xmlCoefs->next;
    }
    nreg=regions.size();
    app_log()<<" Shape: "<<shape<< std::endl;
    app_log()<<" Max Power: "<<pwrs<< std::endl;
    app_log()<<" Grid: "<<ngrid<< std::endl;
    app_log()<<" regions:";
    for(int i(0); i<nreg; i++)
      app_log()<<" "<<regions[i];
    app_log()<< std::endl;
    if(voffsets.size()!=nreg*(grps+1))
      voffsets.resize(nreg*(grps+1),0);
    app_log()<<" voffsets:";
    for(int i(0); i<nreg*(grps+1); i++)
      app_log()<<" "<<voffsets[i];
    app_log()<< std::endl;
    Values.resize( (grps+1)*nreg*pwrs,0 );
    return true;
  }

  bool get(std::ostream& os) const
  {
    return true;
  }

  QMCHamiltonianBase* makeClone(ParticleSet& qp, TrialWaveFunction& psi)
  {
    NumberFluctuations* nf = new NumberFluctuations(qp);
    nf->nreg=nreg;
    nf->pwrs=pwrs;
    nf->grps=grps;
    nf->shape=shape;
    nf->myIndex=myIndex;
    nf->regions=regions;
    nf->voffsets=voffsets;
    nf->ngrid=ngrid;
    nf->Values.resize( (grps+1)*nreg*pwrs,0 );
    return nf;
  }

};
}
#endif


