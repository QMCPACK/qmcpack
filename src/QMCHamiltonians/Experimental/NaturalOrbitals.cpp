//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include <QMCHamiltonians/NaturalOrbitals.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <Numerics/OhmmsBlas.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <Particle/DistanceTableData.h>
#include <Numerics/DeterminantOperators.h>

#include <set>

namespace qmcplusplus
{

NaturalOrbitals::NaturalOrbitals(ParticleSet& elns, TrialWaveFunction& psi)
  :M(1), refPsi(psi), Rprime(elns), overwrite_hdf5(false)
{
  UpdateMode.set(COLLECTABLE,1);
  psi_ratios.resize(elns.getTotalNum());
}

void NaturalOrbitals::resetTargetParticleSet(ParticleSet& P)
{
}

void moveRprime()
{
  PosType newpos;
  for (int i=0; i<OHMMS_DIM; ++i)
    newpos[i]=myRNG();
  for(int i=0; i<steps; i++)
  {
    Phi->evaluate(Rprime, 0, temp_phiV);
    P_r_prime=Dot(temp_phiV,temp_phiV);
    if(P_r_prime>P_r)
    {
      P_r=P_r_prime;
    }
  }
}


NaturalOrbitals::Return_t NaturalOrbitals::evaluate(ParticleSet& P)
{
  const int np=P.getTotalNum();
  const std::vector<DistanceTableData::TempDistType>& temp(P.DistTables[0]->Temp);
  Vector<RealType> tmpn_k(nofK);
  for (int s=0; s<M; ++s)
  {
    PosType newpos;
    for (int i=0; i<OHMMS_DIM; ++i)
      newpos[i]=myRNG();
    //make it cartesian
    newpos=Lattice.toCart(newpos);
    P.makeVirtualMoves(newpos); //updated: temp[i].r1=|newpos-P.R[i]|, temp[i].dr1=newpos-P.R[i]
    refPsi.get_ratios(P,psi_ratios);
//         for (int i=0; i<np; ++i) app_log()<<i<<" "<<psi_ratios[i].real()<<" "<<psi_ratios[i].imag()<< std::endl;
    P.rejectMove(0); //restore P.R[0] to the orginal position
    for (int ik=0; ik < kPoints.size(); ++ik)
    {
      for (int i=0; i<np; ++i)
        kdotp[i]=dot(kPoints[ik],temp[i].dr1_nobox);
      eval_e2iphi(np,kdotp.data(),phases.data());
      RealType nofk_here(std::real(BLAS::dot(np,phases.data(),&psi_ratios[0])));//psi_ratios.data())));
      nofK[ik]+= nofk_here;
      tmpn_k[ik]=nofk_here;
    }
    for (int iq=0; iq < compQ.size(); ++iq)
      for (int i=0; i<mappedQtonofK[iq].size(); ++i)
        compQ[iq] += tmpn_k[mappedQtonofK[iq][i]];
  }
  for (int ik=0; ik<nofK.size(); ++ik)
    nofK[ik] *= norm_nofK;
  for (int iq=0; iq<compQ.size(); ++iq)
    compQ[iq] *= mappedQnorms[iq];
  if (hdf5_out)
  {
    int j=myIndex;
    for (int ik=0; ik<nofK.size(); ++ik,++j)
      P.Collectables[j]+= nofK[ik];
    for (int iq=0; iq<compQ.size(); ++iq,++j)
      P.Collectables[j]+= compQ[iq];
  }
  return 0.0;
}

void NaturalOrbitals::registerCollectables(std::vector<observable_helper*>& h5desc
    , hid_t gid) const
{
  if (hdf5_out)
  {
    //descriptor for the data, 1-D data
    std::vector<int> ng(1);
    //add nofk
    ng[0]=nofK.size();
    observable_helper* h5o=new observable_helper("nofk");
    h5o->set_dimensions(ng,myIndex);
    h5o->open(gid);
    h5o->addProperty(const_cast<std::vector<PosType>&>(kPoints),"kpoints");
    h5o->addProperty(const_cast<std::vector<int>&>(kWeights),"kweights");
    h5desc.push_back(h5o);
    //add compQ
    ng[0]=Q.size();
    h5o=new observable_helper("compQ");
    h5o->set_dimensions(ng,myIndex+nofK.size());
    h5o->open(gid);
    h5o->addProperty(const_cast<std::vector<RealType>&>(Q),"q");
    h5desc.push_back(h5o);
  }
}


void NaturalOrbitals::addObservables(PropertySetType& plist, BufferType& collectables)
{
  if (hdf5_out)
  {
    myIndex=collectables.size();
    collectables.add(nofK.begin(),nofK.end());
    collectables.add(compQ.begin(),compQ.end());
  }
  else
  {
    myIndex=plist.size();
    for (int i=0; i<nofK.size(); i++)
    {
      std::stringstream sstr;
      sstr << "nofk_" <<i;
      int id=plist.add(sstr.str());
    }
    for (int i=0; i<Q.size(); i++)
    {
      std::stringstream sstr;
      sstr << "Q_" <<i;
      int id=plist.add(sstr.str());
    }
  }
}



void NaturalOrbitals::setObservables(PropertySetType& plist)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex);
    copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size());
  }
}

void NaturalOrbitals::setParticlePropertyList(PropertySetType& plist
    , int offset)
{
  if (!hdf5_out)
  {
    copy(nofK.begin(),nofK.end(),plist.begin()+myIndex+offset);
    copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size()+offset);
  }
}

bool NaturalOrbitals::putSpecial(xmlNodePtr cur, ParticleSet& elns, bool rootNode)
{
  OhmmsAttributeSet pAttrib;
  std::string hdf5_flag="yes";
  pAttrib.add(hdf5_flag,"hdf5");
  pAttrib.add(M,"samples");
  pAttrib.put(cur);
  hdf5_out = (hdf5_flag=="yes");
//     app_log()<<" NaturalOrbitals::putSpecial "<< std::endl;
  xmlNodePtr kids=cur->children;
  while (kids!=NULL)
  {
    std::string cname((const char*)(kids->name));
//         app_log()<<" NaturalOrbitals::cname : "<<cname<< std::endl;
    if (cname=="kpoints")
    {
      std::string ctype("manual");
      OhmmsAttributeSet pAttrib;
      pAttrib.add(ctype,"mode");
      pAttrib.add(kgrid,"grid");
      pAttrib.put(kids);
#if OHMMS_DIM==3
      int numqtwists(6*kgrid+3);
      std::vector<int> qk(0);
      mappedQtonofK.resize(numqtwists,qk);
      compQ.resize(numqtwists);
      RealType qn(4.0*M_PI*M_PI*std::pow(Lattice.Volume,-2.0/3.0));
      mappedQnorms.resize(numqtwists,qn*0.5/RealType(M));
      if (twist[0]==0)
        mappedQnorms[kgrid]=qn/RealType(M);
      if (twist[1]==0)
        mappedQnorms[3*kgrid+1]=qn/RealType(M);
      if (twist[2]==0)
        mappedQnorms[5*kgrid+2]=qn/RealType(M);
//             app_log()<<" Jnorm="<<qn<< std::endl;
      Q.resize(numqtwists);
      for (int i=-kgrid; i<(kgrid+1); i++)
      {
        PosType kpt;
        kpt[0]=i-twist[0];
        kpt[1]=i-twist[1];
        kpt[2]=i-twist[2];
        kpt=Lattice.k_cart(kpt);
        Q[i+kgrid]=std::abs(kpt[0]);
        Q[i+kgrid+(2*kgrid+1)]=std::abs(kpt[1]);
        Q[i+kgrid+(4*kgrid+2)]=std::abs(kpt[2]);
      }
      app_log()<<" Using all k-space points with (nx^2+ny^2+nz^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<< std::endl;
      app_log()<<"  My twist is:"<<twist[0]<<"  "<<twist[1]<<"  "<<twist[2]<< std::endl;
      int indx(0);
      int kgrid_squared=kgrid*kgrid;
      for (int i=-kgrid; i<(kgrid+1); i++)
      {
        for (int j=-kgrid; j<(kgrid+1); j++)
        {
          for (int k=-kgrid; k<(kgrid+1); k++)
          {
            if (i*i+j*j+k*k<=kgrid_squared) //if (std::sqrt(i*i+j*j+k*k)<=kgrid)
            {
              PosType kpt;
              kpt[0]=i-twist[0];
              kpt[1]=j-twist[1];
              kpt[2]=k-twist[2];
              //convert to Cartesian: note that 2Pi is multiplied
              kpt=Lattice.k_cart(kpt);
              kPoints.push_back(kpt);
              mappedQtonofK[i+kgrid].push_back(indx);
              mappedQtonofK[j+kgrid+(2*kgrid+1)].push_back(indx);
              mappedQtonofK[k+kgrid+(4*kgrid+2)].push_back(indx);
              indx++;
            }
          }
        }
      }
#endif
#if OHMMS_DIM==2
      int numqtwists(4*kgrid+2);
      std::vector<int> qk(0);
      mappedQtonofK.resize(numqtwists,qk);
      compQ.resize(numqtwists);
      RealType qn(2.0*M_PI/std::sqrt(Lattice.Volume));
      mappedQnorms.resize(numqtwists,qn*0.5/RealType(M));
      if (twist[0]==0)
        mappedQnorms[kgrid]=qn/RealType(M);
      if (twist[1]==0)
        mappedQnorms[3*kgrid+1]=qn/RealType(M);
//             app_log()<<" Jnorm="<<qn<< std::endl;
      Q.resize(numqtwists);
      for (int i=-kgrid; i<(kgrid+1); i++)
      {
        PosType kpt;
        kpt[0]=i-twist[0];
        kpt[1]=i-twist[1];
        kpt=Lattice.k_cart(kpt);
        Q[i+kgrid]=std::abs(kpt[0]);
        Q[i+kgrid+(2*kgrid+1)]=std::abs(kpt[1]);
      }
      app_log()<<" Using all k-space points with (nx^2+ny^2)^0.5 < "<< kgrid <<" for Momentum Distribution."<< std::endl;
      app_log()<<"  My twist is:"<<twist[0]<<"  "<<twist[1]<< std::endl;
      int indx(0);
      int kgrid_squared=kgrid*kgrid;
      for (int i=-kgrid; i<(kgrid+1); i++)
      {
        for (int j=-kgrid; j<(kgrid+1); j++)
        {
          if (i*i+j*j<=kgrid_squared) //if (std::sqrt(i*i+j*j+k*k)<=kgrid)
          {
            PosType kpt;
            kpt[0]=i-twist[0];
            kpt[1]=j-twist[1];
            //convert to Cartesian: note that 2Pi is multiplied
            kpt=Lattice.k_cart(kpt);
            kPoints.push_back(kpt);
            mappedQtonofK[i+kgrid].push_back(indx);
            mappedQtonofK[j+kgrid+(2*kgrid+1)].push_back(indx);
            indx++;
          }
        }
      }
#endif
    }
    kids=kids->next;
  }
  if (rootNode)
  {
    std::stringstream sstr;
    sstr<<"Kpoints";
    for(int i(0); i<OHMMS_DIM; i++)
      sstr<<"_"<<round(100.0*twist[i]);
    sstr<<".dat";
    std::ofstream fout(sstr.str().c_str());
    fout.setf(std::ios::scientific, std::ios::floatfield);
    fout << "# mag_k        ";
    for(int i(0); i<OHMMS_DIM; i++)
      fout << "k_"<<i<<"           ";
    fout << std::endl;
    for (int i=0; i<kPoints.size(); i++)
    {
      float khere(std::sqrt(dot(kPoints[i],kPoints[i])));
      fout<<khere;
      for(int j(0); j<OHMMS_DIM; j++)
        fout<<"   "<<kPoints[i][j];
      fout<< std::endl;
    }
    fout.close();
    sstr.str("");
    sstr<<"Qpoints";
    for(int i(0); i<OHMMS_DIM; i++)
      sstr<<"_"<<round(100.0*twist[i]);
    sstr<<".dat";
    std::ofstream qout(sstr.str().c_str());
    qout.setf(std::ios::scientific, std::ios::floatfield);
    qout << "# mag_q" << std::endl;
    for (int i=0; i<Q.size(); i++)
    {
      qout<<Q[i]<< std::endl;
    }
    qout.close();
  }
  nofK.resize(kPoints.size());
  norm_nofK=1.0/RealType(M);
  return true;
}

bool NaturalOrbitals::get(std::ostream& os) const
{
  return true;
}

QMCHamiltonianBase* NaturalOrbitals::makeClone(ParticleSet& qp
    , TrialWaveFunction& psi)
{
  NaturalOrbitals* myclone=new NaturalOrbitals(qp,psi);
  myclone->resize(kPoints,Q);
  myclone->myIndex=myIndex;
  myclone->kgrid=kgrid;
  myclone->norm_nofK=norm_nofK;
  myclone->mappedQtonofK.resize(mappedQtonofK.size());
  for(int i=0; i<mappedQtonofK.size(); i++)
    myclone->mappedQtonofK[i]=mappedQtonofK[i];
  myclone->hdf5_out=hdf5_out;
  myclone->mappedQnorms=mappedQnorms;
  myclone->M=M;
  return myclone;
}

void NaturalOrbitals::resize(const std::vector<PosType>& kin, const std::vector<RealType>& qin)
{
  //copy kpoints
  kPoints=kin;
  nofK.resize(kin.size());
  //copy q
  Q=qin;
  compQ.resize(qin.size());
}

void NaturalOrbitals::setRandomGenerator(RandomGenerator_t* rng)
{
  //simply copy it
  myRNG=*rng;
}
}

