//////////////////////////////////////////////////////////////////
// (c) Copyright 2008-  by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   National Center for Supercomputing Applications &
//   Materials Computation Center
//   University of Illinois, Urbana-Champaign
//   Urbana, IL 61801
//   e-mail: jnkim@ncsa.uiuc.edu
//
// Supported by 
//   National Center for Supercomputing Applications, UIUC
//   Materials Computation Center, UIUC
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
#include <QMCHamiltonians/MomentumEstimator.h>
#include <QMCWaveFunctions/TrialWaveFunction.h>
#include <Numerics/e2iphi.h>
#include <Numerics/OhmmsBlas.h>
#include <OhmmsData/AttributeSet.h>
#include <Utilities/SimpleParser.h>
#include <Particle/DistanceTableData.h>
#include <Numerics/DeterminantOperators.h>

namespace qmcplusplus 
{

  MomentumEstimator::MomentumEstimator(ParticleSet& elns, TrialWaveFunction& psi)
    :M(4), refPsi(psi), kgrid(4), Lattice(elns.Lattice)
  {
    UpdateMode.set(COLLECTABLE,1);
    psi_ratios.resize(elns.getTotalNum());
    kdotp.resize(elns.getTotalNum());
    phases.resize(elns.getTotalNum());
  }

  void MomentumEstimator::resetTargetParticleSet(ParticleSet& P)
  {
  }

  MomentumEstimator::Return_t MomentumEstimator::evaluate(ParticleSet& P)
  {
    const int np=P.getTotalNum();
    nofK=0.0;
    compQ=0.0;

    //will use temp[i].r1 for the Compton profile
    const vector<DistanceTableData::TempDistType>& temp(P.DistTables[0]->Temp);
    for(int s=0; s<M; ++s)
    {
      PosType newpos;
      for(int i=0; i<OHMMS_DIM;++i) newpos[i]=myRNG();
      //make it cartesian
      newpos=Lattice.toCart(newpos);
      P.makeVirtualMoves(newpos); //updated: temp[i].r1=|newpos-P.R[i]|, temp[i].dr1=newpos-P.R[i]
      refPsi.get_ratios(P,psi_ratios);
      P.rejectMove(0); //restore P.R[0] to the orginal position

      ////debug get_ratios with ratio, use it whenever an OrbitalBase implements get_ratios
      //vector<RealType> r_org(np);
      //for(int i=0; i<np; ++i)
      //{
      //  PosType delta=newpos-P.R[i];
      //  P.makeMove(i,delta);
      //  r_org[i]=refPsi.ratio(P,i);
      //  P.rejectMove(i);
      //  cout << "ratio("<<i<<")=" << r_org[i] << " diff=" << r_org[i]-psi_ratios[i] << endl;
      //}
      //APP_ABORT("Done with test");

      for(int ik=0; ik < kPoints.size(); ++ik)
      {
        RealType kdotp_primed=dot(kPoints[ik],newpos);
        for(int i=0; i<np; ++i) kdotp[i]=kdotp_primed-dot(kPoints[ik],P.R[i]);
	//this is the same
        //for(int i=0; i<np; ++i) kdotp[i]=dot(kPoints[ik],temp[i].dr1);
        eval_e2iphi(np,kdotp.data(),phases.data());
        nofK[ik]+=real(BLAS::dot(np,phases.data(),psi_ratios.data()));
      }
//spherically averaged compton profile
      for(int i=0; i<np; ++i) psi_ratios[i] *= temp[i].rinv1*temp[i].rinv1;
      for(int iq=0; iq < Q.size(); ++iq)
      {
        for(int i=0; i<np; ++i) kdotp[i]=Q[iq]*temp[i].r1;
        eval_e2iphi(np,kdotp.data(),phases.data());
	compQ[iq]+=real(BLAS::dot(np,phases.data(),psi_ratios.data()));
      }
    }
if (hdf5_out)
{
//need normalization factor
     int j=myIndex;
     for(int ik=0; ik<nofK.size(); ++ik,++j) P.Collectables[j]+=norm_nofK*nofK[ik];
     for(int iq=0; iq<compQ.size(); ++iq,++j) P.Collectables[j]+=norm_compQ*compQ[iq];
}
else
{
    for(int ik=0; ik<nofK.size(); ++ik) nofK[ik] *= norm_nofK/RealType(M);
    for(int iq=0; iq<compQ.size(); ++iq) compQ[iq] *= norm_compQ/RealType(M);
}
    return 0.0;
  }

  void MomentumEstimator::registerCollectables(vector<observable_helper*>& h5desc
      , hid_t gid) const
  {
    if(hdf5_out)
    {
    //descriptor for the data, 1-D data
    vector<int> ng(1);

    //add nofk
    ng[0]=nofK.size();
    observable_helper* h5o=new observable_helper("nofk");
    h5o->set_dimensions(ng,myIndex);
    h5o->open(gid);
    h5o->addProperty(const_cast<vector<PosType>&>(kPoints),"kpoints");
    h5o->addProperty(const_cast<vector<int>&>(kWeights),"kweights");
    h5desc.push_back(h5o);

    //add compQ
    ng[0]=Q.size();
    h5o=new observable_helper("compQ");
    h5o->set_dimensions(ng,myIndex+nofK.size());
    h5o->open(gid);
    h5o->addProperty(const_cast<vector<RealType>&>(Q),"q");
    h5desc.push_back(h5o);
    }
  }


  void MomentumEstimator::addObservables(PropertySetType& plist, BufferType& collectables)
  {
    if (!hdf5_out)
    {
      myIndex=plist.size();
      for(int i=0;i<kgrid;i++) for(int j=0;j<kgrid;j++) for(int k=0;k<kgrid;k++)
      {
	  std::stringstream sstr;
	  sstr << "P_" <<i<< "_" << j<<"_"<<k;
	  int id=plist.add(sstr.str());
      }
      for(int i=0;i<Q.size();i++)
      {
	std::stringstream sstr;
	  sstr << "Q_" <<i;
	  int id=plist.add(sstr.str());
      }
    }
    else
    {
      myIndex=collectables.size();
      collectables.add(nofK.begin(),nofK.end());
      collectables.add(compQ.begin(),compQ.end());
    }
  }
  


  void MomentumEstimator::setObservables(PropertySetType& plist)
  {
    //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myIndex);
    if (!hdf5_out)
    {
      std::copy(nofK.begin(),nofK.end(),plist.begin()+myIndex);
      std::copy(compQ.begin(),compQ.end(),plist.begin()+myIndex+nofK.size());
    }
    
  }

  void MomentumEstimator::setParticlePropertyList(PropertySetType& plist
      , int offset)
  {
    //std::copy(gof_r.first_address(),gof_r.last_address(),plist.begin()+myDebugIndex+offset);
  }

  bool MomentumEstimator::putSpecial(xmlNodePtr cur, ParticleSet& elns)
  {
    string ctype("scalar");
    OhmmsAttributeSet pAttrib;
    pAttrib.add(ctype,"mode");
    pAttrib.put(cur);
    if (ctype=="hdf5") hdf5_out=true;
    else hdf5_out=false;
    
    xmlNodePtr kids=cur->children;
    while(kids!=NULL)
    {
      string cname((const char*)(kids->name));
      if (cname=="kpoints")
      {
	string ctype("manual");
	OhmmsAttributeSet pAttrib;
	pAttrib.add(ctype,"mode");
	pAttrib.add(kgrid,"grid");
        pAttrib.put(kids);
	      
// 	if(ctype=="manual")
// 	{
// 	  vector<RealType> kpt_unsorted;
// 	  putContent(kpt_unsorted,kids);
// 	  if(kpt_unsorted.size()%4!=0)
// 	  {
// 	    app_log()<<"Format for K points is \"kx ky kz wgt\". "<<endl;
// 	    APP_ABORT("MomentumEstimator::put");
// 	  }
// 	  int nkpts=kpt_unsorted.size()/4;
// 	  vector<PosType> ktmp(nkpts);
// 	  vector<int> kwgt(nkpts);
// 	  for(int i=0,j=0;i<nkpts;i++)
// 	  {
// 	    ktmp[i][0]=kpt_unsorted[j++];
// 	    ktmp[i][1]=kpt_unsorted[j++];
// 	    ktmp[i][2]=kpt_unsorted[j++];
// 	    kwgt[i]=kpt_unsorted[j++];
// 	  }
// 	  kPoints=ktmp;
// 	  kWeights=kwgt;
// 	}
// 	else if(ctype=="auto")
// 	{
	 Matrix<RealType> BasisMatrix(3);
	 
	 ParticleSet::ParticlePos_t R_cart(1);
         R_cart.setUnit(PosUnit::CartesianUnit);
         ParticleSet::ParticlePos_t R_unit(1);
         R_unit.setUnit(PosUnit::LatticeUnit);
      
	 for (int i=0;i<3;i++)
	  {
	    R_unit[0][0]=0;
	    R_unit[0][1]=0;
	    R_unit[0][2]=0;
	    R_unit[0][i]=1;
	    elns.convert2Cart(R_unit,R_cart);
	    app_log()<<R_cart[0][0]<<" "<<R_cart[0][1]<<" "<<R_cart[0][2]<<endl;
	    for (int j=0;j<3;j++) BasisMatrix(i,j)= R_cart[0][j];
	  }
	  invert_matrix(BasisMatrix, false);
	 app_log()<<" Reciprocal Matrix "<<endl;
	 for (int i=0;i<3;i++) 
	  {
	    for (int j=0;j<3;j++)
	    {
	      app_log()<<BasisMatrix(i,j)<<" ";
	    }
	    app_log()<<endl;
	  }
	 
	  app_log()<<" Using a "<<kgrid<<"x"<< kgrid<<"x"<< kgrid << " cubic grid in k-space for Momentum Distribution."<<endl;
	  for(int i=0;i<kgrid;i++) for(int j=0;j<kgrid;j++) for(int k=0;k<kgrid;k++)
	  {
	    PosType kpt;
	    
	    kpt[0]=RealType(i)*BasisMatrix(0,0)+RealType(j)*BasisMatrix(1,0)+RealType(k)*BasisMatrix(2,0);
	    kpt[1]=RealType(i)*BasisMatrix(0,1)+RealType(j)*BasisMatrix(1,1)+RealType(k)*BasisMatrix(2,1);
	    kpt[2]=RealType(i)*BasisMatrix(0,2)+RealType(j)*BasisMatrix(1,2)+RealType(k)*BasisMatrix(2,2);
	    for(int l=0;l<3;l++) kpt[l]=2.0*M_PI*kpt[l];
	    kPoints.push_back(kpt);
	    kWeights.push_back(1);
	  }
// 	}
      }
      else if(cname=="q")
      {
	string ctype("manual");
	RealType dx(0.1);
	RealType maxX(5.0);
	OhmmsAttributeSet pAttrib;
	pAttrib.add(ctype,"mode");
	pAttrib.add(dx,"dx");
	pAttrib.add(maxX,"max");
        pAttrib.put(kids);
	if(ctype=="manual")
	{
	  vector<RealType> q_unsorted;
	  putContent(q_unsorted,kids);
	  Q=q_unsorted;
	}
	else if(ctype=="auto")
	{
	  RealType X=0;
	  Q.push_back(X);
	  while(X<maxX)
	  {
	    X+=dx;
	    Q.push_back(X);
	  }
	}
      }
//       else if(cname=="qpoints")
//       {
//       }
      kids=kids->next;
    }
    nofK.resize(kPoints.size());
    compQ.resize(Q.size());

    norm_nofK=1.0;// /elns.Lattice.Volume;
    norm_compQ=1.0;//4.0*M_PI*M_PI;// /elns.Lattice.Volume;

    for(int i=0;i<kPoints.size();i++)
    {
      app_log()<<kPoints[i][0]<<" "<<kPoints[i][1]<<" "<<kPoints[i][2]<<" "<<std::sqrt(kPoints[i][0]*kPoints[i][0]+kPoints[i][1]*kPoints[i][1]+kPoints[i][2]*kPoints[i][2])<<endl;
    }

    return true;
  }

  bool MomentumEstimator::get(std::ostream& os) const
  {
    return true;
  }

  QMCHamiltonianBase* MomentumEstimator::makeClone(ParticleSet& qp
      , TrialWaveFunction& psi)
  {
    MomentumEstimator* myclone=new MomentumEstimator(qp,psi);
    myclone->resize(kPoints,Q,kWeights);
    myclone->norm_nofK=norm_nofK;
    myclone->norm_compQ=norm_compQ;
    myclone->myIndex=myIndex;
    myclone->kgrid=kgrid;
    return myclone;
  }

  void MomentumEstimator::resize(const vector<PosType>& kin, const vector<RealType>& qin, const vector<int>& win)
  {
    //copy kpoints
    kPoints=kin;
    nofK.resize(kin.size());
  
    //copy q
    Q=qin;
    compQ.resize(qin.size());
    
    kWeights=win;
  }

  void MomentumEstimator::setRandomGenerator(RandomGenerator_t* rng)
  {
    //simply copy it
    myRNG=*rng;
  }
}

/***************************************************************************
 * $RCSfile$   $Author: jnkim $
 * $Revision: 2945 $   $Date: 2008-08-05 10:21:33 -0500 (Tue, 05 Aug 2008) $
 * $Id: ForceBase.h 2945 2008-08-05 15:21:33Z jnkim $ 
 ***************************************************************************/
