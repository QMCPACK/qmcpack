//////////////////////////////////////////////////////////////////
// (c) Copyright 2003- by Jeongnim Kim
//////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////
//   Jeongnim Kim
//   jnkim@ornl.gov
//////////////////////////////////////////////////////////////////
// -*- C++ -*-
namespace qmcplusplus
{

/** branch engine for VMC using |Phi|^2
 */
struct VMCBranch
{
  inline VMCBranch(BranchBase& bb)
  {}

  /** return true if r*r is zero
   */
  inline bool phaseChanged(double r, double phi) const
  {
    return r*r<numeric_limits<double>::epsilon();
  }

  /** return weight=1
   */
  inline double weight(double enew, double eold) const
  {
    return 1.0;
  }
};

/** branch engine for fixed-node DMC or RMC */
struct FixedNodeCBranch
{
  const BranchBase& brancher;
  inline FixedNodeBranch(BranchBase& bb)
    : brancher(bb)
  {}

  /** return true if r*r is zero
   */
  inline bool phaseChanged(double r, double phi) const
  {
#if defined(QMC_COMPLEX)
    return false;
#else
    return std::cos(psi0) < numeric_limits<RealType>::epsilon();
#endif
  }

  /** return the DMC weight
   */
  inline double weight(double enew, double eold) const
  {
    return brancher.Branch(enew,eold);
  }
};

template<typename BB>
inline void QMCUpdate::advanceWalkerPbyP(Walker_t& thisWalker, int nbranch, int nsub)
{
  BB brancher(*branchEngine);
  Walker_t::Buffer_t& w_buffer(thisWalker.DataSet);
  W.loadWalker(thisWalker,true);
  Psi.copyFromBuffer(W,w_buffer);
  for(int ibranch=0; ibranch<nbranch; ++ibranch)
  {
    int nAcceptTemp(0);
    int nRejectTemp(0);
    RealType eold(thisWalker.Properties(LOCALENERGY));
    RealType enew(eold);
    RealType rr_proposed=0.0;
    RealType rr_accepted=0.0;
    RealType gf_acc=1.0;
    for(int iter=0; iter<nsub; ++iter)
    {
      makeGaussRandomWithEngine(deltaR,RandomGen);
      for(int ig=0; ig<W.groups(); ++ig) //loop over species
      {
        RealType tauovermass = Tau*MassInvS[ig];
        RealType oneover2tau = 0.5/(tauovermass);
        RealType sqrttau = std::sqrt(tauovermass);
        for (int iat=W.first(ig); iat<W.last(ig); ++iat)
        {
          //get the displacement
          GradType grad_iat=Psi.evalGrad(W,iat);
          PosType dr;
          getScaledDrift(tauovermass, grad_iat, dr);
          dr += sqrttau * deltaR[iat];
          RealType rr=tauovermass*dot(deltaR[iat],deltaR[iat]);
          rr_proposed+=rr;

          //singular gradient, reject
          if(rr>m_r2max)
          {
            ++nRejectTemp;
            continue;
          }

          //move is out of bound, reject
          if(!W.makeMoveAndCheck(iat,dr))
          {
            ++nRejectTemp;
            continue;
          }

          PosType newpos(W.R[iat]);
          RealType ratio = Psi.ratioGrad(W,iat,grad_iat);
          bool valid_move=false;
          //first check, if phase changed
          if (brancher.phaseChanged(ratio,Psi.getPhaseDiff()))
          {
            ++nRejectTemp;
            ++nNodeCrossing;
            W.rejectMove(iat);
            Psi.rejectMove(iat);
          }
          else
          {
            RealType logGf = -0.5*dot(deltaR[iat],deltaR[iat]);
            getScaledDrift(tauovermass, grad_iat, dr);
            dr = thisWalker.R[iat] - newpos - dr;
            RealType logGb = -oneover2tau*dot(dr,dr);
            RealType prob = ratio*ratio*std::exp(logGb-logGf);
            if(RandomGen() < prob)
            {
              valid_move=true;
              ++nAcceptTemp;
              W.acceptMove(iat);
              Psi.acceptMove(W,iat);
              rr_accepted+=rr;
              gf_acc *=prob;//accumulate the ratio
            }
            else
            {
              ++nRejectTemp;
              W.rejectMove(iat);
              Psi.rejectMove(iat);
            }
          }
        }//iat
      }//ig
    } //iter=substeps without energy evaluations

    if(UseTMove) nonLocalOps.reset();
    bool advanced=true;
    if(nAcceptTemp>0)
    {
      //need to overwrite the walker properties
      thisWalker.Age=0;
      thisWalker.R = W.R;
      RealType logpsi = Psi.updateBuffer(W,w_buffer,false);
      W.saveWalker(thisWalker);
      if(UseTMove)
        enew= H.evaluate(W,nonLocalOps.Txy);
      else
        enew= H.evaluate(W);
      thisWalker.resetProperty(logpsi,Psi.getPhase(),enew,rr_accepted,rr_proposed,1.0 );
      thisWalker.Weight *= brancher.weight(enew,eold);
      //this is always called at the end
      //H.auxHevaluate(W,thisWalker);
      H.saveProperty(thisWalker.getPropertyBase());
    }
    else
    {
      //all moves are rejected: does not happen normally with reasonable wavefunctions
      advanced=false;
      thisWalker.Age++;
      thisWalker.Properties(R2ACCEPTED)=0.0;
      H.rejectedMove(W,thisWalker);
      ++nAllRejected;
      enew=eold;//copy back old energy
      gf_acc=1.0;
      thisWalker.Weight *= brancher.weight(enew,eold);
    }

    if(UseTMove)
    {
      int ibar = nonLocalOps.selectMove(RandomGen());
      //make a non-local move
      if(ibar)
      {
        int iat=nonLocalOps.id(ibar);
        if(!W.makeMoveAndCheck(iat,nonLocalOps.delta(ibar)))
          continue;
        RealType ratio=Psi.ratio(W,iat,dG,dL);
        W.acceptMove(iat);
        Psi.acceptMove(W,iat);
        W.G += dG;
        W.L += dL;
        RealType logpsi = Psi.evaluateLog(W,w_buffer);
        W.saveWalker(thisWalker);
        //PAOps<RealType,OHMMS_DIM>::copy(W.G,thisWalker.Drift);
        ++NonLocalMoveAccepted;
      }
    }
    nAccept += nAcceptTemp;
    nReject += nRejectTemp;
  }

  //now compute all the other observables per walker
  H.auxHevaluate(W,thisWalker);
}

}

/***************************************************************************
 * $RCSfile: DMCUpdatePbyP.cpp,v $   $Author: mcminis2 $
 * $Revision: 3697 $   $Date: 2009-03-24 19:30:46 -0400 (Tue, 24 Mar 2009) $
 * $Id: DMCUpdatePbyP.cpp 3697 2009-03-24 23:30:46Z mcminis2 $
 ***************************************************************************/
