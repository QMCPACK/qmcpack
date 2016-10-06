//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    



/*
 * Reptile.h
 *
 * This class is a wrapper/utility class for groups of walkers within
 * MCWalkerConfiguration for Reptation Monte Carlo.
 *
 * Walkers are stored as usual in MCWalkerConfiguration.  However, this interface
 * represents the reptile as a circular queue within a segment of the walker list.
 */

#ifndef QMCPLUSPLUS_REPTILE_H
#define QMCPLUSPLUS_REPTILE_H

#include "QMCDrivers/DriftOperators.h"
#include "Configuration.h"
#include "Walker.h"

namespace qmcplusplus
{
  class MCWalkerConfiguration;

  class Reptile:public QMCTraits
  {


  public:
    typedef MCWalkerConfiguration::Walker_t Walker_t;
    //typedef Walker_t::Buffer_t              Buffer_t;
    //    typedef MCWalkerConfiguration::Walker_t Walker_t;
    typedef MCWalkerConfiguration::iterator WalkerIter_t;
    typedef std::vector< Walker_t::ParticlePos_t > ReptileConfig_t;

      std::vector < IndexType > Action;
      std::vector < IndexType > TransProb;

    RealType forwardprob;
    RealType backwardprob;
    RealType forwardaction;
    RealType backwardaction;

    RealType tau;

    IndexType direction, headindex, nbeads;
      MCWalkerConfiguration & w;
    Walker_t *prophead;
    WalkerIter_t repstart, repend;

    inline Reptile (MCWalkerConfiguration & W, WalkerIter_t start, WalkerIter_t end):w (W), repstart (start), repend (end), direction (1), headindex (0), prophead (0)	//, r2prop(0.0), r2accept(0.0),tau(0.0)
    {
      w = W;
      Action.resize (3);
      Action[0] = w.addProperty ("ActionBackward");
      Action[1] = w.addProperty ("ActionForward");
      Action[2] = w.addProperty ("ActionLocal");
      TransProb.resize (2);
      TransProb[0] = w.addProperty ("TransProbBackward");
      TransProb[1] = w.addProperty ("TransProbForward");

      nbeads = repend - repstart;
    }

     ~Reptile ()
    {
    }

    inline IndexType size ()
    {
      return nbeads;
    }

    inline Walker_t & operator[] (IndexType i)
    {
      return getWalker (getBeadIndex (i));
    }

    inline IndexType wrapIndex (IndexType repindex)
    {
      return (repindex % nbeads + nbeads) % nbeads;
    }

    inline Walker_t & getWalker (IndexType i)
    {
      WalkerIter_t bead = repstart + wrapIndex (i);
      return **bead;
    }

    inline IndexType getBeadIndex (IndexType i)
    {
      return wrapIndex (headindex + direction * i);
    }
    inline Walker_t & getBead (IndexType i)
    {
      return getWalker (getBeadIndex (i));
    }
    inline Walker_t & getHead ()
    {
      return getWalker (getBeadIndex (0));
    }
    inline Walker_t & getTail ()
    {
      return getWalker (getBeadIndex (nbeads - 1));
    }
    inline Walker_t & getNext ()
    {
      return getWalker (getBeadIndex (nbeads - 2));
    }
    inline Walker_t & getCenter ()
    {
      return getWalker (getBeadIndex ((nbeads - 1) / 2));
    }
    //inline void setProposedHead(){

    inline void flip ()
    {
      // direction*=-1;
      // headindex = getBeadIndex(nbeads-1);
      headindex = wrapIndex (headindex - direction);
      direction *= -1;
    }

    inline void setDirection (IndexType dir)
    {
      direction = dir;
    }

    inline void setBead (Walker_t & walker, IndexType i)
    {
      IndexType index = getBeadIndex (i);
      Walker_t & newbead (getWalker (index));
      newbead = walker;		//This should be a hard copy
    }

    inline void setHead (Walker_t & overwrite)
    {
      //overwrite last element.
      headindex = getBeadIndex (nbeads - 1);	//sets to position of tail.
      Walker_t & newhead (getBead (0));
      newhead = overwrite;
    }
    //This function does two things:  1.)  Moves the reptile forward 1 step.  2.) Returns the new head.  
    inline Walker_t & getNewHead ()
    {
      //overwrite last element.
      headindex = getBeadIndex (nbeads - 1);	//sets to position of tail.
      return getWalker (headindex);
    }

    void saveAction (Walker_t & walker, IndexType d, RealType val,
		     IndexType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType actionindex = 2;
      if (direction != 0)
	actionindex = (1 - d * direction) / 2;
      walker.Properties (nPsi, Action[actionindex]) = val;
    }

    RealType getDirectionalAction (Walker_t & walker, IndexType d,
				   IndexType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType actionindex = 2;
      if (d != 0)
	actionindex = (1 - direction * d) / 2;

      return walker.Properties (nPsi, Action[actionindex]);
    }

    RealType getLinkAction (Walker_t & new_walker, Walker_t & old_walker,
			    IndexType d, IndexType nPsi = 0)
    {
      RealType af = getDirectionalAction (old_walker, +1, nPsi);
      RealType ab = getDirectionalAction (new_walker, -1, nPsi);
      RealType a0 =
	getDirectionalAction (old_walker, 0,
			      nPsi) + getDirectionalAction (new_walker, 0,
							    nPsi);
      return af + ab + a0;

    }

    void saveTransProb (Walker_t & walker, IndexType d, RealType val,
			IndexType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType transindex = (1 - d * direction) / 2;
      walker.Properties (nPsi, TransProb[transindex]) = val;
    }

    void saveTransProb (ParticleSet & W, IndexType d, RealType val,
			IndexType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType transindex = (1 - d * direction) / 2;
      W.Properties (nPsi, TransProb[transindex]) = val;
    }
    RealType getTransProb (Walker_t & walker, IndexType d, RealType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType transindex = (1 - d * direction) / 2;
      return walker.Properties (nPsi, TransProb[transindex]);
    }
    RealType getTransProb (ParticleSet & W, IndexType d, RealType nPsi = 0)
    {
      //IndexType repdirection=circbuffer.get_direction();
      IndexType transindex = (1 - d * direction) / 2;
      return W.Properties (nPsi, TransProb[transindex]);
    }

    inline void printState ()
    {
      app_log () << "********PRINT REPTILE STATE*********\n";
      app_log () << "Direction=" << direction << "  Headindex=" << headindex
	<< "  tail=" << getBeadIndex (nbeads -
				      1) << "\n  next=" <<
	getBeadIndex (nbeads - 2) << "  nbeads=" << nbeads << std::endl;
      app_log () <<
	"BeadIndex\tWrapIndex\tEnergy\tAction[0]\tAction[1]\tAction[2]\t\n";
      for (int i = 0; i < nbeads; i++)
	{
	  app_log () << i << "\t" << getBeadIndex (i) << "\t" << getBead (i).
	    Properties (LOCALENERGY) << "\t" << getBead (i).
	    Properties (Action[0]) << "\t" << getBead (i).
	    Properties (Action[1]) << "\t" << getBead (i).
	    Properties (Action[2]) << "\n";
	}
      app_log () << "POSITIONS===============:\n";
      for (int i = 0; i < nbeads; i++)
	{
	  //  app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).R[0]<<"\n";
	  //  app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).R[1]<<"\n";
	  app_log () << "BEAD #" << i << " tau = " << tau * i << std::endl;
	  app_log () << getBead (i).R << std::endl;
	}
      app_log () << "GVECS===============:\n";
      for (int i = 0; i < nbeads; i++)
	{
//      app_log()<<i<<"\t1"<<1<<"\t"<<getBead(i).G[0]<<"\n";
//      app_log()<<i<<"\t2"<<2<<"\t"<<getBead(i).G[1]<<"\n";
	  app_log () << "BEAD #" << i << " tau = " << tau * i << std::endl;
	  app_log () << getBead (i).G << std::endl;
	}
      app_log () << "************************************\n";
    }
    inline RealType getTau ()
    {
      return tau;
    }
    inline void setTau (RealType t)
    {
      tau = t;
    }


    //This takes a value of imaginary time "t" and returns a 3N particle position vector, corresponding to a time slice extrapolated
    // from the current reptile.  If t>length of reptile, then return the last bead.  if t<0; return the first bead. 
    inline Walker_t::ParticlePos_t linearInterp (RealType t)
    {
      IndexType nbead = IndexType (t / tau);	//Calculate the lower bound on the timeslice.  t is between binnum*Tau and (binnum+1)Tau
      RealType beadfrac = t / tau - nbead;	//the fractional coordinate between n and n+1 bead
      if (nbead <= 0)
	{
	  ParticleSet::ParticlePos_t result = getHead ().R;
	  return result;
	}
      else if (nbead >= nbeads - 1)
	{
	  ParticleSet::ParticlePos_t result = getTail ().R;
	  return result;
	}

      else
	{
	  Walker_t::ParticlePos_t dR (getBead (nbead + 1).R),
	    interpR (getBead (nbead).R);
	  dR = dR - getBead (nbead).R;

	  interpR = getBead (nbead).R + beadfrac * dR;
	  return interpR;
	}

    }
    inline ReptileConfig_t getReptileSlicePositions (RealType tau,
						     RealType beta)
    {
      IndexType nbeads_new = IndexType (beta / tau);
      ReptileConfig_t new_reptile_coords (0);

      for (IndexType i = 0; i < nbeads_new; i++)
	new_reptile_coords.push_back (linearInterp (tau * i));

      return new_reptile_coords;
    }

    inline void setReptileSlicePositions (ReptileConfig_t & rept)
    {
      if (rept.size () == nbeads)
	{
	  for (int i = 0; i < nbeads; i++)
	    getBead (i).R = rept[i];
	}
      else;
    }

    inline void setReptileSlicePositions (Walker_t::ParticlePos_t R)
    {
      for (int i = 0; i < nbeads; i++)
	getBead (i).R = R;
    }

  };


}
#endif
