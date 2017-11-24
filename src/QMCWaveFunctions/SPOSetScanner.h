//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////
    

#ifndef QMCPLUSPLUS_SPOSET_SCANNER_H
#define QMCPLUSPLUS_SPOSET_SCANNER_H

#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "QMCWaveFunctions/SPOSetBase.h"
#include "OhmmsData/AttributeSet.h"

namespace qmcplusplus
{

  /** a scanner for all the SPO sets.
   */
  class SPOSetScanner
  {
  public:
    typedef std::map<std::string,ParticleSet*>         PtclPoolType;
    typedef std::map<std::string,SPOSetBasePtr>        SPOMapType;
    typedef QMCTraits::RealType                        RealType;
    typedef QMCTraits::ValueType                       ValueType;
    typedef OrbitalSetTraits<ValueType>::ValueVector_t ValueVector_t;
    typedef OrbitalSetTraits<ValueType>::GradVector_t  GradVector_t;
    typedef OrbitalSetTraits<ValueType>::HessVector_t  HessVector_t;
    
    RealType myfabs(RealType s) { return std::fabs(s); }
    template<typename T> std::complex<T> myfabs(std::complex<T>& s) { return std::complex<T>(myfabs(s.real()),myfabs(s.imag())); }
    template<typename T> TinyVector<T, OHMMS_DIM> myfabs(TinyVector<T, OHMMS_DIM>& s) { return TinyVector<T, OHMMS_DIM>(myfabs(s[0]),myfabs(s[1]),myfabs(s[2])); }

    SPOMapType& SPOMap;
    PtclPoolType& PtclPool;
    ParticleSet* ions;
    ParticleSet& target;

    // construction/destruction
    SPOSetScanner(SPOMapType& spomap, ParticleSet& targetPtcl, PtclPoolType& psets): SPOMap(spomap), target(targetPtcl), PtclPool(psets), ions(0){};
    //~SPOSetScanner(){}; 

    // processing scanning
    void put(xmlNodePtr cur)
    {
      app_log() << "Entering the SPO set scanner!" << std::endl;
      // check in the source particle set and search for it in the pool.
      std::string sourcePtcl("ion0");
      OhmmsAttributeSet aAttrib;
      aAttrib.add(sourcePtcl,"source");
      aAttrib.put(cur);
      PtclPoolType::iterator pit(PtclPool.find(sourcePtcl));
      if(pit == PtclPool.end())
        app_log() << "Source particle set not found. Can not be used as reference point." << std::endl;
      else
        ions=(*pit).second;

      // scanning the SPO sets
      xmlNodePtr cur_save=cur;
      SPOSetBasePtr mySPOSet;
      for (SPOMapType::iterator spo_iter=SPOMap.begin(); spo_iter!=SPOMap.end(); spo_iter++)
      {
        app_log() << "  Processing SPO " << spo_iter->first << std::endl;
        mySPOSet = spo_iter->second;
        // scanning the paths
        cur = cur_save->children;
        while (cur != NULL)
        {
          std::string trace_name("no name");
          OhmmsAttributeSet aAttrib;
          aAttrib.add(trace_name,"name");
          aAttrib.put(cur);
          std::string cname;
          getNodeName(cname,cur);
          std::string prefix(spo_iter->first+"_"+cname+"_"+trace_name);
          if (cname == "path")
          {
            app_log() << "    Scanning a " << cname << " called " << trace_name << " and writing to " << prefix+"_v/g/l/report.dat" << std::endl;
            scan_path(cur, mySPOSet, prefix);
          }
          else
          {
            if(cname !="text" && cname !="comment") app_log() << "    Unknown type of scanning " << cname << std::endl;
          }
          cur = cur->next;    
        }
      }
      app_log() << "Exiting the SPO set scanner!" << std::endl << std::endl;
    }

    // scanning a path
    void scan_path(xmlNodePtr cur, SPOSetBasePtr mySPOSet, std::string prefix)
    {
      std::string file_name;
      file_name=prefix+"_v.dat";
      std::ofstream output_v(file_name.c_str());
      file_name=prefix+"_g.dat";
      std::ofstream output_g(file_name.c_str());
      file_name=prefix+"_l.dat";
      std::ofstream output_l(file_name.c_str());
      file_name=prefix+"_report.dat";
      std::ofstream output_report(file_name.c_str());

      int nknots(2);
      int from_atom(-1);
      int to_atom(-1);
      TinyVector<double,OHMMS_DIM> from_pos(0.0, 0.0, 0.0);
      TinyVector<double,OHMMS_DIM> to_pos(0.0, 0.0, 0.0);
      
      OhmmsAttributeSet aAttrib;
      aAttrib.add(nknots,"nknots");
      aAttrib.add(from_atom,"from_atom");
      aAttrib.add(to_atom,"to_atom");
      aAttrib.add(from_pos,"from_pos");
      aAttrib.add(to_pos,"to_pos");
      aAttrib.put(cur);

      // sanity check
      if ( nknots < 2 ) nknots=2;
      // check out the reference atom coordinates
      if (ions)
      {
        if ( from_atom >= 0 && from_atom < ions->R.size() )
          from_pos = ions->R[from_atom];
        if ( to_atom >= 0 && to_atom < ions->R.size() )
          to_pos = ions->R[to_atom];
      }

      // prepare a fake particle set
      ValueVector_t SPO_v, SPO_l, SPO_v_avg, SPO_l_avg;
      GradVector_t SPO_g, SPO_g_avg;
      int OrbitalSize(mySPOSet->size());
      SPO_v.resize(OrbitalSize);
      SPO_g.resize(OrbitalSize);
      SPO_l.resize(OrbitalSize);
      SPO_v_avg.resize(OrbitalSize);
      SPO_g_avg.resize(OrbitalSize);
      SPO_l_avg.resize(OrbitalSize);
      SPO_v_avg=0.0;
      SPO_g_avg=0.0;
      SPO_l_avg=0.0;
      double Delta= 1.0/(nknots-1);
      int elec_count=target.R.size();
      ParticleSet::SingleParticlePos_t zero_pos(0.0,0.0,0.0);
      for(int icount=0, ind=0; icount<nknots; icount++, ind++)
      {
        if ( ind == elec_count ) ind = 0 ;
        target.R[ind][0] = (to_pos[0]-from_pos[0]) * Delta * icount + from_pos[0];
        target.R[ind][1] = (to_pos[1]-from_pos[1]) * Delta * icount + from_pos[1];
        target.R[ind][2] = (to_pos[2]-from_pos[2]) * Delta * icount + from_pos[2];
        target.makeMoveAndCheck(ind, zero_pos);
        mySPOSet->evaluate(target, ind, SPO_v, SPO_g, SPO_l);
        std::ostringstream o;
        o << "x_y_z  " << std::fixed << std::setprecision(7) << target.R[ind][0] << " " << target.R[ind][1] << " " << target.R[ind][2] ;
        output_v << o.str() << " : "  << std::scientific << std::setprecision(12);
        output_g << o.str() << " : "  << std::scientific << std::setprecision(12);
        output_l << o.str() << " : "  << std::scientific << std::setprecision(12);
        for(int iorb=0; iorb<OrbitalSize; iorb++)
        {
          SPO_v_avg[iorb] += myfabs(SPO_v[iorb]);
          SPO_g_avg[iorb] += myfabs(SPO_g[iorb]);
          SPO_l_avg[iorb] += myfabs(SPO_l[iorb]);
          output_v << SPO_v[iorb] << "  ";
          output_g << SPO_g[iorb][0] << "  " << SPO_g[iorb][1] << "  " << SPO_g[iorb][2] << "  ";
          output_l << SPO_l[iorb] << "  ";
        }
        output_v << std::endl;
        output_g << std::endl;
        output_l << std::endl;
      }
#ifdef QMC_COMPLEX
      output_report << "#   Report: Orb   Value_avg I/R  Gradients_avg  Laplacian_avg" << std::endl;
#else
      output_report << "#   Report: Orb   Value_avg   Gradients_avg   Laplacian_avg" << std::endl;
#endif
      for(int iorb=0; iorb<OrbitalSize; iorb++)
        output_report << "\t" << iorb << "    " << std::scientific << SPO_v_avg[iorb]*static_cast<RealType>(1.0/nknots) << "   "
#ifdef QMC_COMPLEX
                  << SPO_v_avg[iorb].imag()/SPO_v_avg[iorb].real() << "   "
#endif
                  << SPO_g_avg[iorb][0]*static_cast<RealType>(1.0/nknots) << "   " << SPO_g_avg[iorb][1]*static_cast<RealType>(1.0/nknots) << "   "
                  << SPO_g_avg[iorb][2]*static_cast<RealType>(1.0/nknots) << "   " << SPO_l_avg[iorb]*static_cast<RealType>(1.0/nknots)
                  << std::fixed << std::endl;
      output_v.close();
      output_g.close();
      output_l.close();
      output_report.close();
    }
  };
}

#endif
