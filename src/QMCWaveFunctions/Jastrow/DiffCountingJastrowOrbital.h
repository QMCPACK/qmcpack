
#ifndef QMCPLUSPLUS_DIFF_COUNTING_JASTROW_ORBITAL_H
#define QMCPLUSPLUS_DIFF_COUNTING_JASTROW_ORBITAL_H

#include "Particle/ParticleSet.h"
#include "QMCWaveFunctions/Jastrow/CountingRegion.h"

#include "QMCWaveFunctions/DiffWaveFunctionComponent.h"

namespace qmcplusplus{

template <class RegionType>
class DiffCountingJastrowOrbital: public DiffWaveFunctionComponent
{
  // variables

public:
  // constructor
  DiffCountingJastrowOrbital(ParticleSet& P)
  {
  }

  // extended checks and array resizing
  void initialize()
  {
  }
  
  void checkOutVariables(const opt_variables_type& optvars)
  {
  }

  void resetTargetParticleSet(ParticleSet& P)
  {
    // 'prepare internal data for a new particle set'
  }

  void resetParameters(const opt_variables_type& optvars)
  {
  }

//  // derivatives of single-particle update for ecp nonlocal quadrature
//  void CountingJastrowOrbital::evaluateTempDerivatives(ParticleSet& P, 
//                           const opt_variables_type& active, 
//                           RealType& ratioval,
//                           std::vector<RealType>& dlogpsi_t,
//                           int iat,
//                           PosType dr)
//  {
//    // assume that current state is determined by having called
//    // evaluateDerivatives(P,...) immediately before this function
//    P.makeMoveAndCheck(iat,dr);
//    ratioval = ratio(P,iat);
//    // all non-temp variables are set to values associated with position P
//    // all temp (_t) variables are set to values for moved position 
//    // evaluate log of F parameter derivs at moved position
//  
//    if(opt_F)
//    {
//      for(int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
//      {
//  
//        std::string id = opt_id[OPT_F][oi];
//        int ia = myVars.getIndex(id);
//        if(ia == -1)
//          continue; // ignore inactive parameters
//        int IJ = opt_index[OPT_F][oi];
//        int I = IJ/num_regions;
//        int J = IJ%num_regions;
//        // coefficient due to symmetry of F: \sum\limits_{I} F_{II} C_I^2 + \sum\limits_{J > I} 2 F_{IJ}*C_I*C_J
//        RealType x = (I==J)?1:2;
//        RealType dJF_val = x*(C->sum_t(I)*C->sum_t(J));
//        dlogpsi_t[ia] += dJF_val;
//      }
//    }
//    // evaluate partial derivatives of G at moved position
//    if(opt_G)
//    {
//      for(int oi = 0; oi < opt_index[OPT_G].size(); ++oi)
//      {
//        std::string id = opt_id[OPT_G][oi];
//        int ia = myVars.getIndex(id);
//        if(ia == -1)
//          continue; // ignore inactive params
//        int I = opt_index[OPT_G][oi];
//        RealType dJG_val = C->sum_t(I);
//        dlogpsi_t[ia] += dJG_val;
//      }
//    }
//  
//    if(opt_C)
//    {
//      // difference; easier to calculate than absolute values
//      static std::vector<RealType> dCdiff;
//      static int max_num_derivs = C->max_num_derivs();
//      dCdiff.resize(max_num_derivs*num_regions);
//      // easy-index functions for evaluateDerivatives calls
//      std::function<RealType&(int,int)> _dCsum  = [&](int I, int p)->RealType&{ return dCsum[p*num_regions + I]; }; 
//      std::function<RealType&(int,int)> _dCdiff = [&](int I, int p)->RealType&{ return dCdiff[p*num_regions + I]; }; 
//      // pointer to C->C[I]->myVars.Index
//      // for pI in { 0 .. C->num_derivs(I) }
//      //   dCindex->[pI]  is the index that corresponds to this parameter in active.
//      //   i.e., active[dCindex->[pI]] <=> C->C[I]->myVars.Index[pI]
//      std::fill(dCdiff.begin(), dCdiff.end(), 0);
//      for(int I = 0; I < num_regions; ++I)
//      {
//        // get the number of active parameters for the Ith counting region
//        opt_variables_type I_vars = C->getVars(I); 
//        int I_num_derivs = I_vars.size();
//        // evaluateTempDerivatives increments difference of derivative to dCdiff 
//        C->evaluateTempDerivatives(P, I, iat, _dCdiff);
//        // loop over parameters for the Ith counting function
//        for(int pI = 0; pI < I_num_derivs; ++pI)
//        {
//          // index for active optimizable variables
//          int ia = I_vars.Index[pI];
//          if(ia == -1)
//            continue; // ignore inactive
//          for(int J = 0; J < num_regions; ++J)
//          {
//            dlogpsi_t[ia] += (_dCsum(J,pI) + _dCdiff(J,pI))*(2*FCsum_t[J] + G[J]);
//          }
//        }
//      }
//  
//    } // end opt_C
//  
//    // move particle back to the original position
//    P.makeMoveAndCheck(iat,-1.0*dr);
//  }


  void evaluateDerivRatios(ParticleSet& VP, const opt_variables_type& optvars, Matrix<ValueType>& dratios)
  {
  }

  void evaluateDerivatives(ParticleSet& P, const opt_variables_type& optvars, 
    std::vector<RealType>& dlogpsi, std::vector<RealType>& dhpsioverpsi)
  {
  }


  //void evaluateDerivatives(ParticleSet& P, 
  //                         const opt_variables_type& active, 
  //                         std::vector<RealType>& dlogpsi, 
  //                         std::vector<RealType>& dhpsioverpsi)
  //{
  //  evaluateExponents(P);
  //  // indices, strings
  //  // evaluate derivatives of F
  //  if(opt_F)
  //  {
  //    for(int oi = 0; oi < opt_index[OPT_F].size(); ++oi)
  //    {
  //
  //      std::string id = opt_id[OPT_F][oi];
  //      int ia = myVars.getIndex(id);
  //      if(ia == -1)
  //        continue; // ignore inactive parameters
  //      int IJ = opt_index[OPT_F][oi];
  //      int I = IJ/num_regions;
  //      int J = IJ%num_regions;
  //      // coefficient due to symmetry of F: \sum\limits_{I} F_{II} C_I^2 + \sum\limits_{J > I} 2 F_{IJ}*C_I*C_J
  //      RealType x = (I==J)?1:2;
  //      RealType dJF_val = C->sum(I)*C->sum(J)*x;
  //      RealType dJF_gg = 0, dJF_lap = 0;
  //      for(int i = 0; i < num_els; ++i)
  //      {
  //         dJF_gg += x*(dot(C->grad(I,i),P.G[i])*C->sum(J) + C->sum(I)*dot(C->grad(J,i),P.G[i]));
  //         dJF_lap += x*(C->lap(I,i)*C->sum(J) + 2*dot(C->grad(I,i),C->grad(J,i)) + C->lap(J,i)*C->sum(I));
  //      }
  //      dlogpsi[ia] += dJF_val;
  //      dhpsioverpsi[ia] += -0.5*dJF_lap - dJF_gg;
  //      
  //    }
  //  }
  //
  //  // evaluate partial derivatives of G
  //  if(opt_G)
  //  {
  //    RealType dJG_val, dJG_gg, dJG_lap;
  //    for(int oi = 0; oi < opt_index[OPT_G].size(); ++oi)
  //    {
  //      std::string id = opt_id[OPT_G][oi];
  //      int ia = myVars.getIndex(id);
  //      if(ia == -1)
  //        continue; // ignore inactive params
  //      int I = opt_index[OPT_G][oi];
  //      RealType dJG_val = C->sum(I);
  //      RealType dJG_gg = dJG_lap = 0;
  //      for(int i = 0; i < num_els; ++i)
  //      {
  //         dJG_gg += dot(C->grad(I,i),P.G[i]);
  //         dJG_lap += C->lap(I,i);
  //      }
  //      dlogpsi[ia] += dJG_val;
  //      dhpsioverpsi[ia] += -0.5*dJG_lap - dJG_gg;
  //    }
  //  }
  //
  //  // evaluate partial derivatives of C
  //  static int deriv_print_index = 0;
  //
  //  if(opt_C)
  //  {
  //    // containers for CountingRegions' evaluateDerivatives calculations
  //    // blocks of dimension n_p x n_C
  //    // ex: dNsum = \sum\limits_k [ [dC1 / dp1] [dC2 / dp1] .. ]
  //    // Where each [dCi / dpj] block is n_p x 1 vector of derivatives of parameters
  //    //   for counting region j as summed over electron coordinate k
  //
  //    // exception: dNFN ggsum is an n_p x 1 vector since it is an evaluation of a quadratic form:
  //    // \sum\limits_{kI} [\nabla_k dC_I/dpj] dot [ (F \nabla_k C)_I ] 
  //    // since we have the premultiplied (F\nabla_k C) vector on hand.
  //    // make a lambda function FCgrad(I,i) which gives the appropriate element of FCgrad[iI]
  //
  //    // clear some vectors
  //    std::fill(FCggsum.begin(),FCggsum.end(),0);
  //    std::fill(FClapsum.begin(),FClapsum.end(),0);
  //
  //    // easy-index functions for evaluateDerivatives calls
  //    std::function<const GradType&(int,int)> _FCgrad = [&](int I, int i)->const GradType&{ return FCgrad[I*num_els + i] ; };
  //    std::function<const RealType&(int,int)> _FClap  = [&](int I, int i)->const RealType&{ return FClap[I*num_els + i] ; };
  //    std::function<RealType&(int,int)> _dCsum        = [&](int I, int p)->RealType&{ return dCsum[p*num_regions + I]; }; 
  //    std::function<RealType&(int,int)> _dCggsum      = [&](int I, int p)->RealType&{ return dCggsum[p*num_regions + I] ; };
  //    std::function<RealType&(int,int)> _dClapsum     = [&](int I, int p)->RealType&{ return dClapsum[p*num_regions + I] ; };
  //    // evaluate FCggsum
  //    for(int I = 0; I < num_regions; ++I)
  //    {
  //      for(int i = 0; i < num_els; ++i)
  //      {
  //        FCggsum[I] += dot(_FCgrad(I,i),P.G[i]);
  //        FClapsum[I] += _FClap(I,i);
  //      }
  //    }
  //    // pointer to C->C[I]->myVars.Index
  //    // for pI in { 0 .. C->num_derivs(I) }
  //    //   dCindex->[pI]  is the index that corresponds to this parameter in active.
  //    //   i.e., active[dCindex->[pI]] <=> C->C[I]->myVars.Index[pI]
  //
  //    // external print block
  //    if(debug && deriv_print_index < debug_seqlen)
  //    {
  //      app_log() << std::endl << "=== evaluateDerivatives ===" << std::endl;
  //      app_log() << "== print current exponent values ==" << std::endl;
  //      evaluateExponents_print(app_log(),P);
  //      app_log() << "== additional counting function terms ==" << std::endl;
  //      app_log() << "P.G: ";
  //      std::copy(P.G.begin(), P.G.end(), std::ostream_iterator<GradType>(app_log(), ", "));
  //      app_log() << std::endl << "FCgrad: ";
  //      std::copy(FCgrad.begin(), FCgrad.end(), std::ostream_iterator<GradType>(app_log(), ", "));
  //      app_log() << std::endl << "FClap: ";
  //      std::copy(FClap.begin(), FClap.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //      app_log() << std::endl << "FCggsum: ";
  //      std::copy(FCggsum.begin(), FCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //      app_log() << std::endl << "FClapsum: ";
  //      std::copy(FClapsum.begin(), FClapsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //      app_log() << std::endl;
  //    }
  //
  //    for(int I = 0; I < num_regions; ++I)
  //    {
  //      // get the number of active parameters for the Ith counting region
  //      opt_variables_type I_vars = C->getVars(I); 
  //      int I_num_derivs = I_vars.size();
  //      // clear arrays before each evaluate
  //      std::fill(dCsum.begin(),dCsum.end(),0);
  //      std::fill(dCggsum.begin(),dCggsum.end(),0);
  //      std::fill(dClapsum.begin(),dClapsum.end(),0);
  //      std::fill(dCFCggsum.begin(),dCFCggsum.end(),0);
  //      // evaluate all derivatives for the Ith counting function
  //      C->evaluateDerivatives(P, I, _FCgrad, _dCsum, _dCggsum, _dClapsum, dCFCggsum);
  //      if(debug && deriv_print_index < debug_seqlen)
  //      {
  //        // print out current index information
  //        app_log() << std::endl;
  //        app_log() << "  == evaluateDerivatives for counting region " << I << ", num_derivs: " << I_num_derivs << " ==" << std::endl;
  //        app_log() << "  Indices: ";
  //        std::copy(I_vars.Index.begin(), I_vars.Index.end(), std::ostream_iterator<int>(app_log(),", "));
  //        app_log() << std::endl << "  Names: ";
  //        for(auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
  //          app_log() << (*it).first << ", ";
  //        app_log() << std::endl << "  Values: ";
  //        for(auto it = I_vars.NameAndValue.begin(); it != I_vars.NameAndValue.end(); ++it)
  //          app_log() << (*it).second << ", ";
  //        // print out values from evaluate derivatives
  //        app_log() << std::endl << "  dCsum: ";
  //        std::copy(dCsum.begin(), dCsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //        app_log() << std::endl << "  dCggsum: ";
  //        std::copy(dCggsum.begin(), dCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //        app_log() << std::endl << "  dClapsum: ";
  //        std::copy(dClapsum.begin(), dClapsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //        app_log() << std::endl << "  dCFCggsum: ";
  //        std::copy(dCFCggsum.begin(), dCFCggsum.end(), std::ostream_iterator<RealType>(app_log(), ", "));
  //        app_log() << std::endl;
  //      }
  //      // loop over parameters for the Ith counting function
  //      for(int pI = 0; pI < I_num_derivs; ++pI)
  //      {
  //        // index for active optimizable variables
  //        int ia = I_vars.Index[pI];
  //        if(ia == -1)
  //          continue; // ignore inactive
  //        // middle laplacian term: 
  //        dhpsioverpsi[ia] += -0.5*(4.0*dCFCggsum[pI]);
  //        if(debug && deriv_print_index < debug_seqlen)
  //        {
  //          app_log() << "    == evaluateDerivatives calculations ==" << std::endl;
  //          app_log() << "    pI: " << pI << ", name: " <<  I_vars.name(pI) <<  ", ia: " << ia << std::endl;
  //          app_log() << "    dCFCggsum: " << dCFCggsum[pI] << std::endl;
  //        }
  //        for(int J = 0; J < num_regions; ++J)
  //        {
  //          dlogpsi[ia] += _dCsum(J,pI)*(2*FCsum[J] + G[J]);
  //          // grad dot grad terms
  //          dhpsioverpsi[ia] += -1.0*( _dCggsum(J,pI)*(2.0*FCsum[J] + G[J]) + _dCsum(J,pI)*2.0*FCggsum[J]  );
  //          // outer laplacian terms
  //          dhpsioverpsi[ia] += -0.5*( 2.0*_dCsum(J,pI)*FClapsum[J] + _dClapsum(J,pI)*(2.0*FCsum[J] + G[J]) ) ;
  //          if(debug && deriv_print_index < debug_seqlen)
  //          {
  //            app_log() << "      J: " << J << std::endl;
  //            app_log() << "      dlogpsi term          : " << _dCsum(J,pI)*(2*FCsum[J] + G[J]) << std::endl;
  //            app_log() << "      dhpsi/psi, graddotgrad: " << -1.0*( _dCggsum(J,pI)*(2.0*FCsum[J] + G[J]) + _dCsum(J,pI)*2.0*FCggsum[J]  ) << std::endl;
  //            app_log() << "      dhpsi/psi, laplacian  : " << -0.5*( 2.0*_dCsum(J,pI)*FClapsum[J] + _dClapsum(J,pI)*(2.0*FCsum[J] + G[J]) ) << std::endl;
  //          }
  //
  //
  //        }
  //      }
  //    }
  //
  //  } // end opt_C
  //  // increment and modulo deriv_print_index
  //  deriv_print_index = deriv_print_index % debug_period;
  //  deriv_print_index++;
  //}

  bool addRegion(RegionType* CR);

  bool addDebug(int seqlen, int period);

};

}
#endif
