//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Miguel Morales, moralessilva2@llnl.gov, Lawrence Livermore National Laboratory
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    


#include "QMCDrivers/QMCFixedSampleLinearOptimize.h"
#include "Particle/HDFWalkerIO.h"
#include "Particle/DistanceTable.h"
#include "OhmmsData/AttributeSet.h"
#include "Message/CommOperators.h"
#include "QMCDrivers/QMCCostFunctionBase.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#if defined(ENABLE_OPENMP)
#include "QMCDrivers/VMC/VMCSingleOMP.h"
#include "QMCDrivers/QMCCostFunctionOMP.h"
#endif
//#include "QMCDrivers/VMC/VMCSingle.h"
//#include "QMCDrivers/QMCCostFunctionSingle.h"
#include "QMCApp/HamiltonianPool.h"
#include "Numerics/Blasf.h"
#include "Numerics/MatrixOperators.h"
#include <cassert>
#if defined(QMC_CUDA)
#include "QMCDrivers/VMC/VMC_CUDA.h"
#include "QMCDrivers/QMCCostFunctionCUDA.h"
#endif
#ifdef HAVE_LMY_ENGINE
#include "formic/utils/matrix.h"
#include "formic/utils/random.h"
#include "formic/utils/lmyengine/var_dependencies.h"
#endif
#include <iostream>
#include <fstream>
#include <stdexcept>

/*#include "Message/Communicate.h"*/

namespace qmcplusplus
{

  using MatrixOperators::product;


QMCFixedSampleLinearOptimize::QMCFixedSampleLinearOptimize(MCWalkerConfiguration& w,
    TrialWaveFunction& psi, QMCHamiltonian& h, HamiltonianPool& hpool, WaveFunctionPool& ppool):
  QMCLinearOptimize(w,psi,h,hpool,ppool), 
#ifdef HAVE_LMY_ENGINE
vdeps(1,std::vector<double>()),
#endif
 Max_iterations(1), exp0(-16), nstabilizers(3),
  stabilizerScale(2.0), bigChange(50), w_beta(0.0),  MinMethod("OneShiftOnly"), GEVtype("mixed"),
  StabilizerMethod("best"), GEVSplit("no"), stepsize(0.25), doAdaptiveThreeShift(false),
  targetExcitedStr("no"), targetExcited(false), block_lmStr("no"), block_lm(false),
  bestShift_i(-1.0), bestShift_s(-1.0), shift_i_input(0.01), shift_s_input(1.00), doOneShiftOnly(false),
  num_shifts(3), nblocks(1), nolds(1), nkept(1), nsamp_comp(0), omega_shift(0.0), max_param_change(0.3),
  max_relative_cost_change(10.0), block_first(true), block_second(false), block_third(false)
{
  IsQMCDriver=false;
  //set the optimization flag
  QMCDriverMode.set(QMC_OPTIMIZE,1);
  //read to use vmc output (just in case)
  RootName = "pot";
  QMCType ="QMCFixedSampleLinearOptimize";
  m_param.add(WarmupBlocks,"warmupBlocks","int");
  m_param.add(Max_iterations,"max_its","int");
  m_param.add(nstabilizers,"nstabilizers","int");
  m_param.add(stabilizerScale,"stabilizerscale","double");
  m_param.add(bigChange,"bigchange","double");
  m_param.add(MinMethod,"MinMethod","string");
  m_param.add(exp0,"exp0","double");
  m_param.add(targetExcitedStr,"targetExcited","string");
  m_param.add(block_lmStr, "block_lm", "string");
  m_param.add(nblocks, "nblocks", "int");
  m_param.add(nolds, "nolds", "int");
  m_param.add(nkept, "nkept", "int");
  m_param.add(nsamp_comp, "nsamp_comp", "int");
  m_param.add(omega_shift,"omega","double");
  m_param.add(max_relative_cost_change,"max_relative_cost_change","double");
  m_param.add(max_param_change,"max_param_change","double");
  m_param.add(shift_i_input, "shift_i", "double");
  m_param.add(shift_s_input, "shift_s", "double");
  m_param.add(num_shifts, "num_shifts", "int");

  #ifdef HAVE_LMY_ENGINE
  //app_log() << "construct QMCFixedSampleLinearOptimize" << endl;
  std::vector<double> shift_scales(3, 1.0);
  EngineObj = new cqmc::engine::LMYEngine(&vdeps, 
                                          false, // exact sampling
                                          true, // ground state?
                                          false, // variance correct,
                                          true, 
                                          false, // print matrices,
                                          true, // build matrices
                                          false, // spam
                                          false, // use var deps?
                                          true, // chase lowest
                                          false, // chase closest
                                          false, // eom
                                          false, 
                                          false, // eom related
                                          false, // eom related
                                          false, // use block?
                                          120000, // number of samples
                                          0,  // number of parameters
                                          60, // max krylov iter
                                          0, // max spam inner iter
                                          1, // spam appro degree
                                          0, // eom relateds 
                                          0, // eom related
                                          0, // eom related
                                          0.0, // omega
                                          0.0, // var weight
                                          1.0e-6, // convergence threshold
                                          0.99, // minimum S singular val
                                          0.0, 
                                          0.0, 
                                          10.0, // max change allowed
                                          1.00, // identity shift
                                          1.00, // overlap shift
                                          0.3, // max parameter change
                                          shift_scales, 
                                          app_log());
  #endif
//   stale parameters
//   m_param.add(eigCG,"eigcg","int");
//   m_param.add(TotalCGSteps,"cgsteps","int");
//   m_param.add(w_beta,"beta","double");
//   quadstep=-1.0;
//   m_param.add(quadstep,"quadstep","double");
//   m_param.add(stepsize,"stepsize","double");
//   m_param.add(exp1,"exp1","double");
//   m_param.add(GEVtype,"GEVMethod","string");
//   m_param.add(GEVSplit,"GEVSplit","string");
//   m_param.add(StabilizerMethod,"StabilizerMethod","string");
//   m_param.add(LambdaMax,"LambdaMax","double");
  //Set parameters for line minimization:
}

/** Clean up the vector */
QMCFixedSampleLinearOptimize::~QMCFixedSampleLinearOptimize()
{
  #ifdef HAVE_LMY_ENGINE
  delete EngineObj;
  #endif
}

QMCFixedSampleLinearOptimize::RealType QMCFixedSampleLinearOptimize::Func(RealType dl)
{
  for (int i=0; i<optparm.size(); i++)
    optTarget->Params(i) = optparm[i] + dl*optdir[i];
  QMCLinearOptimize::RealType c = optTarget->Cost(false);
  //only allow this to go false if it was true. If false, stay false
//    if (validFuncVal)
  validFuncVal = optTarget->IsValid;
  return c;
}

bool QMCFixedSampleLinearOptimize::run()
{

  // if requested, perform the update via the adaptive three-shift or single-shift method
#ifdef HAVE_LMY_ENGINE
  if ( doAdaptiveThreeShift ) return adaptive_three_shift_run();
#endif
  if ( doOneShiftOnly ) return one_shift_run();

  start();
  bool Valid(true);
  int Total_iterations(0);
//size of matrix
  numParams = optTarget->NumParams();
  N = numParams + 1;
//   where we are and where we are pointing
  std::vector<RealType> currentParameterDirections(N,0);
  std::vector<RealType> currentParameters(numParams,0);
  std::vector<RealType> bestParameters(numParams,0);
  for (int i=0; i<numParams; i++)
    bestParameters[i] = currentParameters[i] = optTarget->Params(i);
//   proposed direction and new parameters
  optdir.resize(numParams,0);
  optparm.resize(numParams,0);

  while (Total_iterations < Max_iterations)
  {
    Total_iterations+=1;
    app_log()<<"Iteration: "<<Total_iterations<<"/"<<Max_iterations<< std::endl;
    if (!ValidCostFunction(Valid))
      continue;
//this is the small amount added to the diagonal to stabilize the eigenvalue equation. 10^stabilityBase
    RealType stabilityBase(exp0);
//     reset params if necessary
    for (int i=0; i<numParams; i++)
      optTarget->Params(i) = currentParameters[i];
    myTimers[4]->start();
    RealType lastCost(optTarget->Cost(true));
    myTimers[4]->stop();
//     if cost function is currently invalid continue
    Valid=optTarget->IsValid;
    if (!ValidCostFunction(Valid))
      continue;
    RealType newCost(lastCost);
    RealType startCost(lastCost);
    Matrix<RealType> Left(N,N);
    Matrix<RealType> Right(N,N);
    Matrix<RealType> S(N,N);
//     stick in wrong matrix to reduce the number of matrices we need by 1.( Left is actually stored in Right, & vice-versa)
    optTarget->fillOverlapHamiltonianMatrices(Right,Left);
    S.copy(Left);
    bool apply_inverse(true);
    if(apply_inverse)
    {
      Matrix<RealType> RightT(Left);
      invert_matrix(RightT,false);
      Left=0;
      product(RightT,Right,Left);
//       Now the left matrix is the Hamiltonian with the inverse of the overlap applied ot it.
    }
    //Find largest off-diagonal element compared to diagonal element.
    //This gives us an idea how well conditioned it is, used to stabilize.
    RealType od_largest(0);
    for (int i=0; i<N; i++)
      for (int j=0; j<N; j++)
        od_largest=std::max( std::max(od_largest,std::abs(Left(i,j))-std::abs(Left(i,i))), std::abs(Left(i,j))-std::abs(Left(j,j)));
    app_log()<<"od_largest "<<od_largest<< std::endl;
    //if(od_largest>0)
    //  od_largest = std::log(od_largest);
    //else
    //  od_largest = -1e16;
    //if (od_largest<stabilityBase)
    //  stabilityBase=od_largest;
    //else
    //  stabilizerScale = std::max( 0.2*(od_largest-stabilityBase)/nstabilizers, stabilizerScale);
    app_log()<<"  stabilityBase "<<stabilityBase<< std::endl;
    app_log()<<"  stabilizerScale "<<stabilizerScale<< std::endl;
    int failedTries(0);
    bool acceptedOneMove(false);
    for (int stability=0; stability<nstabilizers; stability++)
    {
      bool goodStep(true);
//       store the Hamiltonian matrix in Right
      for (int i=0; i<N; i++)
        for (int j=0; j<N; j++)
          Right(i,j)= Left(j,i);
      RealType XS(stabilityBase+stabilizerScale*(failedTries+stability));
      for (int i=1; i<N; i++)
        Right(i,i) += std::exp(XS);
      app_log()<<"  Using XS:"<<XS<<" "<<failedTries<<" "<<stability<< std::endl;
      RealType lowestEV(0);
      myTimers[2]->start();
      lowestEV = getLowestEigenvector(Right,currentParameterDirections);
      Lambda = getNonLinearRescale(currentParameterDirections,S);
      myTimers[2]->stop();
//       biggest gradient in the parameter direction vector
      RealType bigVec(0);
      for (int i=0; i<numParams; i++)
        bigVec = std::max(bigVec,std::abs(currentParameterDirections[i+1]));
//       this can be overwritten during the line minimization
      RealType evaluated_cost(startCost);
      if (MinMethod=="rescale")
      {
        if (std::abs(Lambda*bigVec)>bigChange)
        {
          goodStep=false;
          app_log()<<"  Failed Step. Magnitude of largest parameter change: "<<std::abs(Lambda*bigVec)<< std::endl;
          if (stability==0)
          {
            failedTries++;
            stability--;
          }
          else
            stability=nstabilizers;
        }
        for (int i=0; i<numParams; i++)
          optTarget->Params(i) = currentParameters[i] + Lambda*currentParameterDirections[i+1];
        optTarget->IsValid = true;
      }
      else
      {
        for (int i=0; i<numParams; i++)
          optparm[i] = currentParameters[i];
        for (int i=0; i<numParams; i++)
          optdir[i] = currentParameterDirections[i+1];
        TOL = param_tol/bigVec;
        AbsFuncTol=true;
        largeQuarticStep=bigChange/bigVec;
        LambdaMax = 0.5*Lambda;
        myTimers[3]->start();
        if (MinMethod=="quartic")
        {
          int npts(7);
          quadstep = stepsize*Lambda;
          largeQuarticStep=bigChange/bigVec;
          Valid=lineoptimization3(npts,evaluated_cost);
        }
        else
          Valid=lineoptimization2();
        myTimers[3]->stop();
        RealType biggestParameterChange = bigVec*std::abs(Lambda);
        if (biggestParameterChange>bigChange)
        {
          goodStep=false;
          failedTries++;
          app_log()<<"  Failed Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
          if (stability==0)
            stability--;
          else
            stability=nstabilizers;
        }
        else
        {
          for (int i=0; i<numParams; i++)
            optTarget->Params(i) = optparm[i] + Lambda * optdir[i];
          app_log()<<"  Good Step. Largest LM parameter change:"<<biggestParameterChange<< std::endl;
        }
      }

      if (goodStep)
      {
// 	this may have been evaluated allready
// 	newCost=evaluated_cost;
        //get cost at new minimum
        newCost = optTarget->Cost(false);
        app_log()<<" OldCost: "<<lastCost<<" NewCost: "<<newCost<<" Delta Cost:"<<(newCost-lastCost)<< std::endl;
        optTarget->printEstimates();
        //                 quit if newcost is greater than lastcost. E(Xs) looks quadratic (between steepest descent and parabolic)
        // mmorales
        Valid=optTarget->IsValid;
        //if (MinMethod!="rescale" && !ValidCostFunction(Valid))
        if (!ValidCostFunction(Valid))
        {
          goodStep=false;
          app_log()<<"  Good Step, but cost function invalid"<< std::endl;
          failedTries++;
          if(stability>0)
            stability=nstabilizers;
          else
            stability--;
        }
        if (newCost < lastCost && goodStep)
        {
          //Move was acceptable
          for (int i=0; i<numParams; i++)
            bestParameters[i] = optTarget->Params(i);
          lastCost=newCost;
          acceptedOneMove=true;
          if(std::abs(newCost-lastCost)<1e-4)
          {
            failedTries++;
            stability=nstabilizers;
            continue;
          }
        }
        else if (stability>0)
        {
          failedTries++;
          stability=nstabilizers;
          continue;
        }
      }
      app_log().flush();
      app_error().flush();
      if(failedTries>20) break;
        //APP_ABORT("QMCFixedSampleLinearOptimize::run TOO MANY FAILURES");
    }

    if (acceptedOneMove)
    {
      app_log()<<"Setting new Parameters"<<std::endl;
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = bestParameters[i];
    }
    else
    {
      app_log()<<"Revertting to old Parameters"<<std::endl;
      for (int i=0; i<numParams; i++)
        optTarget->Params(i) = currentParameters[i];
    }
    app_log().flush();
    app_error().flush();
  }

  finish();
  return (optTarget->getReportCounter() > 0);
}

/** Parses the xml input file for parameter definitions for the wavefunction optimization.
* @param q current xmlNode
* @return true if successful
*/
bool
QMCFixedSampleLinearOptimize::put(xmlNodePtr q)
{
  std::string useGPU("yes");
  std::string vmcMove("pbyp");
  OhmmsAttributeSet oAttrib;
  oAttrib.add(useGPU,"gpu");
  oAttrib.add(vmcMove,"move");
  oAttrib.put(q);
  m_param.put(q);

  tolower(targetExcitedStr);
  targetExcited = ( targetExcitedStr == "yes" );

  tolower(block_lmStr);
  block_lm = ( block_lmStr == "yes" );

  // get whether to use the adaptive three-shift version of the update
  doAdaptiveThreeShift = ( MinMethod == "adaptive" );
  doOneShiftOnly = ( MinMethod == "OneShiftOnly" );

  // sanity check
  if ( targetExcited && !doAdaptiveThreeShift )
    APP_ABORT("targetExcited = \"yes\" requires that MinMethod = \"adaptive\"");

#ifdef ENABLE_OPENMP
  if ( doAdaptiveThreeShift && (omp_get_max_threads() > 1) ) {
        //throw std::runtime_error("OpenMP threading not enabled with AdaptiveThreeShift optimizer.  Use MPI for parallelism instead, and set OMP_NUM_THREADS to 1.");
        app_log() << "test version of OpenMP threading with AdaptiveThreeShift optimizer" << std::endl;
  }
#endif

  // check parameter change sanity
  if ( max_param_change <= 0.0 ) throw std::runtime_error("max_param_change must be positive in QMCFixedSampleLinearOptimize::put");

  // check cost change sanity
  if ( max_relative_cost_change <= 0.0 ) throw std::runtime_error("max_relative_cost_change must be positive in QMCFixedSampleLinearOptimize::put");

  // check shift sanity
  if ( shift_i_input <= 0.0 ) throw std::runtime_error("shift_i must be positive in QMCFixedSampleLinearOptimize::put");
  if ( shift_s_input <= 0.0 ) throw std::runtime_error("shift_s must be positive in QMCFixedSampleLinearOptimize::put");

  // if this is the first time this function has been called, set the initial shifts
  if ( bestShift_i < 0.0 && doAdaptiveThreeShift ) bestShift_i = shift_i_input;
  if ( doOneShiftOnly ) bestShift_i = shift_i_input;
  if ( bestShift_s < 0.0 ) bestShift_s = shift_s_input;

  xmlNodePtr qsave=q;
  xmlNodePtr cur=qsave->children;
  int pid=OHMMS::Controller->rank();
  while (cur != NULL)
  {
    std::string cname((const char*)(cur->name));
    if (cname == "mcwalkerset")
    {
      mcwalkerNodePtr.push_back(cur);
    }
    cur=cur->next;
  }
  //no walkers exist, add 10
  if (W.getActiveWalkers() == 0)
    addWalkers(omp_get_max_threads());
  NumOfVMCWalkers=W.getActiveWalkers();
  //create VMC engine
  if (vmcEngine ==0)
  {
#if defined (QMC_CUDA)
    if (useGPU == "yes")
      vmcEngine = new VMCcuda(W,Psi,H,psiPool);
    else
#endif
      vmcEngine = new VMCSingleOMP(W,Psi,H,hamPool,psiPool);
    vmcEngine->setUpdateMode(vmcMove[0] == 'p');
    vmcEngine->initCommunicator(myComm);
  }

  vmcEngine->setStatus(RootName,h5FileRoot,AppendRun);
  vmcEngine->process(qsave);

  bool success=true;
  if (optTarget == 0)
  {
#if defined (QMC_CUDA)
    if (useGPU == "yes")
      optTarget = new QMCCostFunctionCUDA(W,Psi,H,hamPool);
    else
#endif
      optTarget = new QMCCostFunctionOMP(W,Psi,H,hamPool);
    optTarget->setStream(&app_log());
    success=optTarget->put(q);
  }
  return success;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  returns a vector of three shift values centered around the provided shift.
///
/// \param[in]      central_shift  the central shift
///
///////////////////////////////////////////////////////////////////////////////////////////////////
std::vector<double> QMCFixedSampleLinearOptimize::prepare_shifts(const double central_shift) const {
  std::vector<double> retval(num_shifts);

  // check to see whether the number of shifts is odd
  if ( num_shifts % 2 == 0 ) throw std::runtime_error("number of shifts must be odd in QMCFixedSampleLinearOptimize::prepare_shifts");

  // decide the central shift index
  int central_index = num_shifts / 2;

  for (int i = 0; i < num_shifts; i++) {
    if ( i < central_index ) 
      retval.at(i) = central_shift / (4.0 * (central_index - i));
    else if ( i > central_index ) 
      retval.at(i) = central_shift * (4.0 * (i - central_index));
    else if ( i == central_index ) 
      retval.at(i) = central_shift;
    //retval.at(i) = central_shift 
    //retval.at(0) = central_shift * 4.0;
    //retval.at(1) = central_shift;
    //retval.at(2) = central_shift / 4.0;
  } 
  return retval;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  prints a header for the summary of each shift's result
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimize::print_cost_summary_header() {
    app_log() << "   " << std::right << std::setw(12) << "shift_i";
    app_log() << "   " << std::right << std::setw(12) << "shift_s";
    app_log() << "   " << std::right << std::setw(20) << "max param change";
    app_log() << "   " << std::right << std::setw(20) << "cost function value";
    app_log() << std::endl;
    app_log() << "   " << std::right << std::setw(12) << "------------";
    app_log() << "   " << std::right << std::setw(12) << "------------";
    app_log() << "   " << std::right << std::setw(20) << "--------------------";
    app_log() << "   " << std::right << std::setw(20) << "--------------------";
    app_log() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  prints a summary of the computed cost for the given shift
///
/// \param[in]      si             the identity shift
/// \param[in]      ss             the overlap shift
/// \param[in]      mc             the maximum parameter change
/// \param[in]      cv             the cost function value
/// \param[in]      ind            the shift index: -1 (for initial state), 0, 1, or 2
/// \param[in]      bi             index of the best shift
/// \param[in]      gu             flag telling whether it was a good update
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimize::print_cost_summary(const double si, const double ss, const RealType mc, const RealType cv, const int ind, const int bi, const bool gu) {
    if ( ind >= 0 ) {
      if ( gu ) {
        app_log() << "   " << std::scientific << std::right << std::setw(12) << std::setprecision(4)  << si;
        app_log() << "   " << std::scientific << std::right << std::setw(12) << std::setprecision(4)  << ss;
        app_log() << "   " << std::scientific << std::right << std::setw(20) << std::setprecision(4)  << mc;
        app_log() << "   " << std::fixed      << std::right << std::setw(20) << std::setprecision(12) << cv;
        //app_log() << "   " << std::right << std::setw(12) << ( ind == 0 ? "big shift" : ( ind == 1 ? "medium shift" : "small shift" ) );
      } else {
        app_log() << "   " << std::right << std::setw(12) << "N/A";
        app_log() << "   " << std::right << std::setw(12) << "N/A";
        app_log() << "   " << std::right << std::setw(20) << "N/A";
        app_log() << "   " << std::right << std::setw(20) << "N/A";
        app_log() << "   " << std::right << std::setw(12) << "bad update";
      }
    } else {
      app_log() << "   " << std::right << std::setw(12) << "N/A";
      app_log() << "   " << std::right << std::setw(12) << "N/A";
      app_log() << "   " << std::right << std::setw(20) << "N/A";
      app_log() << "   " << std::fixed      << std::right << std::setw(20) << std::setprecision(12) << cv;
      app_log() << "   " << std::right << std::setw(12) << "initial";
    }
    if ( ind == bi )
      app_log() << "  <--";
    app_log() << std::endl;
}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Returns whether the proposed new cost is the best compared to the others.
///
/// \param[in]      ii             index of the proposed best cost
/// \param[in]      cv             vector of new costs
/// \param[in]      ic             the initial cost
///
///////////////////////////////////////////////////////////////////////////////////////////////////
bool QMCFixedSampleLinearOptimize::is_best_cost(const int ii, const std::vector<RealType> & cv, const RealType ic) const {
  
  // initialize return value
  bool retval = true;

  //app_log() << "retval = " << retval << std::endl;

  // new cost is only best if it is less than or equal to the other new costs
  for (int i = 0; i < cv.size(); i++) {
    // we only worry about being less than the other value if the other value is within the maximum relative change threshold
    const bool other_is_valid = ( ( ic == 0.0 ? 0.0 : std::abs( ( cv.at(i) - ic ) / ic ) ) < this->max_relative_cost_change );
    retval = ( retval && ( cv.at(ii) <= cv.at(i) || !other_is_valid ) );
    //app_log() << "cv.at(ii)   = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(ii) << " <= "
    //          << "cv.at(i)    = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(i)  << " ?" << std::endl;
    //app_log() << "retval = " << retval << std::endl;
  }

  // new cost is only best if it is less than the initial cost
  retval = ( retval && cv.at(ii) < ic );
  //app_log() << "cv.at(ii)   = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << cv.at(ii) << " <= "
  //          << "ic          = " << std::fixed << std::right << std::setw(20) << std::setprecision(12) << ic        << " ?" << std::endl;
  //app_log() << "retval = " << retval << std::endl;

  // new cost is only best if it's relative change from the initial cost is not too large ( or if the initial cost is exactly zero )
  retval = ( retval && ( ic == 0.0 ? 0.0 : std::abs( ( cv.at(ii) - ic ) / ic ) ) < this->max_relative_cost_change );
  //app_log() << "std::abs( ( cv.at(ii) - ic ) / ic ) = " << std::fixed << std::right << std::setw(20) << std::setprecision(12)
  //          << std::abs( ( cv.at(ii) - ic ) / ic ) << " <= " << this->max_relative_cost_change << " ? " << std::endl;
  //app_log() << "retval = " << retval << std::endl;
  //app_log() << std::endl;

  // return whether the proposed cost is actually the best
  return retval;

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  For each set of shifts, solves the linear method eigenproblem by building and
///         diagonalizing the matrices.
///
/// \param[in]      shfits_i              vector of identity shifts
/// \param[in]      shfits_s              vector of overlap shifts
/// \param[out]     parameterDirections   on exit, the update directions for the different shifts
///
///////////////////////////////////////////////////////////////////////////////////////////////////
void QMCFixedSampleLinearOptimize::solveShiftsWithoutLMYEngine(const std::vector<double> & shifts_i,
                                                               const std::vector<double> & shifts_s,
                                                               std::vector<std::vector<RealType> > & parameterDirections) {

  // get number of shifts to solve
  const int nshifts = shifts_i.size();

  // get number of optimizable parameters
  numParams = optTarget->NumParams();

  // get dimension of the linear method matrices
  N = numParams + 1;

  // prepare vectors to hold the parameter updates
  parameterDirections.resize(nshifts);
  for (int i = 0; i < parameterDirections.size(); i++)
    parameterDirections.at(i).assign(N, 0.0);

  // allocate the matrices we will need
  Matrix<RealType> ovlMat(N,N); ovlMat = 0.0;
  Matrix<RealType> hamMat(N,N); hamMat = 0.0;
  Matrix<RealType> invMat(N,N); invMat = 0.0;
  Matrix<RealType> sftMat(N,N); sftMat = 0.0;
  Matrix<RealType> prdMat(N,N); prdMat = 0.0;

  // build the overlap and hamiltonian matrices
  optTarget->fillOverlapHamiltonianMatrices(hamMat, ovlMat);

  //// print the hamiltonian matrix
  //app_log() << std::endl;
  //app_log() << "printing H matrix:" << std::endl;
  //for (int i = 0; i < hamMat.rows(); i++) {
  //  for (int j = 0; j < hamMat.cols(); j++)
  //    app_log() << " " << std::scientific << std::right << std::setw(14) << std::setprecision(5) << hamMat(i,j);
  //  app_log() << std::endl;
  //}
  //app_log() << std::endl;

  //// print the overlap matrix
  //app_log() << std::endl;
  //app_log() << "printing S matrix:" << std::endl;
  //for (int i = 0; i < ovlMat.rows(); i++) {
  //  for (int j = 0; j < ovlMat.cols(); j++)
  //    app_log() << " " << std::scientific << std::right << std::setw(14) << std::setprecision(5) << ovlMat(i,j);
  //  app_log() << std::endl;
  //}
  //app_log() << std::endl;

  // compute the inverse of the overlap matrix
  invMat.copy(ovlMat);
  invert_matrix(invMat, false);

  // compute the update for each shift
  for (int shift_index = 0; shift_index < nshifts; shift_index++) {

    // prepare to shift the hamiltonain matrix
    sftMat.copy(hamMat);

    // apply the identity shift
    for (int i=1; i<N; i++)
      sftMat(i,i) += shifts_i.at(shift_index);

    // apply the overlap shift
    for (int i=1; i<N; i++)
    for (int j=1; j<N; j++)
      sftMat(i,j) += shifts_s.at(shift_index) * ovlMat(i,j);

    // multiply the shifted hamiltonian matrix by the inverse of the overlap matrix
    qmcplusplus::MatrixOperators::product(invMat, sftMat, prdMat);

    // transpose the result (why?)
    for (int i=0; i<N; i++)
    for (int j=i+1; j<N; j++)
      std::swap(prdMat(i,j), prdMat(j,i));

    // compute the lowest eigenvalue of the product matrix and the corresponding eigenvector
    const RealType lowestEV = this->getLowestEigenvector(prdMat, parameterDirections.at(shift_index));

    // compute the scaling constant to apply to the update
    this->Lambda = this->getNonLinearRescale(parameterDirections.at(shift_index), ovlMat);

    // scale the update by the scaling constant
    for (int i=0; i<numParams; i++)
      parameterDirections.at(shift_index).at(i+1) *= Lambda;

  }

}

///////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Performs one iteration of the linear method using an adaptive scheme that tries three
///         different shift magnitudes and picks the best one.
///         The scheme is adaptive in that it saves the best shift to use as a starting point
///         in the next iteration.
///         Note that the best shift is chosen based on a different sample than that used to
///         construct the linear method matrices in order to avoid over-optimizing on a particular
///         sample.
///
/// \return  ???
///
///////////////////////////////////////////////////////////////////////////////////////////////////
#ifdef HAVE_LMY_ENGINE
bool QMCFixedSampleLinearOptimize::adaptive_three_shift_run() {

  // remember what the cost function grads flag was
  const bool saved_grads_flag = this->optTarget->getneedGrads();

  // remember the initial number of samples
  const int init_num_samp = this->optTarget->getNumSamples();

  // the index of central shift
  const int central_index = num_shifts / 2;

  // get number of optimizable parameters
  numParams = optTarget->NumParams();

  // prepare the shifts that we will try
  const std::vector<double> shifts_i = this->prepare_shifts(this->bestShift_i);
  const std::vector<double> shifts_s = this->prepare_shifts(this->bestShift_s);
  std::vector<double> shift_scales(shifts_i.size(), 1.0);
  for (int i = 0; i < shift_scales.size(); i++) 
    shift_scales.at(i) = shifts_i.at(i) / shift_i_input;

  // ensure the cost function is set to compute derivative vectors
  this->optTarget->setneedGrads(true);

  // prepare previous updates 
  int count = 0;
  while ( block_lm && previous_update.size() < nolds ) {
    previous_update.push_back(formic::ColVec<double>(numParams));
    for (int i = 0; i < numParams; i++) 
      previous_update.at(count).at(i) =  2.0 * ( formic::random_number<double>() - 0.5 );
    count++;
  }

  if ( !EngineObj->full_init() ) {

    // prepare a variable dependency object with no dependencies
    formic::VarDeps real_vdeps(numParams, std::vector<double>());
    vdeps = real_vdeps;
    EngineObj->get_param(&vdeps,
                         false, // exact sampling
                         !targetExcited, 
                         false, // variable deps use?
                         false, // eom
                         false, // ssquare
                         block_lm, 
                         12000, 
                         numParams,
                         omega_shift,
                         max_relative_cost_change,
                         shifts_i.at(central_index), 
                         shifts_s.at(central_index),
                         max_param_change,
                         shift_scales);
  }

  // update shift
  EngineObj->shift_update(shift_scales);

  // turn on wavefunction update mode
  EngineObj->turn_on_update();

  // initialize the engine if we do not use block lm or it's the first part of block lm
  EngineObj->initialize(nblocks, 
                        0, 
                        nkept, 
                        previous_update,
                        false);

  // reset the engine
  EngineObj->reset();

  // generate samples and compute weights, local energies, and derivative vectors
  this->engine_start(EngineObj);

  // get dimension of the linear method matrices
  N = numParams + 1;

  // have the cost function prepare derivative vectors
  EngineObj->energy_target_compute();
  const Return_t starting_cost = EngineObj->target_value();
  const Return_t init_energy = EngineObj->energy_mean();
  
  // print out the initial energy
  app_log() << std::endl
            << "*************************************************************************************************" << std::endl
            << "Solving the linear method equations on the initial sample with initial energy" << std::setw(20) << std::setprecision(12) << init_energy << std::endl
            << "*************************************************************************************************" << std::endl
            << std::endl;
  //const Return_t starting_cost = this->optTarget->LMYEngineCost(true);

  // prepare wavefunction update which does nothing if we do not use block lm
  EngineObj->wfn_update_prep();

  if ( block_lm ) {

    this->optTarget->setneedGrads(true);

    int numOptParams = optTarget->NumParams();

    // reset the engine object 
    EngineObj->reset();

    // finish last sample
    this->finish();

    // take sample
    this->engine_start(EngineObj);
  }

  // say what we are doing
  app_log() << std::endl
            << "*********************************************************" << std::endl
            << "Solving the linear method equations on the initial sample" << std::endl
            << "*********************************************************" << std::endl
            << std::endl;

  // for each set of shifts, solve the linear method equations for the parameter update direction
  std::vector<std::vector<RealType> > parameterDirections;
  #ifdef HAVE_LMY_ENGINE
  // call the engine to perform update
  EngineObj->wfn_update_compute();
  //std::cout << "optimization here 0.5" << std::endl;
  #else
  this->solveShiftsWithoutLMYEngine(shifts_i, shifts_s, parameterDirections);
  #endif

  // size update direction vector correctly
  //for (int i = 0; i < EngineObj->good_solve().size(); i++) 
  //  app_log() << EngineObj->good_solve().at(i) << "  ";
  //app_log() << endl;
  parameterDirections.resize(shifts_i.size());
  for (int i = 0; i < shifts_i.size(); i++) {
    parameterDirections.at(i).assign(N, 0.0);
    if ( true ) {
      for (int j = 0; j < N; j++) 
        parameterDirections.at(i).at(j) = EngineObj->wfn_update().at(i*N+j);
    }
    else 
      parameterDirections.at(i).at(0) = 1.0;
  }

  // now that we are done with them, prevent further computation of derivative vectors
  this->optTarget->setneedGrads(false);

  // prepare vectors to hold the initial and current parameters
  std::vector<RealType> currParams(numParams, 0.0);

  // initialize the initial and current parameter vectors
  for (int i=0; i<numParams; i++)
    currParams.at(i) = this->optTarget->Params(i);

  // create a vector telling which updates are within our constraints
  std::vector<bool> good_update(parameterDirections.size(), true);

  // compute the largest parameter change for each shift, and zero out updates that have too-large changes
  std::vector<RealType> max_change(parameterDirections.size(), 0.0);
  for (int k = 0; k < parameterDirections.size(); k++) {
    for ( int i = 0; i < numParams; i++)
      max_change.at(k) = std::max(max_change.at(k), std::abs(parameterDirections.at(k).at(i+1) / parameterDirections.at(k).at(0)));
    good_update.at(k) = ( good_update.at(k) && max_change.at(k) <= max_param_change );
  }

  // prepare to use the middle shift's update as the guiding function for a new sample
  for (int i=0; i<numParams; i++)
    this->optTarget->Params(i) = currParams.at(i) + parameterDirections.at(central_index).at(i+1);

  // say what we are doing
  app_log() << std::endl
            << "************************************************************" << std::endl
            << "Updating the guiding function with the middle shift's update" << std::endl
            << "************************************************************" << std::endl
            << std::endl;

  // generate the new sample on which we will compare the different shifts
  this->finish();

  // reset the number of samples
  //this->optTarget->setNumSamples(nsamp_comp);
  //nTargetSamples = nsamp_comp;
  //app_log() << "# of sample before correlated sampling is " << nTargetSamples << std::endl;
  //app_log() << "number of samples is" << this->optTarget->getNumSamples() << std::endl;
  app_log() << std::endl
            << "*************************************************************" << std::endl
            << "Generating a new sample based on the updated guiding function" << std::endl
            << "*************************************************************" << std::endl
            << std::endl;
  this->start();
  //app_log() << "number of samples is" << this->optTarget->getNumSamples() << std::endl;

  // say what we are doing
  app_log() << std::endl
            << "******************************************************************" << std::endl
            << "Comparing different shifts' cost function values on updated sample" << std::endl
            << "******************************************************************" << std::endl
            << std::endl;

  // update the current parameters to those of the new guiding function
  for (int i=0; i<numParams; i++)
    currParams.at(i) = optTarget->Params(i);

  // compute cost function for the initial parameters (by subtracting the middle shift's update back off)
  for (int i=0; i<numParams; i++)
    this->optTarget->Params(i) = currParams.at(i) - parameterDirections.at(central_index).at(i+1);
  this->optTarget->IsValid = true;
  const RealType initCost = this->optTarget->LMYEngineCost(false, EngineObj);

  // compute the update directions for the smaller and larger shifts relative to that of the middle shift
  for (int i=0; i<numParams; i++) {
    for (int j = 0; j < parameterDirections.size(); j++) {
      if ( j != central_index ) 
        parameterDirections.at(j).at(i+1) -= parameterDirections.at(1).at(i+1);
    //parameterDirections.at(2).at(i+1) -= parameterDirections.at(1).at(i+1);
    }
  }

  // prepare a vector to hold the cost function value for each different shift
  std::vector<RealType> costValues(num_shifts, 0.0);

  // compute the cost function value for each shift and make sure the change is within our constraints
  for (int k = 0; k < parameterDirections.size(); k++) {
    for (int i=0; i<numParams; i++)
      this->optTarget->Params(i) = currParams.at(i) + ( k == num_shifts ? 0.0 : parameterDirections.at(k).at(i+1) );
    this->optTarget->IsValid = true;
    costValues.at(k) = this->optTarget->LMYEngineCost(false, EngineObj);
    good_update.at(k) = ( good_update.at(k) && std::abs( (initCost - costValues.at(k)) / initCost ) < max_relative_cost_change );
    //app_log() << std::abs( (starting_cost - costValues.at(k)) / starting_cost ) << "  ";
    if (!good_update.at(k))
      costValues.at(k) = std::abs(1.5*initCost) + 1.0;
  }
  //app_log() << endl;

  //for (int i = 0; i < good_update.size(); i++)
  //  app_log() << good_update.at(i) << "  ";
  //app_log() << endl;

  // find the best shift and the corresponding update direction
  const std::vector<RealType> * bestDirection = 0;
  int best_shift = -1;
  for (int k = 0; k < costValues.size() && std::abs( (initCost - initCost) / initCost ) < max_relative_cost_change; k++)
    if ( this->is_best_cost(k, costValues, initCost) && good_update.at(k) ) {
      best_shift = k;
      bestDirection = &parameterDirections.at(k);
    }

  // print the results for each shift
  app_log() << std::endl;
  this->print_cost_summary_header();
  this->print_cost_summary(0.0, 0.0, 0.0, initCost, -1, best_shift, true);
  for (int k = 0; k < good_update.size(); k++)
    this->print_cost_summary(shifts_i.at(k), shifts_s.at(k), max_change.at(k), costValues.at(k), k, best_shift, good_update.at(k));

  // if any of the shifts produced a good update, apply the best such update and remember those shifts for next time
  if ( bestDirection ) {

    this->bestShift_i = shifts_i.at(best_shift);
    this->bestShift_s = shifts_s.at(best_shift);
    for (int i=0; i<numParams; i++)
      optTarget->Params(i) = currParams.at(i) + ( best_shift == central_index ? 0.0 : bestDirection->at(i+1) );
    app_log() << std::endl
              << "*****************************************************************************" << std::endl
              << "Applying the update for shift_i = "
              << std::scientific << std::right << std::setw(12) << std::setprecision(4) << this->bestShift_i
              << "     and shift_s = "
              << std::scientific << std::right << std::setw(12) << std::setprecision(4) << this->bestShift_s
              << std::endl
              << "*****************************************************************************" << std::endl
              << std::endl;

  // otherwise revert to the old parameters and set the next shift to be larger
  } else {
    this->bestShift_i *= 10.0;
    this->bestShift_s *= 10.0;
    for (int i=0; i<numParams; i++)
      optTarget->Params(i) = currParams.at(i) - parameterDirections.at(central_index).at(i+1);
    app_log() << std::endl
              << "***********************************************************" << std::endl
              << "Reverting to old parameters and increasing shift magnitudes" << std::endl
              << "***********************************************************" << std::endl
              << std::endl;
  }

  // save the update for future linear method iterations
  if ( block_lm && bestDirection ) {
    
    // save the difference between the updated and old variables
    formic::ColVec<double> update_dirs(numParams, 0.0);
    for (int i = 0; i < numParams; i++) 
      update_dirs.at(i) = bestDirection->at(i+1) + parameterDirections.at(central_index).at(i+1);
    previous_update.insert(previous_update.begin(), update_dirs);
    
    // eliminate the oldest saved update if necessary
    while ( previous_update.size() > nolds )
      previous_update.pop_back();
  }

  // return the cost function grads flag to what it was
  this->optTarget->setneedGrads(saved_grads_flag);

  // perform some finishing touches for this linear method iteration
  this->finish();

  // set the number samples to be initial one
  this->optTarget->setNumSamples(init_num_samp);
  nTargetSamples = init_num_samp;

  //app_log() << "block first second third end " << block_first << block_second << block_third << endl; 
  // return whether the cost function's report counter is positive
  return (optTarget->getReportCounter() > 0);

}
#endif

bool QMCFixedSampleLinearOptimize::one_shift_run() {

  // ensure the cost function is set to compute derivative vectors
  optTarget->setneedGrads(true);

  // generate samples and compute weights, local energies, and derivative vectors
  start();

  // get number of optimizable parameters
  numParams = optTarget->NumParams();

  // get dimension of the linear method matrices
  N = numParams + 1;

  // prepare vectors to hold the initial and current parameters
  std::vector<RealType> currentParameters(numParams, 0.0);

  // initialize the initial and current parameter vectors
  for (int i=0; i<numParams; i++)
    currentParameters.at(i) = optTarget->Params(i);

  // prepare vectors to hold the parameter update directions for each shift
  std::vector<RealType> parameterDirections;
  parameterDirections.assign(N, 0.0);

  // compute the initial cost
#ifdef QMC_CUDA
  // Ye : can't call computedCost directly, internal data was not correct for ham,ovl matrices.
  // more investiation is needed.
  const RealType initCost = optTarget->Cost(true);
#else
  const RealType initCost = optTarget->computedCost();
#endif

  // say what we are doing
  app_log() << std::endl
            << "*****************************************" << std::endl
            << "Building overlap and Hamiltonian matrices" << std::endl
            << "*****************************************" << std::endl;

  // allocate the matrices we will need
  Matrix<RealType> ovlMat(N,N); ovlMat = 0.0;
  Matrix<RealType> hamMat(N,N); hamMat = 0.0;
  Matrix<RealType> invMat(N,N); invMat = 0.0;
  Matrix<RealType> prdMat(N,N); prdMat = 0.0;

  // build the overlap and hamiltonian matrices
  optTarget->fillOverlapHamiltonianMatrices(hamMat, ovlMat);
  invMat.copy(ovlMat);

  // prepare vector to hold largest parameter change for each shift
  RealType max_change(0.0);

  // apply the identity shift
  for (int i=1; i<N; i++)
  {
    hamMat(i,i) += bestShift_i;
    if(invMat(i,i)==0) invMat(i,i) = bestShift_i*bestShift_s;
  }

  // compute the inverse of the overlap matrix
  invert_matrix(invMat, false);

  // apply the overlap shift
  for (int i=1; i<N; i++)
    for (int j=1; j<N; j++)
      hamMat(i,j) += bestShift_s * ovlMat(i,j);

  // multiply the shifted hamiltonian matrix by the inverse of the overlap matrix
  qmcplusplus::MatrixOperators::product(invMat, hamMat, prdMat);

  // transpose the result (why?)
  for (int i=0; i<N; i++)
    for (int j=i+1; j<N; j++)
      std::swap(prdMat(i,j), prdMat(j,i));

  // compute the lowest eigenvalue of the product matrix and the corresponding eigenvector
  const RealType lowestEV = getLowestEigenvector(prdMat, parameterDirections);

  // compute the scaling constant to apply to the update
  Lambda = getNonLinearRescale(parameterDirections, ovlMat);

  // scale the update by the scaling constant
  for (int i=0; i<numParams; i++)
    parameterDirections.at(i+1) *= Lambda;

  // now that we are done building the matrices, prevent further computation of derivative vectors
  optTarget->setneedGrads(false);

  // prepare to use the middle shift's update as the guiding function for a new sample
  for (int i=0; i<numParams; i++)
    optTarget->Params(i) = currentParameters.at(i) + parameterDirections.at(i+1);

  RealType bigVec(0);
  for (int i=0; i<numParams; i++)
    bigVec = std::max(bigVec,std::abs(parameterDirections.at(i+1)));
  app_log() << std::endl << "Among totally " << numParams << " optimized parameters, "
            << "largest LM parameter change : "
            << bigVec << std::endl;

  // compute the new cost
  optTarget->IsValid = true;
  const RealType newCost = optTarget->Cost(false);

  app_log() << std::endl
            << "******************************************************************************" << std::endl
            << "Init Cost = "
            << std::scientific << std::right << std::setw(12) << std::setprecision(4) << initCost
            << "    New Cost = "
            << std::scientific << std::right << std::setw(12) << std::setprecision(4) << newCost
            << "  Delta Cost = "
            << std::scientific << std::right << std::setw(12) << std::setprecision(4) << newCost - initCost
            << std::endl
            << "******************************************************************************" << std::endl;

  if ( !optTarget->IsValid || std::isnan(newCost)) {
    app_log() << std::endl << "The new set of parameters is not valid. Revert to the old set!" << std::endl;
    for (int i=0; i<numParams; i++)
      optTarget->Params(i) = currentParameters.at(i);
    bestShift_s=bestShift_s*4.0;
  } else {
    if ( bestShift_s > 1.0e-2 ) bestShift_s=bestShift_s/4.0;
    // say what we are doing
    app_log() << std::endl << "The new set of parameters is valid. Updating the trial wave function!" << std::endl;
  }

  app_log() << std::endl
            << "*****************************************************************************" << std::endl
            << "Applying the update for shift_i = "
            << std::scientific << std::right << std::setw(12) << std::setprecision(4) << bestShift_i
            << "     and shift_s = "
            << std::scientific << std::right << std::setw(12) << std::setprecision(4) << bestShift_s
            << std::endl
            << "*****************************************************************************" << std::endl;

  // perform some finishing touches for this linear method iteration
  finish();

  // return whether the cost function's report counter is positive
  return (optTarget->getReportCounter() > 0);

}

}
