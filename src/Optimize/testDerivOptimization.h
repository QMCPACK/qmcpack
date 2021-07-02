//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Mark Dewing, mdewing@anl.gov, Argonne National Laboratory
//
// File created by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_TESTDERIVOPTIMIZATION_H
#define QMCPLUSPLUS_TESTDERIVOPTIMIZATION_H

#include "Optimize/OptimizeBase.h"
#include "OhmmsData/ParameterSet.h"
#include <fstream>

using qmcplusplus::app_log;

/**
 * @file testDerivOptimization.h
 * Tests for variational parameter derivatives
 *
 * This optimization type will compare finite difference derivatives with the analytic derivatives.
 * It can also output a file that can be used with 'qmca' to get error bars.
 *
 * The input is specified with the 'optimize' tag and the 'test' method.
 * For example:
 * \code
 *  <loop max="4">
 *   <qmc method="opt" move="pbyp">
 *     <optimize method="test" output_param_file="yes">
 *     </optimize>
 *     ...
 *   </qmc>
 *  </loop>
 * \endcode
 *
 * The output_param_file is optional and defaults to 'no'.
 * If 'yes', a file named "<project id>.param.s000.scalar.dat" is created.  Each iteration of the optimizer
 * loop outputs one line in the file. (The above example will produce 4 entries in the file).
 * This is a hack to enable computing error bars on the parameter gradients.
 *
 */

template<class T>
class testDerivOptimization : public MinimizerBase<T>
{
public:
  typedef T Return_t;
  typedef typename MinimizerBase<T>::ObjectFuncType ObjectFuncType;

  ObjectFuncType* TargetFunc;

  testDerivOptimization(const std::string& RootName, ObjectFuncType* atarget = 0)
      : TargetFunc(atarget), first_(true), output_param_file_(false), param_deriv_index_(0), RootName_(RootName)
  {
    if (atarget)
      setTarget(atarget);
  }

  ~testDerivOptimization() override {}

  void setTarget(ObjectFuncType* fn)
  {
    TargetFunc = fn;
    NumParams_ = TargetFunc->getNumParams();
    resizeAllArray(NumParams_);
    for (int i = 0; i < NumParams_; i++)
      Parms_[i] = std::real(TargetFunc->Params(i));

    if (output_param_file_ && first_)
    {
      int namelen = RootName_.length();
      // Assume that the RootName has the series suffix (".s000").
      // Remove the series suffix to get the project id
      std::string fname = RootName_.substr(0, namelen - 5) + ".param.s000.scalar.dat";
      param_deriv_file_.open(fname);
      param_deriv_file_ << "# Index ";
      for (int i = 0; i < NumParams_; i++)
      {
        param_deriv_file_ << " " << TargetFunc->getParamName(i);
      }
      param_deriv_file_ << std::endl;
      first_ = false;
    }
  }

  void resizeAllArray(int newSize)
  {
    a_xi_.resize(newSize, 0.0);
    Parms_.resize(newSize, 0.0);
  }

  bool optimize(ObjectFuncType* fn) override
  {
    setTarget(fn);
    return optimize();
  }

  bool optimize()
  {
    dfunc(Parms_, a_xi_);
    //make sure wave function has the right parameters for next runs
    for (int i = 0; i < NumParams_; i++)
      TargetFunc->Params(i) = Parms_[i];
    TargetFunc->Report();
    return true;
  }

  bool get(std::ostream&) const;

  bool put(std::istream&);

  /**  Parse the xml file for parameters
   *
   */

  bool put(xmlNodePtr cur) override
  {
    std::string output_file("no");
    OhmmsAttributeSet attrib;
    attrib.add(output_file, "output_param_file");
    attrib.put(cur);

    output_param_file_ = (output_file != "no");

    return true;
  }

  Return_t func(std::vector<Return_t> _p)
  {
    for (int i = 0; i < NumParams_; ++i)
      TargetFunc->Params(i) = _p[i];
    return TargetFunc->Cost();
  }


  void dfunc(const std::vector<Return_t>& RT, std::vector<Return_t>& FG)
  {
    ///To test we simply output the analytic and numeric gradients of the cost function. Make sure they agree.
    std::vector<Return_t> Dummy(FG);
    TargetFunc->GradCost(Dummy, RT, 1e-5);
    TargetFunc->GradCost(FG, RT, 0.0);

    if (output_param_file_)
    {
      param_deriv_file_ << param_deriv_index_ << " ";
      param_deriv_index_++;
    }

    app_log() << "Param_Name  Value    Numeric            Analytic       Percent" << std::endl;
    for (int k = 0; k < NumParams_; k++)
    {
      std::string vname = TargetFunc->getParamName(k);
      if (Dummy[k] != 0)
        app_log() << vname << " " << RT[k] << "  " << Dummy[k] << "  " << FG[k] << "  "
                  << 100 * (Dummy[k] - FG[k]) / Dummy[k] << std::endl;
      else
        app_log() << vname << " " << RT[k] << "  " << Dummy[k] << "  " << FG[k] << "   inf" << std::endl;
      if (output_param_file_)
        param_deriv_file_ << std::setprecision(10) << FG[k] << " ";
    }
    if (output_param_file_)
      param_deriv_file_ << std::endl;
    app_log() << std::endl;
  }

private:
  int NumParams_;
  std::vector<Return_t> a_xi_;
  std::vector<Return_t> Parms_;

  std::ofstream param_deriv_file_;
  bool first_;
  bool output_param_file_;
  int param_deriv_index_;
  std::string RootName_;
};

#endif
