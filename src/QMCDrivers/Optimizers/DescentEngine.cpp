
//Prototype code for a descent engine
//Some functions modeled on those in linear method engine, but trying to avoid use of formic wrappers for things like mpi.
//Long-term(?) amibition is to eventually get rid of most of formic except for parts essential for the adaptive three shift LM and BLM.


#include <cmath>
#include <vector>
#include <string>
#include "QMCDrivers/Optimizers/DescentEngine.h"
#include "Message/CommOperators.h"

namespace qmcplusplus
{
DescentEngine::DescentEngine(Communicate* comm, const xmlNodePtr cur)
    : engineTargetExcited(false),
      myComm(comm),
      excited("no"),
      flavor("RMSprop"),
      TJF_2Body_eta(.01),
      TJF_1Body_eta(.01),
      F_eta(.001),
      Gauss_eta(.001),
      CI_eta(.01),
      Orb_eta(.001),
      ramp_eta(false),
      ramp_num(30)
{

  m_param.add(excited, "targetExcited", "string");
  
  //Type of descent method being used
  m_param.add(flavor, "flavor", "string");
  m_param.add(TJF_2Body_eta, "TJF_2Body_eta", "double");
   
  
  m_param.add(TJF_1Body_eta, "TJF_1Body_eta", "double");
  m_param.add(F_eta, "F_eta", "double");
 m_param.add(CI_eta, "CI_eta", "double");
  m_param.add(Gauss_eta, "Gauss_eta", "double");
  m_param.add(Orb_eta, "Orb_eta", "double");
  m_param.add(ramp_etaStr,"Ramp_eta","string");
  m_param.add(ramp_num,"Ramp_num","int");

  processXML(cur);
}


bool DescentEngine::processXML(const xmlNodePtr cur)
{
 
 m_param.put(cur);

  engineTargetExcited = (excited == "yes");
 
 
  ramp_eta = (ramp_etaStr == "yes");

  return true;
}

void DescentEngine::clear_samples(const size_t numOptimizables)
{
  avg_le_der_samp.resize(numOptimizables);
  avg_der_rat_samp.resize(numOptimizables);
  LDerivs.resize(numOptimizables);

  numParams = numOptimizables;

  std::fill(avg_le_der_samp.begin(), avg_le_der_samp.end(), 0.0);
  std::fill(avg_der_rat_samp.begin(), avg_der_rat_samp.end(), 0.0);

  w_sum       = 0;
  e_avg       = 0;
  e_sum       = 0;
  eSquare_sum = 0;
  eSquare_avg = 0;
}

void DescentEngine::setEtemp(std::vector<double> etemp)
{
  e_sum       = etemp[0];
  w_sum       = etemp[1];
  eSquare_sum = etemp[2];
  e_avg       = e_sum / w_sum;
  eSquare_avg = eSquare_sum / w_sum;

  int my_rank = myComm->rank();
  if (my_rank == 0)
  {
    std::cout << "e_sum: " << e_sum << std::endl;
    std::cout << "w_sum: " << w_sum << std::endl;
    std::cout << "e_avg: " << e_avg << std::endl;
    std::cout << "eSquare_sum: " << eSquare_sum << std::endl;
    std::cout << "eSquare_avg: " << eSquare_avg << std::endl;
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
///
/// \param[in]  der_rat_samp   <n|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  le_der_samp    <n|H|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  ls_der_samp    <|S^2|Psi_i>/<n|Psi> (i = 0 (|Psi>), 1, ... N_var )
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::take_sample(std::vector<double>& der_rat_samp,
                                std::vector<double>& le_der_samp,
                                std::vector<double>& ls_der_samp,
                                double vgs_samp,
                                double weight_samp)
{
  const size_t numOptimizables = der_rat_samp.size() - 1;


  for (int i = 0; i < numOptimizables; i++)
  {
    avg_le_der_samp.at(i) += le_der_samp.at(i + 1);
    avg_der_rat_samp.at(i) += der_rat_samp.at(i + 1);
  }
}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that Take Sample Data from the Host Code
///
/// \param[in]  local_en       local energy
/// \param[in]  vgs_samp       |<n|value_fn>/<n|guiding_fn>|^2
/// \param[in]  weight_samp    weight for this sample
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::take_sample(double local_en, double vgs_samp, double weight_samp) {}

////////////////////////////////////////////////////////////////////////////////////////////////////////
/// \brief  Function that reduces all vector information from all processors to the root
///         processor
///
////////////////////////////////////////////////////////////////////////////////////////////////////////
void DescentEngine::sample_finish()
{
  int my_rank = myComm->rank();

  myComm->allreduce(avg_le_der_samp);
  myComm->allreduce(avg_der_rat_samp);

  for (int i = 0; i < LDerivs.size(); i++)
  {
    avg_le_der_samp.at(i)  = avg_le_der_samp.at(i) / w_sum;
    avg_der_rat_samp.at(i) = avg_der_rat_samp.at(i) / w_sum;

    if (my_rank == 0)
    {
      std::cout << "Parameter # " << i << " Hamiltonian term: " << avg_le_der_samp.at(i) << std::endl;
      std::cout << "Parameter # " << i << " Overlap term: " << avg_der_rat_samp.at(i) << std::endl;
    }

    if (!engineTargetExcited)
    {
      LDerivs.at(i) = 2 * avg_le_der_samp.at(i) - e_avg * (2 * avg_der_rat_samp.at(i));
    }
  }
}


//Function for updating parameters during descent optimization
void DescentEngine::updateParameters(int stepNum, int descentNum)
{
  int my_rank = myComm->rank();
  if (my_rank == 0)
  {
    std::cout << "Number of Parameters: " << numParams << std::endl;

    std::cout << "Parameter Type step sizes: "
              << " TJF_2Body_eta=" << TJF_2Body_eta << " TJF_1Body_eta=" << TJF_1Body_eta << " F_eta=" << F_eta
              << " CI_eta=" << CI_eta << " Orb_eta=" << Orb_eta << std::endl;
  }

  // Get set of derivatives for current (kth) optimization step
  std::vector<double> curDerivSet = derivRecords.at(derivRecords.size() - 1);
  std::vector<double> prevDerivSet;

  if (!taus.empty())
  {
    // Get set of derivatives for previous (k-1th) optimization step
    prevDerivSet = derivRecords.at(derivRecords.size() - 2);
  }

  double denom;
  double numer;
  double v;
  double corNumer;
  double corV;

  double epsilon = 1e-8;
  double type_Eta;

  double tau;
  // Update parameters according to specified flavor of gradient descent method

  // RMSprop corresponds to the method used by Booth and co-workers
  if (flavor.compare("RMSprop") == 0)
  {
    if (my_rank == 0)
      std::cout << "Using RMSprop" << std::endl;

    // To match up with Booth group paper notation, prevLambda is lambda_k-1,
    // curLambda is lambda_k, nextLambda is lambda_k+1
    double curLambda  = .5 + .5 * std::sqrt(1 + 4 * lambda * lambda);
    double nextLambda = .5 + .5 * std::sqrt(1 + 4 * curLambda * curLambda);
    double gamma      = (1 - curLambda) / nextLambda;

    // Define damping factor that turns off acceleration of the algorithm
    // small value of d corresponds to quick damping and effectively using
    // steepest descent
    double d           = 100;
    double decayFactor = std::exp(-(1 / d) * (stepNum));
    gamma              = gamma * decayFactor;

    double rho = .9;

    for (int i = 0; i < numParams; i++)
    {
    
      double curSquare = std::pow(curDerivSet.at(i), 2);

      // Need to calculate step size tau for each parameter inside loop
      // In RMSprop, the denominator of the step size depends on a a running average of past squares of the parameter derivative
      if (derivsSquared.size() < numParams)
      {
        curSquare = std::pow(curDerivSet.at(i), 2);
      }
      else if (derivsSquared.size() >= numParams)
      {
        curSquare = rho * derivsSquared.at(i) + (1 - rho) * std::pow(curDerivSet.at(i), 2);
      }

      denom = std::sqrt(curSquare + epsilon);

      //The numerator of the step size is set according to parameter type based on input choices
      type_Eta = this->setStepSize(i);
      tau      = type_Eta / denom;

      // Include an additional factor to cause step size to eventually decrease to 0 as number of steps taken increases
      double stepLambda = .1;

      double stepDecayDenom = 1 + stepLambda * stepNum;
      tau                   = tau / stepDecayDenom;


      //Update parameter values
      //If case corresponds to being after the first descent step
      if (taus.size() >= numParams)
      {
        double oldTau = taus.at(i);

        currentParams.at(i) = (1 - gamma) * (currentParams.at(i) - tau * curDerivSet.at(i)) +
            gamma * (paramsCopy.at(i) - oldTau * prevDerivSet.at(i));
      }
      else
      {
        tau = type_Eta;

        currentParams.at(i) = currentParams.at(i) - tau * curDerivSet.at(i);
      }

      if (taus.size() < numParams)
      {
        // For the first optimization step, need to add to the vectors
        taus.push_back(tau);
        derivsSquared.push_back(curSquare);
      }
      else
      {
        // When not on the first step, can overwrite the previous stored values
        taus[i]          = tau;
        derivsSquared[i] = curSquare;
      }

      paramsCopy[i] = currentParams[i];
    }

    // Store current (kth) lambda value for next optimization step
    lambda = curLambda;
  }
  // Random uses only the sign of the parameter derivatives and takes a step of random size within a range.
  else if (flavor.compare("Random") == 0)
  {
    if (my_rank == 0)
      std::cout << "Using Random" << std::endl;

    for (int i = 0; i < numParams; i++)
    {
      denom        = 1;
      double alpha = ((double)rand() / RAND_MAX);
      double sign  = std::abs(curDerivSet[i]) / curDerivSet[i];
      if (std::isnan(sign) && my_rank == 0)
      {
        std::cout << "Got a nan, choosing sign randomly with 50-50 probability" << std::endl;

        double t = ((double)rand() / RAND_MAX);
        if (t > .5)
        {
          sign = 1;
        }
        else
        {
          sign = -1;
        }
      }
      if (my_rank == 0)
        std::cout << "This is random alpha: " << alpha << " with sign: " << sign << std::endl;

      currentParams.at(i) = currentParams.at(i) - tau * alpha * sign;
    }
  }

  else
  {
    // ADAM method
    if (flavor.compare("ADAM") == 0)
    {
      if (my_rank == 0)
        std::cout << "Using ADAM" << std::endl;

      for (int i = 0; i < numParams; i++)
      {
        double curSquare = std::pow(curDerivSet.at(i), 2);
        double beta1     = .9;
        double beta2     = .99;
        if (descentNum == 0)
        {
          numerRecords.push_back(0);
          denomRecords.push_back(0);
        }
        numer = beta1 * numerRecords[i] + (1 - beta1) * curDerivSet[i];
        v     = beta2 * denomRecords[i] + (1 - beta2) * curSquare;

        corNumer = numer / (1 - std::pow(beta1, descentNum + 1));
        corV     = v / (1 - std::pow(beta2, descentNum + 1));

        denom = std::sqrt(corV) + epsilon;

        type_Eta = this->setStepSize(i);
        tau      = type_Eta / denom;

        currentParams.at(i) = currentParams.at(i) - tau * corNumer;

        if (taus.size() < numParams)
        {
          // For the first optimization step, need to add to the vectors
          taus.push_back(tau);
          derivsSquared.push_back(curSquare);
          denomRecords[i] = v;
          numerRecords[i] = numer;
        }
        else
        {
          // When not on the first step, can overwrite the previous stored values
          taus[i]          = tau;
          derivsSquared[i] = curSquare;
          denomRecords[i]  = v;
          numerRecords[i]  = numer;
        }

        paramsCopy[i] = currentParams.at(i);
      }
    }
    // AMSGrad method, similar to ADAM except for form of the step size denominator
    else if (flavor.compare("AMSGrad") == 0)
    {
      if (my_rank == 0)
        std::cout << "Using AMSGrad" << std::endl;


      for (int i = 0; i < numParams; i++)
      {
        double curSquare = std::pow(curDerivSet.at(i), 2);
        double beta1     = .9;
        double beta2     = .99;
        if (descentNum == 0)
        {
          numerRecords.push_back(0);
          denomRecords.push_back(0);
        }

        numer = beta1 * numerRecords[i] + (1 - beta1) * curDerivSet[i];
        v     = beta2 * denomRecords[i] + (1 - beta2) * curSquare;
        v     = std::max(denomRecords[i], v);

        denom    = std::sqrt(v) + epsilon;
        type_Eta = this->setStepSize(i);
        tau      = type_Eta / denom;

        currentParams.at(i) = currentParams.at(i) - tau * numer;

        if (taus.size() < numParams)
        {
          // For the first optimization step, need to add to the vectors
          taus.push_back(tau);
          derivsSquared.push_back(curSquare);
          denomRecords[i] = v;
          numerRecords[i] = numer;
        }
        else
        {
          // When not on the first step, can overwrite the previous stored values
          taus[i]          = tau;
          derivsSquared[i] = curSquare;
          denomRecords[i]  = v;
          numerRecords[i]  = numer;
        }

        paramsCopy[i] = currentParams.at(i);
      }
    }
  }


  /*
  //During the hybrid method,store 5 vectors of parameter differences over the course of a descent section
  if (doHybrid && ((descentNum + 1) % (descent_len / 5) == 0)) {
      app_log() << "Step number in macro-iteration is " << stepNum % descent_len
                << " out of expected total of " << descent_len
                << " descent steps." << std::endl;
    storeVectors(paramsForDiff);
  }
  */
}

// Helper method for setting step size according parameter type.
double DescentEngine::setStepSize(int i)
{
  double type_eta;


  std::string name = engineParamNames[i];
  

  int type = engineParamTypes[i];

  //Step sizes are assigned according to parameter type identified from the variable name.
  //Other parameter types could be added to this section as other wave function ansatzes are developed.
  if ((name.find("uu") != std::string::npos) || (name.find("ud") != std::string::npos))
  {
    type_eta = TJF_2Body_eta;
  }
  //If parameter name doesn't have "uu" or "ud" in it and is of type 1, assume it is a 1 body Jastrow parameter.
  else if (type == 1)
  {
    type_eta = TJF_1Body_eta;
  }
  else if (name.find("F_") != std::string::npos)
  {
    type_eta = F_eta;
  }
  else if (name.find("CIcoeff_") != std::string::npos || name.find("CSFcoeff_") != std::string::npos)
  {
    type_eta = CI_eta;
  }
  else if (name.find("orb_rot_") != std::string::npos)
  {
    type_eta = Orb_eta;
  }
  else if (name.find("g") != std::string::npos)
  {
    //Gaussian parameters are rarely optimized in practice but the descent code allows for it.
    type_eta = Gauss_eta;
  }
  else
  {
    //If there is some other parameter type that isn't in one of the categories with a default/input, use a conservative default step size.
    type_eta = .001;
  }

  if (ramp_eta && descentNum < ramp_num)
  {
    type_eta = type_eta * (descentNum + 1) / ramp_num;
  }

  return type_eta;
}

void DescentEngine::setupUpdate(int& paramNum,
                                std::vector<std::string>& paramNames,
                                std::vector<int>& paramTypes,
                                std::vector<double>& initialParams)
{
  numParams = paramNum;


  for (int i = 0; i < numParams; i++)
  {
    engineParamNames.push_back(paramNames[i]);
    engineParamTypes.push_back(paramTypes[i]);
    paramsCopy.push_back(initialParams[i]);
    currentParams.push_back(initialParams[i]);
  }
}


} // namespace qmcplusplus
