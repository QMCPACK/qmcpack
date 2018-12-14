
#include "Configuration.h"
#include "Numerics/OhmmsBlas.h" 
#include <Utilities/FairDivide.h>
#include "AFQMC/Utilities/Utils.h"
#include "AFQMC/config.h" 
#include "AFQMC/Propagators/AFQMCSharedPropagator.h"
#include "AFQMC/Walkers/WalkerConfig.hpp"

namespace qmcplusplus 
{

namespace afqmc
{

void AFQMCSharedPropagator::parse(xmlNodePtr cur)
{
  using qmcplusplus::app_log;

  if(cur == NULL)
    APP_ABORT("Error in AFQMCSharedPropagator::parse(): Null xml ptr.\n");

  // defaults
  vbias_bound=50.0;
  free_projection=false;
  hybrid=true;
  importance_sampling=true;
  apply_constrain=true;

  xmlNodePtr curRoot=cur;

  std::string sub("yes");
  std::string constrain("yes");
  std::string hyb("yes");
  std::string freep("no");
  std::string impsam("yes");
  ParameterSet m_param;
  m_param.add(vbias_bound,"vbias_bound","double");
  m_param.add(constrain,"apply_constrain","std::string");
  m_param.add(impsam,"importance_sampling","std::string");
  m_param.add(hyb,"hybrid","std::string");
  m_param.add(freep,"free_projection","std::string");

  //m_param.add(sz_pin_field_file,"sz_pinning_field_file","std::string");
  //m_param.add(sz_pin_field_mag,"sz_pinning_field","double");

  m_param.put(cur);

  std::transform(constrain.begin(),constrain.end(),constrain.begin(),(int (*)(int)) tolower);
  if(constrain == "no" || constrain == "false") apply_constrain = false;
  std::transform(impsam.begin(),impsam.end(),impsam.begin(),(int (*)(int)) tolower);
  if(impsam == "false" || impsam == "no") importance_sampling = false;
  std::transform(hyb.begin(),hyb.end(),hyb.begin(),(int (*)(int)) tolower);
  if(hyb == "no" || hyb == "false") hybrid = false;
  std::transform(freep.begin(),freep.end(),freep.begin(),(int (*)(int)) tolower);
  if(freep == "yes" || freep == "true") free_projection=true; 

    app_log()<<"\n\n --------------- Parsing Propagator input ------------------ \n\n";

  if(free_projection) {

    //if (tau_release > 0) {
      //app_log()<<" AFQMC with free projection after tau = " << tau_release << "\n"
             //<<" WARNING: Will Set: \n"
             //<<"               -apply_constrain     = no \n"
             //<<"               -importance_sampling = no \n"
             //<<" after this time.\n"
             //<<" Setting: \n"
             //<<"               -hybrid       = yes \n"
             //<<"               -importance_sampling = yes\n"
             //<<"               -apply_constrain     = yes\n"
             //<<"               -free_projection     = no\n"
             //<<" now.\n";
      //importance_sampling = true;
      //apply_constrain = true;
      //hybrid = true;
      //free_projection = false;
    //} else {
      app_log()<<" AFQMC with free projection. \n"
             <<" WARNING: Setting: \n"
             <<"               -apply_constrain     = no \n"
             <<"               -importance_sampling = no \n"
             <<"               -hybrid       = yes \n";
      importance_sampling=false;
      hybrid=true;
      apply_constrain=false;
    //}

  } else { 

    if(!importance_sampling) 
      app_log()<<" WARNING: importance_sampling=no without free projection does not make sense. \n"; 

    if(hybrid)
      app_log()<<" Using hybrid method to calculate the weights during the propagation." <<"\n";
    else
      app_log()<<" Using local energy method to calculate the weights during the propagation." <<"\n";

  }

  app_log()<<std::endl;

}

void AFQMCSharedPropagator::reset_nextra(int nextra) 
{
  if(nextra==0) return;
  if(last_nextra != nextra ) {
    last_nextra = nextra;
    for(int n=0; n<nextra; n++) {
      int n0,n1;
      std::tie(n0,n1) = FairDivideBoundary(n,TG.getNCoresPerTG(),nextra);
      if(TG.getLocalTGRank()>=n0 && TG.getLocalTGRank()<n1) {
        last_task_index = n;
        break;
      }
    }
    local_group_comm = std::move(shared_communicator(TG.TG_local().split(last_task_index)));
  }
  if(last_task_index < 0 || last_task_index >= nextra)
    APP_ABORT("Error: Problems in AFQMCSharedPropagator::reset_nextra()\n");
}

void AFQMCSharedPropagator::assemble_X(size_t nsteps, size_t nwalk, RealType sqrtdt, 
                          CMatrix_ref& X, CMatrix_ref & vbias, CMatrix_ref& MF, 
                          CMatrix_ref& HWs, bool addRAND) 
{
  // remember to call vbi = apply_bound_vbias(*vb);
  // X[m,ni,iw] = rand[m,ni,iw] + im * ( vbias[m,iw] - vMF[m]  )
  // HW[ni,iw] = sum_m [ im * ( vMF[m] - vbias[m,iw] ) * 
  //                     ( rand[m,ni,iw] + halfim * ( vbias[m,iw] - vMF[m] ) ) ] 
  //           = sum_m [ im * ( vMF[m] - vbias[m,iw] ) * 
  //                     ( X[m,ni,iw] - halfim * ( vbias[m,iw] - vMF[m] ) ) ] 
  // MF[ni,iw] = sum_m ( im * X[m,ni,iw] * vMF[m] )   
 
  TG.local_barrier();
  ComplexType im(0.0,1.0);
  ComplexType halfim(0.0,0.5);
  int nCV = int(X.shape()[0]);
  boost::multi_array_ref<ComplexType,3> X3D(X.origin(),extents[X.shape()[0]][nsteps][nwalk]);
  // generate random numbers
  if(addRAND)
  {
    int i0,iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),int(X.num_elements()),
                                         TG.TG_local().size());  
    sampleGaussianFields_n(X.origin()+i0,iN-i0,*rng);
  }

  // construct X
  std::fill_n(HWs.origin(),HWs.num_elements(),ComplexType(0));  
  std::fill_n(MF.origin(),MF.num_elements(),ComplexType(0));  
  int m0,mN;
  std::tie(m0,mN) = FairDivideBoundary(TG.TG_local().rank(),nCV,TG.TG_local().size());   
  TG.local_barrier();
  for(int m=m0; m<mN; ++m) { 
    auto X_m = X3D[m];
    auto vb_ = vbias[m].origin();
    auto vmf_t = sqrtdt*apply_bound_vbias(vMF[m],1.0);
    auto vmf_ = sqrtdt*vMF[m];
    // apply bounds to vbias
    for(int iw=0; iw<nwalk; iw++) 
      vb_[iw] = apply_bound_vbias(vb_[iw],sqrtdt);
    for(int ni=0; ni<nsteps; ni++) {
      auto X_ = X3D[m][ni].origin();
      auto hws_ = HWs[ni].origin();
      auto mf_ = MF[ni].origin();
      for(int iw=0; iw<nwalk; iw++) {
        // No force bias term when doing free projection.
        ComplexType vdiff = free_projection?ComplexType(0.0, 0.0):(im*(vb_[iw]-vmf_t));
        X_[iw] += vdiff; 
        hws_[iw] -= vdiff * ( X_[iw] - 0.5*vdiff );
        mf_[iw] += im * X_[iw] * vmf_;
      }
    }
  }
  TG.local_barrier();
}

}


}


