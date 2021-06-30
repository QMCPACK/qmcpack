
#include "Configuration.h"
#include "Utilities/FairDivide.h"
#include "AFQMC/Memory/utilities.hpp"
#include "AFQMC/Utilities/Utils.hpp"
#include "AFQMC/config.h"
#include "AFQMCBasePropagator.h"
#include "AFQMC/Walkers/WalkerConfig.hpp"

// TODO: Remove this
#if defined(ENABLE_CUDA)
#include "AFQMC/Numerics/detail/CUDA/Kernels/construct_X.cuh"
#endif

namespace qmcplusplus
{
namespace afqmc
{
void AFQMCBasePropagator::parse(xmlNodePtr cur)
{
  using qmcplusplus::app_log;

  if (cur == NULL)
    APP_ABORT("Error in AFQMCBasePropagator::parse(): Null xml ptr.\n");

  // defaults
  vbias_bound         = 50.0;
  free_projection     = false;
  hybrid              = true;
  importance_sampling = true;
  apply_constrain     = true;
  // this is wrong!!! must get batched capability from SDet, not from input
  nbatched_propagation = 0;
  nbatched_qr          = 0;
  if (number_of_devices() > 0)
    nbatched_propagation = -1;
  if (number_of_devices() > 0)
  {
    // get better bounds later on
    if (NMO < 1024 && NAEA < 512)
      nbatched_qr = -1;
  }

  xmlNodePtr curRoot = cur;

  std::string sub("yes");
  std::string constrain("yes");
  std::string hyb("yes");
  std::string freep("no");
  std::string impsam("yes");
  std::string external_field("");
  std::string P1ev("no");
  double extfield_scale(1.0);
  ParameterSet m_param;
  m_param.add(vbias_bound, "vbias_bound");
  m_param.add(constrain, "apply_constrain");
  m_param.add(impsam, "importance_sampling");
  m_param.add(hyb, "hybrid");
  m_param.add(external_field, "external_field");
  m_param.add(extfield_scale, "external_field_scale");
  m_param.add(P1ev, "printP1eigval");
  if (TG.TG_local().size() == 1)
    m_param.add(nbatched_propagation, "nbatch");
  if (TG.TG_local().size() == 1)
    m_param.add(nbatched_qr, "nbatch_qr");
  m_param.add(freep, "free_projection");

  //m_param.add(sz_pin_field_file,"sz_pinning_field_file");
  //m_param.add(sz_pin_field_mag,"sz_pinning_field");

  m_param.put(cur);

  std::transform(constrain.begin(), constrain.end(), constrain.begin(), (int (*)(int))tolower);
  if (constrain == "no" || constrain == "false")
    apply_constrain = false;
  std::transform(impsam.begin(), impsam.end(), impsam.begin(), (int (*)(int))tolower);
  if (impsam == "false" || impsam == "no")
    importance_sampling = false;
  std::transform(hyb.begin(), hyb.end(), hyb.begin(), (int (*)(int))tolower);
  if (hyb == "no" || hyb == "false")
    hybrid = false;
  std::transform(freep.begin(), freep.end(), freep.begin(), (int (*)(int))tolower);
  if (freep == "yes" || freep == "true")
    free_projection = true;
  std::transform(P1ev.begin(), P1ev.end(), P1ev.begin(), (int (*)(int))tolower);
  if (P1ev == "yes" || P1ev == "true")
    printP1eV = true;

  app_log() << "\n\n --------------- Parsing Propagator input ------------------ \n\n";

  if (nbatched_propagation != 0)
    app_log() << " Using batched propagation with a batch size: " << nbatched_propagation << "\n";
  else
    app_log() << " Using sequential propagation. \n";
  if (nbatched_qr != 0)
    app_log() << " Using batched orthogonalization in back propagation with a batch size: " << nbatched_qr << "\n";
  else
    app_log() << " Using sequential orthogonalization in back propagation. \n";
  app_log() << " vbias_bound: " << vbias_bound << std::endl;

  if (free_projection)
  {
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
    app_log() << " AFQMC with free projection. \n"
              << " WARNING: Setting: \n"
              << "               -apply_constrain     = no \n"
              << "               -importance_sampling = no \n"
              << "               -hybrid       = yes \n";
    importance_sampling = false;
    hybrid              = true;
    apply_constrain     = false;
    //}
  }
  else
  {
    if (!importance_sampling)
      app_log() << " WARNING: importance_sampling=no without free projection does not make sense. \n";

    if (hybrid)
      app_log() << " Using hybrid method to calculate the weights during the propagation."
                << "\n";
    else
      app_log() << " Using local energy method to calculate the weights during the propagation."
                << "\n";
  }

  // MAM: make this more generic, needs changes for noncollinear
  if (external_field != std::string(""))
  {
    //    read_external_field(H1ext);

    // add placeholder for beta component
    P1.emplace_back(P1Type(tp_ul_ul{0, 0}, tp_ul_ul{0, 0}, 0, aux_alloc_));

    spin_dependent_P1 = true;
    H1ext.reextent({2, NMO, NMO});
    TG.Node().barrier();
    if (TG.Node().root())
    {
      std::ifstream in(external_field.c_str());
      for (int i = 0; i < NMO; i++)
        for (int j = 0; j < NMO; j++)
          in >> H1ext[0][i][j];
      for (int i = 0; i < NMO; i++)
        for (int j = 0; j < NMO; j++)
          in >> H1ext[1][i][j];
      if (in.fail())
        APP_ABORT(" Error: Problems with external field.");
      ma::scal(ComplexType(extfield_scale), H1ext[0]);
      ma::scal(ComplexType(extfield_scale), H1ext[1]);
    }
    TG.Node().barrier();
  }

  app_log() << std::endl;
}

void AFQMCBasePropagator::reset_nextra(int nextra)
{
  if (nextra <= 0)
    return;
  if (last_nextra != nextra)
  {
    last_nextra = nextra;
    for (int n = 0; n < nextra; n++)
    {
      int n0, n1;
      std::tie(n0, n1) = FairDivideBoundary(n, TG.getNCoresPerTG(), nextra);
      if (TG.getLocalTGRank() >= n0 && TG.getLocalTGRank() < n1)
      {
        last_task_index = n;
        break;
      }
    }
    local_group_comm = shared_communicator(TG.TG_local().split(last_task_index, TG.TG_local().rank()));
  }
  if (last_task_index < 0 || last_task_index >= nextra)
    APP_ABORT("Error: Problems in AFQMCBasePropagator::reset_nextra()\n");
}

/*
void AFQMCBasePropagator::assemble_X(size_t nsteps, size_t nwalk, RealType sqrtdt, 
                          StaticMatrix& X, StaticMatrix & vbias, StaticMatrix& MF, 
                          StaticMatrix& HWs, bool addRAND) 
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
  int nCV = int(X.size(0));
  // generate random numbers
  if(addRAND)
  {
    int i0,iN;
    std::tie(i0,iN) = FairDivideBoundary(TG.TG_local().rank(),int(X.num_elements()),
                                         TG.TG_local().size());  
    sampleGaussianFields_n(make_device_ptr(X.origin())+i0,iN-i0,*rng);
  }

  // construct X
  fill_n(make_device_ptr(HWs.origin()),HWs.num_elements(),ComplexType(0));  
  fill_n(make_device_ptr(MF.origin()),MF.num_elements(),ComplexType(0));  

// leaving compiler switch until I decide how to do this better
// basically hide this decision somewhere based on the value of pointer!!!
#ifdef ENABLE_CUDA
  kernels::construct_X(nCV,nsteps,nwalk,free_projection,sqrtdt,vbias_bound,
                       to_address(vMF.origin()),
                       to_address(vbias.origin()),
                       to_address(HWs.origin()),
                       to_address(MF.origin()),
                       to_address(X.origin())
                      );
#else
  boost::multi::array_ref<ComplexType,3> X3D(to_address(X.origin()),
                        {long(X.size(0)),long(nsteps),long(nwalk)});
  int m0,mN;
  std::tie(m0,mN) = FairDivideBoundary(TG.TG_local().rank(),nCV,TG.TG_local().size());   
  TG.local_barrier();
  for(int m=m0; m<mN; ++m) { 
    auto X_m = X3D[m];
    auto vb_ = to_address(vbias[m].origin());
    auto vmf_t = sqrtdt*apply_bound_vbias(vMF[m],1.0);
    auto vmf_ = sqrtdt*vMF[m];
    // apply bounds to vbias
    for(int iw=0; iw<nwalk; iw++) 
      vb_[iw] = apply_bound_vbias(vb_[iw],sqrtdt);
    for(int ni=0; ni<nsteps; ni++) {
      auto X_ = X3D[m][ni].origin();
      auto hws_ = to_address(HWs[ni].origin());
      auto mf_ = to_address(MF[ni].origin());
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
#endif
}
*/

} // namespace afqmc


} // namespace qmcplusplus
