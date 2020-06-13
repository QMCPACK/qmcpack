#ifndef AFQMC_TAB_TO_KL_H
#define AFQMC_TAB_TO_KL_H

namespace kernels
{

template<typename T, typename Q>
void Tab_to_Kl(int nwalk, int nocc, int nchol, Q const* Tab, T*  Kl)
{
  for(int w=0; w<nwalk; ++w) {
    for(int a=0; a<nocc; a++) {
      auto Tab_(Tab + ((w*nocc+a)*nocc + a)*nchol);
      auto Kl_(Kl + w*nchol);
      for(int c=0; c<nchol; ++c)
        Kl_[c] += static_cast<T>(Tab_[c]);
    }
  }
}

}


#endif
