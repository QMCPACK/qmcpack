#include <cmath>

namespace qmcplusplus 
{
  template<typename T>
    inline T diff(const Tensor<T,3>& a, const Tensor<T,3>& b)
    {
      T res=0.0;
      for(int i=0; i<9; ++i) res+=(a(i)-b(i))*(a(i)-b(i));
      return res/9.0;
    }

  template<typename T>
    struct TestFunc {

      T k0,k1,k2;
      T d2factor;
      Tensor<T,3> kk;

      TestFunc(int nk0=1, int nk1=1, int nk2=1) {
        const T twopi = 8.0*std::atan(1.0);
        k0=twopi*static_cast<double>(nk0);
        k1=twopi*static_cast<double>(nk1);
        k2=twopi*static_cast<double>(nk2);
        d2factor = -(k0*k0+k1*k1+k2*k2);
        kk(0,0)=-k0*k0; kk(0,1)=-k0*k1; kk(0,2)=-k0*k2;
        kk(1,0)=-k1*k0; kk(1,1)=-k1*k1; kk(1,2)=-k1*k2;
        kk(2,0)=-k2*k0; kk(2,1)=-k2*k1; kk(2,2)=-k2*k2;
      }

      //inline double operator()(const TinyVector<double,3>& pos) {
      //  return sin(k0*pos[0])*sin(k1*pos[1])*sin(k2*pos[2]);
      //}
      //inline double operator()(double x, double y, double z) {
      //  return sin(k0*x)*sin(k1*y)*sin(k2*z);
      //}
      //
      template<class PV>
        inline T f(const PV& pos) {
          return std::sin(k0*pos[0])*std::sin(k1*pos[1])*std::sin(k2*pos[2]);
        }

      template<class PV>
        inline TinyVector<T,3> df(const PV& pos) {
          return TinyVector<T,3>(k0*std::cos(k0*pos[0])*std::sin(k1*pos[1])*std::sin(k2*pos[2]),
              k1*std::sin(k0*pos[0])*std::cos(k1*pos[1])*std::sin(k2*pos[2]),
              k2*std::sin(k0*pos[0])*std::sin(k1*pos[1])*std::cos(k2*pos[2]));
        }

      template<class PV>
        inline T d2f(const PV& pos) {
          return d2factor*f(pos);
        }

      template<class PV>
        inline void d2f(const PV& pos, Tensor<T,3>& hess) 
        {
          hess=f(pos)*kk;
        }

      //inline T d2f(double x, double y, double z) {
      //  return d2factor*f(x,y,z);
      //}

    };
  template<typename T>
    struct ComboFunc {

      std::vector<T> C;
      std::vector<TestFunc<T>*> F;

      ComboFunc() {}
      ~ComboFunc() 
      {
        for(int i=0; i<F.size(); i++) delete F[i];
      }

      void push_back(T c, TestFunc<T>* fn) { 
        C.push_back(c); 
        F.push_back(fn);
      }

      template<typename PV>
        inline T f(const PV& pos)
        {
          T res=0;
          for(int i=0; i<C.size(); i++) res += C[i]*F[i]->f(pos);
          return res;
        }

      template<typename PV>
        inline TinyVector<T,3> df(const PV& pos) {
          TinyVector<T,3> res(0.0,0.0,0.0);
          for(int i=0; i<C.size(); i++) res += C[i]*F[i]->df(pos);
          return res;
        }

      template<class PV>
        inline T d2f(const PV& pos) {
          T res=0;
          for(int i=0; i<C.size(); i++) res += C[i]*F[i]->d2f(pos);
          return res;
        }

      template<class PV>
        inline void d2f(const PV& pos, Tensor<T,3>& hess) 
        {
          hess=0.0;
          Tensor<T,3> h;
          for(int i=0; i<C.size(); i++) 
          {
            F[i]->d2f(pos,h);
            hess+=C[i]*h;
          }
        }

    };
}

