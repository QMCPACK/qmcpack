#ifndef OHMMS_QMC__SPHERICAL_CARTESIAN_TENSOR_H
#define OHMMS_QMC__SPHERICAL_CARTESIAN_TENSOR_H

/** SphericalTensor
 * @author Jeongnim Kim
 * @brief evaluates the Real Spherical Harmonics for Gaussian Orbitals
 \f[
 r^l \Re (Y_l^m(\theta,\phi)). 
 \f]
 *
 The template parameter T is the value_type, e.g. double, and the 
 template parameter Point_t is a vector type which must have the
 operator[] defined.
 */
template<class T, class Point_t>
class SphericalTensor {
public : 

  typedef T value_type;
  typedef Point_t pos_type;
  typedef SphericalTensor<T,Point_t> This_t;

  ///constructor
  explicit SphericalTensor(const int lmax);

  ///makes a table of \f$ r^{l} \Re (Y_l^m) \f$ and their gradients up to lmax.
  void evaluate(const Point_t& p);     

  ///returns the index \f$ l(l+1)+m \f$
  inline int index(int l, int m) const {return (l*(l+1))+m;}

  ///returns the value of \f$ r^{l} \Re (Y_l^m) \f$ given l,m
  inline value_type getYlm(int l, int m) const
   {return Ylm[index(l,m)];}

  ///returns the gradient of \f$ r^{l} \Re (Y_l^m) \f$ given l,m
  inline Point_t getGradYlm(int l, int m) const
  {return gradYlm[index(l,m)];}

  ///returns the value of \f$ r^{l} \Re (Y_l^m) \f$ given index lm
  inline value_type getYlm(int lm) const {return Ylm[lm];}

  ///returns the gradient of \f$ r^{l} \Re (Y_l^m) \f$ given index lm
  inline Point_t getGradYlm(int lm) const {return gradYlm[lm];}

  inline int size() const { return Ylm.size();}

  inline int lmax() const { return Lmax;}

  int Lmax;
  std::vector<value_type> Ylm;
  std::vector<Point_t> gradYlm;

};

template<class SCT, unsigned L> 
struct SCTFunctor { 
  typedef typename SCT::value_type value_type;
  typedef typename SCT::pos_type pos_type;
  static inline void apply(std::vector<value_type>& Ylm, std::vector<pos_type>& gYlm, const pos_type& p) {
    SCTFunctor<SCT,L-1>::apply(Ylm,gYlm,p);
  }
};

template<class SCT> 
struct SCTFunctor<SCT,1> { 
  typedef typename SCT::value_type value_type;
  typedef typename SCT::pos_type pos_type;
  static inline void apply(std::vector<value_type>& Ylm, std::vector<pos_type>& gYlm, const pos_type& p) {
    const value_type L1 = sqrt(3.0);

    Ylm[1]=L1*p[1];
    Ylm[2]=L1*p[2];
    Ylm[3]=L1*p[0];

    gYlm[1]=pos_type(0.0,L1,0.0);
    gYlm[2]=pos_type(0.0,0.0,L1);
    gYlm[3]=pos_type(L1,0.0,0.0);
  }
};

template<class SCT> 
struct SCTFunctor<SCT,2> { 
  typedef typename SCT::value_type value_type;
  typedef typename SCT::pos_type pos_type;
  static inline void apply(std::vector<value_type>& Ylm, std::vector<pos_type>& gYlm, const pos_type& p) {
    SCTFunctor<SCT,1>::apply(Ylm,gYlm,p);
    value_type x=p[0], y=p[1], z=p[2];
    value_type x2=x*x, y2=y*y, z2=z*z;
    value_type xy=x*y, xz=x*z, yz=y*z; 

    const value_type L0 = sqrt(1.25);
    const value_type L1 = sqrt(15.0);
    const value_type L2 = sqrt(3.75);

    Ylm[4]=L1*xy;
    Ylm[5]=L1*yz;
    Ylm[6]=L0*(2.0*z2-x2-y2);
    Ylm[7]=L1*xz;
    Ylm[8]=L2*(x2-y2);

    gYlm[4]=pos_type(L1*y,L1*x,0.0);
    gYlm[5]=pos_type(0.0,L1*z,L1*y);
    gYlm[6]=pos_type(-2.0*L0*x,-2.0*L0*y,4.0*L0*z);
    gYlm[7]=pos_type(L1*z,0.0,L1*x);
    gYlm[8]=pos_type(2.0*L2*x,-2.0*L2*y,0.0);
  }
};

template<class SCT> 
struct SCTFunctor<SCT,3> { 
  typedef typename SCT::value_type value_type;
  typedef typename SCT::pos_type pos_type;
  static inline void apply(std::vector<value_type>& Ylm, std::vector<pos_type>& gYlm, const pos_type& p) {
    SCTFunctor<SCT,2>::apply(Ylm,gYlm,p);
  }
};


template<class T, class Point_t>
SphericalTensor<T, Point_t>::SphericalTensor(const int lmax) : Lmax(lmax){ 
  Ylm.resize((lmax+1)*(lmax+1));
  gradYlm.resize((lmax+1)*(lmax+1));
}

template<class T, class Point_t>
void SphericalTensor<T,Point_t>::evaluate(const Point_t& p) {

  const value_type norm = 1.0/sqrt(16.0*atan(1.0));

  Ylm[0]=1.0;
  gradYlm[0]=0.0;

  switch (Lmax)
  {
    case(0): break;
    case(1): SCTFunctor<This_t,1>::apply(Ylm,gradYlm,p); break;
    case(2): SCTFunctor<This_t,2>::apply(Ylm,gradYlm,p); break;
    //case(3): SCTFunctor<This_t,3>::apply(Ylm,gradYlm,p); break;
    //case(4): SCTFunctor<This_t,4>::apply(Ylm,gradYlm,p); break;
    //case(5): SCTFunctor<This_t,5>::apply(Ylm,gradYlm,p); break;
    defaults: 
             std::cerr << "Lmax>2 is not valid." << std::endl;
             break;
  }

  for(int i=0; i<Ylm.size(); i++) Ylm[i] *= norm;
  for(int i=0; i<Ylm.size(); i++) gradYlm[i] *= norm;
}

#endif
