//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: 
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-

#ifndef QMCPLUSPLUS_MINIAPPS_GRAPHITE_HPP
#define QMCPLUSPLUS_MINIAPPS_GRAPHITE_HPP
namespace qmcplusplus
{
  inline int count_electrons(const ParticleSet &ions)
  {
    return ions.getTotalNum() * 4;
  }

  /** expand 4-atom graphite 
   * @param ions particle set
   * @param tmat tiling matrix
   * @param scale scaling factor
   */
  template<typename T>
    Tensor<T,3>  tile_cell(ParticleSet& ions, Tensor<int,3>& tmat, T scale)
    {
      Tensor<T,3> graphite={4.65099, 0.0, 0.0, 
                           -2.3255, 4.02788,0.0, 
                           0.0,0.0, 12.67609393};
      ions.Lattice.BoxBConds=1; 
      ions.Lattice.set(graphite);
      ions.create(4);
      ions.R.InUnit=0;
      ions.R[0]={0.0,0.0,0.0}; 
      ions.R[1]={0.0,2.68525,0.0}; 
      ions.R[2]={0.0,0.0,6.33805}; 
      ions.R[3]={2.3255, 1.34263, 6.33805}; 

      SpeciesSet& species(ions.getSpeciesSet());
      int icharge= species.addAttribute("charge");//charge_tag);
      species.addSpecies("C");

      expandSuperCell(ions,tmat);

      ions.resetGroups();

      return graphite;
    }

  template<typename PT>
    void graphite_4x4(PT& R)
    {
      typedef typename PT::Type_t postype;

      R[ 0]=  postype(  0.000000,    0.000000,    0.000000);
      R[ 1]=  postype( -0.000000,    2.685253,    0.000000);
      R[ 2]=  postype(  0.000000,    0.000000,    6.338047);
      R[ 3]=  postype(  2.325497,    1.342626,    6.338047);
      R[ 4]=  postype( -2.325497,    4.027879,    0.000000);
      R[ 5]=  postype( -2.325497,    6.713132,    0.000000);
      R[ 6]=  postype( -2.325497,    4.027879,    6.338047);
      R[ 7]=  postype(  0.000000,    5.370505,    6.338047);
      R[ 8]=  postype( -4.650994,    8.055758,    0.000000);
      R[ 9]=  postype( -4.650994,   10.741010,    0.000000);
      R[10]=  postype( -4.650994,    8.055758,    6.338047);
      R[11]=  postype( -2.325497,    9.398384,    6.338047);
      R[12]=  postype( -6.976491,   12.083637,    0.000000);
      R[13]=  postype( -6.976491,   14.768889,    0.000000);
      R[14]=  postype( -6.976491,   12.083637,    6.338047);
      R[15]=  postype( -4.650994,   13.426263,    6.338047);
      R[16]=  postype(  4.650994,    0.000000,    0.000000);
      R[17]=  postype(  4.650994,    2.685253,    0.000000);
      R[18]=  postype(  4.650994,    0.000000,    6.338047);
      R[19]=  postype(  6.976491,    1.342626,    6.338047);
      R[20]=  postype(  2.325497,    4.027879,    0.000000);
      R[21]=  postype(  2.325497,    6.713132,    0.000000);
      R[22]=  postype(  2.325497,    4.027879,    6.338047);
      R[23]=  postype(  4.650994,    5.370505,    6.338047);
      R[24]=  postype(  0.000000,    8.055758,    0.000000);
      R[25]=  postype( -0.000000,   10.741010,    0.000000);
      R[26]=  postype(  0.000000,    8.055758,    6.338047);
      R[27]=  postype(  2.325497,    9.398384,    6.338047);
      R[28]=  postype( -2.325497,   12.083637,    0.000000);
      R[29]=  postype( -2.325497,   14.768889,    0.000000);
      R[30]=  postype( -2.325497,   12.083637,    6.338047);
      R[31]=  postype(  0.000000,   13.426263,    6.338047);
      R[32]=  postype(  9.301988,    0.000000,    0.000000);
      R[33]=  postype(  9.301988,    2.685253,    0.000000);
      R[34]=  postype(  9.301988,    0.000000,    6.338047);
      R[35]=  postype( 11.627485,    1.342626,    6.338047);
      R[36]=  postype(  6.976491,    4.027879,    0.000000);
      R[37]=  postype(  6.976491,    6.713132,    0.000000);
      R[38]=  postype(  6.976491,    4.027879,    6.338047);
      R[39]=  postype(  9.301988,    5.370505,    6.338047);
      R[40]=  postype(  4.650994,    8.055758,    0.000000);
      R[41]=  postype(  4.650994,   10.741010,    0.000000);
      R[42]=  postype(  4.650994,    8.055758,    6.338047);
      R[43]=  postype(  6.976491,    9.398384,    6.338047);
      R[44]=  postype(  2.325497,   12.083637,    0.000000);
      R[45]=  postype(  2.325497,   14.768889,    0.000000);
      R[46]=  postype(  2.325497,   12.083637,    6.338047);
      R[47]=  postype(  4.650994,   13.426263,    6.338047);
      R[48]=  postype( 13.952982,    0.000000,    0.000000);
      R[49]=  postype( 13.952982,    2.685253,    0.000000);
      R[50]=  postype( 13.952982,    0.000000,    6.338047);
      R[51]=  postype( 16.278479,    1.342626,    6.338047);
      R[52]=  postype( 11.627485,    4.027879,    0.000000);
      R[53]=  postype( 11.627485,    6.713132,    0.000000);
      R[54]=  postype( 11.627485,    4.027879,    6.338047);
      R[55]=  postype( 13.952982,    5.370505,    6.338047);
      R[56]=  postype(  9.301988,    8.055758,    0.000000);
      R[57]=  postype(  9.301988,   10.741010,    0.000000);
      R[58]=  postype(  9.301988,    8.055758,    6.338047);
      R[59]=  postype( 11.627485,    9.398384,    6.338047);
      R[60]=  postype(  6.976491,   12.083637,    0.000000);
      R[61]=  postype(  6.976491,   14.768889,    0.000000);
      R[62]=  postype(  6.976491,   12.083637,    6.338047);
      R[63]=  postype(  9.301988,   13.426263,    6.338047);
    }


  template<typename JeeIType>
    void buildJeeI(JeeIType&  JeeI, double rcut)
    {
      using Func=typename JeeIType::FuncType;
      using RealType=typename Func::real_type;
      {
        std::vector<RealType> params=
        {8.227710241e-06,2.480817653e-06,-5.354068112e-06,-1.112644787e-05,
         -2.208006078e-06,5.213121933e-06,-1.537865869e-05,8.899030233e-06,
         6.257255156e-06,3.214580988e-06,-7.716743107e-06,-5.275682077e-06,
         -1.778457637e-06,7.926231121e-06,1.767406868e-06,5.451359059e-08,
         2.801423724e-06,4.577282736e-06,7.634608083e-06,-9.510673173e-07,
         -2.344131575e-06,-1.878777219e-06,3.937363358e-07,5.065353773e-07,
         5.086724869e-07,-1.358768154e-07 };
        Func* functor=new Func;
        functor->cutoff_radius = rcut;
        functor->resize(3,3);
        functor->Parameters=params;
        functor->reset_gamma();

        JeeI.addFunc(0,0,0,functor);
      }
      {
        std::vector<RealType> params=
        {-6.939530224e-06,2.634169299e-05,4.046077477e-05,-8.002682388e-06,
         -5.396795988e-06,6.697370507e-06,5.433953051e-05,-6.336849668e-06,
         3.680471431e-05,-2.996059772e-05,1.99365828e-06,-3.222705626e-05,
         -8.091669063e-06,4.15738535e-06,4.843939112e-06,3.563650208e-07,
         3.786332474e-05,-1.418336941e-05,2.282691374e-05,1.29239286e-06,
         -4.93580873e-06,-3.052539228e-06,9.870288001e-08,1.844286407e-06,
         2.970561871e-07,-4.364303677e-08};
        Func* functor=new Func;
        functor->cutoff_radius = rcut;
        functor->resize(3,3);
        functor->Parameters=params;
        functor->reset_gamma();

        JeeI.addFunc(0,0,1,functor);
      }
      JeeI.check_complete();
    }

  template<typename J2Type>
    void buildJ2(J2Type&  J2, double rcut)
    {
      using Func=typename J2Type::FuncType;
      using RealType=typename Func::real_type;
      const int npts=10;
      std::string optimize("no");
      //RealType rcut=6.4;
      RealType dr=rcut/static_cast<RealType>(npts);
      std::vector<RealType> X(npts+1);
      for(int i=0; i<npts; ++i) X[i]=static_cast<RealType>(i)*dr;

      { //add uu/dd
        std::vector<RealType> Y=
        {0.4711f, 0.3478f, 0.2445f, 0.1677f, 0.1118f,
         0.0733f, 0.0462f, 0.0273f, 0.0145f, 0.0063f, 0.0f};
        // 0.0733f, 0.0462f, 0.0273f, 0.0145f, 0.0063f, 0.0f};
        std::string suu("uu");
        Func* f=new Func;
        f->initialize(npts,X,Y,-0.25,rcut,suu,optimize);
        J2.addFunc(0,0,f);
      }
      { //add ud/du
        std::vector<RealType> Y=
        {0.6715f, 0.4433f, 0.2901f, 0.1889f, 0.1227f,
         0.0793f, 0.0496f, 0.0292f, 0.0152f, 0.0061f, 0.0f};

        std::string suu("ud");
        Func* f=new Func;
        f->initialize(npts,X,Y,-0.5,rcut,suu,optimize);
        J2.addFunc(0,1,f);
      }
    }

  template<typename J1Type>
    void buildJ1(J1Type&  J1, double rcut)
    {
      using Func=typename J1Type::FuncType;
      using RealType=typename Func::real_type;
      const int npts=10;
      std::string optimize("no");
      //RealType rcut=6.4;
      RealType dr=rcut/static_cast<RealType>(npts);
      std::vector<RealType> X(npts+1);
      for(int i=0; i<npts; ++i) X[i]=static_cast<RealType>(i)*dr;

      std::vector<RealType> Y=
      {0.4711f, 0.3478f, 0.2445f, 0.1677f, 0.1118f,
        0.0733f, 0.0462f, 0.0273f, 0.0145f, 0.0063f, 0.0f};
      std::string suu("uu");
      Func* f=new Func;
      f->initialize(npts,X,Y,-0.25,rcut,suu,optimize);
      J1.addFunc(0,f);
    }
}
#endif
