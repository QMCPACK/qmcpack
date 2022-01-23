////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//
// File created by:
////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-

namespace qmcplusplus
{
inline int count_electrons(const ParticleSet& ions) { return ions.getTotalNum() * 12; }

inline auto create_prim_lattice()
{
  Lattice::Tensor_t nio_cell = {7.8811, 7.8811, 0.0, -7.8811, 7.8811, 0.0, 0.0, 0.0, 15.7622};
  ParticleSet::ParticleLayout prim_lat;
  // set PBC in x,y,z directions
  prim_lat.BoxBConds = 1;
  // set the lattice
  prim_lat.set(nio_cell);
  return prim_lat;
}

void tile_cell(ParticleSet& ions, const Tensor<int, 3>& tmat)
{
  auto prim_lat(create_prim_lattice());
  // create Ni and O by group
  std::vector<int> nio_group{16, 16};
  ParticleSet prim_ions(prim_lat);
  // create particles by groups
  prim_ions.create(nio_group);
  // using lattice coordinates
  prim_ions.R.InUnit = PosUnit::Lattice;

  // O
  prim_ions.R[0]  = {0.5, 0.0, 0.25};
  prim_ions.R[1]  = {0.0, 0.0, 0.75};
  prim_ions.R[2]  = {0.25, 0.25, 0.5};
  prim_ions.R[3]  = {0.5, 0.5, 0.25};
  prim_ions.R[4]  = {0.75, 0.25, 0.0};
  prim_ions.R[5]  = {0.0, 0.5, 0.75};
  prim_ions.R[6]  = {0.25, 0.75, 0.5};
  prim_ions.R[7]  = {0.75, 0.75, 0.0};
  prim_ions.R[8]  = {0.0, 0.0, 0.25};
  prim_ions.R[9]  = {0.25, 0.25, 0.0};
  prim_ions.R[10] = {0.5, 0.0, 0.75};
  prim_ions.R[11] = {0.75, 0.25, 0.5};
  prim_ions.R[12] = {0.0, 0.5, 0.25};
  prim_ions.R[13] = {0.25, 0.75, 0.0};
  prim_ions.R[14] = {0.5, 0.5, 0.75};
  prim_ions.R[15] = {0.75, 0.75, 0.5};
  // Ni
  prim_ions.R[16] = {0.0, 0.0, 0.0};
  prim_ions.R[17] = {0.0, 0.5, 0.0};
  prim_ions.R[18] = {0.25, 0.25, 0.75};
  prim_ions.R[19] = {0.5, 0.0, 0.5};
  prim_ions.R[20] = {0.5, 0.5, 0.5};
  prim_ions.R[21] = {0.75, 0.25, 0.25};
  prim_ions.R[22] = {0.25, 0.75, 0.75};
  prim_ions.R[23] = {0.75, 0.75, 0.25};
  prim_ions.R[24] = {0.5, 0.0, 0.0};
  prim_ions.R[25] = {0.0, 0.0, 0.5};
  prim_ions.R[26] = {0.25, 0.25, 0.25};
  prim_ions.R[27] = {0.5, 0.5, 0.0};
  prim_ions.R[28] = {0.75, 0.25, 0.75};
  prim_ions.R[29] = {0.0, 0.5, 0.5};
  prim_ions.R[30] = {0.25, 0.75, 0.25};
  prim_ions.R[31] = {0.75, 0.75, 0.75};

  SpeciesSet& species(prim_ions.getSpeciesSet());
  int icharge = species.addAttribute("charge"); // charge_tag);
  species.addSpecies("O");
  species.addSpecies("Ni");

  ions.getSpeciesSet() = prim_ions.getSpeciesSet();
  expandSuperCell(prim_ions, tmat, ions);
}

template<typename J2Type>
void buildJ2(J2Type& J2, double rcut)
{
  using Func     = typename J2Type::FuncType;
  using RealType = typename Func::real_type;
  const int npts = 10;
  std::string optimize("no");
  //RealType rcut=6.4;
  rcut        = std::min(rcut, 5.5727792532);
  RealType dr = rcut / static_cast<RealType>(npts);
  std::vector<RealType> X(npts + 1);
  for (int i = 0; i < npts; ++i)
    X[i] = static_cast<RealType>(i) * dr;

  { //add uu/dd
    std::vector<RealType> Y = {0.28622356,     0.1947736865,  0.1319544873, 0.08893394669,
                               0.05695575776,  0.03565958405, 0.0220695026, 0.01296086466,
                               0.006601006996, 0.00278714433, 0.0};
    std::string suu("uu");
    Func* f = new Func;
    f->initialize(npts, X, Y, -0.25, rcut, suu, optimize);
    J2.addFunc(0, 0, f);
  }
  { //add ud/du
    std::vector<RealType> Y =
        {0.3689309537,  0.2226722029,  0.1484296802,   0.09617039126,  0.0591878654, 0.03660855878,
         0.02262411664, 0.01322279598, 0.006736329049, 0.002871931038, 0.0};

    std::string suu("ud");
    Func* f = new Func;
    f->initialize(npts, X, Y, -0.5, rcut, suu, optimize);
    J2.addFunc(0, 1, f);
  }
}


#if 0
template <typename J2Type> void buildJ2(J2Type &J2, double rcut)
{
  using Func     = typename J2Type::FuncType;
  using RealType = typename Func::real_type;
  const int npts = 10;
  std::string optimize("no");
  rcut        = std::min(rcut, 5.5727792532);

  { // add uu/dd
    std::vector<RealType> Y = {
        0.28622356,     0.1947736865,  0.1319544873, 0.08893394669,
        0.05695575776,  0.03565958405, 0.0220695026, 0.01296086466,
        0.006601006996, 0.00278714433, 0.0};
    std::string suu("uu");
    Func *f = new Func;
    f->setupParameters(npts, rcut, 0.0, Y);
    J2.addFunc(0, 0, f);
  }
  { // add ud/du
    std::vector<RealType> Y = {0.3689309537,
                               0.2226722029,
                               0.1484296802,
                               0.09617039126,
                               0.0591878654,
                               0.03660855878,
                               0.02262411664,
                               0.01322279598,
                               0.006736329049,
                               0.002871931038,
                               0.0};
    std::string suu("ud");
    Func *f = new Func;
    f->setupParameters(npts, rcut, 0.0, Y);
    J2.addFunc(0, 1, f);
  }
}
#endif

template<typename J1Type>
void buildJ1(J1Type& J1, double rcut)
{
  using Func     = typename J1Type::FuncType;
  using RealType = typename Func::real_type;
  const int npts = 10;
  std::string optimize("no");
  rcut = std::min(rcut, 4.8261684030);
  //RealType rcut=6.4;
  RealType dr = rcut / static_cast<RealType>(npts);
  std::vector<RealType> X(npts + 1);
  for (int i = 0; i < npts; ++i)
    X[i] = static_cast<RealType>(i) * dr;

  { // oxygen
    std::vector<RealType> Y = {-0.2249112633,
                               -0.1847494689,
                               -0.115481408,
                               -0.04000122947,
                               0.01731711068,
                               0.05360131926,
                               0.05983040879,
                               0.03955999983,
                               0.0173998007,
                               0.005162164083,
                               0.0};
    std::string oxygen("O");
    Func* f = new Func;
    f->initialize(npts, X, Y, 0, rcut, oxygen, optimize);
    J1.addFunc(0, f);
  }
  { // nio
    std::vector<RealType> Y = {-1.64485534,
                               -1.470658909,
                               -1.078893976,
                               -0.6878964509,
                               -0.3907004509,
                               -0.1962103494,
                               -0.08512755539,
                               -0.02752356864,
                               -0.00401798318,
                               0.0007665934444,
                               0.0};
    std::string nickel("Ni");
    Func* f = new Func;
    f->initialize(npts, X, Y, 0, rcut, nickel, optimize);
    J1.addFunc(1, f);
  }
}
#if 0
template <typename J1Type> void buildJ1(J1Type &J1, double rcut)
{
  using Func     = typename J1Type::FuncType;
  using RealType = typename Func::real_type;
  const int npts = 10;
  std::string optimize("no");
  rcut        = std::min(rcut, 4.8261684030);

  // oxygen
  std::vector<RealType> Y = {-0.2249112633,
                             -0.1847494689,
                             -0.115481408,
                             -0.04000122947,
                             0.01731711068,
                             0.05360131926,
                             0.05983040879,
                             0.03955999983,
                             0.0173998007,
                             0.005162164083,
                             0.0};
  std::string suu("O");
  Func *f = new Func;
  f->setupParameters(npts, rcut, -0.25, Y);
  J1.addFunc(0, f);

  // nickel
  Y = {-1.64485534,
       -1.470658909,
       -1.078893976,
       -0.6878964509,
       -0.3907004509,
       -0.1962103494,
       -0.08512755539,
       -0.02752356864,
       -0.00401798318,
       0.0007665934444,
       0.0};
  suu = "Ni";
  f   = new Func;
  f->setupParameters(npts, rcut, -0.25, Y);
  J1.addFunc(1, f);
}
#endif
template<typename JeeIType>
void buildJeeI(JeeIType& JeeI, double rcut)
{
  // not added yet
}
} // namespace qmcplusplus
