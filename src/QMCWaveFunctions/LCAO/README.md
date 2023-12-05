LCAOrbitalSet : public SPOSet
    molecular orbitals

SoaLocalizedBasisSet<COT> : public BasisSetBase
    COT(center orbital type)

COT = SoaAtomicBasisSet.h<ROT, SH>
    ROT(radial orbital type)
    builder: AOBasisBuilder

ROT : MultiFunctorAdapter MultiQuinticSpline1D
    builder: RadialOrbitalSetBuilder

SH : SoaSphericalTensor SoaCartesianTensor
