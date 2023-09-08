#ifndef QMCPLUSPLUS_PARTICLESETTRAITS_H
#define QMCPLUSPLUS_PARTICLESETTRAITS_H

#include <config.h>

#include "OhmmsData/RecordProperty.h"
#include "OhmmsPETE/Tensor.h"
#include "OhmmsPETE/TinyVector.h"
#include "Particle/Lattice/CrystalLattice.h"
#include "Particle/ParticleBase/ParticleAttrib.h"
#include "type_traits/complex_help.hpp"

namespace qmcplusplus
{
template <typename T>
struct ParticleSetTraits
{
    enum
    {
        DIM = OHMMS_DIM
    };
    using RealType = RealAlias<T>;
    using ComplexType = std::complex<RealType>;
    using ValueType = T;
    using IndexType = int;
    using PosType = TinyVector<RealType, DIM>;
    using GradType = TinyVector<ValueType, DIM>;
    // using HessType = Tensor<ValueType, DIM>;
    using TensorType = Tensor<ValueType, DIM>;
    // using GradHessType = TinyVector<Tensor<ValueType, DIM>, DIM>;
    // using IndexVector = Vector<IndexType>;
    // using ValueVector = Vector<ValueType>;
    // using ValueMatrix = Matrix<ValueType>;
    // using GradVector = Vector<GradType>;
    // using GradMatrix = Matrix<GradType>;
    // using HessVector = Vector<HessType>;
    // using HessMatrix = Matrix<HessType>;
    // using GradHessVector = Vector<GradHessType>;
    // using GradHessMatrix = Matrix<GradHessType>;
    // using VGLVector = VectorSoaContainer<ValueType, DIM + 2>;

    using FullPrecRealType = double;
    using FullPrecComplexType = std::complex<double>;
    using FullPrecValueType = std::conditional_t<IsComplex_t<T>::value,
        FullPrecComplexType, FullPrecRealType>;

    using PropertySetType = RecordNamedProperty<FullPrecRealType>;
};

template <typename T>
struct LatticeParticleTraits
{
    enum
    {
        DIM = OHMMS_DIM
    };
    using RealType = typename ParticleSetTraits<T>::RealType;

    using ParticleLayout = CrystalLattice<RealType, DIM>;
    using SingleParticleIndex = typename ParticleLayout::SingleParticleIndex;
    using SingleParticlePos = typename ParticleLayout::SingleParticlePos;
    using ParticleTensorType = typename ParticleLayout::Tensor_t;

    using FullPrecRealType = typename ParticleSetTraits<T>::FullPrecRealType;
    using FullPrecComplexType =
        typename ParticleSetTraits<T>::FullPrecComplexType;
    using FullPrecValueType = typename ParticleSetTraits<T>::FullPrecValueType;

    using FullPrecGradType = TinyVector<FullPrecValueType, DIM>;

    using Index_t = int;
    using Scalar_t = FullPrecRealType;
    using Complex_t = FullPrecComplexType;
    using Tensor_t  = Tensor<RealType, OHMMS_DIM>;

    using ParticleIndex = ParticleAttrib<Index_t>;
    using ParticleScalar = ParticleAttrib<Scalar_t>;
    using ParticlePos = ParticleAttrib<SingleParticlePos>;
    using ParticleTensor = ParticleAttrib<ParticleTensorType>;

    using ParticleGradient = ParticleAttrib<FullPrecGradType>;
    using ParticleLaplacian = ParticleAttrib<FullPrecValueType>;
    using SingleParticleValue = FullPrecValueType;
};
} // namespace qmcplusplus
#endif
