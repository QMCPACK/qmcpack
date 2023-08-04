//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National
// Laboratory
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of
//                    Illinois at Urbana-Champaign Mark A. Berrill,
//                    berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois
// at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_COMPOSITE_SPOSETT_H
#define QMCPLUSPLUS_COMPOSITE_SPOSETT_H

#include "QMCWaveFunctions/BasisSetBase.h"
#include "QMCWaveFunctions/SPOSetBuilder.h"
#include "QMCWaveFunctions/SPOSetBuilderFactory.h"
#include "QMCWaveFunctions/SPOSetT.h"

namespace qmcplusplus
{
template <typename T>
class CompositeSPOSetT : public SPOSetT<T>
{
public:
	using ValueVector = typename SPOSetT<T>::ValueVector;
	using GradVector = typename SPOSetT<T>::GradVector;
	using ValueMatrix = typename SPOSetT<T>::ValueMatrix;
	using GradMatrix = typename SPOSetT<T>::GradMatrix;
	using HessMatrix = typename SPOSetT<T>::HessMatrix;
	using GGGMatrix = typename SPOSetT<T>::GGGMatrix;

	/// component SPOSets
	std::vector<std::unique_ptr<SPOSetT<T>>> components;
	/// temporary storage for values
	std::vector<ValueVector> component_values;
	/// temporary storage for gradients
	std::vector<GradVector> component_gradients;
	/// temporary storage for laplacians
	std::vector<ValueVector> component_laplacians;
	/// store the precomputed offsets
	std::vector<int> component_offsets;

	CompositeSPOSetT(const std::string& my_name);
	/**
	 * @TODO: do we want template copy constructor
	 * (i.e., copy from other with different type argument)?
	 */
	CompositeSPOSetT(const CompositeSPOSetT& other);
	~CompositeSPOSetT() override;

	std::string
	getClassName() const override
	{
		return "CompositeSPOSetT";
	}

	/// add a sposet component to this composite sposet
	void
	add(std::unique_ptr<SPOSetT<T>> component);

	/// print out component info
	void
	report();

	// SPOSet interface methods
	/// size is determined by component sposets and nothing else
	inline void
	setOrbitalSetSize(int norbs) override
	{
	}

	std::unique_ptr<SPOSetT<T>>
	makeClone() const override;

	void
	evaluateValue(const ParticleSet& P, int iat, ValueVector& psi) override;

	void
	evaluateVGL(const ParticleSet& P, int iat, ValueVector& psi,
		GradVector& dpsi, ValueVector& d2psi) override;

	/// unimplemented functions call this to abort
	inline void
	not_implemented(const std::string& method)
	{
		APP_ABORT("CompositeSPOSetT::" + method + " has not been implemented");
	}

	// methods to be implemented in the future (possibly)
	void
	evaluate_notranspose(const ParticleSet& P, int first, int last,
		ValueMatrix& logdet, GradMatrix& dlogdet,
		ValueMatrix& d2logdet) override;
	void
	evaluate_notranspose(const ParticleSet& P, int first, int last,
		ValueMatrix& logdet, GradMatrix& dlogdet,
		HessMatrix& ddlogdet) override;
	void
	evaluate_notranspose(const ParticleSet& P, int first, int last,
		ValueMatrix& logdet, GradMatrix& dlogdet, HessMatrix& ddlogdet,
		GGGMatrix& dddlogdet) override;
};

} // namespace qmcplusplus

#endif
