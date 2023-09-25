//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License. See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_CREATE_BSPLINE_READERT_H
#define QMCPLUSPLUS_CREATE_BSPLINE_READERT_H

#include <memory>
#include <string>

namespace qmcplusplus
{
/// forward declaration
template <typename T>
class BsplineReaderBaseT;
template <typename T>
class EinsplineSetBuilderT;

/** create a reader which handles complex (double size real) splines, C2R or C2C
 * case spline storage and computation precision is double
 */
template <typename T>
std::unique_ptr<BsplineReaderBaseT<T>>
createBsplineComplexDoubleT(
    EinsplineSetBuilderT<T>* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles complex (double size real) splines, C2R or C2C
 * case spline storage and computation precision is float
 */
template <typename T>
std::unique_ptr<BsplineReaderBaseT<T>>
createBsplineComplexSingleT(
    EinsplineSetBuilderT<T>* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles real splines, R2R case
 *  spline storage and computation precision is double
 */
template <typename T>
std::unique_ptr<BsplineReaderBaseT<T>>
createBsplineRealDoubleT(
    EinsplineSetBuilderT<T>* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles real splines, R2R case
 *  spline storage and computation precision is float
 */
template <typename T>
std::unique_ptr<BsplineReaderBaseT<T>>
createBsplineRealSingleT(
    EinsplineSetBuilderT<T>* e, bool hybrid_rep, const std::string& useGPU);

} // namespace qmcplusplus
#endif
