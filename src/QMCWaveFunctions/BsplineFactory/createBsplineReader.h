//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2019 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_CREATE_BSPLINE_READER_H
#define QMCPLUSPLUS_CREATE_BSPLINE_READER_H

#include <string>

namespace qmcplusplus
{
///forward declaration
struct BsplineReaderBase;
class EinsplineSetBuilder;

/** create a reader which handles complex (double size real) splines, C2R or C2C case
 *  spline storage and computation precision is double
 */
BsplineReaderBase* createBsplineComplexDouble(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles complex (double size real) splines, C2R or C2C case
 *  spline storage and computation precision is float
 */
BsplineReaderBase* createBsplineComplexSingle(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles real splines, R2R case
 *  spline storage and computation precision is double
 */
BsplineReaderBase* createBsplineRealDouble(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU);

/** create a reader which handles real splines, R2R case
 *  spline storage and computation precision is float
 */
BsplineReaderBase* createBsplineRealSingle(EinsplineSetBuilder* e, bool hybrid_rep, const std::string& useGPU);

} // namespace qmcplusplus
#endif
