#//////////////////////////////////////////////////////////////////////////////////////
#// This file is distributed under the University of Illinois/NCSA Open Source License.
#// See LICENSE file in top directory for details.
#//
#// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
#//
#// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
#//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
#//
#// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
#//////////////////////////////////////////////////////////////////////////////////////
    
    


PROJECT(benchmark)

SET(BENCH numerics)

#ADD_EXECUTABLE( fft1d  fft1d.cpp)
#TARGET_LINK_LIBRARIES(fft1d qmcutil)
#
#ADD_EXECUTABLE( fft1d_bench  fft1d_bench.cpp)
#TARGET_LINK_LIBRARIES(fft1d_bench qmcutil)

FOREACH(p ${BENCH})
  ADD_EXECUTABLE( ${p}  ${p}.cpp)
  TARGET_LINK_LIBRARIES(${p} qmcutil)
  IF(MPI_LIBRARY)
    TARGET_LINK_LIBRARIES(${p} ${MPI_LIBRARY})
  ENDIF(MPI_LIBRARY)
ENDFOREACH(p ${BENCH})
