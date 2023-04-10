set(COMMON_SRCS
    GKIntegration.cc
    CoulombPot.cc
    NLPP.cc
    DFTAtom.cc
    ScreenedPot.cc
    Fitting.cc
    RadialWF.cc
    SplinePot.cc
    CubicSplineCommon.cc
    GeneralPot.cc
    Potential.cc
    IO.cc
    MatrixOps.cc
    Grid.cc
    PotentialBase.cc
    IOASCII.cc
    Functionals.cc)

add_library(common EXCLUDE_FROM_ALL ${COMMON_SRCS})

target_compile_options(common PUBLIC
    $<$<NOT:$<CXX_COMPILER_ID:NVHPC>>:
        -fno-fast-math
    >
    $<$<CXX_COMPILER_ID:NVHPC>:
        -Kieee
    >
>)

target_link_libraries(common PUBLIC qmc_external_multi)
