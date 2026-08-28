# for findpackage FindcuTENSOR.cmake
find_path(cuTENSOR_INCLUDE_DIR
    NAMES cutensor.h
    PATHS ${CUDAToolkit_INCLUDE_DIRS}
)

find_library(cuTENSOR_LIBRARY
    NAMES cutensor
    PATHS ${CUDAToolkit_LIBRARY_DIR}
)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(cuTENSOR
    REQUIRED_VARS cuTENSOR_LIBRARY cuTENSOR_INCLUDE_DIR
)

if (cuTENSOR_FOUND)
    add_library(cuTENSOR::cuTENSOR UNKNOWN IMPORTED)
    set_target_properties(cuTENSOR::cuTENSOR PROPERTIES
        IMPORTED_LOCATION ${cuTENSOR_LIBRARY}
        INTERFACE_INCLUDE_DIRECTORIES ${cuTENSOR_INCLUDE_DIR}
    )
endif()
