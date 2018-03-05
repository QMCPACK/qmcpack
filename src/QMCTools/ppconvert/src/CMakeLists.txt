
# include(src/common/CMakeLists.txt)

SUBDIRS(common)

add_executable(ppconvert EXCLUDE_FROM_ALL "")

target_sources(ppconvert
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParseCommand.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/XMLWriterClass2.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/NLPPClass.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParserClass.cc"
)

target_link_libraries(ppconvert common ${LAPACK_LIBRARIES})
set_property(TARGET ppconvert PROPERTY INSTALL_RPATH "${LAPACK_LINK_DIRECTORIES}")
set_property(TARGET ppconvert PROPERTY BUILD_WITH_INSTALL_RPATH ON)
set_property(TARGET ppconvert PROPERTY LINK_FLAGS "${LAPACK_LINKER_FLAGS}")
target_include_directories(ppconvert PUBLIC 
  ${MKL_INCLUDE_DIRECTORIES}
  )
