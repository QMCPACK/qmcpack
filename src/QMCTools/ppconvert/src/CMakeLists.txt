add_executable(ppconvert EXCLUDE_FROM_ALL "")

add_subdirectory(common)

target_sources(ppconvert
  PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParseCommand.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/XMLWriterClass2.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/NLPPClass.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParserClass.cc"
)

target_link_libraries(ppconvert common Math::BLAS_LAPACK)
