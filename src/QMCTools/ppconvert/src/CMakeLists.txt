
SUBDIRS(common)

add_executable(ppconvert
    "${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParseCommand.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/XMLWriterClass2.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/NLPPClass.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParserClass.cc"
)

target_link_libraries(ppconvert common mkl_intel_lp64 mkl_sequential mkl_core)

target_include_directories(ppconvert PUBLIC
#  "${CMAKE_CURRENT_SOURCE_DIR}/common"
  )
