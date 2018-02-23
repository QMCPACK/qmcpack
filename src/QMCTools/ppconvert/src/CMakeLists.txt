
SUBDIRS(common)

add_executable(ppconvert
    "${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParseCommand.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/XMLWriterClass2.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/NLPPClass.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/ParserClass.cc"
)

target_link_libraries(ppconvert common ${LAPACK_LIBRARIES})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${LAPACK_LINKER_FLAGS}")
target_include_directories(ppconvert PUBLIC 
  ${MKL_INCLUDE_DIRECTORIES}
  )
