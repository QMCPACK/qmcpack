add_executable(ppconvert) # Build ppconvert with standard make invocation

add_subdirectory(common)

target_sources(
  ppconvert
  PUBLIC "${CMAKE_CURRENT_SOURCE_DIR}/CubicSpline.cc" "${CMAKE_CURRENT_SOURCE_DIR}/ParseCommand.cc"
         "${CMAKE_CURRENT_SOURCE_DIR}/XMLWriterClass2.cc" "${CMAKE_CURRENT_SOURCE_DIR}/NLPPClass.cc"
         "${CMAKE_CURRENT_SOURCE_DIR}/ParserClass.cc")

target_link_libraries(ppconvert common platform_LA)

add_test_target_in_output_location(ppconvert bin)
install(TARGETS ppconvert RUNTIME DESTINATION bin)
