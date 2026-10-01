if(BUILD_UNIT_TESTS)
	set(UTEST_EXE test_ppconvert_parser)
	set(UTEST_NAME deterministic-unit_${UTEST_EXE})

	add_executable(${UTEST_EXE} test_parser.cpp ../src/ParserClass.cc ../src/ParseCommand.cc)
	target_include_directories(${UTEST_EXE} PRIVATE ../src)
	target_link_libraries(${UTEST_EXE} catch_main)

	add_unit_test(${UTEST_NAME} 1 1 $<TARGET_FILE:${UTEST_EXE}>)
	set_tests_properties(${UTEST_NAME} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
endif()

add_subdirectory(Oxygen)
