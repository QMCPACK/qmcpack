INCLUDE(${CMAKE_ROOT}/Modules/FindPythonLibs.cmake)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

SET (BOOST_PYTHON $ENV{BOOST_HOME}/libs/python/src)
SET (BOOST_PYTHON_SRCS
    ${BOOST_PYTHON}/numeric.cpp
    ${BOOST_PYTHON}/list.cpp
    ${BOOST_PYTHON}/long.cpp
    ${BOOST_PYTHON}/dict.cpp
    ${BOOST_PYTHON}/tuple.cpp
    ${BOOST_PYTHON}/str.cpp
    ${BOOST_PYTHON}/slice.cpp
    ${BOOST_PYTHON}/aix_init_module.cpp
    ${BOOST_PYTHON}/converter/from_python.cpp
    ${BOOST_PYTHON}/converter/registry.cpp
    ${BOOST_PYTHON}/converter/type_id.cpp
    ${BOOST_PYTHON}/object/enum.cpp
    ${BOOST_PYTHON}/object/class.cpp
    ${BOOST_PYTHON}/object/function.cpp
    ${BOOST_PYTHON}/object/inheritance.cpp
    ${BOOST_PYTHON}/object/life_support.cpp
    ${BOOST_PYTHON}/object/pickle_support.cpp
    ${BOOST_PYTHON}/errors.cpp
    ${BOOST_PYTHON}/module.cpp
    ${BOOST_PYTHON}/converter/builtin_converters.cpp
    ${BOOST_PYTHON}/converter/arg_to_python_base.cpp
    ${BOOST_PYTHON}/object/iterator.cpp
    ${BOOST_PYTHON}/object_protocol.cpp
    ${BOOST_PYTHON}/object_operators.cpp
    ${BOOST_PYTHON}/wrapper.cpp
    )

ADD_LIBRARY(boost_python SHARED ${BOOST_PYTHON_SRCS})
