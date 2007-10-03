SET (PROGRAM_OPTIONS $ENV{BOOST_HOME}/libs/program_options/src)

SET (BOOST_PROGRAM_OPTIONS_SRCS 
${PROGRAM_OPTIONS}/cmdline.cpp      
${PROGRAM_OPTIONS}/options_description.cpp  
${PROGRAM_OPTIONS}/utf8_codecvt_facet.cpp
${PROGRAM_OPTIONS}/winmain.cpp
${PROGRAM_OPTIONS}/config_file.cpp  
${PROGRAM_OPTIONS}/parsers.cpp              
${PROGRAM_OPTIONS}/value_semantic.cpp
${PROGRAM_OPTIONS}/convert.cpp      
${PROGRAM_OPTIONS}/positional_options.cpp   
${PROGRAM_OPTIONS}/variables_map.cpp
        )

ADD_LIBRARY(boost_program_options ${BOOST_PROGRAM_OPTIONS_SRCS})
