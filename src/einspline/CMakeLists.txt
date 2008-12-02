#SET (EINSPLINE_DIR $ENV{EINSPLINE_HOME})

SET(ESLIBSRC
   ${EINSPLINE_DIR}/src/bspline_create.c      
   ${EINSPLINE_DIR}/src/multi_bspline_create.c  
   ${EINSPLINE_DIR}/src/bspline_data.c        
   ${EINSPLINE_DIR}/src/nubspline_create.c    
   ${EINSPLINE_DIR}/src/nubasis.c             
   ${EINSPLINE_DIR}/src/nugrid.c              
   ${EINSPLINE_DIR}/src/multi_nubspline_create.c
   )

 IF(HAVE_SSE)
   SET(ESLIBSRC ${ESLIBSRC}
     ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_s.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_s_cpp.cc       
     ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_c.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_c_cpp.cc       
     )
 ELSE(HAVE_SSE)
   SET(ESLIBSRC ${ESLIBSRC}
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_s.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_s_cpp.cc       
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_c.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_c_cpp.cc       
     )
 ENDIF(HAVE_SSE)

 IF(HAVE_SSE2)
   SET(ESLIBSRC ${ESLIBSRC}
   ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_d.c 
   ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_d_cpp.cc 
   ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_z.c 
   ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_z_cpp.cc 
   ${EINSPLINE_DIR}/src/multi_nubspline_eval_sse_z.c 
   ${EINSPLINE_DIR}/src/multi_nubspline_eval_sse_z_cpp.cc
   )
 ELSE(HAVE_SSE2)
   SET(ESLIBSRC ${ESLIBSRC}
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_d.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_d_cpp.cc 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_z.c 
     ${EINSPLINE_DIR}/src/multi_bspline_eval_std_z_cpp.cc 
     ${EINSPLINE_DIR}/src/multi_nubspline_eval_std_z.c 
     ${EINSPLINE_DIR}/src/multi_nubspline_eval_std_z_cpp.cc  
   )
 ENDIF(HAVE_SSE2)

 SET(ESLIBH
 nugrid.h
 nubspline_structs.h
 nubspline_eval_std_s.h
 nubspline_eval_std_d.h
 nubspline_eval_std_c.h
 nubspline_eval_sse_z.h
 nubspline_eval_sse_s.h
 nubspline_eval_sse_d.h
 nubspline_create.h
 nubspline_base.h
 nubspline.h
 nubasis.h
 multi_bspline_structs.h
 multi_bspline_eval_std_s.h
 multi_bspline_eval_std_c.h
 multi_bspline_eval_sse_z_unrolled.h
 multi_bspline_eval_sse_s.h
 multi_bspline_eval_sse_d.h
 multi_bspline_eval_sse_c.h
 multi_nubspline.h
 multi_bspline_create.h
 multi_nubspline_create.h
 multi_nubspline_structs.h
 bspline_structs.h
 bspline_eval_std_z.h
 bspline_eval_std_s.h
 bspline_eval_std_d.h
 bspline_eval_std_c_ansi.h
 bspline_eval_sse_z.h
 bspline_eval_sse_s.h
 bspline_eval_sse_d.h
 bspline_eval_altivec_s.h
 bspline_create.h
 bspline_base.h
 bspline.h
 blip_create.h
 aligned_alloc.h
 nubspline_eval_std_z.h
 nubspline_eval_sse_c.h
 bspline_eval_std_c.h
 multi_bspline_eval_std_d.h
 multi_bspline_eval_std_s_impl.h
 multi_bspline_eval_sse_z_impl.h
 multi_bspline_eval_z.h
 multi_bspline_eval_std_z_impl.h
 multi_bspline_eval_sse_d_impl.h
 multi_bspline_eval_c.h
 multi_bspline_eval_s.h
 multi_bspline_eval_sse_c_impl.h
 multi_bspline_eval_std_d_impl.h
 bspline_eval_sse_c.h
 multi_bspline_eval_sse_s_impl.h
 multi_bspline.h
 multi_bspline_eval_std_c_impl.h
 multi_bspline_eval_sse_z.h
 multi_bspline_eval_d.h
 multi_nubspline_eval_s.h
 multi_nubspline_eval_c.h
 multi_nubspline_eval_d.h
 multi_nubspline_eval_z.h
 multi_nubspline_eval_sse_z_impl.h
   )

FOREACH(eh ${ESLIBH})
  #CONFIGURE_FILE(${EINSPLINE_DIR}/src/${eh} ${CMAKE_CURRENT_SOURCE_DIR}/${eh})
  CONFIGURE_FILE(${EINSPLINE_DIR}/src/${eh} ${CMAKE_CURRENT_BINARY_DIR}/${eh})
ENDFOREACH(eh)

ADD_LIBRARY(einspline ${ESLIBSRC})
