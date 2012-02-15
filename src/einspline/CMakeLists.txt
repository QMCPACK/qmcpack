SET(SRCS
  ${EINSPLINE_DIR}/src/bspline_create.c        
  ${EINSPLINE_DIR}/src/bspline_data.c          
  ${EINSPLINE_DIR}/src/multi_bspline_create.c  
  ${EINSPLINE_DIR}/src/multi_nubspline_create.c
  ${EINSPLINE_DIR}/src/nubspline_create.c      
  ${EINSPLINE_DIR}/src/nubasis.c               
  ${EINSPLINE_DIR}/src/nugrid.c                
)

IF(HAVE_SSE)
  SET(SRCS ${SRCS} 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_s.c   
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_s_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_c.c   
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_c_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_d.c   
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_d_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_z.c   
    ${EINSPLINE_DIR}/src/multi_bspline_eval_sse_z_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_nubspline_eval_sse_z.c 
    ${EINSPLINE_DIR}/src/multi_nubspline_eval_sse_z_cpp.cc   
    )
else()
  SET(SRCS ${SRCS} 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_s.c 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_s_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_c.c 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_c_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_d.c 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_d_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_z.c 
    ${EINSPLINE_DIR}/src/multi_bspline_eval_std_z_cpp.cc       
    ${EINSPLINE_DIR}/src/multi_nubspline_eval_std_z.c 
    ${EINSPLINE_DIR}/src/multi_nubspline_eval_std_z_cpp.cc   
    )
endif()

if(HAVE_CUDA)
  SET(SRCS  ${SRCS}
    ${EINSPLINE_DIR}/src/multi_bspline_create_cuda.cu  
    ${EINSPLINE_DIR}/src/bspline_create_cuda.cu 
    )
  CUDA_ADD_LIBRARY(einspline ${SRCS})
else()
  ADD_LIBRARY(einspline ${SRCS})
endif()

ADD_EXECUTABLE(time_multi ${EINSPLINE_DIR}/src/time_multi_new.c)
target_link_libraries(time_multi einspline)

SET(HFILES
 bspline_base.h bspline_create.h bspline_eval_sse_c.h        
 bspline_eval_sse_d.h bspline_eval_sse_s.h                   
 bspline_eval_sse_z.h bspline_eval_std_c.h                   
 bspline_eval_std_d.h bspline_eval_std_s.h                   
 bspline_eval_std_z.h bspline.h bspline_structs.h            
 fbspline.h fmulti_bspline.h  fnubspline.h                   
 multi_bspline.h                                             
 multi_bspline_create.h    multi_bspline_structs.h           
 multi_bspline_eval_c.h    multi_bspline_eval_d.h            
 multi_bspline_eval_s.h    multi_bspline_eval_z.h            
 multi_nubspline.h                                           
 multi_nubspline_create.h    multi_nubspline_structs.h       
 multi_nubspline_eval_c.h    multi_nubspline_eval_d.h        
 multi_nubspline_eval_s.h    multi_nubspline_eval_z.h        
 nubspline_base.h nubspline_create.h                         
 nubspline_eval_sse_s.h nubspline_eval_std_s.h               
 nubspline_eval_sse_c.h nubspline_eval_std_c.h               
 nubspline_eval_sse_d.h nubspline_eval_std_d.h               
 nubspline_eval_sse_z.h nubspline_eval_std_z.h               
 nubspline_structs.h nubasis.h nubspline.h nugrid.h          
 multi_bspline_structs_cuda.h multi_bspline_create_cuda.h    
 multi_bspline_eval_cuda.h bspline_eval_cuda.h               
 bspline_structs_cuda.h bspline_create_cuda.h                
 bspline_base_cuda.h)

FOREACH(a ${HFILES})
  INSTALL(FILES ${EINSPLINE_DIR}/src/${a} DESTINATION ${PROJECT_BINARY_DIR}/include/einspline)
ENDFOREACH()

