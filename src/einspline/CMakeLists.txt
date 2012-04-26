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
 bspline_base_cuda.h
 )

#FOREACH(a ${HFILES})
#  #INSTALL(FILES ${EINSPLINE_HOME}/src/${a} DESTINATION ${PROJECT_BINARY_DIR}/include/einspline)
#  configure_file(${EINSPLINE_HOME}/src/${a} ${PROJECT_BINARY_DIR}/include/einspline/${a} COPYONLY)
#ENDFOREACH()

set(SRCS )

SET(SRCS ${SRCS}
  bspline_create.c        
  bspline_data.c          
  multi_bspline_create.c  
  multi_nubspline_create.c
  nubspline_create.c      
  nubasis.c               
  nugrid.c                
)

IF(HAVE_SSE)
  SET(SRCS ${SRCS} 
    multi_bspline_eval_sse_s.c   
    multi_bspline_eval_sse_c.c   
    multi_bspline_eval_sse_s_cpp.cc       
    multi_bspline_eval_sse_c_cpp.cc       
    )
else()
  SET(SRCS ${SRCS} 
    multi_bspline_eval_std_s.c 
    multi_bspline_eval_std_c.c 
    multi_bspline_eval_std_s_cpp.cc       
    multi_bspline_eval_std_c_cpp.cc       
    )
endif()


IF(HAVE_SSE2)
  if(EINSPLINE_SSE_BUG)
    set(SRCS ${SRCS}  
      multi_bspline_eval_sse2_d.c    
      multi_bspline_eval_sse2_d_cpp.cc 
      multi_bspline_eval_sse2_z.c   
      multi_bspline_eval_sse2_z_cpp.cc       
      )
  else()
    set(SRCS ${SRCS} 
      multi_bspline_eval_sse_d.c    
      multi_bspline_eval_sse_d_cpp.cc 
      multi_bspline_eval_sse_z.c   
      multi_bspline_eval_sse_z_cpp.cc       
      )
  endif()

  SET(SRCS ${SRCS} 
    multi_nubspline_eval_sse_z.c 
    multi_nubspline_eval_sse_z_cpp.cc   
    )

else()
  SET(SRCS ${SRCS} 
    multi_bspline_eval_std2_d.c 
    multi_bspline_eval_std2_z.c 
    multi_bspline_eval_std2_d_cpp.cc       
    multi_bspline_eval_std2_z_cpp.cc       
    multi_nubspline_eval_std_z.c 
    multi_nubspline_eval_std_z_cpp.cc   
    )
endif()


if(HAVE_CUDA)
  SET(SRCS  ${SRCS}
    multi_bspline_create_cuda.cu  
    bspline_create_cuda.cu 
    )
  CUDA_ADD_LIBRARY(einspline ${SRCS})
else()
  ADD_LIBRARY(einspline ${SRCS})
endif()

#ADD_EXECUTABLE(time_multi time_multi_new.c)
#target_link_libraries(time_multi einspline)
#add_dependencies(time_multi ${PROJECT_BINARY_HOME}/include/einspline/bspline.h)
