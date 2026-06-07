add_executable(examples
  examples.cpp
) 

target_link_libraries(examples
  PRIVATE fornfdm 
)

if(USE_OPENMP)
  target_link_libraries(examples
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()