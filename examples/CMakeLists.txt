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

add_executable(basic_api
  basic_api.cpp
) 

target_link_libraries(basic_api
  PRIVATE fornfdm 
)

if(USE_OPENMP)
  target_link_libraries(basic_api
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()

add_executable(periodic
  periodic.cpp
) 

target_link_libraries(periodic
  PRIVATE fornfdm 
)

if(USE_OPENMP)
  target_link_libraries(periodic
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()

add_executable(diffusion
  diffusion.cpp
) 

target_link_libraries(diffusion
  PRIVATE fornfdm 
)

if(USE_OPENMP)
  target_link_libraries(diffusion
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()

add_executable(internals
  internals.cpp
) 

target_link_libraries(internals
  PRIVATE fornfdm 
)

if(USE_OPENMP)
  target_link_libraries(internals
    PRIVATE OpenMP::OpenMP_CXX
  )
endif()