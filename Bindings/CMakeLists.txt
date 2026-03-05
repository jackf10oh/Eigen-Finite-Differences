# allow pybind11 to find python 
set(PYBIND11_FINDPYTHON ON)

# Find pybind11
find_package(pybind11 REQUIRED)

# Link the shared library to your pybind11 module
pybind11_add_module(PyFds
  PyFds.cpp  
)

# Link the Python module with your shared library
target_link_libraries(PyFds 
  PRIVATE FdmStencilExprs
)

# Apply static linking flags as compile options
target_link_libraries(PyFds 
  PRIVATE 
  -static 
  -lpthread 
  -static-libgcc 
  -static-libstdc++
)


