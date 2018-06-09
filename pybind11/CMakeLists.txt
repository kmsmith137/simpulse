link_directories(${CMAKE_SOURCE_DIR}/cpp)

# Reference: http://pybind11.readthedocs.io/en/stable/compiling.html#pybind11-add-module
pybind11_add_module(simpulse_pybind11 
   simpulse_pybind11.cpp 
   wrap_single_pulse.cpp
   wrap_phase_models.cpp
   wrap_von_mises_profile.cpp
)

# FIXME do I want PUBLIC or PRIVATE here?
target_link_libraries(simpulse_pybind11 PUBLIC simpulse)

install(TARGETS simpulse_pybind11
  COMPONENT python
  LIBRARY DESTINATION ${PYTHON_MODULE_INSTALL_DIR}/simpulse)

