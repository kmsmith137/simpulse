cmake_minimum_required(VERSION 2.8.12)
project(simpulse)
set(CMAKE_MACOSX_RPATH 0)   # needed to avoid warning from 'cmake'

message("")
message("Reminder: here are some influential variables, and their values (to change, do 'cmake -DVAR=VALUE ...' or 'export VAR=VALUE')")
message("   CMAKE_INSTALL_PREFIX: ${CMAKE_INSTALL_PREFIX}")
message("   CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
message("   CMAKE_INCLUDE_PATH: ${CMAKE_INCLUDE_PATH}")
message("   CMAKE_LIBRARY_PATH: ${CMAKE_LIBRARY_PATH}")
message("   CMAKE_MODULE_PATH: ${CMAKE_MODULE_PATH}")
message("   CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread -Wall -O3 -march=native")

set(PYBIND11_CPP_STANDARD -std=c++11)
set(PYBIND11_PYTHON_VERSION 2.7)


####################################################################################################
#
# Find external prerequisites.
#
# FIXME there is some ugly boilerplate here that I will improve some day.


set(FFTW3_DOC_URL "https://github.com/kmsmith137/proton/blob/master/doc/installing_prerequisites.md#user-content-fftw3")
set(PYBIND11_DOC_URL "https://github.com/kmsmith137/proton/blob/master/doc/installing_prerequisites.md#user-content-pybind11")
set(PREREQUISITES_MISSING "FALSE")

message("")
find_package(FFTW3)
find_package(PYBIND11_WRAPPER)


if (FFTW3_FOUND)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DHAVE_FFTW3")
  set(SIMPULSE_LIBS ${SIMPULSE_LIBS} ${FFTW3_LIBRARY})
  include_directories(${FFTW3_INCLUDE_DIR})
else()
   message("-- Fatal: fftw3 is required")
   set(PREREQUISITES_MISSING "TRUE")
endif ()

if (NOT pybind11_FOUND)
   message("-- Fatal: must have pybind11")
   set(PREREQUISITES_MISSING "TRUE")
endif()

if (PREREQUISITES_MISSING)
   message("")
   message(FATAL_ERROR "Required prerequisites are missing, see above")
endif()


####################################################################################################


# Determine PYTHON_MODULE_INSTALL_DIR.
# Needs to come after find_package(pybind11).
#
# FIXME to do:
#   - figure out how to error-check the python comand
#   - figure out how to display helpful output if the command failed (command line + stderr)
#   - move this to its own .cmake file

execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import distutils.sysconfig as sc; print sc.get_python_lib(prefix='${CMAKE_INSTALL_PREFIX}')"
    OUTPUT_VARIABLE PYTHON_MODULE_INSTALL_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Is python module install dir in sys.path?
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import sys; print int('${PYTHON_MODULE_INSTALL_DIR}' in sys.path)"
    OUTPUT_VARIABLE PYTHON_MODULE_INSTALL_DIR_IN_PYTHONPATH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Does toplevel install dir exist?
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import os; print int(os.path.isdir('${CMAKE_INSTALL_PREFIX}'))"
    OUTPUT_VARIABLE INSTALL_DIR_EXISTS
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Is toplevel install dir user-writeable?
execute_process(
    COMMAND "${PYTHON_EXECUTABLE}" -c "import os; print int(os.access('${CMAKE_INSTALL_PREFIX}', os.W_OK))"
    OUTPUT_VARIABLE INSTALL_DIR_WRITEABLE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)


####################################################################################################


# 'make uninstall' needs to come after determination of PYTHON_MODULE_INSTALL_DIR
add_custom_target(uninstall COMMAND ${PROJECT_SOURCE_DIR}/cmake/uninstall.py)


add_subdirectory(cpp)
add_subdirectory(include)
add_subdirectory(pybind11)
add_subdirectory(python)
add_subdirectory(tests)


# Symlink python modules into toplevel directory.
#
# Note that we symlink 'python/simpulse' in add_custom_target(), but symlink simpulse_pybind11.so in add_custom_command().
# It would be better to use add_custom_command() for both, but due to artificial cmake constraints, we can't have an output
# named 'simpulse' and a library named 'libsimpulse.so' in the same project.
#
# FIXME: this has the disadvantage that the 'python/simpulse' symlink doesn't get cleaned up by 'make clean'.  Can this be fixed?
#
# Note: when cutting-and-pasting this boilerplate into other CMakeLists.txt files,
# don't forget to change 'toplevel_symlinks' to 'SUBDIR_symlinks'.

add_custom_target(toplevel_symlinks ALL DEPENDS ${simpulse_pybind11_basename}
                  COMMAND ${PROJECT_SOURCE_DIR}/cmake/symlink.py python/simpulse simpulse)

add_custom_command(OUTPUT ${simpulse_pybind11_basename}
                   COMMAND ${PROJECT_SOURCE_DIR}/cmake/symlink.py $<TARGET_FILE:simpulse_pybind11> ${simpulse_pybind11_basename})


####################################################################################################


message("")
message("This configuration looks good!")
message("   CMAKE_CXX_FLAGS: ${CMAKE_CXX_FLAGS}")
message("   INCLUDE_PATH: ${CMAKE_INCLUDE_PATH}")
message("   PYTHON_EXECUTABLE: ${PYTHON_EXECUTABLE}")
message("   PYTHON_MODULE_INSTALL_DIR: ${PYTHON_MODULE_INSTALL_DIR}")
message("   LIBS: ${SIMPULSE_LIBS}")
message("")


if (NOT PYTHON_MODULE_INSTALL_DIR_IN_PYTHONPATH)
   message("-- Warning: python module install directory ${PYTHON_MODULE_INSTALL_DIR} is not in your PYTHONPATH")
   message("   You may want to add 'export PYTHONPATH=\$PYTHONPATH:${PYTHON_MODULE_INSTALL_DIR}' to your ~/.bashrc and log in again.")
   message("")
endif()

if (NOT INSTALL_DIR_EXISTS)
   message("-- Warning: cmake install directory ${CMAKE_INSTALL_PREFIX} does not exist (or is not a directory).  You may need to create it.")
   message("")
elseif(NOT INSTALL_DIR_WRITEABLE)
   message("-- Warning: cmake install directory ${CMAKE_INSTALL_PREFIX} is not user-writeable.")
   message("   You may need 'sudo' to install (or rerun with 'cmake -DCMAKE_INSTALL_PREFIX=\$HOME ...')")
   message("")
endif()
