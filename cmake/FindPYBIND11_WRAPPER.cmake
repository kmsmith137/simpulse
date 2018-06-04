# Initializes variables
#
#   pybind11_FOUND - true if pybind11 and all required components found on the system
#   pybind11_VERSION - pybind11 version in format Major.Minor.Release
#   pybind11_INCLUDE_DIRS - Directories where pybind11 and python headers are located.
#   pybind11_INCLUDE_DIR - Directory where pybind11 headers are located.
#   pybind11_DEFINITIONS - Definitions necessary to use pybind11, namely USING_pybind11.
#   pybind11_LIBRARIES - compile flags and python libraries (as needed) to link against.
#   pybind11_LIBRARY - empty.
#
# So far, I haven't had to use these variables anywhere (except pybind11_FOUND) because
# they get "hidden" by calls to pybind11_add_module().
#
# This is mostly a thin wrapper around find_package(pybind11), but gives more diagnostics.


find_path(PYBIND11_QUICKLOOK "pybind11/pybind11.h")

if ("${PYBIND11_DOC_URL}" STREQUAL "")
   set(PYBIND11_DOC_URL "https://github.com/pybind/pybind11")
endif()

if (NOT PYBIND11_QUICKLOOK)
   set(pybind11_FOUND, "FALSE")

   message("-- Couldn't find pybind11 library.")
   message("   The following file was not found: pybind11/pybind11.h")
   message("   If pybind11 is installed, maybe you need -DCMAKE_PREFIX_PATH=xx on your cmake command line?")
   message("      (Note: -DCMAKE_INCLUDE_PATH=xx probably won't work, since we need pybind11 .cmake files too.)")
   message("   If you need to install pybind11, start here: ${PYBIND11_DOC_URL}")
   message("")
else()
   find_package(pybind11)
   if (NOT pybind11_FOUND)
      message("-- Something strange has happened: pybind11/pybind11.h was found, but find_package(pybind11) failed.")
      message("")
      message("   If you got a multi-paragraph cmake error above beginning 'By not providing \"Findpybind11.cmake\" in CMAKE_MODULE_PATH...'")
      message("   then the problem is that cmake is not finding the pybind11 .cmake files.  This can happen if CMAKE_INCLUDE_PATH points")
      message("   to the pybind11 .h files, but CMAKE_MODULE_PATH does not point to the pybind11 .cmake files.  In this case, you should try")
      message("   setting CMAKE_MODULE_PATH (or CMAKE_PREFIX_PATH).")
      message("")
      message("   Otherwise, you may be in uncharted debugging territory!  Feel free to email me (KMS).")
      message("")
   endif()
endif()
