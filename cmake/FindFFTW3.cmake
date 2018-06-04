# Initializes the following variables:
#
#   FFTW3_FOUND
#   FFTW3_INCLUDE_DIR
#   FFTW3_LIBRARY


set(FFTW3_FOUND "TRUE")
set(FFTW3_MISSING_FILES "")

# Get hints from pkg-config.
# I determined that 'fftw3' is the correct pkg-config name with 'pkg-config --list-all | grep fftw'.
find_package(PkgConfig)
pkg_check_modules(PC_FFTW3 QUIET fftw3)

find_path(FFTW3_INCLUDE_DIR fftw3.h
  HINTS ${PC_FFTW3_INCLUDEDIR} ${PC_FFTW3_INCLUDE_DIRS})  # watch out for underscore: INCLUDEDIR vs INCLUDE_DIRS

find_library(FFTW3_LIBRARY fftw3
  HINTS ${PC_FFTW3_LIBDIR} ${PC_FFTW3_LIBRARY_DIRS})  # watch out for underscore: LIBDIR vs LIBRARY_DIRS

if (NOT FFTW3_INCLUDE_DIR)
    set(FFTW3_MISSING_FILES "${FFTW3_MISSING_FILES} fftw3.h")
    set(FFTW3_FOUND "FALSE")
endif()

if (NOT FFTW3_LIBRARY)
    set(FFTW3_MISSING_FILES "${FFTW3_MISSING_FILES} libfftw3")
    set(FFTW3_FOUND "FALSE")
endif()


if (NOT FFTW3_FOUND)
   if ("${FFTW3_DOC_URL}" STREQUAL "")
      set(FFTW3_DOC_URL "http://www.fftw.org/")
   endif()

   message("-- Couldn't find fftw3 library.")
   message("   The following file(s) were not found: ${FFTW3_MISSING_FILES}")
   message("   If fftw3 is installed, maybe you need -DCMAKE_PREFIX_PATH=xx, -DCMAKE_INCLUDE_PATH=xx, or -DCMAKE_LIBRARY_PATH=xx on your cmake command line?")
   message("   If you need to install fftw3, start here: ${FFTW3_DOC_URL}")
   message("")
endif()
