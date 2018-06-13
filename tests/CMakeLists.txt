link_directories(${CMAKE_SOURCE_DIR}/cpp)
link_directories(${CMAKE_SOURCE_DIR}/pybind11)

# List of all unit tests which will be run by 'make test'.
# When adding a new unit test, don't forget to add it here!

add_custom_target(
   check DEPENDS 
   run-test-pulsars
)


####################################################################################################
#
#
# Symlink python modules into toplevel directory.
#
# Note that we symlink 'python/simpulse' in add_custom_target(), but symlink simpulse_pybind11.so in add_custom_command().
# It would be better to use add_custom_command() for both, but due to artificial cmake constraints, we can't have an output
# named 'simpulse' and a library named 'libsimpulse.so' in the same project.
#
# FIXME: this has the disadvantage that the 'python/simpulse' symlink doesn't get cleaned up by 'make clean'.  Can this be fixed?
#
# Note: when cutting-and-pasting this boilerplate into other CMakeLists.txt files,
# don't forget to change 'tests_symlinks' to 'OTHERDIR_symlinks'.


add_custom_target(tests_symlinks ALL DEPENDS ${simpulse_pybind11_basename}
                  COMMAND ${PROJECT_SOURCE_DIR}/cmake/symlink.py ../python/simpulse simpulse)

add_custom_command(OUTPUT ${simpulse_pybind11_basename}
                   COMMAND ${PROJECT_SOURCE_DIR}/cmake/symlink.py $<TARGET_FILE:simpulse_pybind11> ${simpulse_pybind11_basename})


####################################################################################################


add_custom_command(
    OUTPUT test-pulsars-dontedit.py
    COMMAND cp -f ${PROJECT_SOURCE_DIR}/tests/test-pulsars.py test-pulsars-dontedit.py
    DEPENDS test-pulsars.py
)

# FIXME when more unit tests are added, how to avoid repeating dependencies?
# FIXME is there a way to avoid listing all the python files here??  (Just putting 'copy_python_files' as a dependency isn't enough)
#
# FIXME here and below, the "env DYLD_LIBRARY_PATH=..." is an osx hack to ensure that libraries in the build
# directory have higher priority than libraries installed with 'make install'.  Is there a better way??!

add_custom_command(
   OUTPUT run-test-pulsars
   COMMAND env DYLD_LIBRARY_PATH=../cpp:$$DYLD_LIBRARY_PATH ./test-pulsars-dontedit.py && touch run-test-pulsars
   DEPENDS test-pulsars-dontedit.py tests_symlinks copy_python_files simpulse_pybind11 simpulse ../python/simpulse/__init__.py ../python/simpulse/utils.py
)
