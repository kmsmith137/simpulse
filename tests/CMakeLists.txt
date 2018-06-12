configure_file(test-pulse.py test-pulse.py COPYONLY)


####################################################################################################
#
# Symlink python modules into 'tests' directory.
#
# Note: when cutting-and-pasting this boilerplate into other CMakeLists.txt files,
# don't forget to change 'tests_symlinks' to 'OTHERDIR_symlinks'.


add_custom_target(tests_symlinks ALL DEPENDS ${simpulse_pybind11_basename} simpulse)

add_custom_command(OUTPUT ${simpulse_pybind11_basename}
                   COMMAND ${PROJECT_SOURCE_DIR}/cmake/symlink.py $<TARGET_FILE:simpulse_pybind11> ${simpulse_pybind11_basename}
                   DEPENDS simpulse_pybind11)

add_custom_command(OUTPUT simpulse
                   COMMAND ln -s ../python/simpulse simpulse)


####################################################################################################
