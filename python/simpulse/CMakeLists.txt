install(FILES
  __init__.py
  utils.py
  DESTINATION ${PYTHON_MODULE_INSTALL_DIR}/simpulse)

# Copy python files into 'build/python/simpulse', so they can be imported
# from the build directory without 'make install'.
#
# FIXME what happens for an in-source build?  Does this still work?

add_custom_target(copy_python_files ALL DEPENDS __init__.py utils.py)

add_custom_command(
    OUTPUT __init__.py
    COMMAND cp -f ${PROJECT_SOURCE_DIR}/python/simpulse/__init__.py __init__.py
    DEPENDS __init__.py
)

add_custom_command(
    OUTPUT utils.py
    COMMAND cp -f ${PROJECT_SOURCE_DIR}/python/simpulse/utils.py utils.py
    DEPENDS utils.py
)
