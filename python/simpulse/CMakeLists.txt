install(FILES
  __init__.py
  utils.py
  DESTINATION ${PYTHON_MODULE_INSTALL_DIR}/simpulse)

# Copy python files into 'build/python/proton', so they can be imported
# from the build directory without 'make install'.

configure_file(__init__.py __init__.py COPYONLY)
configure_file(utils.py utils.py COPYONLY)
