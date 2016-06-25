# Makefile.local must define the following variables
#   LIBDIR      install dir for C++ libraries
#   INCDIR      install dir for C++ headers
#   PYDIR       install dir for python/cython modules
#   CPP         C++ compiler command line
#
# Some optional variables which I only use for osx/clang:
#   CPP_LFLAGS      extra linker flags when creating a .so or executable file from .o files
#   LIBS_PYMODULE   any extra libraries needed to link a python extension module (osx needs -lPython)
#
# See site/Makefile.local.* for examples.

include Makefile.local

all: libsimpulse.so cython/simpulse.so

install: libsimpulse.so cython/simpulse.so
	cp -f simpulse.hpp $(INCDIR)/
	cp -f libsimpulse.so $(LIBDIR)/
	cp -f cython/simpulse.so $(PYDIR)/

uninstall:
	rm -f $(INCDIR)/simpulse.hpp $(LIBDIR)/libsimpulse.so $(PYDIR)/simpulse.so

clean:
	rm -f *~ *.o *.so cython/*~ cython/*.so cython/simpulse.cpp visual_check/*~ visual_check/plot*.png

%.o: %.cpp simpulse.hpp
	$(CPP) -c -o $@ $<

libsimpulse.so: single_pulse.o
	$(CPP) $(CPP_LFLAGS) -o $@ -shared $^ -lfftw3

cython/simpulse.cpp: cython/simpulse.pyx cython/simpulse_pxd.pxd cython/simpulse_cython.hpp simpulse.hpp
	cython --cplus $<

cython/simpulse.so: cython/simpulse.cpp libsimpulse.so
	$(CPP) $(CPP_LFLAGS) -Wno-unused-function -shared -o $@ $< -lsimpulse $(LIBS_PYMODULE)
