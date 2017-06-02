# Makefile.local must define the following variables
#   LIBDIR      install dir for C++ libraries
#   INCDIR      install dir for C++ headers
#   PYDIR       install dir for python/cython modules
#   CPP         C++ compiler command line, including all flags
#
# Some optional variables which I only use for osx/clang:
#   CPP_LFLAGS      extra linker flags when creating a .so or executable file from .o files
#   LIBS_PYMODULE   any extra libraries needed to link a python extension module (osx needs -lPython)
#
# Note that CPP must include -I<dir> flags which point to the python and numpy include
# directories.  These directories can be determined from python as follows:
#   print distutils.sysconfig.get_python_inc()    # python includes
#   print numpy.get_include()                     # numpy includes
#
# See site/Makefile.local.* for examples.

include Makefile.local

ifndef CPP
$(error Fatal: Makefile.local must define CPP variable)
endif

ifndef INCDIR
$(error Fatal: Makefile.local must define INCDIR variable)
endif

ifndef LIBDIR
$(error Fatal: Makefile.local must define LIBDIR variable)
endif

ifndef PYDIR
$(error Fatal: Makefile.local must define PYDIR variable)
endif


####################################################################################################


all: libsimpulse.so cython/simpulse.so

install: libsimpulse.so cython/simpulse.so
	mkdir -p $(INCDIR) $(LIBDIR) $(PYDIR)
	cp -f simpulse.hpp simpulse_internals.hpp $(INCDIR)/
	cp -f libsimpulse.so $(LIBDIR)/
	cp -f cython/simpulse.so $(PYDIR)/

uninstall:
	rm -f $(INCDIR)/simpulse.hpp $(INCDIR)/simpulse_internals.hpp $(LIBDIR)/libsimpulse.so $(PYDIR)/simpulse.so

clean:
	rm -f *~ *.o *.so cython/*~ cython/*.so cython/simpulse.cpp visual_check/*~ visual_check/plot*.png site/*~

%.o: %.cpp simpulse.hpp simpulse_internals.hpp
	$(CPP) -c -o $@ $<

libsimpulse.so: constant_acceleration_phase_model.o single_pulse.o von_mises_profile.o
	$(CPP) $(CPP_LFLAGS) -o $@ -shared $^ -lfftw3

cython/simpulse.cpp: cython/simpulse.pyx cython/simpulse_pxd.pxd cython/simpulse_cython.hpp simpulse.hpp
	cython --cplus $<

cython/simpulse.so: cython/simpulse.cpp libsimpulse.so
	$(CPP) $(CPP_LFLAGS) -Wno-unused-function -shared -o $@ $< -lsimpulse $(LIBS_PYMODULE)
