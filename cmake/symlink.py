#!/usr/bin/env python
#
# Called in various places in CMakeLists.txt, to create symlinks.
#
# We symlink pybind11/python modules into the build directory and some of its
# subdirectories (build/tests, build/sphinx, etc.) so that they can be imported
# from python, without first doing 'make install'.

import os
import sys

if len(sys.argv) != 3:
    print >>sys.stderr, 'usage: symlink.py <src_abspath> <dst_basename>'
    sys.exit(1)

src_abspath = sys.argv[1]
dst_basename = sys.argv[2]
src_relpath = os.path.relpath(src_abspath, os.getcwd())

if dst_basename != os.path.basename(src_abspath):
    # FIXME improve this somehow
    print >>sys.stderr, '\n!!! Internal error: basename mismatch in cmake/symlink.py !!!'
    print >>sys.stderr, '!!! (email KMS to debug, this is probably a pybind11 configuration issue) !!!\n'
    sys.exit(1)

os.symlink(src_relpath, dst_basename)
