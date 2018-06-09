#!/usr/bin/env python
#
# CMAKE_INSTALL_MANIFEST_FILES

import os
import sys

print 'Command-line args:', sys.argv

flag = False

if not os.path.exists('install_manifest.txt'):
    print >>sys.stderr, "No install_manifest.txt"

for f in open('install_manifest.txt'):
    f = f.rstrip('\r\n')
    if os.path.exists(f):
        print 'Uninstalling', f
        os.remove(f)
        flag = True

if not flag:
    print "Nothing to do for 'make uninstall'"
