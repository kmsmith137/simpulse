#!/usr/bin/env python
#
# This script is called by 'make uninstall'

import os
import sys

flag = False

if not os.path.exists('install_manifest.txt'):
    print >>sys.stderr, "No install_manifest.txt"

def generate_installed_files():
    for f in open('install_manifest.txt'):
        f = f.rstrip('\r\n')
        yield f

        # If .py is being uninstalled, then uninstall .pyc too
        if f.endswith('.py'):
            yield f + 'c'

for f in generate_installed_files():
    if os.path.exists(f):
        print 'Uninstalling', f
        os.remove(f)
        flag = True

    # If directory is "emptied out", remove the empty directory.
    # (This way of doing it isn't optimally efficient, but that's OK!)

    d = os.path.dirname(f)

    while d != f:
        try:
            os.rmdir(d)
        except:
            break

        print 'Removing empty directory', d
        (d, f) = (os.path.dirname(d), d)
        flag = True

if not flag:
    print "Nothing to do for 'make uninstall'"
