#!/usr/bin/env python

import build_helpers

# If called recursively in superbuild, a global persistent HeavyHandedUninstaller will be returned.
u = build_helpers.get_global_heavy_handed_uninstaller()

u.uninstall_headers('simpulse.hpp')
u.uninstall_headers('simpulse_internals.hpp')
u.uninstall_headers('simpulse/*.hpp')
u.uninstall_headers('simpulse/')
u.uninstall_libraries('libsimpulse.*')
u.uninstall_python_package('simpulse')

# If called recursively in superbuild, run() will not be called here.
if __name__ == '__main__':
    u.run()
