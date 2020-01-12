#!/usr/bin/env python

import build_helpers

# If called recursively in superbuild, a global persistent LegacyUninstaller will be returned.
lu = build_helpers.get_global_legacy_uninstaller()

lu.uninstall_headers('simpulse.hpp')
lu.uninstall_headers('simpulse_internals.hpp')
lu.uninstall_headers('simpulse/*.hpp')
lu.uninstall_headers('simpulse/')
lu.uninstall_libraries('libsimpulse.*')
lu.uninstall_pyfiles('simpulse*.so')
lu.uninstall_pyfiles('simpulse/*.py')
lu.uninstall_pyfiles('simpulse/*.pyc')
lu.uninstall_pyfiles('simpulse/')

# If called recursively in superbuild, run() will not be called here.
if __name__ == '__main__':
    lu.run()
