"""simpulse: simulating pulses in radio astronomy"""

from . import utils

try:
    from .simpulse_pybind11 import *
except:
    raise ImportError("Toplevel 'simpulse' python module couldn't import the 'simpulse_pybind11' submodule.  Maybe you need 'make install' in the simpulse build directory?")
