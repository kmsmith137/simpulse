simpulse: C++/python library for simulating pulses in radio astronomy

Right now there isn't much here, just a class representing a single channelized, dispersed, 
scattered pulse. In the future I may add more features (e.g. pulsars and scintillated pulses).

### INSTALLATION

- You'll need the following prerequisites.
    - C++ compiler which supports C++11
    - Cython (installation hint: `pip install Cython`)
    - FFTW3 (installation hint: `sudo yum install fftw-devel` or `brew install fftw`)

- Create a file `Makefile.local` defining a few Makefile variables (see brief 
  description in `Makefile`, or one of the examples in `site/`)
  
- Do `make all install`

- You can make a few example plots with
```
cd visual_check
./plot-pulses.py
```
