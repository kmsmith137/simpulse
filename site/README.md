The site/ directory contains Makefile.local files for a few
machines which I use frequently:

  - Makefile.local.kms_laptop16     osx 10.11 (El Capitan)
  - Makefile.local.orangutan        ubuntu 16.04
  - Makefile.local.frb1             centos 7.3.1611

plus a few others which are unlikely to be of interest.
To use one of these files, do:

  ln -s site/Makefile.local.<MACHINE> Makefile.local

from the toplevel simpulse directory.
