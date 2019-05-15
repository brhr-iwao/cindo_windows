
## CNDO2/INDO for Windows

### Overview
CNDO/2 (denotes "Complete Neglect Differential Overlap/2") and INDO (denotes "Intermediate Neglect of Differential Overlap") are semi-empirical molecular orbital methods developed by [J.A.Pople](https://www.nobelprize.org/prizes/chemistry/1998/pople/facts/) and co-workers.
These methods assume the complete neglect of differential overlap and the neglect of differential overlap in electron interaction integrals except those involving one center, respectively.   
The Windows program is released by the courtesy of Prof. S. Sahu, provides Fortran 90 implementation of CNDO/2 and INDO models.

### How to use
1. Download [cindo.zip](https://github.com/brhr-iwao/cindo_windows/releases/download/190515/cindo.zip) and unzip it.
2. Launch the Windows command prompt.
3. Change the current working directory to the cindo.exe contained directory.
4. Type "cindo" and press enter key to launch "cindo.exe".
5. Follow the command line instructions and enter the proper values or characters.
   You can abort the program anytime during the execution by pressing CTRL + C.

### How to compile from the source
#### Prerequisites
- GNU Fortran (gfortran) ([MinGW](http://www.mingw.org) build)
- GNU Make ([MinGW](http://www.mingw.org) build)
- [MinGW BLAS and LAPACK](https://sourceforge.net/projects/mingw-cross/files/mingw32-lapack-3.2.1-1.zip)
- [The original source code](http://www.cpc.cs.qub.ac.uk/summaries/AECN_v1_0.htm) which is available from [Computer Physics Communication Program library, Queen's University of Belfast](http://www.cpc.cs.qub.ac.uk).

  GNU Fortran and Make ([MinGW](http://www.mingw.org) build) can be installed on Windows PC with [MinGW Installer](https://sourceforge.net/projects/mingw/files/Installer).

#### Steps
1. Replace the original cindo.f90 and Makefile with the provided files in this repository to build the program on Windows with GNU Fortran.
2. Put blas.dll and lapak.dll on the source directory. Blas.dll and lapack.dll are packed in [mingw32-lapack-3.2.1-1.zip](https://sourceforge.net/projects/mingw-cross/files/mingw32-lapack-3.2.1-1.zip).
3. Launch the Windows command prompt and change the current working directory to the source directory.
4. Type "make" and press the enter key to make.

### Disclaimer and Non-profit use Licence Agreement
[CNDO2/INDO for Windows](https://github.com/brhr-iwao/cindo_windows) is subject to [the CRC Program Library Licence](http://cpc.cs.qub.ac.uk/licence/licence.html). Non-profit use and redistribution without modification are approved. [CNDO2/INDO for Windows](https://github.com/brhr-iwao/cindo_windows) is provided "as is" without any warranty. There is no guarantee of correctness of calculation results.
Cite the program as [S. Sahu and A. Shukla, Fortran 90 implementation of the Hartree-Fock approach within the CNDO/2 and INDO models, Computer Physics Communications, 180, 724, 2009 ](https://doi.org/10.1016/j.cpc.2008.11.004). The preprint is available from [arXiv.org](http://arxiv.org/pdf/0812.3690).
