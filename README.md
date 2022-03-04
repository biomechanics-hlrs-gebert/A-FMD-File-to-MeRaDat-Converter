# File Format to Meta Converter
[![DOI](https://zenodo.org/badge/405168166.svg)](https://zenodo.org/badge/latestdoi/405168166)
[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
![Architecture](https://img.shields.io/badge/Architecture-x86_64-blue)
![OS](https://img.shields.io/badge/OS-Linux-blue)
![version](https://img.shields.io/badge/version-1.0.0-blue)

Program to convert \*.vtk files with structured points to the MeRaDat format and vice versa. MeRaDat expands to **Me**ta/**Ra**w **Dat**a-format. Part of my PhD research project. 

The program essentially introduces the meta file to the process chain.
##  MeRaDat
The restart handling is not adopted for this program as the need for compute time and the probability of failure are negligible.
## Usage:
The program currently only accepts *.vtk files with a proper meta basename.
```./fmf_v1.0.0_x86_64 <basename>.vtk```

## Requirements
* x86 64bit Hardware
* Linux x86 64Bit Installation with Bash or Zsh
* GNU Compiler Collection (GCC), especially with gfortran
* An installation of Open-MPI
* Geberts libraries. Managed by: ```./manage_geb-lib.sh```

The program must be compiled with:
* Global integer kind=64Bit, signed
* Meta-format integer kind=64Bit, signed
* MPI integer kind=32Bit

The installation of Open MPI is simplified with the install script of the repository "Overview" of the biomechanics-hlrs-gebert organization @GitHub.
### Optional: Gnu debugging
* [gdb](https://www.gnu.org/software/gdb/)
* [tmpi](https://github.com/Azrael3000/tmpi)
* [tmux](https://github.com/tmux/tmux/wiki)
## Build
It's tested and therefore recommended to build and run the program as follows.
### Set up the Environment
```vim ./auxiliaries/system_environments/<system>.sh```
```source ./environment.source <system>``` 

* Set an architecture/a system
  * Give the absolute base path of your mpi-installation
  * Alternatively give the proper module names of your compute cluster

### Run make:
Build the program:    ```make```
Create documentation: ```make docs```

### Uninstall:
```make clean && rm -r <your program directory>```
## Acknowledgements 
Plain text parsed via [strings module](https://gbenthien.net/strings/index.html) by George Benthien from San Diego.
## Arbitrary
Use this program at your own risk.

