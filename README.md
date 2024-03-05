
![logo](images/logo_ndm.png)


## NDM

![image](images/cea_small.png)
[DES/ISAS/DMN/SRMP](https://www.universite-paris-saclay.fr/laboratoires/service-de-recherches-de-metallurgie-physique-des/isas/dmn)

### Overview

``NDM`` is fortran code computing empirical potential molecular dynamics (dynamique moléculaire en potentiels empiriques)
* etc TODO

See [NDM TODO better introduction](https://inis.iaea.org/collection/NCLCollectionStore/_Public/37/064/37064766.pdf)


NDM is built around
- [Fortran 2008+](https://fr.wikipedia.org/wiki/Fortran)
- [CMake](https://cmake.org/)
- [GNU compilers products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).
- [OneApi Intel products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).

All the input data are prescribed non-interactively,
using command lines and/or ASCII files.

NDM can be compiled and run on any Unix-like system (and TODO windows).


### Installation

```
git clone --branch ndm2021_cv ssh://gitolite@ssh-codev-tuleap.intra.cea.fr:2044/ndm/NDM.git NDM
```

### Compilation

### Compiling NDM with Cmake

#### Prerequisities
It is recommended to use Intel's oneAPI suite including Intel's mpi implementation, the mpiifort wrapper and MKL.
To use lammps it must be compiled as a shared library. The shared library Lammps compilation should create the required liblammps.so dynamic library. It is recommended to compile Lammps with the same compiler suite as NDM. see lammps doc

#### Basic build, in the NDM directory :

Steps to build NDM with the new CMakeLists.txt (requires cmake 3.20) :

```
mkdir build; cd build
cmake ..
cmake --build . --parallel 4
```

#### Advanced build :
Configuration parameters should be set in a preset file, examples can be found in cmake_files/. It is possible to define preprocessor directives, set a specific compiler, define release/debug flags, etc...

```
cmake -C <absolute_path_to_preset_file> -B <absolute_path_to_build_directory> -S <absolute_path_to_CMakeLists.txt>
cmake --build <absolute_path_to_build_directory> --parallel <n_proc>
```

the bash script compile_ndm.sh can be found in the scripts/ directory for convenience.
When Cmake fails to find the relevant dependencies it is often required to define additionnal variables to indiquate the libraries locations such as MPI_HOME, MPI_ROOT, etc... They are listed in the bash script.





Compilation instructions are provided in the documentation included in the
distribution repository
- file [README_ndm_compilation.md](READMES/README_ndm_compilation.md)

Also see TODO
- github documentation [ndm-docs TODO](https://jpc.github.io/ndm-docs/contents/installation.html).


### Integration Tests TODO

- Data are small,
  located at [NDM/examples directories](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM) branch `ndm2021_cv`.

- Data are big,
  located in a [NDM_TESTS separate git repository](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM_TESTS.git) TODO.
