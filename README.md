
![logo](images/logo_ndm.png)


## NDM

![image](images/cea_small.png)
<!--- TODO: following is a dead link -->
[DES/ISAS/DMN/SRMP](https://www.universite-paris-saclay.fr/laboratoires/service-de-recherches-de-metallurgie-physique-des/isas/dmn)

### Overview

``NDM`` is fortran code computing empirical potential molecular dynamics (dynamique moléculaire en potentiels empiriques)


NDM is built around
- [Fortran 2008+](https://fr.wikipedia.org/wiki/Fortran)
- [CMake](https://cmake.org/)
- [GNU compilers products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).
- [OneApi Intel products](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html).

All the input data are prescribed non-interactively,
using command lines and/or ASCII files.

NDM can be compiled and run on any Unix-like system .

### Compilation

NDM can be compiled using cmake (version 3.20).

#### Prerequisities
It is recommended to use Intel's oneAPI suite including Intel's mpi implementation, the mpiifort wrapper and MKL but GNU compilers are also available.

Following are environment setup examples for some known clusters:

|   Cluster         |        NDM           | NDM + LAMMPS   |
|--------------------|----------------------|----------------|
| Topaze/Irene INTEL | <pre>module load fortran/inteloneapi/24.0.0<br>module load mpi/intelmpi/24.0.0</pre> | x |
| Topaze/Irene MIX | <pre>module load fortran/inteloneapi/24.0.0<br>module load mpi/openmpi/4.1.6.4</pre> | <pre>module load fortran/inteloneapi/24.0.0<br>module load mpi/openmpi/4.1.6.4<br>module load lammps/2Aug2023</pre> |
| Topaze/Irene GNU   | <pre>module load fortran/gcc/13.2.0<br>module load mpi/openmpi/4.1.4</pre> | x |
| Gatsby GNU   | <pre>module load openmpi/gcc/4.0.2<br>module load gcc/11.2.0<br>module load mkl</pre>| x |
| Gatsby INTEL   | <pre>module load mpi<br>module load compiler<br>module load mkl</pre> | <pre>module load lammps/29Sep2021-u2</pre> |


#### Basic cmake build, in the NDM directory :

Steps to build NDM using the `CMakeLists.txt` file (requires cmake 3.20) :

```
mkdir build; cd build
cmake ..
make -j
```

This should detect the required compilers and dependencies if they are available and find a suitable configuration.
If not, make sure the prerequisities are met and consider defining advanced options manually.

#### NDM+LAMMPS build

```
cmake .. -D NDM_PACKAGE_LIST=LAMMPS
```

#### NDM+MILADY build

```
export MILADY_ROOT=PATH/TO/ML/DIRECTORY
cmake .. -D NDM_PACKAGE_LIST=MILADY
```
<!---
#### Advanced options :

Options can be passed through the command line :
```
cmake .. -D variable=value
```

Or using a cmake preset file :
```
cmake .. -C <preset_file.cmake>
```

To set the environment and specify source or build directories you can use the `scripts/compile_ndm.sh` bash script.
-->

#### Available options

- Compilers : `CMAKE_Fortran_COMPILER`, `CMAKE_CXX_COMPILER`
- preprocessor definitions : `NDM_COMPILE_DEFINITION`
- compilation configuration : `CMAKE_BUILD_TYPE` (`RELEASE` or `DEBUG`)
- compilation flags : `CMAKE_<lang>_FLAGS_<config>`, for example `CMAKE_Fortran_FLAGS_RELEASE`
- MPI : `MPI_HOME`
- LAMMPS library directory : `LAMMPS_HOME`
- LAMMPS library name : `LAMMPS_LIBRARY_NAME` (for example `lammps_serial`, `lammps_mpi` or `lammps`)
- Additional dependencies for LAMMPS: `LAMMPS_EXTRA_LIBRARIES` (for example `gomp` to use open mpi LAMMPS)

examples:

  cmake .. -D NDM_PACKAGE_LIST="LAMMPS"
  
  cmake .. -D NDM_PACKAGE_LIST=MILADY
  
  cmake .. -DCMAKE_BUILD_TYPE=DEBUG -DCMAKE_Fortran_FLAGS_DEBUG="-O0 -g -C -fpe-all=0 -traceback"

in this last example the debug options are passed explicitly (for ifort) 

To specify multiple values, for example for preprocessor definitions or packages, use semicolon separator:
```
cmake .. -D NDM_PACKAGE_LIST="MPI;LAMMPS"
```

<!---
#### Compile NDM with LAMMPS

You should go to the READMES directory for detailled instrictions !
To compile NDM with LAMMPS a version of LAMMPS must be compiled as a library with all the necessary header files.
On Gatsby the following command can be used:
```
module load lammps/29Sep2021-u2
```
To compile a different version of LAMMPS from the sources, see the relevant README file in the NDM `READMES` directory.

Once LAMMPS is installed it can be added to the compilation using the NDM_PACKAGE_LIST variable:
```
cmake .. -D NDM_PACKAGE_LIST="LAMMPS"
```
If CMake doesn't find LAMMPS, the following variables can be set :
- LAMMPS_LIBRARY_NAME : LAMMPS is often compiled with a custom library name depending on the configuration using LAMMPS LAMMPS_MACHINE option. 
For example `lammps/29Sep2021-u2` has both `lammps_serial` and `lammps_mpi` installed. Since LAMMPS 2021 the executable and library have a similar name (lmp_serial and lammps_serial for example) but not for older versions. To use a specific library or when CMake doesn't find any, the library name must be set explicitly.

- LAMMPS_HOME : In case of a local installation it may be necessary to specify the directory where the library (liblammps.a or liblammps.so) or .pc file is located.

- LAMMPS_EXTRA_LIBRARIES : Depending on LAMMPS configuration and packages it may have further dependencies. 
CMake will find some of these dependencies automatically but not all of them since they are not listed in the .pc file.
When such dependencies are missing the linker will show `undefined references` at link stage which usually hint toward which library is missing and what LAMMPS package is causing the error.
These libraries can be identified at the end of a successful LAMMPS compilation in the list of link libraries or in the executable dependencies.
Then they can be added manually to the NDM compilation using the LAMMPS_EXTRA_LIBRARIES variable with either just the name or the full path to the library.
This is for example necessary when compiling LAMMPS with Open MPI or with the PYTHON package.
-->


### Integration Tests

After compilation with tests (`NDM_OPT_TEST=ON`, by default) you can run tests from the build directory using the `ctest` command:
```
# Run all tests
ctest

# Run ndm-only tests:
ctest -R ndm

# Run milady tests:
ctest -R mld
```

<!---
TODO: introduce non-regression tests:
For now, tests only check that NDM runs without failure.

- Data are small,
  located at [NDM/examples directories](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM) branch `ndm2021_cv`.

- Data are big,
  located in a [NDM_TESTS separate git repository](https://codev-tuleap.intra.cea.fr/plugins/git/ndm/NDM_TESTS.git) TODO.
-->

### Code coverage

Tests coverage gives information on which blocks and functions of the code are used when running tests and can be useful to set relevant tests.

Only works with intel ifort, compile NDM with option `NDM_OPT_COVERAGE=ON` and tests. After compilation, run tests with `ctest` and read the output html file with any web browser, for example:

```
ctest
firefox CODE_COVERAGE.HTML
```
