## Debug CMake problems

As an example

### Sometimes CMake have problems

```
cmake -S/home/vanwambeke/GitHub/MLD/MILADY -B/home/vanwambeke/GitHub/MLD/build -DMLD_OPT_TRACE=OFF
-- The Fortran compiler identification is Intel 2021.2.0.20210228
-- The C compiler identification is Intel 2021.2.0.20210228
-- Detecting Fortran compiler ABI info
-- Detecting Fortran compiler ABI info - done
-- Check for working Fortran compiler: /usr/local/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/ifort - skipped
-- Detecting C compiler ABI info
-- Detecting C compiler ABI info - done
-- Check for working C compiler: /usr/local/intel/oneapi/compiler/2021.2.0/linux/bin/intel64/icc - skipped
-- Detecting C compile features
-- Detecting C compile features - done
-- CMAKE_MODULE_PATH : /home/vanwambeke/GitHub/MLD/MILADY/cmake_files
-- Could NOT find MPI_C (missing: MPI_C_WORKS)
-- Found MPI_Fortran: /usr/local/iopenmpi/lib64/libmpi_usempif08.so (found version "3.1")
CMake Error at /usr/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:230 (message):
  Could NOT find MPI (missing: MPI_C_FOUND) (found version "3.1")
Call Stack (most recent call first):
  /usr/share/cmake/Modules/FindPackageHandleStandardArgs.cmake:594 (_FPHSA_FAILURE_MESSAGE)
  /usr/share/cmake/Modules/FindMPI.cmake:1830 (find_package_handle_standard_args)
  CMakeLists.txt:50 (find_package)


-- Configuring incomplete, errors occurred!
See also "/home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeOutput.log".
See also "/home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeError.log".
problem cmake command
```

### Fix them

RTFM(1) and do

(1) *Read The Fucking Message*

```
See also "/home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeOutput.log".
See also "/home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeError.log".
```

And so on

```
pluma /home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeOutput.log &
pluma /home/vanwambeke/GitHub/MLD/build/CMakeFiles/CMakeError.log &
```

And searching string `warning` and/or `error` 
or `something else`.

To find

```
/usr/lib64/gcc/x86_64-suse-linux/11/../../../../x86_64-suse-linux/bin/ld: 
warning: libutil.so, 
needed by /usr/local/iopenmpi/lib64/libmpi_usempif08.so, not found 
(try using -rpath or -rpath-link)
```

To conclude that on Sasha SUSE computer library */usr/local/iopenmpi* 
have been compiled on another SUSE computer (with not same versioned system)

May be.






