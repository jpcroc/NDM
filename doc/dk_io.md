---
project: dk_io library
version: 0.1
author: Paul Fossati
email: paul.fossati@cea.fr
docmark_alt: @
md_extensions: markdown.extensions.tables
graph: true
incl_src: false
source: false
---

# Fortran library for file i/o

This library contains APIs to read and write files from Fortran programs based on the [[FileObject(type)]] derived type. Files compressed in supported formats are processed transparently and do not require any special treatment on the caller part.
It also contains subroutines to read and write a variety of structure files formats used in atomic-scale modeling codes.

## Supported compression formats

The following compression formats are supported.

 Format    | Extension | Library
-----------|-----------|-----------
 gzip      | .gz       | zlib
 bzip2     | .bz2      | libbzip2
 zstandard | .zst      | libzstd


Compressed file i/o use external libraries, which may have to be enabled at compile time. To check the available format, run the [[cfginfo(program)]] program with the `--info` option.

## Supported structure files formats

This table summarises the different features supported when reading the supported file formats.

 Format  | Positions | Velocities | Forces | Spin    | Masses | Charges | Auxiliaries | Multiple frames | Symmetry
---------|-----------|------------|--------|---------|--------|---------|-------------|-----------------|----------
 Abinit  | 🕒        | ❌         | ❌     | 🕒     | ❌    | ❌      | ❌         |❌              | ❌
 Atomeye | ✔        | ✔          | ❌ [0] | ❌ [0] | 🕒    | ❌      | ✔          |❌              | ❌
 CASTEP  | ✔        | ✔          | ❌     | 🕒     | ❌    | ❌      | ❌          |❌              | 🕒 [1]
 CIF     | ✔        | ❌         | ❌     | ❌     | ❌    | ❌      | ❌          |🕒              | ✔
 DL_POLY | ✔        | ✔          | ✔     | ❌     | ❌    | ❌      | ❌          |🕒              | ❌
 GULP    | ✔        | ❌         | ❌     | ❌     | ❌    | ❌      | ❌          |🕒              | ✔
 LAMMPS  | ✔        | 🕒         | ❌     | ❌     | ❌    | ❌      | ❌          |❌              | ❌
 VASP    | ✔        | 🕒         | ❌     | 🕒     | ❌    | ❌      | ❌          |❌              | ❌
 XYZ     | ✔        | 🕒         | 🕒 [0] | 🕒 [0] | 🕒    | 🕒      | ✔          |🕒              | ❌


[0]: Forces and magnetic moments can be read from the auxiliary fields.

[1]: The full structure can be build from the irreductible crystal sites if the SYMMETRY_OPS block is present.

# Utilities to handle structure files

This package includes also the [[cfginfo(program)]] and [[cfgconvert(program)]] programs.
