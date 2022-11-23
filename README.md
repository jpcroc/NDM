
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
