
## Prerequisites

Prerequisite are *only* a PYTHON-3 environment.

- Using [Miniconda](http://conda.pydata.org/miniconda.html). *It is easy*.
- Mandatory prerequisites python3, PyQt5, numpy, pandas, jupyter etc.
- See List python environment configuration

Example of activate a miniconda environment on is246206:

```bash

# set conda in $PATH etc.
source ${HOME}/commons/conda_initialize_is246206

conda env list
  # conda environments:
  #
  base                  *  /export/home/catA/miniconda
  py3qt5                   /export/home/catA/miniconda/envs/py3qt5
  pypinns                  /export/home/catA/miniconda/envs/pypinns

conda activate py3qt5

which python
  /export/home/catA/miniconda/envs/py3qt5/bin/python

envs -e PATH
  PATH                           = /export/home/catA/miniconda/envs/py3qt5/bin
                                   /export/home/catA/miniconda/condabin
                                   .
                                   /home/catA/wambeke
                                   /home/catA/wambeke/commons                             

```


### Execute an howto example

- Do that after reading below chapter for best comprehension.
- Directly (**outside jupyter**) CLI execution example, **for simple test**.

```bash
cd .../PACKAGESPY/doc

python
>> import PACKAGESPY_howto_level_2 as HT2
>> HT2.test_2()
```


### Some MINICONDA tricks

Here for memory.


#### List python configurarion

For example.

```
conda info

     active environment : py3qt5
    active env location : /export/home/catA/miniconda/envs/py3qt5
            shell level : 2
       user config file : /home/catA/wambeke/.condarc
 populated config files :
          conda version : 22.9.0
    conda-build version : not installed
         python version : 3.9.12.final.0
       virtual packages : __cuda=11.6=0
                          __linux=5.17.12=0
                          __glibc=2.33=0
                          __unix=0=0
                          __archspec=1=x86_64
       base environment : /export/home/catA/miniconda  (writable)
      conda av data dir : /export/home/catA/miniconda/etc/conda
  conda av metadata url : None
           channel URLs : https://repo.anaconda.com/pkgs/main/linux-64
                          https://repo.anaconda.com/pkgs/main/noarch
                          https://repo.anaconda.com/pkgs/r/linux-64
                          https://repo.anaconda.com/pkgs/r/noarch
          package cache : /export/home/catA/miniconda/pkgs
                          /home/catA/wambeke/.conda/pkgs
       envs directories : /export/home/catA/miniconda/envs
                          /home/catA/wambeke/.conda/envs
               platform : linux-64
             user-agent : conda/22.9.0 requests/2.28.1 CPython/3.9.12 Linux/5.17.12-100.fc34.x86_64 fedora/34 glibc/2.33
                UID:GID : 94740:470572
             netrc file : None
           offline mode : False

```

#### List python environment configuration

For example.

```
conda activate py3qt5
conda env export > ./doc/conda_environment_py3qt5.yaml
cat ./doc/conda_environment_py3qt5.yaml
  name: py3qt5
  channels:
    - conda-forge
    - defaults
  dependencies:
    - _libgcc_mutex=0.1=conda_forge
    - _openmp_mutex=4.5=2_kmp_llvm
    - argon2-cffi=21.3.0=pyhd3eb1b0_0
    - argon2-cffi-bindings=21.2.0=py39h7f8727e_0
    - asttokens=2.0.5=pyhd3eb1b0_0
    - attrs=21.4.0=pyhd3eb1b0_0
    - backcall=0.2.0=pyhd3eb1b0_0
    - beautifulsoup4=4.11.1=py39h06a4308_0
    - blas=1.0=mkl
    - bleach=4.1.0=pyhd3eb1b0_0
    - bottleneck=1.3.4=py39hce1f21e_0
    - brotli=1.0.9=he6710b0_2
    - brotlipy=0.7.0=py39h27cfd23_1003
    - ca-certificates=2022.07.19=h06a4308_0
    - certifi=2022.9.24=py39h06a4308_0
    - cffi=1.15.0=py39hd667e15_1
    - charset-normalizer=2.0.4=pyhd3eb1b0_0
    - cryptography=37.0.1=py39h9ce1e76_0
    - cycler=0.11.0=pyhd3eb1b0_0
    - cython=0.29.28=py39h295c915_0
    - dbus=1.13.18=hb2f20db_0
    - debugpy=1.5.1=py39h295c915_0
    - decorator=5.1.1=pyhd3eb1b0_0
    - defusedxml=0.7.1=pyhd3eb1b0_0
    - distro=1.5.0=pyhd3eb1b0_1
    - entrypoints=0.4=py39h06a4308_0
    - executing=0.8.3=pyhd3eb1b0_0
    - expat=2.4.4=h295c915_0
    - fontconfig=2.13.1=h6c09931_0
    - fonttools=4.25.0=pyhd3eb1b0_0
    - freetype=2.11.0=h70c0345_0
    - giflib=5.2.1=h7b6447c_0
    - glib=2.69.1=h4ff587b_1
    - gst-plugins-base=1.14.0=h8213a91_2
    - gstreamer=1.14.0=h28cd5cc_2
    - icu=58.2=he6710b0_3
    - idna=3.3=pyhd3eb1b0_0
    - intel-openmp=2021.4.0=h06a4308_3561
    - ipykernel=6.9.1=py39h06a4308_0
    - ipython=8.3.0=py39h06a4308_0
    - ipython_genutils=0.2.0=pyhd3eb1b0_1
    - ipywidgets=7.6.5=pyhd3eb1b0_1
    - jedi=0.18.1=py39h06a4308_1
    - jinja2=3.0.3=pyhd3eb1b0_0
    - jpeg=9e=h7f8727e_0
    - jupyter=1.0.0=py39h06a4308_7
    - jupyter_client=7.2.2=py39h06a4308_0
    - jupyter_console=6.4.3=pyhd3eb1b0_0
    - jupyter_core=4.10.0=py39h06a4308_0
    - jupyterlab_pygments=0.1.2=py_0
    - jupyterlab_widgets=1.0.0=pyhd3eb1b0_1
    - kiwisolver=1.4.2=py39h295c915_0
    - lcms2=2.12=h3be6417_0
    - ld_impl_linux-64=2.38=h1181459_1
    - libffi=3.3=he6710b0_2
    - libgcc-ng=12.1.0=h8d9b700_16
    - libpng=1.6.37=hbc83047_0
    - libsodium=1.0.18=h7b6447c_0
    - libstdcxx-ng=12.1.0=ha89aaad_16
    - libtiff=4.2.0=h2818925_1
    - libuuid=1.0.3=h7f8727e_2
    - libwebp=1.2.2=h55f646e_0
    - libwebp-base=1.2.2=h7f8727e_0
    - libxcb=1.15=h7f8727e_0
    - libxml2=2.9.14=h74e7548_0
    - libxslt=1.1.35=h4e12654_0
    - libzlib=1.2.12=h166bdaf_2
    - llvm-openmp=14.0.4=he0ac6c6_0
    - lxml=4.8.0=py39h1f438cf_0
    - lz4-c=1.9.3=h295c915_1
    - markupsafe=2.1.1=py39h7f8727e_0
    - matplotlib=3.5.1=py39h06a4308_1
    - matplotlib-base=3.5.1=py39ha18d171_1
    - matplotlib-inline=0.1.2=pyhd3eb1b0_2
    - mistune=0.8.4=py39h27cfd23_1000
    - mkl=2021.4.0=h06a4308_640
    - mkl-service=2.4.0=py39h7f8727e_0
    - mkl_fft=1.3.1=py39hd3c417c_0
    - mkl_random=1.2.2=py39h51133e4_0
    - munkres=1.1.4=py_0
    - nbclient=0.5.13=py39h06a4308_0
    - nbconvert=6.4.4=py39h06a4308_0
    - nbformat=5.3.0=py39h06a4308_0
    - ncurses=6.3=h7f8727e_2
    - nest-asyncio=1.5.5=py39h06a4308_0
    - notebook=6.4.11=py39h06a4308_0
    - numexpr=2.8.1=py39h807cd23_2
    - numpy=1.22.3=py39he7a7128_0
    - numpy-base=1.22.3=py39hf524024_0
    - openssl=1.1.1q=h7f8727e_0
    - packaging=21.3=pyhd3eb1b0_0
    - pandas=1.4.2=py39h295c915_0
    - pandas-datareader=0.10.0=pyhd3eb1b0_0
    - pandocfilters=1.5.0=pyhd3eb1b0_0
    - parso=0.8.3=pyhd3eb1b0_0
    - pcre=8.45=h295c915_0
    - pexpect=4.8.0=pyhd3eb1b0_3
    - pickleshare=0.7.5=pyhd3eb1b0_1003
    - pillow=9.0.1=py39h22f2fdc_0
    - pip=21.2.4=py39h06a4308_0
    - ply=3.11=py_1
    - prometheus_client=0.13.1=pyhd3eb1b0_0
    - prompt-toolkit=3.0.20=pyhd3eb1b0_0
    - prompt_toolkit=3.0.20=hd3eb1b0_0
    - ptyprocess=0.7.0=pyhd3eb1b0_2
    - pure_eval=0.2.2=pyhd3eb1b0_0
    - pycparser=2.21=pyhd3eb1b0_0
    - pygments=2.11.2=pyhd3eb1b0_0
    - pyomo=6.4.2=py39h5a03fae_0
    - pyopenssl=22.0.0=pyhd3eb1b0_0
    - pyparsing=3.0.4=pyhd3eb1b0_0
    - pyqt=5.9.2=py39h2531618_6
    - pyrsistent=0.18.0=py39heee7806_0
    - pysocks=1.7.1=py39h06a4308_0
    - python=3.9.12=h12debd9_1
    - python-dateutil=2.8.2=pyhd3eb1b0_0
    - python-fastjsonschema=2.15.1=pyhd3eb1b0_0
    - python_abi=3.9=2_cp39
    - pytz=2022.1=py39h06a4308_0
    - pyyaml=6.0=py39h7f8727e_1
    - pyzmq=22.3.0=py39h295c915_2
    - qt=5.9.7=h5867ecd_1
    - qtconsole=5.3.0=pyhd3eb1b0_0
    - qtpy=2.0.1=pyhd3eb1b0_0
    - readline=8.1.2=h7f8727e_1
    - requests=2.28.0=py39h06a4308_0
    - send2trash=1.8.0=pyhd3eb1b0_1
    - setuptools=61.2.0=py39h06a4308_0
    - sip=4.19.13=py39h295c915_0
    - six=1.16.0=pyhd3eb1b0_1
    - soupsieve=2.3.1=pyhd3eb1b0_0
    - sqlite=3.38.5=hc218d9a_0
    - stack_data=0.2.0=pyhd3eb1b0_0
    - terminado=0.13.1=py39h06a4308_0
    - testpath=0.6.0=py39h06a4308_0
    - tk=8.6.12=h1ccaba5_0
    - tornado=6.1=py39h27cfd23_0
    - traitlets=5.1.1=pyhd3eb1b0_0
    - typing-extensions=4.1.1=hd3eb1b0_0
    - typing_extensions=4.1.1=pyh06a4308_0
    - tzdata=2022a=hda174b7_0
    - urllib3=1.26.9=py39h06a4308_0
    - wcwidth=0.2.5=pyhd3eb1b0_0
    - webencodings=0.5.1=py39h06a4308_1
    - wheel=0.37.1=pyhd3eb1b0_0
    - widgetsnbextension=3.5.2=py39h06a4308_0
    - xz=5.2.5=h7f8727e_1
    - yaml=0.2.5=h7b6447c_0
    - zeromq=4.3.4=h2531618_0
    - zlib=1.2.12=h7f8727e_2
    - zstd=1.5.2=ha4553b6_0
    - pip:
      - alabaster==0.7.12
      - anyio==3.6.1
      - babel==2.10.3
      - click==8.1.3
      - colorama==0.4.5
      - docutils==0.17.1
      - gitdb==4.0.9
      - gitpython==3.1.27
      - greenlet==1.1.3
      - imagesize==1.4.1
      - importlib-metadata==4.12.0
      - jsonschema==3.2.0
      - jupyter-book==0.13.1
      - jupyter-cache==0.4.3
      - jupyter-server==1.18.1
      - jupyter-server-mathjax==0.2.6
      - jupyter-sphinx==0.3.2
      - latexcodec==2.0.1
      - linkify-it-py==1.0.3
      - markdown-it-py==1.1.0
      - mdit-py-plugins==0.2.8
      - myst-nb==0.13.2
      - myst-parser==0.15.2
      - nbdime==3.1.1
      - psutil==5.9.1
      - pybtex==0.24.0
      - pybtex-docutils==1.0.2
      - pydata-sphinx-theme==0.8.1
      - smmap==5.0.0
      - sniffio==1.2.0
      - snowballstemmer==2.2.0
      - sphinx==4.5.0
      - sphinx-book-theme==0.3.3
      - sphinx-comments==0.0.3
      - sphinx-copybutton==0.5.0
      - sphinx-design==0.1.0
      - sphinx-external-toc==0.2.4
      - sphinx-jupyterbook-latex==0.4.6
      - sphinx-multitoc-numbering==0.1.3
      - sphinx-thebe==0.1.2
      - sphinx-togglebutton==0.3.2
      - sphinxcontrib-applehelp==1.0.2
      - sphinxcontrib-bibtex==2.5.0
      - sphinxcontrib-devhelp==1.0.2
      - sphinxcontrib-htmlhelp==2.0.0
      - sphinxcontrib-jsmath==1.0.1
      - sphinxcontrib-qthelp==1.0.3
      - sphinxcontrib-serializinghtml==1.1.5
      - sqlalchemy==1.4.40
      - uc-micro-py==1.0.1
      - websocket-client==1.4.0
      - zipp==3.8.1
  prefix: /export/home/catA/miniconda/envs/py3qt5

```
