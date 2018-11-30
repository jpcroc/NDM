#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2018  CEA/DEN
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
# 
# See http://www.salome-platform.org or email : webmaster.salome@opencascade.com
# %% LICENSE_END

import sys

"""
***********************************************
* prerequisites environment variables 
* to use iradinaGUI
***********************************************

#example of set env to use iradinaGUI
bash

export CONDADIR=/volatile/wambeke/miniconda2
export PATH=$CONDADIR/bin:$PATH
source activate pyscientific # as python3 PyQt5
  
# useful environ var
IRADINAGUI_ROOT_DIR=.../iradinaGUI
IRADINA_ROOT_DIR=.../iradina  # where .de code is located


${IRADINAGUI_ROOT_DIR}/iradinaGUI -h    # help
"""



#####################################
_importErrorMessage = r"""
***********************************************
* iradinaGUI configuration is incorrect.
***********************************************

in 'classical iradinaGUI installation' configuration.
user usually have:

  >> bash
  >> IRADINA_ROOT_DIR=...       # where .de code is located
  >> IRADINAGUI_ROOT_DIR=...    # useful environ var
  >> which iradina_linux64.exe  # iradina executable (from C.Borschel)
  >> ${IRADINAGUI_ROOT_DIR}/iradinaGUI

Minimum python configuration is:
  python 3.5
  PyQt5, PyQt5.QtWebEngine
  xml
  numpy
  matplotlib
  pandas

Tricks:
  Easily get local python (as no superuser Linux or Windows) see:
  https://conda.io/miniconda.html
  https://conda.io/docs/install/quick.html
  and:
  >> conda create --name pyscientific \
           python=3.5 pyqt=5 \
           pip jupyter matplotlib numpy pandas pandas-datareader \
           scipy sympy jsonschema pyyaml libxml2

  >> export CONDADIR=.../miniconda3
  >> export PATH=$CONDADIR/bin:$PATH
  >> source activate pyscientific 
"""

#####################################
_rootFeaturesErrorMessage = r"""
***********************************************
* iradina code configuration is incorrect.
***********************************************

in 'classical iradina code compilation/installation' configuration.
user usually have:

  >>> TODO
"""

import subprocess as SP

def error(message):
  print("ERROR: %s" % message)
  return

def TestImports():
  """
  test of prerequisites import python for IradinaGui
  message for problem(s), aborting immediatly. 
  """
  # TODO imports = "PyQt5 PyQt5.QtWebEngineWidgets salomepy xyzpy numpy pandas".split()
  imports = "PyQt5 salomepy xyzpy zlib numpy  pandas".split()

  res = "OK"
  for ii in imports:
    try:
      exec("import %s" % ii)
    except:
      error("problem with 'import %s'" % ii)
      res = "KO"

  if res == "KO":
    print(_importErrorMessage)
    # sys.exit("------> Fix it. Aborting.\n")
    # sys.exit("------> Fix it. Now it's at your owns risks.\n")
  return res

def TestIradinaFeatures():
  """
  test of iradina code features
  TODO
  """
  cmd = "iradina.exe -h"
  res = SP.Popen(cmd, shell=True, stdout=SP.PIPE, stderr=SP.PIPE).communicate()
  if res[1] != "":
    error("problem with '%s': %s" % (cmd, res[1])) #stderr
  stdout = res[0].split()
  prereq = "qt python".split()
  res = "OK"
  for ii in prereq:
    if ii not in stdout:
      error("problem with necessary ROOT compilation feature '%s'" % ii)
      res = "KO"

  if res == "KO":
    print(_rootFeaturesErrorMessage)
    print(_importErrorMessage)
    #print("------> Fix it. Now it's at your owns risks.\n")
    sys.exit("------> Fix it. Aborting.\n")
  return res


# TestImports()
# TestIradinaFeatures()

