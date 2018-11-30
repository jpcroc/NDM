#!/usr/bin/env python
#-*- coding:utf-8 -*-

#  Copyright (C) 2010-2018  CEA/DEN
#
#  This library is free software; you can redistribute it and/or
#  modify it under the terms of the GNU Lesser General Public
#  License as published by the Free Software Foundation; either
#  version 2.1 of the License.
#
#  This library is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public
#  License along with this library; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA

"""
This file is the main API for config configparser for iradinaGUI
"""

import os
import platform
import shutil
import fnmatch
import sys

import iradinapy.debug as DBG
import iradinapy.returnCode as RCO
import iradinapy.utilsIra as UTS
import iradinapy.dateTime as DATT

import configparserpy.configParserUtils as CPAU

verbose = False


_configUserStr = """\
# User config file iradinaGUI_user.cfg
# user may modify here configuration at his convenience...
# original configuration values are in config file iradinaGUI_default.cfg
# see: https://wiki.python.org/moin/ConfigParserExamples

# ... at your own risks ...
# WARNING: names of values are in lowercase evaluation

# example...

[MainWindow]
title = Iradina GUI
sizex = 800
sizey = 550

"""

_configDefaultStr = """\
# User config file iradinaGUI_default.cfg
# user find here default configuration values
# see: https://wiki.python.org/moin/ConfigParserExamples

# ALL USER MODIFICATIONS HERE WILL BE REMOVED
# WARNING: names of values are in lowercase evaluation

[General]
# usermode implies what data to show: advanced, simple (for now)
usermode = simple

[MainWindow]
title = Iradina GUI
sizex = 800
sizey = 550

"""

"""
#TODO future

[DefaultEnvironment]
# TODO not used yet
# iradina code root dir (.../iradina)
iradina_root_dir = default
# iradina GUI working dir (.../IRADINAGUI_WORKDIR)
iradinagui_workdir = default
# iradina GUI LOGS dir  (.../IRADINAGUI_WORKDIR/LOGS)
iradinagui_logdir = default

[Logger]
logdir = LOG

"""

_mainConfig = [None] # EZ use global unique main config (mutable)

def getMainConfig():
  return _mainConfig[0] # have to be set in ConfigManager.getMainConfig()

class ConfigManager(object):
  """
  Manages the read/write of config files of iradinaGUI,
  and merges if useful.

  | file iradinaGUI_user.cfg
  | file iradinaGUI_default.cfg
  """

  def __init__(self, runner):
    self.runner = runner
    self.logger = runner.getLogger()
    self.options = runner.getOptions()
    self.fileUser = "iradinaGUI_user.cfg"
    self.fileDefault = "iradinaGUI_default.cfg"
    # better store string as ready to create new config instance(s)
    # and readable for debug
    self.configUserStr = None # done one time
    self.configDefaultStr = None # done one time
    self._configUserStr = _configUserStr
    self._configDefaultStr = _configDefaultStr

  def getWorkdir(self):
    return self.options.workdir

  def getRealPath(self, name):
    res = os.path.join(self.getWorkdir(), name)
    return os.path.realpath(res)

  def checkFileExist(self, filename):
    """filename as name relative to workdir"""
    return os.path.isfile(self.getRealPath(filename))

  def assertUserDefaultFiles(self):
    """
    if inexisting, create config files user and default.
    relative to workdir
    """
    aFile = self.getRealPath(self.fileUser)
    if not os.path.isfile(aFile):
      # inexisting, create it
      open(aFile, "w").write(self._configUserStr)

    # unconditionnaly override it
    aFile = self.getRealPath(self.fileDefault)
    with open(aFile, "w") as f:
      f.write(self._configDefaultStr)
      # add elements densities for user setting
      # warning: names of elements are in lowercase evaluation !shit!
      import mendeleievpy.localMendeleiev as MDLV
      f.write("# volumetric mass density (g/cm3)\n")
      f.write("[Densities]\n")
      for ele in MDLV.getSymbolElementsList():
        # capitalize not in reading config, so write lower to be unambigous
        f.write("%s = %s\n" % (ele.lower(), MDLV.getDensity_wiki(ele)))

      f.write("\n[AtomicWeight]\n")
      for ele in MDLV.getSymbolElementsList():
        # capitalize not in reading config, so write lower to be unambigous
        f.write("%s = %s\n" % (ele.lower(), MDLV.getAtomicWeight_wiki(ele)))

  def _getConfig(self, fileName):
    """create config from a file name, reading file, new instance config"""
    cfg = CPAU.getConfigFromStr(open(aFile).read())
    return cfg

  def getUserConfig(self):
    """new instance user config"""
    if self.configUserStr == None:
      aFile = self.getRealPath(self.fileUser)
      self.configUserStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromStr(self.configUserStr)
    return RCO.ReturnCode("OK", "get user config", cfg)

  def getDefaultConfig(self):
    """new instance default config"""
    if self.configDefaultStr == None:
      aFile = self.getRealPath(fileDefault)
      self.configDefaultStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromStr(self.configDefaultStr)
    return RCO.ReturnCode("OK", "get default config", cfg)

  def getMainConfig(self):
    """
    main as merged config of default plus overrides of user
    new instance main config
    """
    self.assertUserDefaultFiles()
    aFile = self.getRealPath(self.fileUser)
    self.configDefaultStr = open(aFile, "r").read()
    aFile = self.getRealPath(self.fileDefault)
    self.configUserStr = open(aFile, "r").read()
    cfg = CPAU.getConfigFromDefaultAndUserStr(self.configDefaultStr, self.configUserStr)
    DBG.write("getMainConfig", cfg.toDict())
    return RCO.ReturnCode("OK", "get main config", cfg)

  def setMainConfig(self, cfg):
    """set a user main config as global, is user choice"""
    _mainConfig[0] = cfg


#########################################################################
# user mode (advanced or else) to show or hide some data in treeview
#########################################################################

# list of hidden items by association with value of item.getTreePyName()

_advanced = [] # [] as nothing hidden

_simple = """
Controls
dataInformations.release
ion_distribution
enter_y
enter_z
beam_spread
straggling_model
cell_count_y
cell_count_z
cell_size_y
cell_size_z
periodic_boundary_x
periodic_boundary_y
periodic_boundary_z
CompositionFileType
UseDensityMultiplicator
special_geometry
display_interval
storage_interval
status_update_interval
store_transmitted_ions
store_exiting_recoils
store_exiting_limit
#store_energy_deposit
store_ion_paths
store_recoil_cascades
store_path_limit
flight_length_type
flight_length_constant
detailed_sputtering
min_energy
seed1
seed2
normalize_output
transport_type
multiple_collisions 
scattering_calculation
do_not_store_damage
store_range3d
store_info_file 
pythonManager
IonSurfEnergy
Components*.Density
Materials*.TargetConcentration
""".split()
# fnmatch pattern [seq] matches any character in seq

_modesNames = "simple advanced".split()

_modes = {
  "simple": _simple,
  "advanced": _advanced,
}

_currentMode = ["simple"] # list as global mutable
if verbose: print("set _currentMode", _currentMode[0])

"""
def getUserName():
  res = os.getenv('USERNAME')
  if res == None:
    res = os.getenv('USER')
  if res is None:
    raise Exception("can't get user name in env var 'USER' or 'USERNAME'")
  return res

# set as default if not wambeke or developper
if getUserName() in "wambeke christian".split():
  _currentMode = ["advanced"] # list as global mutable
else:
  _currentMode = ["simple"] # list as global mutable
"""

def getExistingModes():
  # DBG.write("congigIra.isHidden modes", _modes, DBG.isDeveloper())
  return _modesNames

def setCurrentMode(modeName):
  if not modeName in getExistingModes():
    raise Exception("unknown mode '%s'" % modeName)
  _currentMode[0] = modeName

def getCurrentMode():
  res = _currentMode[0]
  if verbose: print("getCurrentMode", res)
  return res

def isHidden(item, nameAttr=None, modeName=None):
  """
  avoid Components[*].Density because fnmatch pattern [seq] matches any character in seq.
  use Components*.Density instead
  """
  if modeName is None:
    mode = _currentMode[0]
  else:
    mode = modeName
  name = item.getTreePyName()
  if nameAttr is not None:
    name += "." + nameAttr
  hidden = _modes[mode]
  res = False
  for i in hidden:
    if fnmatch.fnmatch(name, "*" + i):
      # if verbose: print("fnmatch hidden ", name, "*" + i)
      res = True
      break
    # else:
    #  if "Density" in name:print("not fnmatch hidden ", name, "*" + i)
  # DBG.write("congigIra.isHidden '%s' %s" % (name, res), "")
  return res

