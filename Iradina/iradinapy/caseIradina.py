#!/usr/bin/env python
# -*- coding: utf-8 -*-

# %% LICENSE_SALOME_CEA_BEGIN
# Copyright (C) 2008-2015  CEA/DEN
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


import os
import platform
import math
from PyQt5 import QtGui, QtCore, QtWidgets as QTW
import subprocess as SP
import pprint as PP

from xyzpy.baseXyz import _XyzConstrainBase, ListOfBaseXyz

import xyzpy.utilsXyz as UXYZ
import xyzpy.classFactoryXyz as CLFX
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

# set classes Xyz in factory xyzpy.utilsXyz
import xyzpy.intFloatListXyz as IFLX
import mendeleievpy.projectileDrt
import mendeleievpy.targetDrt
import mendeleievpy.elementDrt

import iradinapy.configIra as CFGIRA
import iradinapy.debug as DBG

_Mole = 6.02e23

verbose = False
verboseEvent = False

###############################################################
# common methods
###############################################################

# ...

###############################################################
# classes
###############################################################

###############################################################
class BoolFalseIra(IFLX.BoolXyz):
  """
  ['False', 'True'] with write config file iradina strCfg [0, 1]
  default False
  """
  _defaultValue = "False"
  _toCfg = {'False': 0, 'True': 1}

  def strCfg(self):
    res = self._toCfg[str(self)]
    return res

###############################################################
class BoolTrueIra(BoolFalseIra):
  """
  ['False', 'True'] with write config file iradina strCfg [0, 1]
  default True
  """
  _defaultValue = "True"
  pass

###############################################################
class BeamSpreadIra(IFLX.FloatPosXyz):
  """
  float positive with default 1.
  """
  _defaultValue = 1.
  pass

###############################################################
class IonE0Ira(IFLX.FloatPosXyz):
  """
  float positive with default 10e3 keV
  """
  _defaultValue = 1e3 # keV
  pass

###############################################################
class IonMIra(IFLX.FloatPosXyz):
  """
  Ion mass (g/mol) float positive with default 1.
  """
  _defaultValue = 1.
  pass

###############################################################
class IonAngleIra(IFLX.IntRangeXyz):
  """
  angle incidence xy (degree)
  integer range -90 to 90
  default 0.
  """
  _allowedRange = [-90, 90]
  _defaultValue = 0.

  def toXml(self, **kwargs):
    import mendeleievpy.localMendeleiev as MDLV
    """set tooltip_1 attibute xml for tooltip as long name element (value column 1)"""
    res = super(IonAngleIra, self).toXml(**kwargs)
    alpha = math.radians(self)
    sin = math.sin(alpha)
    cos = math.cos(alpha)
    tooltip_1 =  "ion_vx = %.3f\nion_vy = %.3f" % (cos, sin) # sin,cos as display tooltip
    # avoid accents in xml as not well-formed (invalid token)
    try:
      tooltip_1 += " - " + MDLV.nameToFr[tooltip_1]
    except:
      pass
    res.attrib["tooltip_1"] = MDLV.removeAccents(tooltip_1)
    return res

  def toVx(self):
    alpha = math.radians(self)
    cos = math.cos(alpha)
    return cos

  def toVy(self):
    alpha = math.radians(self)
    sin = math.sin(alpha)
    return sin


###############################################################
class IonVxIra(IFLX.FloatPosXyz):
  """
  vector incidence x float positive with default 1.
  """
  _defaultValue = 1.
  pass

###############################################################
class IonVyIra(IFLX.FloatPosXyz):
  """
  vector incidence y float positive with default 0.
  """
  _defaultValue = 0.
  pass

###############################################################
class IonVzIra(IFLX.FloatPosXyz):
  """
  vector incidence z float positive with default 0
  """
  _defaultValue = 0.
  pass

###############################################################
class MinEnergyIra(IFLX.FloatPosXyz):
  """
  float positive with default 5.
  """
  _defaultValue = 5.
  pass

###############################################################
class MaxNoIonIra(IFLX.IntPosXyz):
  """
  int positive with default 20000
  """
  _defaultValue = 20000
  pass

###############################################################
class StorageIntervalIra(IFLX.IntPosXyz):
  """
  int positive with default 2000
  """
  _defaultValue = 2000
  pass

###############################################################
class Seed1Ira(IFLX.IntPosXyz):
  """
  int positive with default 123
  """
  _defaultValue = 123
  pass

###############################################################
class Seed2Ira(IFLX.IntPosXyz):
  """
  int positive with default 456
  """
  _defaultValue = 456
  pass

###############################################################
class StatusUpdateIntervalIra(IFLX.IntPosXyz):
  """
  int positive with default 10000
  """
  _defaultValue = 10000
  pass

###############################################################
class CellCountxIra(IFLX.IntSupEq1Xyz):
  """
  int positive with default 100
  """
  _defaultValue = 100
  pass

###############################################################
class CellCountyIra(IFLX.IntSupEq1Xyz):
  """
  int positive with default 1
  """
  _defaultValue = 1
  pass

###############################################################
class CellCountzIra(IFLX.IntSupEq1Xyz):
  """
  int positive with default 1
  """
  _defaultValue = 1
  pass

###############################################################
class CellDepthxIra(IFLX.FloatPosXyz):
  """
  float positive with default 1000. (nm)
  """
  _defaultValue = 1000.

  def toXml(self, **kwargs):
    """set tooltip_1 attibute xml for tooltip as long name element (value column 1)"""
    res = super(CellDepthxIra, self).toXml(**kwargs)
    value = self/self.parentAsAttribute().cell_count_x
    tooltip_1 =  "implies cell_size_x = %.3f" % value # as display tooltip
    res.attrib["tooltip_1"] = tooltip_1
    return res

  pass

###############################################################
class CellSizexIra(IFLX.FloatPosXyz):
  """
  float positive with default 10.
  """
  _defaultValue = 10.
  pass

###############################################################
class CellSizeyIra(IFLX.FloatPosXyz):
  """
  float positive with default 100.
  """
  _defaultValue = 100.
  pass

###############################################################
class CellSizezIra(IFLX.FloatPosXyz):
  """
  float positive with default 100.
  """
  _defaultValue = 100.
  pass

###############################################################
class StorePathLimitIra(IFLX.IntSupEq1Xyz):
  """
  int positive with default 50
  """
  _defaultValue = 50
  pass

###############################################################
class DisplayIntervalIra(IFLX.IntSupEq1Xyz):
  """
  int positive with default 20
  """
  _defaultValue = 20
  pass

###############################################################
class StragglingModelIra(IFLX.StrInListXyz):
  _allowedList = [
    "(0) No straggling",
    "(1) Bohr straggling",
    "(2) Chu correction [PRA 13 (1976) 2057]",
    "(3) Chu + Yang correction [NIMB 61 (1991) 149]",
  ]
  _defaultValue = "(3) Chu + Yang correction [NIMB 61 (1991) 149]"

  def strCfg(self):
    """get index from  '(i) blah blah' """
    res = str(self)[1:2]
    return res

###############################################################
class SimulationTypeIra(IFLX.StrInListXyz):
  _allowedList = [
    "(0) Full cascade",
    "(3) Ion profile only",
    "(4) Quick calculation of damage",
  ]
  _defaultValue = "(0) Full cascade"

  def strCfg(self):
    """get index from  '(i) blah blah' """
    res = str(self)[1:2]
    return res

###############################################################
class CompositionFileTypeIra(IFLX.StrInListXyz):
  """
  If the file is just one column of values then FileType should be set to 1.
  If the file contains 4 columns (x, y, z, value) then set it to 0.
  """
  _allowedList = [
    "(0) 4 columns (x, y, z, value)",
    "(1) 1 column of values only",
  ]
  _defaultValue = "(0) 4 columns (x, y, z, value)"

  def strCfg(self):
    """get index from  '(i) blah blah' """
    res = str(self)[1:2]
    return res

###############################################################
class FlightLengthTypeIra(IFLX.StrInListXyz):
  _allowedList = [
    "(0) Poisson distributed flight length and impact pars",
    "(1) Atomic spacing",
    "(2) Constant flightpath in nm",
    "(3) algorithm like-SRIM",
  ]
  _defaultValue = "(0) Poisson distributed flight length and impact pars"
  pass

  def strCfg(self):
    """get index from  '(i) blah blah' """
    res = str(self)[1:2]
    parent = self.parentAsAttribute()
    if parent is None:
      return res
    if parent.simulation_type.strCfg() == "4":
      # logger.debug("FlightLengthTypeIra mandatory 3 with simulation_type 4")
      return "3"
    return res


###############################################################
class IonDistributionIra(IFLX.StrInListXyz):
  _allowedList = [
    "(0) Random ion entry positions",
    "(1) Centered",
    "(2) specified position (see enter_y, enter_z)",
    "(3) random square around position (see enter_y, enter_z, beam_spread)",
  ]
  _defaultValue = "(1) Centered"
  pass

  def strCfg(self):
    """get index from  '(i) blah blah' """
    res = str(self)[1:2]
    return res

###############################################################
class ElementCountIra(IFLX.IntRangeXyz):
  """
  as define MAX_NO_MATERIALS 20       /* Maximum number of different materials */
  """
  _allowedRange = [1, 20]
  pass

###############################################################
class FloatListIra(IFLX.StrXyz):
  """accept [val1, val2, val3]"""
  pass

###############################################################
class MaterialNameIra(IFLX.StrXyz):
  """
  string no more 24 characters
  target.c: if(strlen(MaterialName)>=25){MaterialName[24]='\0';}
  """
  _lenMax = 23 # precaution ...0 to 24 including '\0'
  _defaultValue = "NotSet"
  pass


###############################################################
class DensityIra(IFLX.FloatPosXyz):
  """
  float positive with default 0.
  supposed g/cm3, have to convert to at/cm3
  """
  _defaultValue = 0.

  def getCalculatedValue(self):
    """
    prorata ElementConc(s) and Components densities
    a trivial approximation
    """
    parent = self.parentAsAttribute()
    concs = [float(i.ElementConc) for i in parent.Components]
    densities = [float(i.Density) for i in parent.Components]
    # DBG.write("DensityIra.getCalculatedValue", [concs, densities])
    try:
      concs = self.normalize(concs)
      for i, c in enumerate(concs):
        densities[i] = c * densities[i]
      return sum(densities)
    except:
      logger.warning("problem DensityIra.getCalculatedValue from %s" % [concs, densities])
      return None

  def toAtomCm3(self):
    """prorata ElementConc(s) and Components atomic weights"""
    parent = self.parentAsAttribute()
    concs = [float(i.ElementConc) for i in parent.Components]
    weights = [float(i.Component.atomicWeight) for i in parent.Components] # weight per mole
    # DBG.write("DensityIra.toAtomCm3", [concs, weights])
    try:
      concs = self.normalize(concs)
      for i, c in enumerate(concs):
        weights[i] = c * weights[i]
      avgWeight = sum(weights) # average mole weight of "average atom"
      res = float(self)/avgWeight*_Mole # g/cm3 to at/cm3
      return res
    except:
      logger.warning("problem DensityIra.toAtomCm3 from %s" % [concs, weights])
      return None

  def toXml(self, **kwargs):
    """set tooltip_1 attibute xml for tooltip as long name element (value column 1)"""
    res = super(DensityIra, self).toXml(**kwargs)
    atomCm3 = self.toAtomCm3()
    if atomCm3 is not None:
      tooltip_1 =  "%.2e (at/cm3)" % atomCm3
    else:
      tooltip_1 =  "%s (at/cm3)" % atomCm3
    res.attrib["tooltip_1"] = tooltip_1
    return res

  def getActionsContextMenu(self):
    """append action 'Set defaut value'"""
    actions = super(DensityIra, self).getActionsContextMenu()
    actions.append(self._createAction(
      'Set value (%.3f) (a trivial approximation from components densities)' % self.getCalculatedValue(), None,
      'Set approximate calculated value', self.setCalculatedValue, "setdefaut"))
    return actions

  def setCalculatedValue(self):
    new = self.getCalculatedValue()
    if new is not None:
      self.setValueByControllerSignal(new)
    return

  def normalize(self, aList):
    total = float(sum(aList))
    if total == 0.:
      logger.warning("DensityIra problem total=0 in normalize from %s" % aList)
      return [i for i in aList] # do nothing
    return [i/total for i in aList]

  def strCfg(self):
    """iradina needs density of the material in atoms/cm3"""
    res = "%.3g" % self.toAtomCm3()
    return res


########################################################################################
class DensityTargetComponentIra(IFLX.FloatPosXyz):
  """TODO could set default density from wiki localMendeleiev or user config"""
  _defaultValue = 0.

  def getWikipediaValue(self):
    symbol = self.parentAsAttribute().Component.symbolElement
    import mendeleievpy.localMendeleiev as MDLV
    res = MDLV.getDensity_wiki(symbol)
    return res

  def getUserConfigValue(self):
    symbol = self.parentAsAttribute().Component.symbolElement
    import iradinapy.configIra as CFGIRA
    cfg = CFGIRA.getMainConfig() # supposed global setMainConfig() done
    res = cfg['Densities'][symbol]
    return res

  def getActionsContextMenu(self):
    """append action 'Set defaut value'"""
    actions = super(DensityTargetComponentIra, self).getActionsContextMenu()
    symbol = self.parentAsAttribute().Component.symbolElement
    actions.append(self._createAction(
      'Set wikipedia default value (%s)' % self.getWikipediaValue(), None,
      'from wikipedia default values', self.setDefaultDensityWiki, "setdefaut"))
    actions.append(self._createAction(
      'Set user default value (%s)' %  self.getUserConfigValue(), None,
      'from user default values', self.setDefaultDensityUser, "setdefaut"))
    return actions

  def setDefaultDensityWiki(self):
    new = self.getWikipediaValue()
    self.setValueByControllerSignal(new)
    return

  def setDefaultDensityUser(self):
    new = self.getUserConfigValue()
    self.setValueByControllerSignal(new)
    return

###############################################################
class CaseIra(_XyzConstrainBase):
  """
  general informations about case and for launch iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("IonBeam", "IonBeamIra"),
    ("Target", "TargetIra"),
    ("Simulation", "SimulationIra"),
 ]
  _icon = "caseira"

  _helpDict = {
    "IonBeam": (u"Define incident ion beam parameters", u""),
    "Target": (u"Define target materials and geometry parameters", u""),
    "Simulation": (u"Define general simulation parameters", u""),
  }

  _defaultVersion = "1.0.0"

  def __init__(self):
    super(CaseIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def setDefaultValues(self):
    self.IonBeam.setDefaultValues()
    self.Simulation.setDefaultValues()
    self.Target.setDefaultValues()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGIRA.isHidden(self, nameAttr)


###############################################################
class IonBeamIra(_XyzConstrainBase):
  """
  general informations about beam iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("ion", "IsotopeIra"),
    ("ionE0", "IonE0Ira"),
    ("ion_angle_xy", "IonAngleIra"),
    #("ion_vx", "IonVxIra"),
    #("ion_vy", "IonVyIra"),
    #("ion_vz", "IonVzIra"),
    ("ion_distribution", "IonDistributionIra"),
    ("enter_y", "FloatPosXyz"),
    ("enter_z", "FloatPosXyz"),
    ("beam_spread", "BeamSpreadIra"),
  ]
  _icon = "ionbeamira"

  _ionHelp = u"""\
Ion entry distribution model
0 = Completely random distribution of entry positions of the ion on the x = 0 plane. 
    The other three parameters are ignored.
1 = All ions enter at the center of the target (like in TRIM). 
    The other three parameters are ignored .
2 = All ions enter at the defined (y, z) position in the x = 0 plane.
    The position must be specified by enter_y and enter_z (nm).
3 = Like 2, but the ions are spread randomly around the entry point.
    beam_spread specifies how much the y and z position
    is spread around the defined entry point (nm).
"""

  _helpDict = {
    "ionZ": (u"atomic number of incident ion (number of protons)", u""),
    "ionM": (u"Mass of the ion (g/mol)", u""),
    "ionE0": (u"Impinging energy (keV)\nnote: will be converted as (eV) for iradina code", u""),
    "ion_angle_xy": (u"Ion beam incident angle plane x-y (degree 0 to 90)", u""),
    "ion_vx": (u"Ion beam orientation vector x", u""),
    "ion_vy": (u"Ion beam orientation vector y", u""),
    "ion_vz": (u"Ion beam orientation vector z", u""),
    "ion_distribution": (_ionHelp, u""),
    "enter_y": (u"Entry point y (nm)", u""),
    "enter_z": (u"Entry point z (nm)", u""),
    "beam_spread": (u"""how much the y and z position is spread 
around the defined entry point (nm).
relevant if ion_distribution is 3""", u""),
  }

  def __init__(self):
    super(IonBeamIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = CFGIRA.isHidden(self, nameAttr)
    if res: return res # hidden user mode

    # res is False, other hidden option(s)
    idistrib = int(str(self.ion_distribution)[1:2])
    if nameAttr in ["enter_y", "enter_z"]:
      if idistrib not in [2, 3]: return True
    if nameAttr in ["beam_spread"]:
      if idistrib not in [3]: return True # only relevant for option 3
    return res


###############################################################
class SimulationIra(_XyzConstrainBase):
  """
  general informations about simulation iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("launch_in_background", "BoolFalseIra"),
    ("max_no_ions", "MaxNoIonIra"),
    ("display_interval", "DisplayIntervalIra"),
    ("storage_interval", "StorageIntervalIra"),
    ("status_update_interval", "StatusUpdateIntervalIra"),
    ("store_transmitted_ions", "BoolFalseIra"),
    ("store_exiting_recoils", "BoolFalseIra"),
    ("store_exiting_limit", "IntPosXyz"),
    ("store_energy_deposit", "BoolFalseIra"),
    ("store_ion_paths", "BoolFalseIra"),
    ("store_recoil_cascades", "BoolFalseIra"),
    ("store_path_limit", "StorePathLimitIra"),
    ("simulation_type", "SimulationTypeIra"),
    ("flight_length_type", "FlightLengthTypeIra"),
    ("flight_length_constant", "FloatPosXyz"),
    ("detailed_sputtering", "BoolFalseIra"),
    ("min_energy", "MinEnergyIra"),
    ("seed1", "IntSupEq1Xyz"),
    ("seed2", "IntSupEq1Xyz"),
    # ("OutputFileBaseName", "StrXyz"), # ./output/test
    ("normalize_output", "BoolTrueIra"),
    ("transport_type", "BoolTrueIra"),
    ("multiple_collisions", "BoolFalseIra"),
    ("scattering_calculation", "BoolFalseIra"),
    ("do_not_store_damage", "BoolFalseIra"),
    ("store_range3d", "BoolFalseIra"),
    ("store_info_file", "BoolTrueIra"),
  ]
  _icon = "simulationira"

  _seedHelp = u"""Initial seeds for the pseudo random number generator.
    If you choose the same seeds for the same simulation, 
    the program will always generate exactly the same results.
    If you want to run the program again and add up results from multiple runs
    to obtain better statistics, you must change these seeds"""

  _helpDict = {
    "launch_in_background": (u"Iradina calculus in background", u""),
    "max_no_ions": (u"Number of ions that should be simulated impinging on the target", u""),

    "display_interval": (u"display the progress every n ions (as log)", ""),
    "storage_interval": (u"""store in file intermediate results every n ions,
allows killing the running program and see (at least)
the results up to the last storing point""", ""),

    "status_update_interval": (u"""\
Instruct iradina how often it should generate a status file,
which can be monitored by other programs.
To synchronizing GUI by example, if useful.
(when the -g option is supplied via the command line)""", ""),

    "store_transmitted_ions": (u"""creates a file where it stores
the exit position, direction and energy 
for each ion that leaves the target""", ""),

    "store_exiting_recoils": (u"""creates a file where it stores
the exit positions, direction and energy 
of recoils leaving the simulation volume""", ""),

    "store_exiting_limit": (u"""Since storing above takes up memory and cpu time,
set the max number of leaving recoils to be stored.""", ""),

    "store_energy_deposit": (u"""Instructs iradina to create and store two arrays,
that sum up electronic and nuclear energy loss in each cell.""", ""),

    "store_ion_paths": (u"""Instruct iradina to store in file the exact flying paths of ions.
this options generates amounts of data and slows down the program significantly.""", ""),

    "store_recoil_cascades": (u"""Instruct iradina to store in file the exact recoils.
this options generates amounts of data and slows down the program significantly.""", ""),

    "store_path_limit": (u"""Instruct iradina to stop storing ion paths and recoils after a number of ions.
This way you get a few paths for visualization,
but most of the simulation will be still be fast.""", ""),

    "simulation_type": (u"""\
0 = Full cascade (flight_length_type 0)
    Full simulation including detailed following of all recoils in
    collision cascades and deatiled calculation of damage.
    (Similar to TRIM’s 'Detailed calculation with full damage cascades').
3 = Ion profile distribution only (flight_length_type 0)
    Recoils are not followed, damage is not calculated. 
    Faster. 
    This does only work when transport_type is set 1. 
    Otherwise iradina always considers recoils.
4 = Quick calculation of damage (implies mandatory flight_length_type 3)
    (i.e. Modified NRT(Kinchin-Pease like) flight path is SRIM like
""", u""),

    "flight_length_type": (u"""\
Between two collisions, the ions travels on a straight free flight path 
and only loses energy to electrons. 
this is how this free flight length is calculated.
0 = The flight length is varied randomly according to a Poisson
    distribution with an average flight length according to the
    mean interatomic distance.
1 = The flight length is constant and always corresponds to the
    mean interatomic distance.
2 = The flight length is constant and always exactly the number 
    specified by the flight_length_constant parameter (nm).
3 = Mandatory value when simulation_type is 'Quick calculation of damage'
    (i.e. Modified NRT(Kinchin-Pease like) flight path is SRIM like
""", u""),
    "min_energy": (u"Minimum energy below which all projectiles are stopped", u""),

    "detailed_sputtering": (u"""Calculating sputter yields from MC simulations is not trivial.
Read the documentation please.""", u""),

    "seed1": (_seedHelp, u""),
    "seed2": (_seedHelp, u""),

    "normalize_output": (u"""Instruct iradina to store output results in units of (1/cm3) per (ions/cm2)
This allows the direct calculation of the concentration of implanted ions or 
some defect type for a specific fluence. 
Note: this is a little different from TRIM: 
      iradina calculates ions/cm2 by assuming that the plane in cm2 is perpendicular to the ion beam.
      In TRIM, cm2 corresponds to sample surface.""", u""),

    "transport_type": (u"""There are two different versions implemented.
0 = Much more accurate (especially important for sputter calculations).
1 = Much faster (more similar to the computing flow in Corteo),
    yields similar ion and damage distributions but much worse sputter yields.""", u""),

    "multiple_collisions": (u"""Let projectiles do more than one collision per flightpath (in annular cylinders).
This number defines the number of extra collisions. 
So, 'False' means just the one 'normal' collision.""", u""),

    "scattering_calculation": (u"""Makes iradina use the fast database method for calculation of the scattering angle.
'True' makes iradina use MAGIC.
Note: using MAGIC only works when transport_type is set to '0'.""", u""),

    "do_not_store_damage": (u"""do not store the files describing the various distribution of damage. 
The damage is still calculated, just not stored.""", u""),
    "store_range3d": (u"""Creates a file that lists the final position of 
of each implanted ion, exact coordinates, not cells.
Similar to the Range3D file created by TRIM.""", u""),
    "store_info_file": (u"""Create a file with additional information on the 
simulation each time results are stored.
This information includes version of iradina, time of the simulation etc.""", u""),
    #"": (u"", u""),
  }

  def __init__(self):
    super(SimulationIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    if nameAttr in ["launch_in_background"]:
      return True # unconditionnaly hidden for the moment

    res = CFGIRA.isHidden(self, nameAttr)  # as advanced mode, or not
    if res: return res

    if self.simulation_type.strCfg() == "4":
      if nameAttr in ["flight_length_type", "flight_length_constant"]:
        res = True

    if self.flight_length_type.strCfg() != "2":
      if nameAttr == "flight_length_constant":
        res = True

    return res

########################################################################################
class IsotopeIra(_XyzConstrainBase):
  """
  general informations about IsotopeIra
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files Drt
    ("symbolElement", "StrSymbolIsotopeDrt"),
    ("isotopeNumero", "StrXyz"),
    ("atomicNumber", "AtomicNumberDrt"),
    ("atomicWeight", "FloatPosXyz"),
  ]
  _helpDict = {
    "symbolElement": (u"element symbol", u""),
    "isotopeNumero": (u"isotope numero, number of nucleons (could be unknown for big Z)", u""),
    "atomicNumber": (u"atomic number, number of protons", u""),
    "atomicWeight": (u"isotope, molar mass (g)", u""),
  }
  _icon = "isotopedrt"

  def __init__(self):
    super(IsotopeIra, self).__init__()
    self._defautNameAsRoot = "Isotope"
    self.setIsCast(True)
    # self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    super(IsotopeIra, self).setDefaultValues()
    self.symbolElement = "Fe"

  def __setattr__(self, name, value):
    # print("%s.__setattr__" % self.__class__.__name__)
    super(IsotopeIra, self).__setattr__(name, value)
    if name == "symbolElement":  # only one warning to loop of setattr
      self.on_attributesChange(True)
    return

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = False  # default

    if nameAttr == "displacementEnergy":
      if self.useLibraryDisplacementEnergy != "False": res = True
    # print("%s.isHidden %s" % (nameAttr, res))
    return res

  def on_attributesChange(self, verbose=False):
    if verboseEvent: print("%s.on_attributesChange" % self._className)
    self.checkValues(verbose)
    controller = self.getRoot().getController()
    # aStr = self.toStrXml()
    if controller != None:
      if verboseEvent:
        print("refresh views signal of controller %s views %s" % \
              (controller.objectName(), str([str(v.objectName()) for v in controller.getViews()])))
      controller.refreshViewsSignal.emit()
    else:
      if verboseEvent: print("%s new values: controller None, no refresh views" % self.__class__.__name__)  # ,"\n",aStr

  def checkValues(self, verbose=True):
    import mendeleievpy.elementDrt as EDRT
    try:
      symbol = str(self.symbolElement)
    except:
      symbol = "Fe"
    ele = EDRT.ElementDrt()
    if "_" in symbol:
      s, an, aw = symbol.split("_")  # (non interruptible) atomic affectation like Li_6_5.9634
      ele.setSymbol(s)
      self.atomicNumber = int(ele.atomicNumber)
      self.isotopeNumero = an
      self.atomicWeight = float(aw)
      super(IsotopeIra, self).__setattr__("symbolElement", s)  # inherited to not loop
    else:
      ele.setSymbol(symbol)
      self.atomicNumber = int(ele.atomicNumber)
      self.isotopeNumero = ele.getIsotopeNumero()
      self.atomicWeight = float(ele.atomicWeight)
    if verboseEvent: print("checkValues", symbol, int(ele.atomicNumber), float(ele.atomicWeight))
    return

  def getActionsContextMenu(self):
    """append action 'Append file projectile'"""
    actions = super(IsotopeIra, self).getActionsContextMenu()
    actions.append(self._createAction(
      'Browse element', None, 'Mendeleiv table browser', self.browseElement, "browseelement"))
    return actions

  def browseElement(self):
    self.symbolElement.browseElementsDialog()
    return

########################################################################################
class ConcentrationIra(IFLX.FloatRangeXyz):
  """initial value as 1."""
  _defaultValue = 1.
  # allows 'direct' chemical stoechoimetry inputs (as saccharose C12H22O11)
  _allowedRange = [0., 30.]
  pass

########################################################################################
class TargetComponentIra(_XyzConstrainBase):
  """
  general informations about TargetComponentIra
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files Drt
    ("Component", "IsotopeIra"),
    ("Density", "DensityTargetComponentIra"),
    ("ElementConc", "ConcentrationIra"),
    ("ElementDispEnergy", "FloatPosXyz"),
    ("ElementLattEnergy", "FloatPosXyz"),
    ("ElementSurfEnergy", "FloatPosXyz"),
    #("IonSurfEnergy", "FloatPosXyz"),
  ]

  _Note = """
Note: You can extract these values from TRIM
      but you should know that these may be different for compounds.
"""

  _helpDict = {
    "Component": (u"element component", u""),
    "Density": (u"""\
Volumetric mass density of element (g/cm3)
This value is optional,
only used for trivial approximation of material density 
(from  components densities and concentrations).
you should find a fair documented value by yourself.
""", u"mass per unit volume (rho)"),
    "ElementConc": (u"""Concentration of element in this material (value in [0, 1]).
Could be set as stoechometric number.
Sum of all ElementConc in this material may be > 1, a normalization will be done.""", u""),
    "ElementDispEnergy": (u"displacement energy barrier for element (eV)" + _Note, u""),
    "ElementLattEnergy": (u"lattice binding energy barrier for element (eV)" + _Note, u""),
    "ElementSurfEnergy": (u"surface binding energy barrier for element (eV)" + _Note, u""),
 }
  _icon = "targetcomponentdrt"

  def __init__(self):
    super(TargetComponentIra, self).__init__()
    self._defautNameAsRoot = "TargetComponent"
    self.setIsCast(True)
    # self._setAllAttributesList()
    self.setDefaultValues()

  def setDefaultValues(self):
    super(TargetComponentIra, self).setDefaultValues()
    self.ElementDispEnergy = 25
    self.ElementLattEnergy = 3
    self.ElementSurfEnergy = 3

    # self.symbolElement = "Fe"

  def __setattr__(self, name, value):
    # print("%s.__setattr__" % self.__class__.__name__)
    super(TargetComponentIra, self).__setattr__(name, value)
    #if name == "symbolElement":  # only one warning to loop of setattr
    if name == "Component":  # only one warning to loop of setattr
      self.on_attributesChange(True)
    return

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = CFGIRA.isHidden(self, nameAttr) # mode adanced or simple etc...
    if res: return True # hidden

    res = False  # default
    if nameAttr == "displacementEnergy":
      if self.useLibraryDisplacementEnergy != "False": res = True
    # print("%s.isHidden %s" % (nameAttr, res))
    return res

  def on_attributesChange(self, verbose=False):
    # TODO logger.todo("TODO %s.on_attributesChange Density" % self._className)
    # TODO self.checkValues(verbose)
    controller = self.getRoot().getController()
    # aStr = self.toStrXml()
    if controller != None:
      if verboseEvent:
        print("refresh views signal of controller %s views %s" % \
              (controller.objectName(), str([str(v.objectName()) for v in controller.getViews()])))
      controller.refreshViewsSignal.emit()
    else:
      if verboseEvent: print("%s new values: controller None, no refresh views" % self.__class__.__name__)  # ,"\n",aStr

  def checkValues(self, verbose=True):
    import mendeleievpy.elementDrt as EDRT
    try:
      symbol = str(self.symbolElement)
    except:
      symbol = "Fe"
    ele = EDRT.ElementDrt()
    if "_" in symbol:
      s, an, aw = symbol.split("_")  # (non interruptible) atomic affectation like Li_6_5.9634
      ele.setSymbol(s)
      self.atomicNumber = int(ele.atomicNumber)
      self.isotopeNumero = an
      self.atomicWeight = float(aw)
      super(TargetComponentIra, self).__setattr__("symbolElement", s)  # inherited to not loop
    else:
      ele.setSymbol(symbol)
      self.atomicNumber = int(ele.atomicNumber)
      self.isotopeNumero = ele.getIsotopeNumero()
      self.atomicWeight = float(ele.atomicWeight)
    if verboseEvent: print("checkValues", symbol, int(ele.atomicNumber), float(ele.atomicWeight))
    return

  def getActionsContextMenu(self):
    """append action 'Append file projectile'"""
    actions = super(TargetComponentIra, self).getActionsContextMenu()
    actions.append(self._createAction(
      'Browse element', None, 'Mendeleiv table browser', self.browseElement, "browseelement"))
    return actions

  def browseElement(self):
    self.symbolElement.browseElementsDialog()
    return


###############################################################
class ListOfTargetComponentsIra(ListOfBaseXyz):
  _allowedClasses = [TargetComponentIra]
  _icon = "listtargetcomponentsdrt"

  '''
  def getActionsContextMenu(self):
    """append action 'Append file Component'"""
    actions = super(ListOfTargetComponentsIra, self).getActionsContextMenu()
    actions.append(self._createAction(
      'Append Component', None, 'Browse looking for new component', self.browseDialog, "targetcomponentdrt"))
    return actions

  def browseDialog(self):
    self.addItem(self._allowedClasses[0])
    # f = self[-1].fileComponent
    # return f.browseDialog()
  '''

###############################################################
class TargetIra(_XyzConstrainBase):
  """
  general informations about target iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("straggling_model", "StragglingModelIra"),
    ("Materials", "ListOfMaterialIra"),
    ("Structure", "StructureIra"),
  ]
  _icon = "targetira"

  _helpDict = {
    #"": (u"", u""),
    "straggling_model": (u"""for 'material' version of iradina code
0 = No straggling.
1 = Bohr straggling [Boh48].
2 = With Chu corrections [Chu76].
3 = Chu and Yang’s corrections [YOW91]

Use number 3, it is not perfect but the best implemented.
All models need the same simulation time, as values are precalculated and tabulated anyway.
""", u""),
    "Materials": (u"""
It is strongly recommended
just include the materials you really need for the current simulation
iradina will create big matrices for every combinations of two elements in the target""", u""),
  }

  def __init__(self):
    super(TargetIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGIRA.isHidden(self, nameAttr)

###############################################################
class MaterialIra(_XyzConstrainBase):
  """
  general informations about material of target iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("name", "MaterialNameIra"),
    ("IsVacuum", "BoolFalseIra"),
    # ("Components", "ListOfTargetComponentsIra"),
    # ("ElementCount", "ElementCountIra"),
    ("TargetConcentration", "ConcentrationIra"),
    ("Density", "DensityIra"),
    #("ElementsZ", "FloatListIra"),
    #("ElementsM", "FloatListIra"),
    #("ElementsConc", "FloatListIra"),
    #("ElementsDispEnergy", "FloatListIra"),
    #("ElementsLattEnergy", "FloatListIra"),
    #("ElementsSurfEnergy", "FloatListIra"),
    ("IonSurfEnergy", "FloatPosXyz"),
    ("Components", "ListOfTargetComponentsIra"),
  ]
  _icon = "materialira"

  _helpDict = {
    "name": (u"name of the material ('FeO' for example, without square brackets)", u""),
    "IsVacuum": (u"vacuum material (only one, last material)", u""),
    "TargetConcentration": (u"""\
Volumic concentration of this material in target (value in [0, 1]).
Sum of all TargetConcentration of materials may be > 1, a normalization will be done.
Useful if more than one material in target.""", u""),
    "Density": (u"""Material density (g/cm3)
you should find a fair documented value by yourself.
Note: iradina needs density in atoms/cm3,
      value will be automatically converted""", u""),

    #"ElementsConc": (u"Concentration (normalized sum will be computed as 1)", u""),
    #"ElementsDispEnergy": (u"Displacement energy (eV)", u""),
    #"ElementsLattEnergy": (u"Lattice energy (eV)", u""),
    #"ElementsSurfEnergy": (u"Surface energy (eV)", u""),

    "IonSurfEnergy": (u"""Surface binding energy of the ion.
Note: The ions gains this, when entering the material,
      and looses it when exiting.
      Under oblique angles, the ion is refracted.
      This parameter may be omitted as it is
      usually not of great importance.""", u""),
    #"": (u"", u""),
  }

  def __init__(self):
    super(MaterialIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = CFGIRA.isHidden(self, nameAttr)
    if res: return res # hidden user mode

    # res is False, other hidden option(s)
    if bool(self.IsVacuum):
      if nameAttr not in ["name", "IsVacuum", "TargetConcentration"]: return True
    return res

  def getElementCount(self):
    return len(self.Components)

  def getElementsSymbol(self):
    return [i.Component.symbolElement for i in self.Components]

  def getElementsZ(self):
    return [i.Component.atomicNumber for i in self.Components]

  def getElementsM(self):
    return [i.Component.atomicWeight for i in self.Components]

  def getElementsConc(self):
    res = [float(i.ElementConc) for i in self.Components]
    return self.normalize(res)

  def getElementsDispEnergy(self):
    res = [float(i.ElementDispEnergy) for i in self.Components]
    return res

  def getElementsLattEnergy(self):
    res = [float(i.ElementLattEnergy) for i in self.Components]
    return res

  def getElementsSurfEnergy(self):
    res = [float(i.ElementSurfEnergy) for i in self.Components]
    return res

  def normalize(self, aList):
    total = float(sum(aList))
    return [i/total for i in aList]

###############################################################
class ListOfMaterialIra(ListOfBaseXyz):
  """
  Important note: it is strongly recommended NOT to create one file with all materials that you
  know for all of your simulations! You should just include the materials you really need for the
  current simulation, because iradina will create 2.6 MByte scattering matrices for every possible
  combination of two elements in the target! So the memory usage increases with the square of
  the number of different elements! If your materials contain 92 different elements, iradina needs
  22 GByte of memory in the 4-MSB version or 352 GByte in the 6-MSB version
  """
  _allowedClasses = [MaterialIra]
  _icon = "listmaterialira"

  def getAllAttributesName(self):
    res = [i.getAttributeName() for i in self]
    return res

###############################################################
class TypeCompositionIra(IFLX.StrInListXyz):
  _allowedList = [
    "homogenous one material",
    "heterogenous random multiple materials",
  ]
  _defaultValue = "homogenous one material"

  def strCfg(self):
    """get info"""
    res = "# typeComposition generated as '%s'" + str(self)
    return res

###############################################################
class SeedIra(IFLX.IntPosXyz):
  """
  python random.seed(a=None)
  to initialize internal state of the random number generator
  as None or no argument seeds from current time or
  from an operating system specific randomness source if available
  """
  _defaultValue = 0 # as None
  pass

###############################################################
class StructureIra(_XyzConstrainBase):
  """
  general informations about target iradina
  """
  _attributesList = [  # list, not a dict because sequential order list is used in files
    ("typeComposition", "TypeCompositionIra"),
    ("seedComposition", "SeedIra"),
    ("cell_count_x", "CellCountxIra"),
    ("cell_count_y", "CellCountyIra"),
    ("cell_count_z", "CellCountzIra"),
    ("total_depth_x", "CellDepthxIra"),
    # ("cell_size_x", "CellSizexIra"),
    ("cell_size_y", "CellSizeyIra"),
    ("cell_size_z", "CellSizezIra"),
    ("periodic_boundary_x", "BoolFalseIra"),
    ("periodic_boundary_y", "BoolTrueIra"),
    ("periodic_boundary_z", "BoolTrueIra"),
    ("CompositionFileType", "CompositionFileTypeIra"),
    # ("CompositionFileName", "StrXyz"), # ./Composition.in
    ("UseDensityMultiplicator", "BoolFalseIra"), # 0
    # ("DensityMultiplicatorFileName", "StrXyz"), # ./DensityMultiplicator.in or none
    ("special_geometry", "BoolFalseIra"),
  ]
  _icon = "structureira"

  _helpDict = {
    "cell_count_x": (u"Number of slices in the depth of the material", u""),
    "total_depth_x": (u"""Sample x total depth (nm)"
implies cell_size_x=total_depth_x/cell_count_x
Cell sizes below 1 nm are not recommended""", u""),
    # "cell_size_x": (u"calculated cell_depth_x/cell_count_x (nm)", u""),
    "cell_size_y": (u"""Sample y elementary cell depth (nm)
Cell sizes below 1 nm are not recommended""", u""),
    "cell_size_z": (u"""Sample z elementary cell depth (nm)
Cell sizes below 1 nm are not recommended""", u""),

    "CompositionFileType": (u"""\
0 = the composition file must contain four columns.
    The first three columns are the x–, y–, and z–coordinate of a cell 
    and the fourth column is the material number for this cell.
1 = the composition file must contain one column.
    which is the material number for this cell.
Note: Set it at 0 for security.
""", u""),

    "UseDensityMultiplicator": (u"""\
If various parts of your target consist of the same material,
but have different densities, you do not need to define different materials.
Instead you can activate the UseDensityMultiplicator option by setting it to 'True'.
Then you may specify an additional composition file (see DensityMultiplicator.in), which
contains a multiplicator for the density of each cell in the target.""", u""),

    "special_geometry": (u"""\
If iradina has been compiled with a non-standard target geometry,
this special geometry can be switched off by setting special_geometry to 'False'.
depend of version of iradina ('nanowire' or 'NP-nanoparticle')
False = use standard grid to determine material
True  = the current material and surface are not checked using the rectangular cell grid, 
        but other parameters instead.
""", u""),
    "typeComposition": (u"""\
type of target structure:
homogenous one material (as first item of Materials)
heterogenous random multiple materials
etc. (TODO)
""", u""),
    "seedComposition": (u"""\
used if random is used in Composition.in file generation
to initialize internal state of the random number generator
(0) seeds from current time
(>0) seeds with value (as fixed)
""", u""),
  }

  def __init__(self):
    super(StructureIra, self).__init__()
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    res = CFGIRA.isHidden(self, nameAttr)
    if nameAttr == "CompositionFileType":
      # init and hidden, as always at "(0) 4 columns (x, y, z, value)"
      # for typeComposition choices and methods
      # see 'def toFileCompositionIn' and 'def fn_heterogenous_random_multiple_materials' etc.
      res = True
    if nameAttr == "UseDensityMultiplicator":
      # init and hidden, as always at "(0)"
      res = True
    if nameAttr == "special_geometry":
      # init and hidden, as always at "(0)"
      res = True
    if nameAttr == "seedComposition":
      if "random" in str(self.typeComposition):
        return False
      else:
        return True
    return res


  def getCellSizeX(self):
    res = self.total_depth_x/self.cell_count_x
    return res


# factory pattern using xyzpy.utilsXyz._dictOfXyzClass

CLFX.appendAllXyzClasses([
  BoolFalseIra, BoolTrueIra,
  ConcentrationIra, SeedIra,

  BeamSpreadIra,
  IonE0Ira, IonMIra, IonAngleIra, IonVxIra, IonVyIra, IonVzIra,
  IonDistributionIra,
  MinEnergyIra, MaxNoIonIra, StorageIntervalIra,
  Seed1Ira, Seed2Ira,
  StatusUpdateIntervalIra,
  CellCountxIra, CellCountyIra, CellCountzIra,
  CellDepthxIra, CellSizexIra, CellSizeyIra, CellSizezIra,

  TypeCompositionIra,
  CompositionFileTypeIra,
  SimulationTypeIra, StorePathLimitIra, DisplayIntervalIra,
  FlightLengthTypeIra,

  StructureIra, StragglingModelIra,

  DensityIra, FloatListIra,
  MaterialNameIra, ElementCountIra,

  IsotopeIra, DensityTargetComponentIra,
  TargetComponentIra, ListOfTargetComponentsIra,

  MaterialIra, ListOfMaterialIra,
  IonBeamIra, SimulationIra, TargetIra,
  CaseIra,
])
