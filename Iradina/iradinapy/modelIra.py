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
import glob
import pprint as PP

import xyzpy.loggingXyz as LOG
import traceback
from xyzpy.baseXyz import _XyzConstrainBase, BaseXyz
#from xyzpy.intFloatListXyz import _XyzImmBase #, ListOfBaseXyz
import xyzpy.utilsXyz as UXYZ
from xyzpy.utilsXyz import newElement as ETELE
from xyzpy.guiXyz.dialogXmlXyz import DialogXyz, DialogXmlXyz
import xml.etree.ElementTree as ET
import salomepy.utilsWorkdir as UTW
import xyzpy.classFactoryXyz as CLFX
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()

#set classes Xyz in factory xyzpy.utilsXyz
import xyzpy.intFloatListXyz as IFLX

#set classes Ura in factory xyzpy.utilsXyz
import iradinapy.analysisIra
import iradinapy.controlSimulationIra

import mendeleievpy.projectileDrt
import mendeleievpy.targetDrt
import mendeleievpy.elementDrt

import iradinapy.configIra as CFGIRA
import iradinapy.debug as DBG

verbose = False
verboseEvent = True
debug = False


########################################################################################
class ModelIra(_XyzConstrainBase):
  """
  general instance to group all iradina data files in a directory
  input file iradina release ?
  """
  _attributesList = [ #list, not a dict because sequential order list is used in files
    # just for test & debug
    # ("Isotope", "IsotopeDrt"),
    # ("Projectile", "ProjectileDrt"),
    # ("Element", "ElementDrt"),
    # ("TargetComponent", "TargetComponentDrt"),
    # ("Target", "TargetDrt"),
    ("Controls", "ControlSimulationIra"),
    ("Analysis", "AnalysisIra"),
  ]
  
  _helpDict = {
    "Controls": (u"Known settings for calculus", u""),
    "Analysis": (u"iradina analysis data parameters", u""),
  }
  _icon = "iradinapy.resources.iradinagui"

  
  def __init__(self):
    super(ModelIra, self).__init__()
    self._defautNameAsRoot = "Iradina"
    self.setIsCast(True)
    self._setAllAttributesList()

  def isHidden(self, nameAttr):
    """to know if attribute is currently displayed in treeView and other dialog widget"""
    return CFGIRA.isHidden(self, nameAttr)

  def setDefaultValues(self):
    self.Controls.setDefaultValues()
    self.Analysis.setDefaultValues()

  def setFromFileIra(self, fileName, verbose=False):
    """override inherited method"""
    raise Exception("Cannot set from file iradinaGUI Analysis.xml, use method setFromDirectoryIra(aDirWithIraFiles) instead")

  def toStrIra(self):
    """override inherited method"""
    raise Exception("Cannot write to file iradinaGUI Analysis.xml, use method toDirectoryIra(aDirWithIraFiles) instead")

  def getEtudeWorkdirExpanded(self):
    res = self.Analysis.dataInformations.getEtudeWorkdir()
    res = os.path.realpath(os.path.expandvars(res))
    if "$" in res:
      logger.error("impossible to expand etude directory: '%s'" % res)
      return None
    return res

  def getEtudeWorkdirBrut(self, expanded=True):
    return self.Analysis.dataInformations.getEtudeWorkdirBrut()

  def getHistoryFile(self):
    return self.Analysis.historyFileManager

  def getActionsContextMenu(self):
    actions = super(ModelIra, self).getActionsContextMenu()
    actions.append( self._createAction('change user mode', None, 'change user mode of data display', self.userMode, "usermode") )
    actions.append( self._createAction('user expanding', None, 'custom reset expanding iradina tree', self.userExpand, "expand") )
    return actions

  def userExpand(self):
    """
    filename patterns '*,??'
    warning not for '[' ']' as 'alist[*]'
    no found way to quote meta-character
    https://docs.python.org/2/library/fnmatch.html
    """
    expanded = r"""
Iradina
Analysis
dataInformations
dataManager
macroManager
functions
libraries
macros
userFileManager
""".split()
    #print("userExpand\n%s" % PP.pformat(expanded))
    self .getController().UserExpandSignal.emit(expanded)

  def userMode(self):
    """
    change user mode
    """
    cur = CFGIRA.getCurrentMode()
    mode = CLFX.getXyzClassFromName("UserModeIra")(cur)
    aDialog = DialogXyz(parent=self.getDesktop()) # open in center desktop
    aWidget = mode.createEditor(self.getDesktop())
    aDialog.setUpWidgetLayout([aWidget])
    aDialog.setMinimumSize(200, 100)
    aDialog.setWindowTitle("change user mode")
    aDialog.hideButtons("Reset Help".split())
    aDialog.exec_()
    res = aDialog.choice
    if aDialog.choice == "Cancel": return

    res = aWidget.getValue()
    # DBG.write("change userMode %s - > %s" % (cur, res), " ")
    self.getController().UserModeSignal.emit(res)
    return True

#factory pattern using xyzpy.utilsXyz._dictOfXyzClass
CLFX.appendAllXyzClasses( [ModelIra] )



