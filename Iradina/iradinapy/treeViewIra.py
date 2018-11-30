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
import sys
import pprint as PP
from PyQt5 import QtGui, QtCore, QtWidgets as QTW

from xyzpy.guiXyz.treeXmlXyz import TreeXmlXyz, TreeXmlXyzItem
import xyzpy.loggingXyz as LOG

logger = LOG.getLogger()
verbose = False
verboseEvent = False

"""
cosmetic stuff for treeView Iradina
"""


####################################################
class TreeXmlXyzItemIra(TreeXmlXyzItem):

  def setFromXml(self, data):
    """
    assume display role tooltip role of expanded
    basename of file or directory in column 1 instead column 2
    """
    # if verbose: print("TreeXmlXyzItemIra.setFromXml")
    super(TreeXmlXyzItemIra, self).setFromXml(data)

    along = self._along  # set by super(TreeXmlXyzItemIra...
    if "$" in along:  # may be environ variable... expand for tooltip
      filename = os.path.expandvars(along)
    else:
      filename = along

    ###############################################
    # VariablePythonXyz in VariablesInterpreterPythonXyz: display classic
    if data.attrib["typeClass"] == 'VariablePythonXyz':
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, filename)
      return

    ###############################################
    # files/directories for Iradina...
    if os.path.isfile(filename):
      # print("setData %s" % filename)
      super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, filename)
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, along)
      super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, os.path.basename(filename))
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.EditRole, along)  # editrole override displayrole
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.DisplayRole, "")

    if os.path.isdir(filename):
      # print("setData %s" % filename)
      super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, filename)

    ###############################################
    # .dat files Parameters for Iradina
    search = "#COLUMN:"
    nb = len(search)
    aStr = along.strip()
    if search in aStr[0:nb]:  # is .dat Parameters for Iradina
      # print("aStr '%s'" % aStr)
      if True:  # TODO try:
        name, index, title, unit = [i.strip() for i in aStr[nb:].split("|")]
        tooltip = "index: %s\nname: %s\ntitle: %s\nunit: %s" % (index, name, title, unit)
        display = " : ".join([name, index, title, unit])
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, tooltip)
        super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, tooltip)
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, name)
        super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.EditRole, along)  # editrole override displayrole
        super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.DisplayRole, display)
      else:  # except:
        logger.warning("fix problem format in .dat file: '%s'" % along)
        pass

    ###############################################
    # DataIra in ListOf as filename
    try:
      if data.attrib["typeClass"] == 'DataIra':
        try:
          index = "(%s)" % data.attrib["index"]
        except:
          index = ""
        # print("attributes '%s'" % PP.pformat(data.attrib))
        # print "data[0]", data[0].text
        nameIni = str(data[0].text)  # data[0] = name of file is first attribute
        nameReal = os.path.expandvars(nameIni)
        nameBase = os.path.basename(str(data[0].text))
        # nameSimple = os.path.splitext(nameBase)[0]
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, nameBase + index)
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, nameIni)
        super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, nameReal)
    except:
      pass

    '''###############################################
    # DataIra in ListOf as filename
    try:
      if data.attrib["typeClass"] == 'SimulationIra':
        super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, "Analysis")
        #super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, nameIni)
        #super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.ToolTipRole, nameReal)
    except:
      pass'''

    ###############################################
    # ExpressionIra in ListOf as filename
    if True:  # TODO try:
      if data.attrib["typeClass"] == 'ExpressionIra':
        try:
          index = "(%s)" % data.attrib["index"]
        except:
          index = ""
        # print("attributes '%s'" % PP.pformat(data.attrib))
        # print "data", data.text
        nameIni = data.text.strip()
        if nameIni == "": nameIni = "="
        try:
          name, expression = [i.strip() for i in nameIni.split("=")]
          super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.DisplayRole, name)
          super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, nameIni)
          super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.DisplayRole, expression)
          super(TreeXmlXyzItem, self).setData(0, QtCore.Qt.ToolTipRole, nameIni)
          super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.EditRole, nameIni)  # editrole
          super(TreeXmlXyzItem, self).setData(1, QtCore.Qt.DisplayRole, expression)
        except:
          logger.warning("fix problem ExpressionIra '%s'" % data.text)
          pass
    else:  # except:
      pass

    return


####################################################
class TreeViewIra(TreeXmlXyz):
  class COLS:
    labels = ['Name', 'Value', 'Attributes']
    Tag = 0
    Text = 1
    Attributes = 2

  def __init__(self, parent=None):
    super(TreeViewIra, self).__init__(parent)

    self.setHeaderLabels(self.COLS.labels)
    """
    #set in app.setFont
    self.setFont(QtGui.QFont("Monospace", 9))
    font = QtGui.QFont(self.font()) #copy: self.font() is const
    font.setFamily("Monospace")
    font.setPointSize(9)
    self.setFont(font)
    """
    self.setAlternatingRowColors(True)
    pal = self.palette()
    pal.setColor(pal.Base, QtGui.QColor(200, 230, 200))
    pal.setColor(pal.Text, QtGui.QColor(0, 0, 0))
    self.setPalette(pal)
    self._TreeXmlXyzItemClass = TreeXmlXyzItemIra

