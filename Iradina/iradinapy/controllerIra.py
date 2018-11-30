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
import time
import shutil
import pprint as PP
import subprocess as SP

from PyQt5 import QtGui, QtCore, QtWidgets as QTW

import xyzpy.loggingXyz as LOG
from xyzpy.controllerXyz import ControllerXyz
import xyzpy.actionsFactoryXyz as ACFX

# append factory classes
import xyzpy.intFloatListXyz as IFLX  # append factory classes

import iradinapy.modelIra
from iradinapy.modelIra import ModelIra

from salomepy.qTabMultipleTextCentral import QTabMultipleTextCentral
# from qtRootCanvasWidget import QtRootCanvasWidget

import xyzpy.utilsXyz as UXYZ
import salomepy.utilsWorkdir as UTW
import iradinapy.iradinaFilePatterns as UFPA
import iradinapy.iradinaSettings as USET

logger = LOG.getLogger()
verbose = True
verboseEvent = True

_MyDir = os.path.split(os.path.realpath(__file__))[0]

########################################################################################
class QTabMultipleTextCentralIra(QTabMultipleTextCentral):
  # -1 as useless as not created
  # user use inheritage for modify TAB_... order or appearence tabs
  TAB_LOG_CMD = 0
  TAB_OTHERTEXTEDIT = -1
  TAB_FILESYSTEM = 1

  # additional inherited
  # TAB_ROOTCANVAS = 3

  def _init_tabs_inherited(self):
    return
    """
    if self.TAB_ROOTCANVAS != -1:
      wid = QtRootCanvasWidget()
      wid.saveFileExt = ".iradina"
      wid.tabName = "Canvas"  # as foreign clear short label name
      # wid.attName = "ROOTCanvasWidget"  # as immmutable name
      # self.activeTabs[self.TAB_ROOTCANVAS] = wid
    """


########################################################################################
class ControllerIra(ControllerXyz):
  """
  class for manage request and action to/from views and model of iradinaGUI
  as MVC pattern
  one controller for one model and some views
  """

  # http://pyqt.sourceforge.net/Docs/PyQt4/new_style_signals_slots.html
  # http://qt-project.org/doc/qt-5.1/qtwidgets/qaction.html ... Detailed Description

  LaunchIradinaSignal = QtCore.pyqtSignal(object)
  LaunchAllTestsSignal = QtCore.pyqtSignal(object)
  LaunchIradinaCodeTestsSignal = QtCore.pyqtSignal(object)
  LoadIradinaDefaultModelSignal = QtCore.pyqtSignal(object)
  LoadIradinaModelXmlSignal = QtCore.pyqtSignal(object)
  SaveIradinaModelXmlSignal = QtCore.pyqtSignal(object)
  RefreshIradinaModelSignal = QtCore.pyqtSignal(object)
  LoadDataFileInModelSignal = QtCore.pyqtSignal(object)
  UpdateEtudeSignal = QtCore.pyqtSignal(object)
  UpdateEtudeDataSignal = QtCore.pyqtSignal(object)
  UpdateEtudeRootlogonSignal = QtCore.pyqtSignal(object)
  CreateCpackEtudeSignal = QtCore.pyqtSignal(object)
  CreateDocSignal = QtCore.pyqtSignal(object)
  GitCommitEtudeSignal = QtCore.pyqtSignal(object)
  ClearModelSignal = QtCore.pyqtSignal(object)
  ExecPythonCodeSignal = QtCore.pyqtSignal(object)
  ExecRootCodeSignal = QtCore.pyqtSignal(object)
  ExecRootProcessFileSignal = QtCore.pyqtSignal(object)
  ExecSearchMethodInIradinaSignal = QtCore.pyqtSignal(object)
  ExecPrintROOTContextSignal = QtCore.pyqtSignal(object)
  ExecPostTreatmentsSignal = QtCore.pyqtSignal(object)
  UserExpandSignal = QtCore.pyqtSignal(list)
  UserModeSignal = QtCore.pyqtSignal(str)
  IradinaGuiHelpSignal = QtCore.pyqtSignal(object)
  IradinaCodeHelpSignal = QtCore.pyqtSignal(object)

  def __init__(self, *args, **kwargs):
    super(ControllerIra, self).__init__(*args, **kwargs)
    self.initializeDone = False
    self.workDir = None  # will not change
    self.currentDir = None  # could change, subdirectory of workDir (a priori)
    self.centralLogView = None
    self.iradinaNameFileSaveXml = None
    self.iradinaNameDirSaveIra = None
    self.actions = []
    self.treeViews = []
    self.docks = []
    self.toolBars = []
    self.__initialize()
    self._iradinaXmlDefaultName = "iradinaGui.xml"
    # self.setIpcController(ControllerIpcIra())

    # connect
    self.UpdateEtudeSignal.connect(self._updateEtude)
    self.UpdateEtudeDataSignal.connect(self._updateDataEtude)
    #self.UpdateEtudeRootlogonSignal.connect(self._updateRootlogonEtude)
    #self.GitCommitEtudeSignal.connect(self._gitCommitEtude)
    #self.CreateCpackEtudeSignal.connect(self._createCpackEtude)
    #self.CreateDocSignal.connect(self._createDocEtude)
    #self.ExecPythonCodeSignal.connect(self._ExecPythonCode)
    #self.ExecRootCodeSignal.connect(self._ExecRootCode)
    #self.ExecRootProcessFileSignal.connect(self._ExecRootProcessFile)
    #self.ExecSearchMethodInIradinaSignal.connect(self._ExecSearchMethodInIradina)
    # self.ExecPrintROOTContextSignal.connect(self._ExecPrintROOTContext)
    self.ExecPostTreatmentsSignal.connect(self._ExecPostTreatments)
    self.UserExpandSignal.connect(self._UserExpand)
    self.UserModeSignal.connect(self._UserMode)
    self.LoadDataFileInModelSignal.connect(self._loadDataFileInModel)
    #self.modelChangeSignal.connect(self.SetIradinaGUIDIR)

    self.isController = True

  def __initialize(self):
    if self.initializeDone == True:
      logger.error("initialize only once for controller %s" % (self.objectName()))
      return
    self.__initializeWorkdir()
    self.__createActions()
    self.__addToolBars()
    self.__addCentral()
    self.__addDocks()
    if self._desktop:  # salome?
      logger.debug("ControllerIra: desktop have set dock and central widget")
      for dock in self.docks:
        self._desktop.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
        dock.show()
      self._desktop.setCentralWidget(self.centralLogView)
      self._desktop.centralWidget().show()
      for tb in self.toolBars:
        self._desktop.addToolBar(tb)

  def __addCentral(self):
    import salomepy.qMainWindowForLog as QMFL

    self.centralLogView = QMFL.QMainWindowForLog(centralWidget=QTabMultipleTextCentralIra())
    if self._desktop:  # salome? no need dock of QMainWindowForLog
      for dock in self.centralLogView.docks: dock.hide()  # no need
      for tb in self.centralLogView.toolBars: tb.hide()

    logger.step("create central view for log iradina")
    self.setView(self.centralLogView)

    tabs = self.centralLogView.getTabs()
    filetab = self.centralLogView.getTabByName("Explore Dir")
    rootPath = os.path.join("$IRADINAGUI_WORKDIR")
    filetab.setDirRootPath(rootPath, filters=[])
    """
    self.setCentralWidget(central)
    self.centralWidget().resize(self.centralWidget().size())
    self.centralWidget().show
    """

  def __addDocks(self):
    import iradinapy.treeViewIra as TVIRA
    # self.docks = []
    dock = QTW.QDockWidget("IradinaObjects")
    treeView = TVIRA.TreeViewIra()
    self.treeViews.append(treeView)
    logger.step("create treeView for iradina")
    dock.setWidget(treeView)
    dock.setAllowedAreas(QtCore.Qt.LeftDockWidgetArea)
    self.docks.append(dock)
    self.setView(treeView)  # contains treeViews
    """
    for dock in self.docks:
      self.addDockWidget(QtCore.Qt.LeftDockWidgetArea, dock)
      dock.show()
    """

  def __addToolBars(self):
    # self.toolBars = []
    tb = QTW.QToolBar("EditIra")  # self.addToolBar("Edit")
    for action in self.actions:
      tb.addAction(action)
    # act = ACFX.getCommonActionByName("GeneralHelp")
    # if act != None: tb.addAction(act)
    self.toolBars.append(tb)

  def __createActions(self):
    """create actions for self widget AND other widgets through ACFX.addInCommonActions"""
    logger.debug("create actions %s" % (self.objectName()))
    # self.actions = []
    action = ACFX.QActionXyz(name="LoadIradinaDefaultModel", text="Load Default Iradina data")
    ok = action.setAction(slot=self.LoadIradinaDefaultModelAction, signal=self.LoadIradinaDefaultModelSignal,
                          shortcut=None, tooltip=u"New Iradina data", icon="iradinapy.resources.iradinagui")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LoadIradinaModelXml", text="Load Iradina data xml")
    ok = action.setAction(slot=self.LoadIradinaModelXmlAction, signal=self.LoadIradinaModelXmlSignal,
                          shortcut=None, tooltip=u"Load Iradina data from file xml", icon="openxml")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="SaveIradinaModelXml", text="Save Iradina data XmlXyz format")
    ok = action.setAction(slot=self.SaveIradinaModelXmlAction, signal=self.SaveIradinaModelXmlSignal,
                          shortcut=None, tooltip=u"Save Iradina data to file xml", icon="savexml")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    action = ACFX.QActionXyz(name="LaunchIradina", text="Launch Iradina calculus")
    ok = action.setAction(slot=self.LaunchIradinaAction, signal=self.LaunchIradinaSignal,
                          shortcut=None, tooltip=u"Launch Iradina calculus", icon="run")
    if ok:
      ACFX.addInCommonActions(action)
      self.actions.append(action)

    import iradinapy.debug as DBG
    if DBG.isDeveloper():
      action = ACFX.QActionXyz(name="LaunchAllTests", text="Launch tests")
      ok = action.setAction(slot=self.LaunchAllTestsAction, signal=self.LaunchAllTestsSignal,
                            shortcut=None, tooltip=u"Launch all tests", icon="test")
      if ok:
        ACFX.addInCommonActions(action)
        self.actions.append(action)

      action = ACFX.QActionXyz(name="LaunchIradinaCodeTests", text="Launch IradinaCode tests")
      ok = action.setAction(slot=self.LaunchIradinaCodeTestsAction, signal=self.LaunchIradinaCodeTestsSignal,
                            shortcut=None, tooltip=u"Launch all IradinaCode tests", icon="testIradinaCode")
      if ok:
        ACFX.addInCommonActions(action)
        self.actions.append(action)

    action = ACFX.QActionXyz(name="RefreshIradinaModel", text="Refresh views")
    ok = action.setAction(slot=self.RefreshIradinaModelAction, signal=self.RefreshIradinaModelSignal,
                          shortcut=None, tooltip=u"Refresh IradinaObjects tree view", icon="refresh")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="ClearModel", text="Clear Model")
    ok = action.setAction(slot=self.ClearModelAction, signal=self.ClearModelSignal,
                          shortcut=None, tooltip=u"Clear Iradina model", icon="clearModel")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="GuiHelp", text="iradina GUI help")
    ok = action.setAction(slot=self.IradinaGuiHelpAction, signal=self.IradinaGuiHelpSignal,
                          shortcut=None, tooltip=u"iradina GUI help", icon="helpIraGui")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

    action = ACFX.QActionXyz(name="IradinaHelp", text="iradina CODE help")
    ok = action.setAction(slot=self.IradinaCodeHelpAction, signal=self.IradinaCodeHelpSignal,
                          shortcut=None, tooltip=u"iradina CODE help", icon="helpIraCode")
    if ok:
      ACFX.addInCommonActions(action)  # may be not in common, only for controller?
      self.actions.append(action)

  def IradinaGuiHelpAction(self):
    nameBrowser = UXYZ.getBrowser()
    tmp = "$IRADINAGUI_ROOT_DIR/doc/build/html/index.html".split("/")
    tmp = os.path.join(*tmp)
    nameUrlHelp = os.path.expandvars(tmp)
    if os.path.exists(nameUrlHelp):
      cmd = "%s %s &" % (nameBrowser, nameUrlHelp)
      proc = SP.Popen(cmd, shell=True)
    else:
      logger.error("inexisting name Url Help: '%s'" % nameUrlHelp )

  def IradinaCodeHelpAction(self):
    nameBrowser = UXYZ.getBrowser()
    tmp = "$IRADINAGUI_ROOT_DIR/doc/src/iradinaDocuments/20140804_iradina_manual.pdf".split("/")
    tmp = os.path.join(*tmp)
    nameUrlHelp = os.path.expandvars(tmp)
    if os.path.exists(nameUrlHelp):
      cmd = "%s %s &" % (nameBrowser, nameUrlHelp)
      proc = SP.Popen(cmd, shell=True)
    else:
      logger.error("inexisting name Url Help: '%s'" % nameUrlHelp )


  def __initializeWorkdir(self):
    """
    initialize user working directory $IRADINAGUI_WORKDIR if not existing
    """
    nameVar = "IRADINAGUI_WORKDIR"
    workDir = os.getenv(nameVar)
    if workDir == None:
      homeDir = os.getenv("HOME")
      workDir = os.path.join(homeDir, "IRADINAGUI_WORKDIR")
    workDir = os.path.realpath(workDir)
    os.environ[nameVar] = workDir
    if not os.path.exists(workDir):
      os.makedirs(workDir)
    self.workDir = workDir
    self.currentDir = workDir

  def getEtudeWorkdirExpanded(self):
    if self._model == None:
      return None
    else:
      return self._model.getEtudeWorkdirExpanded()

  ############################################################################
  """
  slots:
  all name beginning with uppercase are slot without argument
  they are also postfixed ... xxxxAction
  return status is only to significate method is ok, but nothing about action
  """

  ############################################################################
  def ShowAllViewsAction(self):
    logger.debug("ShowAllViewsAction")
    for view in self.getViews():
      view.show()
    return True

  def ClearModelAction(self):
    logger.debug("ClearModelAction")
    self.clearModel()
    return True

  def RefreshIradinaModelAction(self):
    """refresh Views"""
    logger.debug("RefreshIradinaModelAction")
    # aDataXml =  self._model.toXml() #as a copy
    # self.SetIradinaGUIDIR() #IradinaGUIDIR change?
    for treeView in self.treeViews:
      # treeView.update()
      treeView.refreshModel()
    return True

  def ExpandAllAction(self):
    logger.debug("ExpandAllAction")
    for treeView in self.treeViews:
      treeView.expandAll()
    return True

  def LaunchAllTestsAction(self):
    # cmd = "AllTestLauncher.py 2>&1" #unittests strerr to stdout
    cmd = "cd $IRADINAGUI_ROOT_DIR; pwd ; AllTestLauncher.sh"
    logger.debug("LaunchAllTestsAction '%s'" % cmd)
    self.centralLogView.launchCmdIntoPopen(cmd)
    return True

  def LaunchIradinaCodeTestsAction(self):
    cmd = "AllIradinaCodeTestLauncher.sh"
    logger.debug("LaunchIradinaCodeTestsAction '%s'" % cmd)
    self.centralLogView.launchCmdIntoPopen(cmd)
    return True

  def CreateEtudeWorkdir(self):
    etudeDir = self._model.getEtudeWorkdirExpanded()
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    if os.path.isfile(etudeDir):
      QTW.QMessageBox.warning(self._desktop, "warning",
            "etude directory existing yet as file, fix it:\n'%s'" % etudeDir)
      return False

    if os.path.isdir(etudeDir):
      # QTW.QMessageBox.warning(self._desktop, "warning", "etude directory existing yet:\n'%s'" % etudeDir)
      pass
    else:
      UTW.makeDir(etudeDir)
    # for subDir in "data macros doc".split():
    for subDir in "output".split():
      aDir = os.path.join(etudeDir, subDir)
      if not os.path.isdir(aDir):
        UTW.makeDir(aDir)
        logger.info("create directory '%s'" % aDir)
    QTW.QMessageBox.information(self._desktop, "info", "current etude directory:\n'%s'" % etudeDir)
    #self._createGitignore(etudeDir)
    #self._createREADMEfile(etudeDirBrut)
    #self._createDoxyfile(etudeDirBrut)
    #self._createRootLogon(etudeDirBrut)
    #self._gitInit(etudeDir)
    #self._gitCommit("initial commit", etudeDir)
    return True

  def createFilePattern(self, nameFile, patternName):
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      QTW.QMessageBox.warning(self._desktop, "warning", "etude directory not existing, create it.")
      self.CreateEtudeWorkdir()

    replaces = [("@FILE@", nameFile)]
    existingPatterns = UFPA.getPatternKeys()
    if not patternName in existingPatterns:
      _, ext = os.path.splitext(nameFile)
      patternNameByDefault = "aUserFile%s" % ext
      contents = UFPA.getFilePatterns(patternNameByDefault, replaces=replaces)
    else:
      contents = UFPA.getFilePatterns(patternName, replaces=replaces)
    with open(nameFile, "w") as f:
      f.write(contents)
    return True

  def _insertFileInListOf(self, index, aFile, theListOf):
    """no insert if existing yet"""
    for i in theListOf:
      if str(i) == aFile: return True  # existing yet TODO compare true expandvars name file
    try:
      newItem = theListOf._allowedClasses[0](aFile)  # immutable
    except:
      newItem = theListOf._allowedClasses[0]()  # mutable as data List Of
      newItem.name = aFile
    if index >= 0:
      theListOf.insert(index, newItem)
    else:  # index -1
      theListOf.append(newItem)
    return True

  def _createREADMEfile(self, etudeDirBrut):
    """
    create file README.txt
    """
    name = "README.txt"
    nameFile = os.path.join(etudeDirBrut, name)
    nameFileExp = os.path.expandvars(nameFile)
    res = self.createFilePattern(nameFileExp, name)
    ok = self._insertFileInListOf(0, nameFile, self._model.Analysis.userFileManager)
    self.RefreshIradinaModelSignal.emit(None)  # changed model
    return True

  def _updateEtude(self, toUpdate="all"):
    """
    override files in etude, with message, and synchronize controller model
    """
    return False
    # TODO obsolete ?
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      self.CreateEtudeWorkdir()

    logger.warning("TODO update etude %s create other directories yes or no" % etudeDirBrut)
    return True


    logger.debug("update etude %s" % etudeDirBrut)
    ###### dir functions
    functions = self._model.Analysis.macroManager.functions
    # new functions instance for ultimate replace in model
    arg1 = functions.__class__()
    for i in functions:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.functions = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])

    ###### dir macros
    macros = self._model.Analysis.macroManager.macros
    # new macros instance for ultimate replace in model
    arg1 = macros.__class__()
    for i in macros:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.macros = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])

    ###### dir libraries
    libs = self._model.Analysis.macroManager.libraries
    # new libraries instance for ultimate replace in model
    arg1 = libs.__class__()
    for i in libs:
      newFile = self._smartCopyFileInEtude(i, "macros")
      arg1.append(newFile)
    cmd = ".Analysis.macroManager.libraries = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])
    self.UpdateEtudeRootlogonSignal.emit(None)  # automatic
    return True

  def _loadDataFileInModel(self, aFile):
    return False
    # TODO obsolete ?
    name = aFile
    nameFileExp = os.path.expandvars(name)
    nameFile = os.path.realpath(nameFileExp)
    if not os.path.isfile(nameFile):
      ret = mbox.warning(self._desktop, "error",
                         "inexisting data file:\n%s" % nameFile,
                         mbox.Ok)
      if ret == mbox.Ok: return False  # nothing to do
    theListOf = self._model.Analysis.dataManager
    for datai in theListOf:  # mutable class
      ii = datai.getNameExpanded()  # compare true expandvars name file
      if ii == nameFile:
        self.RefreshIradinaModelSignal.emit(None)  # no changed model, may be changed file contents?
        self.RefreshIradinaModelSignal.emit(None)  # changed contents precaution
        return True  # existing yet
    ok = self._insertFileInListOf(-1, name, theListOf)
    # logger.debug("_loadDataFileInModel:\n  %s %s %s" % (aFile, nameFile, PP.pformat(theListOf)))
    self.RefreshIradinaModelSignal.emit(None)  # changed model
    return True

  def _ExecPostTreatments(self):
    logger.warning("TODO _ExecPostTreatments")
    #from matplotlibpy.matplotlibWindowToolbar import MatplotlibWindowToolbar
    #self.MatplotlibWindow = MatplotlibWindowToolbar()
    #self.MatplotlibWindow.show()

    #from pandaspy.pandasOscarMainWidget import PandasTabWidget, PandasMainWidget
    #from pandaspy.pandasMainWidgetXyz import PandasMainDialogXyz
    import pandaspy.pandasMainWidgetXyz as PDMW
    self.PandasWindow = PDMW.PandasMainDialogXyz("IradinaGUI Plot Viewer")
    self.PandasWindow.show()
    return True

  def _updateDataEtude(self, toUpdate="all"):
    """
    override files in etude/data, with message, and synchronize controller model
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()
    if not os.path.isdir(etudeDirExp):
      self.CreateEtudeWorkdir()

    logger.debug("updateDataEtude %s" % etudeDirBrut)

    ###### data files
    dataManager = self._model.Analysis.dataManager
    # new dataManager instance for ultimate replace in model
    arg1 = dataManager.__class__()
    for i in dataManager:
      if toUpdate in ["all", i.name]:
        newFile = self._smartCopyFileDataInEtude(i, "data")
      else:
        newFile = i.__class__();
        newFile.name = str(i.name)
      arg1.append(newFile)

    cmd = ".Analysis.dataManager = args[1]"
    self.setModelItemValueSignalList.emit([cmd, arg1])
    return True

  def _splitFileName(self, aFile):
    """TODO future small tests avoid bug"""
    f = str(aFile)
    aDir, aName = os.path.split(f)  # as ('$HOME', 'aa.txt')
    aDirExp, aNameExp = os.path.split(os.path.realpath(os.path.expandvars(f)))
    return aDir, aName, aDirExp, aNameExp

  def _copyFile(self, originFile, targetFile):
    """copy expanded files and question if override"""
    logger.debug("copy file:\noriginFile %s\ntargetFile %s\n" % (originFile, targetFile))
    if originFile == targetFile: return False  # nothing to do
    if os.path.isfile(targetFile):
      mbox = QTW.QMessageBox
      ret = mbox.warning(self._desktop, "warning",
                         "override file:\n'%s'?" % targetFile,
                         mbox.Yes | mbox.No)
      if ret == mbox.No: return False  # nothing to do
    # copy file...
    shutil.copyfile(originFile, targetFile)
    # message...
    text = "\ncopyFile:\n  originFile %s\n  targetFile %s\n" % (originFile, targetFile)
    self.centralLogView.insertText(text)
    hist = self._model.getHistoryFile()
    hist.appendHistoryCopyOf(originFile, targetFile)
    return True

  def _smartCopyFileInEtude(self, aFile, targetDir):
    """
    have to be smart to copy and override if necessary
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()

    aDir, aName, aDirExp, aNameExp = self._splitFileName(aFile)

    originFile = os.path.join(aDirExp, aNameExp)
    targetFile = os.path.join(etudeDirExp, targetDir, aNameExp)
    targetFileBrut = os.path.join(etudeDirBrut, targetDir, aNameExp)

    ok = self._copyFile(originFile, targetFile)

    if ok:
      newFile = aFile.__class__(targetFileBrut)
    else:
      newFile = aFile.__class__(str(aFile))  # duplicate reference without copy
    return newFile

  def _smartCopyFileDataInEtude(self, aFile, targetDir):
    """
    have to be smart to copy and override if necessary
    assume dataManager items
    """
    etudeDirBrut = self._model.getEtudeWorkdirBrut()
    etudeDirExp = self._model.getEtudeWorkdirExpanded()

    aDir, aName, aDirExp, aNameExp = self._splitFileName(aFile.name)

    originFile = os.path.join(aDirExp, aNameExp)
    targetFile = os.path.join(etudeDirExp, targetDir, aNameExp)
    targetFileBrut = os.path.join(etudeDirBrut, targetDir, aNameExp)

    ok = self._copyFile(originFile, targetFile)

    if ok:
      newFile = aFile.__class__();
      newFile.name = targetFileBrut
    else:
      newFile = aFile.__class__();
      newFile.name = targetFileBrut  # duplicate reference without copy
    return newFile

  def _launchEtude(self):
    """
    TODO
    """
    #mbox = QTW.QMessageBox
    #ret = mbox.question(self._desktop, "question", " TODO !!! launch background ?", mbox.Yes | mbox.No)
    #logger.warning("ControllerIra._launchEtude() TODO...")
    self._model.Analysis.toFileIra()
    return

  def assertDataDirectory(self, etudeDir):
    dataDir = os.path.realpath(os.path.join(etudeDir, "..", "data"))
    if os.path.isdir(dataDir):
      logger.info("use existing corteo data directory:\n%s" % dataDir)
    else:
      name = os.path.join(*"${IRADINAGUI_ROOT_DIR} iradinaCode data_4bit".split())
      origDir = os.path.realpath(os.path.expandvars(name))
      logger.warning("create user data corteo directory:\n%s ->\n%s" % (origDir, dataDir))
      UTW.copyDir(origDir, dataDir)

  def LaunchIradinaAction(self):
    if self._model is None:
      logger.warning("Inexisting iradina data tree")
      mbox = QTW.QMessageBox
      ret = mbox.question(self._desktop, "question", "Inexisting iradina data tree,\nCreate one ?", mbox.Yes | mbox.No)
      if ret == mbox.Yes:
        # logger.debug("Create iradina data tree")
        self.LoadIradinaDefaultModelAction()
      return False
    etudeDir = self._model.getEtudeWorkdirExpanded()
    logger.info("Launch iradina in %s" % etudeDir)
    self.CreateEtudeWorkdir()
    self._updateEtude()
    self._launchEtude()
    self.assertDataDirectory(etudeDir) # assert etudeDir/../data directory exists (corteo database)

    logger.trace("save iradina model xml")
    res = self.SaveIradinaModelXmlAction(withMessage=False)
    if res == False:
      # no save as abort
      logger.warning("no save xml as abort:\n%s" % etudeDir)
      return True
    # if self.iradinaNameDirSaveIra == None:  # cancel save... cancel run
    #   return True
    if platform.system() == "Windows":
      name = "launchIRADINA.bat"
      nameFile = os.path.join(etudeDir, name)
      strLaunch = self.getStrLaunchBat(etudeDir)
      with open(nameFile, "w") as f:
        f.write(strLaunch)
      # useless .bat windows self.chmodarwx(nameFile)
      cmd = "START /B %s" % nameFile # windows popen accept ONE line
      # cmd = strLaunch
      msg = "Launch iradina bat script:\n%s" % cmd

    else: #"Linux etc."
      name = "launchIRADINA.sh"
      nameFile = os.path.join(etudeDir, name)
      strLaunch = self.getStrLaunchSh(etudeDir)
      with open(nameFile, "w") as f:
        f.write(strLaunch)
      self.chmodarwx(nameFile)
      cmd = strLaunch # linux popen accept multiples lines
      msg = "Launch iradina bash script:\n%s" % cmd

    logger.info(msg)
    mbox = QTW.QMessageBox
    ret = mbox.question(self._desktop, "question", "Iradina data files created\nLaunch Iradina code ?", mbox.Yes | mbox.No)
    if ret == mbox.Yes:
      self.centralLogView.launchCmdIntoPopen(cmd)
    return True


  ##################################################################
  def getStrLaunchBat(self, etudeDir, args="4bit etc."):
    # data Corteo in etudeDir/../data
    # logger.debug("current etude directory %s" % etudeDir)

    code_exe = os.path.join(*"${IRADINAGUI_ROOT_DIR} iradinaCode iradina_mingw64.exe".split())
    code_exe = os.path.realpath(os.path.expandvars(code_exe))
    data_dir = os.path.join(etudeDir, "..", "data")
    data_dir = os.path.realpath(os.path.expandvars(data_dir))

    logger.info("iradina code file %s" % code_exe)
    logger.info("iradina data dir  %s" % data_dir)

    if not os.path.isdir(data_dir):
      msg = "rem inexisting corteo data directory: %s\nrem fix it." % data_dir
      logger.critical(msg)
      return msg

    mater_file = "Materials.in"
    struc_file = "Structure.in"
    compo_file = "Composition.in"

    material_file = os.path.join(etudeDir, mater_file)
    structure_file = os.path.join(etudeDir, struc_file)
    composition_file = os.path.join(etudeDir, compo_file)
    cmd = """\
@echo off
SET studyDir={0}
SET iradinaExe={1}
CD %studyDir%

REM tree /A
REM dir

REM to see iradina progress (if log as file)
REM type {0}/iradina.log

REM get iradina log as file
REM %iradinaExe% -p 9 -data ../data -c ./Configuration.in > ./iradina.log

REM get iradina log as stdout
%iradinaExe% -p 9 -data ../data -c ./Configuration.in
SET iradinaExitCode=%errorlevel%

ECHO END of iradina, exit code is %iradinaExitCode%
EXIT %iradinaExitCode%
""".format(etudeDir, code_exe)
    return cmd

  ##################################################################
  def getStrLaunchSh(self, etudeDir, args="4bit etc."):
    # data Corteo in etudeDir/../data
    # logger.debug("current etude directory %s" % etudeDir)

    code_exe = os.path.join(*"${IRADINAGUI_ROOT_DIR} iradinaCode iradina_linux64.exe".split())
    code_exe = os.path.realpath(os.path.expandvars(code_exe))
    data_dir = os.path.join(etudeDir, "..", "data")
    data_dir = os.path.realpath(os.path.expandvars(data_dir))

    logger.info("iradina code file %s" % code_exe)
    logger.info("iradina data dir  %s" % data_dir)

    if not os.path.isdir(data_dir):
      msg = "# inexisting corteo data directory: %s\n# fix it." % data_dir
      logger.critical(msg)
      return msg

    mater_file = "Materials.in"
    struc_file = "Structure.in"
    compo_file = "Composition.in"

    material_file = os.path.join(etudeDir, mater_file)
    structure_file = os.path.join(etudeDir, struc_file)
    composition_file = os.path.join(etudeDir, compo_file)
    cmd = """\
set -x
studyDir={0}
iradinaExe={1}
cd $studyDir
tree
ls -alt
( $iradinaExe -p 9 -data ../data -c ./Configuration.in  | tee ./iradina.log )&
# to see iradina progress
# tail -f {0}/iradina.log
""".format(etudeDir, code_exe)
    return cmd

  def SaveIradinaModelXmlAction(self, withMessage=True):
    if self._model == None:
      QTW.QMessageBox.warning(self._desktop, "warning", "No Iradina data to save")
      return False
    etudeDir = self._model.getEtudeWorkdirExpanded()
    logger.debug("SaveIradinaModelXmlAction %s" % etudeDir)
    if not os.path.isdir(etudeDir):  # create without question if inexisting
      self.CreateEtudeWorkdir()
    nameFile = os.path.join(etudeDir, self._iradinaXmlDefaultName)
    self._model.toFileXml(nameFile)
    if withMessage:
      QTW.QMessageBox.information(self._desktop, "info", "saved file:\n'%s'" % nameFile)
    return True

  def LoadIradinaDefaultModelAction(self):
    logger.debug("LoadIradinaDefaultModel")
    if True: #try:
      aData = ModelIra()
      aData.setDefaultValues()
      self.clearModel()
      self.setModel(aData)
      self.RefreshIradinaModelSignal.emit(None)
      # self.ExpandAllAction()
      self._model.userExpand()
      return True
    else: #except:
      import traceback
      trace = traceback.format_exc()
      # traceback.print_exc() #better explicit verbose problem
      QTW.QMessageBox.warning(self._desktop, "warning",
                                "Load Iradina Default data: houston! we have a problem \n\n%s" % trace)
      return False

  def LoadIradinaModelXmlAction(self):
    from widgetpy.salomeQFileDialog import SalomeQFileDialog
    aDialog = SalomeQFileDialog(parent=self._desktop)
    aDir = USET.getExpandedVar("_IRADINAGUI_WORKDIR")
    nameFile = aDialog.browseFileDialog('Load Iradina file xml', aDir, "(*.xml *.XML)", [])
    logger.debug("LoadIradinaModelXml %s" % nameFile)
    if nameFile == "": return True  # cancel
    realPath = os.path.realpath(nameFile)
    if not os.path.isfile(realPath):
      QTW.QMessageBox.warning(self._desktop, "warning",
                                "Load Iradina xml data: not a file \n'%s'" % realPath)
      return False
    try:
      aData = UXYZ.fromFileXml(realPath)
      logger.debug("LoadIradinaModelXml: data loaded from file")
      if aData.__class__.__name__ != "ModelIra":
        QTW.QMessageBox.warning(self._desktop, "warning", \
                                  "Load Iradina data: '%s' is not a ModelIra instance from \n'%s'" % \
                                  (aData.__class__.__name__, realPath))
        return False
      self.clearModel()
      self.setModel(aData)
      # logger.debug("LoadIradinaModelXml: set model data done")
      self.RefreshIradinaModelSignal.emit(None)
      # self.ExpandAllAction()
      self._model.userExpand()
      return True
    except Exception as e:
      # traceback.print_exc()
      QTW.QMessageBox.warning(self._desktop, "warning",
                                "Load Iradina xml data: problem loading\n%s" % e)
      return False

  def loadXmlFile(self, aFileOrDir):
    realPath = os.path.expandvars(aFileOrDir)
    if os.path.isdir(realPath):
      realPath = os.path.join(realPath, "iradinaGui.xml")  # default name
    if not os.path.isfile(realPath):
      logger.error("inexisting file: %s" % realPath)
      return

    etudeDirToLoad, _ = os.path.split(realPath)
    USET.setEnvVar("IRADINAGUIDIR", etudeDirToLoad)
    aData = UXYZ.fromFileXml(realPath)

    self.clearModel()
    self.setModel(aData)

    self.checkFromCpack()

    self.RefreshIradinaModelSignal.emit(None)
    self.ExpandAllAction()

  def isPresentInFile(self, aFile, strOrListOfStr):
    """return True if str or one str of [str1, str2,...] is in file"""
    try:
      with open(aFile, "r") as f:
        contents = f.read()
      if type(strOrListOfStr) == str:
        if strOrListOfStr in contents: return True
      elif type(strOrListOfStr) == list:
        for aStr in strOrListOfStr:
          if aStr in contents: return True
    except:
      pass
    return False

  def _UserExpand(self, aList):
    logger.debug("User expand\n%s" % PP.pformat(aList))
    for treeView in self.treeViews:
      treeView.userExpandModel(aList)
    return True

  def _UserMode(self, aStr):
    import iradinapy.configIra as CFGIRA
    logger.debug("User change mode '%s'" % aStr)
    CFGIRA.setCurrentMode(aStr)
    self.RefreshIradinaModelSignal.emit(None)
    return True


def launchFromSalomePyConsole():
  from PyQt4 import QtGui
  from iradinapy.controllerIra import ControllerIra
  desktop = QTW.QMainWindow()
  desktop.resize(1000, 700)
  ctrl = ControllerIra(desktop=desktop)
  desktop.show()
  return desktop


if __name__ == '__main__':
  from salomepy.onceQApplication import OnceQApplication

  app = OnceQApplication()
  aWindow = launchFromSalomePyConsole()
  app.exec_()
