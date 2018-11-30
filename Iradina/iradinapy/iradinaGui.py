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
This file is the main API file for iradinaGUI

| Warning: NO '__main__ ' call allowed, 
|          Use '../iradinaGUI' (in parent directory)
| 
| Usage: see file ../iradinaGUI
"""

import sys

_KOSYS = 1 # avoid import src

# Compatibility python 2/3 for input function
# input stays input for python 3 and input = raw_input for python 2
try: 
    input = raw_input
except NameError: 
    pass
  
########################################################################
# NO __main__ entry allowed, use '../iradinaGUI' (in parent directory)
########################################################################
if __name__ == "__main__":
    msg = """
ERROR: 'iradinaGui.py' is not main command entry (CLI) for iradinaGUI.
       Use '../iradinaGUI' instead.\n\n"""
    sys.stderr.write(msg)
    sys.exit(_KOSYS)

import os
import sys
import pprint as PP
import argparse as AP
import gettext

import iradinapy # for __version__
import iradinapy.debug as DBG # Easy print stderr (for DEBUG only)
import iradinapy.returnCode as RCO # Easy (ok/ko, why) return methods code
import iradinapy.utilsIra as UTS

# get path to src
rootdir = os.path.realpath( os.path.join(os.path.dirname(__file__), "..") )
DBG.write("iradinaGUI rootdir", rootdir)
srcdir = os.path.join(rootdir, "iradinapy")
cmdsdir = os.path.join(rootdir, "commands")

'''
if DBG.isDeveloper():
  workdirdefault = os.path.realpath(os.path.join(rootdir, "..", "IRADINAGUI_WORKDIR"))
else:
  workdirdefault = os.path.expandvars(os.path.join("$HOME", "IRADINAGUI_WORKDIR"))
'''

# load resources for internationalization
gettext.install("iradinaGUI", os.path.join(srcdir, "i18n"))

_LANG = os.getenv("LANG") # original locale


########################################################################
# utility methods
########################################################################

def launchIra(command):
  """
  launch iradinaGUI as subprocess.Popen
  command as string ('iradinaGUI --help' for example)
  used for unittest, or else...

  :return: RCO.ReturnCode
  """
  if "iradinaGUI" not in command.split()[0]:
    raise Exception(_("Not a valid command for launchIra: '%s'") % command)
  env = dict(os.environ) # copy
  # theorically useless, in user environ $PATH,
  # on ne sait jamais
  # https://docs.python.org/2/library/os.html
  # On some platforms, including FreeBSD and Mac OS X,
  # setting environ may cause memory leaks.
  # see test/initializeTest.py
  if rootdir not in env["PATH"].split(":"):
    env["PATH"] = rootdir + ":" + env["PATH"]
  # TODO setLocale not 'fr' on subprocesses, why not?
  # env["LANG"] == ''
  res = UTS.Popen(command, env=env, logger=None) # no logger! unittest!
  return res
    
def getVersion():
  """get version number as string"""
  return iradinapy.__version__
 
def assumeAsList(strOrList):
  """return a list as sys.argv if string"""
  if type(strOrList) is list:
    return list(strOrList) # copy
  else:
    res = strOrList.split(" ")
    return [r for r in res if r != ""] # supposed string to split for convenience


########################################################################
# parser arguments iradinaGUI class
########################################################################
class ArgumentParserNoExit(AP.ArgumentParser):
  """change Exiting method as no exit"""
  def exit(self, status=0, message=None):
    if message:
      self._print_message(message, AP._sys.stderr)
    return

  def getLogLevels(self):
    import iradinapy.loggingIra as LOG
    return LOG._knownLevels

  def getLogLevelsStr(self):
    import iradinapy.loggingIra as LOG
    return LOG._knownLevelsStr

  #################################################
  # filter methods
  #################################################

  def filter_logLevel(self, aStr):
    import iradinapy.loggingIra as LOG
    DBG.write("filter_logLevel(self, string) %s" % (aStr), True)
    value = LOG.filterLevel(aStr)
    DBG.write("filter_logLevel(self, string) %s -> %s" % (aStr, value), True)
    return value

  def filter_list(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    try:
      exec("value = [%s]" % string)
    except:
      msg = "%r is not a list" % string
      raise Exception(msg)
    return value

  def filter_list_float(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    try:
      exec("value = [%s]" % string)
      value = [float(v) for v in value]
    except:
      msg = "%r is not a list of float" % string
      raise Exception(msg)
    return value

  def filter_range(self, string):
    """
    parser filter from string 'vmin,vmax'
    returns list (if not error with python exec(value=[xx,yy]))
    """
    try:
      exec("value = [%s]" % string)
      value = [float(v) for v in value]
    except:
      msg = "%r is not a list of float" % string
      raise Exception(msg)
    if len(value) != 2:
      msg = "%r is not a range of 2 float" % string
      raise Exception(msg)
    return value

  def filter_list_int(self, string):
    """
    parser filter from string 'xx,yy,zz,...'
    returns list (if not error with python exec(value=[xx,yy,zz,...]))
    """
    try:
      exec("value = [%s]" % string)
      value = [int(v) for v in value]
    except:
      msg = "%r is not a list of int" % string
      raise Exception(msg)
    return value

  def filter_square(self, string):
    value = int(string)
    sqrtv = value**(.5)
    if sqrtv != int(sqrtv):
      msg = "%r is not a perfect square" % string
      raise Exception(msg)
    return value

  def filter_int_positive(self, string):
    try:
      value = int(string)
    except:
      msg = "%r is not a correct positive integer number" % string
      raise Exception(msg)
    if value < 0:
      msg = "%r is not a positive integer number" % string
      raise Exception(msg)
    return value

  def filter_float_positive(self, string):
    try:
      value = float(string)
    except:
      msg = "%r is not a correct positive float number" % string
      raise Exception(msg)
    if value < 0:
      msg = "%r is not a positive float number" % string
      raise Exception(msg)
    return value

  def filter_existing_file(self, string):
    try:
      ok = os.path.isfile(string)
    except:
      msg = "%r is not a correct existing file name" % string
      raise Exception(msg)
    if not ok:
      msg = "%r is not an existing file" % string
      raise Exception(msg)
    return value

  def filter_workdir(self, string):
    try:
      wdir = os.path.realpath(os.path.expandvars(string))
      # DBG.write("arg workdir", (wdir, rootdir), True)
      ok = os.path.isdir(string)
      value = wdir
    except:
      msg = "%r is not a correct existing directory" % string
      raise Exception(msg)
    if rootdir in wdir:
      msg = "working directory: %r\nhave to be outside directory: %r" % (wdir, rootdir)
      raise Exception(msg)
    # not existing directory: have to be created by user, or automatic.
    # if not ok:
    #  msg = "%r is not an existing directory" % wdir
    #  raise Exception(msg)
    return value


########################################################################
# Ira class
########################################################################
class Ira(object):
  """
  The main class that stores all the commands of iradinaGui
  (usually known as 'runner' argument in Command classes)
  """
  def __init__(self, logger):
    """Initialization

    :param logger: (Logger) The logger to use
    """

    # Read the iradinaGui
    # iradinaGUI <options>
    # (the list of possible options is at the beginning of this file)

    self.logger = logger # the logger that will be use
    self.arguments = None # args are postfixes options
    self.options = None # the main options passed to iradinaGui
    self.config = None # the config values from file workdir/iradinaGUI_*.cfg

    self.parser = self._getParser()
    self._confirmMode = True

  def __repr__(self):
    aDict = {
      "arguments": self.arguments,
      "options": self.options,
    }
    tmp = PP.pformat(aDict)
    res = "Ira(\n %s\n)\n" % tmp[1:-1]
    return res

  def getId(self):
    """assimiled as integer incremented on _idCommandHandlers"""
    return -1 # assimiled as root of main first command

  def getLogger(self):
    if self.logger is None: # could use owner Ira instance logger
      import iradinapy.loggingIra as LOG
      self.logger = LOG.getDefaultLogger()
      self.logger.critical("Ira logger not set, unexpected situation, fixed as default")
      return self.logger
    else:                   # could use local logger
      return self.logger

  def getOptions(self):
    return self.options

  def getConfig(self):
    return self.config

  def assumeAsList(self, strOrList):
    # DBG.write("Ira assumeAsList", strOrList, True)
    return assumeAsList(strOrList)

  def _getParser(self):
    """
    Define all possible <options> for iradinaGui CLI: 'iradinaGUI <options>'
    (internal use only)
    """
    workdirdefault = os.getenv("IRADINAGUI_WORKDIR")
    if workdirdefault is None:
      workdirdefault = ""
    workdirdefault = os.path.realpath(os.path.expandvars(workdirdefault))

    # parser = AP.ArgumentParser(description='launch iradinaGUI', argument_default=None)
    parser = ArgumentParserNoExit(description='Graphical user interface to launch iradina code', argument_default=None, add_help=False)

    parser.add_argument(
      '-h', '--help',
      help='show this help message (and do no exit)',
      default=False,
      action='store_true'
    )
    parser.add_argument(
      '-d', '--doc',
      help='show html documentation',
      default=False,
      action='store_true'
    )
    parser.add_argument(
      '-g', '--gui',
      help='show GUI',
      default=False,
      action='store_true'
    )
    parser.add_argument(
      '-v', '--verbose',
      help='set log level verbosity: %s default=%s' % (parser.getLogLevelsStr(), "INFO"),
      type=parser.filter_logLevel,
      default= "INFO",
      metavar='logLevel'
    )
    parser.add_argument(
      '-w', '--workdir',
      help="current working directory default=%s" % workdirdefault,
      type=parser.filter_workdir,
      default=workdirdefault,
      metavar='dirName'
    )
    """# other examples
    parser.add_argument(
      '-p', '--pattern',
      help="file pattern for unittest files ['test_*.py'|'*Test.py'...]",
      default="test_???_*.py",  # as alphabetical ordered test site
      metavar='filePattern'
    )
    parser.add_argument(
      '-t', '--type',
      help="type of output: ['std'(standart ascii)|'xml'|'html']",
      default="std",
      choices=['std', 'xml', 'html'],
      metavar='outputType'
    )"""

    # DBG.write('dir(parser)', dir(parser))
    return parser

  def parseArguments(self, arguments):
    args = self.assumeAsList(arguments)
    options = self.parser.parse_args(args)
    DBG.write("IradinaGUI options", options)
    return options

  def execute_cli(self, cli_arguments):
    """select first argument as a command in directory 'commands', and launch on arguments

    :param cli_arguments: (str or list) The iradinaGUI CLI arguments (as sys.argv)
    """
    logger = self.getLogger() # shortcut

    args = self.assumeAsList(cli_arguments)
    # no arguments do print help and returns
    if len(args) == 0:
      self.print_help()
      return RCO.ReturnCode("OK", "No arguments, as '--help'")

    self.options = self.parseArguments(args)
    options = self.getOptions() # shortcut

    if options.doc:
      self.show_doc()
      return RCO.ReturnCode("OK", "exit if arguments contains '--doc'")

    if not os.path.isdir(options.workdir):
      msg = "inexisting working directory: %s\ncreate it." % options.workdir
      return RCO.ReturnCode("KO", msg)

    # set main handler level
    logger.setLevelMainHandler(self.options.verbose)
    logger.trace("arguments options:\n%s" % PP.pformat(self.options.__dict__))

    # if the help option has been called, print command help
    if self.options.help:
      self.print_help()
      # and not continue
      return RCO.ReturnCode("OK", "exit if arguments contains '--help'")

    import iradinapy.configIra as CFGIRA
    cfgMgr = CFGIRA.ConfigManager(self)
    # as main config, current config from option.workdir/iradinaGUI_*.cfg
    rc = cfgMgr.getMainConfig() # current config from options.workdir/xxx.cfg
    if not rc.isOk():
      return rc
    self.config = rc.getValue()
    cfgMgr.setMainConfig(self.config) # set as global with CFGIRA.getMainConfig()

    # cvw TODO logger.debug("Configuration values:\n%s", self.config.writeToStr())

    # cfg = self.config.toCatchAll()
    # TODO [DefaultEnvironment]
    # DBG.write("DefaultEnvironment", self.config.writeToStr(), True)
    # TODO set environ var if not set? not yet
    # if cfg.DefaultEnvironment.iradina_root_dir == "default":
    # etc.

    # logger.setFileHandlerForCommand(self, cmdInstance)

    # Run the main command using options and config
    returnCode = RCO.ReturnCode("OK", "seems finished on (%s)" % self.options)

    if self.options.gui:
      msg = "BEGIN GUI on options:\n%s" % self.options
      logger.step(msg)
      returnCode = self.runGUI()

      if type(returnCode) != RCO.ReturnCode:
        logger.critical(msg)
        raise Exception("ReturnCode GUI type unexpected %s" % type(returnCode))

      msg = "END GUI on options\n%s\n%s" % (self.options, str(returnCode))
      logger.step(msg)

      """# no xml log file, (yet)
      import iradinapy.linksXml as LKXML
      LKXML.setAttribLinkForCommand(cmdInstance, "full_launched_cmd", strArgs)
      LKXML.setAttribLinkForCommand(cmdInstance, "cmd_res", returnCode.toXmlPassed())
      """

    # close logger/main handler have to be caller stuff...
    return returnCode

  def get_help(self):
    """get general help colored string"""
    msg = self.getColoredVersion() + "\n\n"
    msg += self.parser.format_help()
    return msg

  def print_help(self):
    """prints iradinaGui general help"""
    self.logger.info(self.get_help())

  def show_doc(self):
    """show iradinaGui general documentation"""
    import xyzpy.utilsXyz as UXYZ
    nameBrowser = UXYZ.getBrowser()  #
    fileIndex = os.path.join(rootdir, "doc", "build", "html", "index.html")
    cmd = "%s %s &" % (nameBrowser, fileIndex)
    self.logger.warning("active show doc command:\n%s" % cmd)
    UTS.Popen(cmd)

  def getColoredVersion(self):
    """get colored iradinaGui version message"""
    version = getVersion()
    msg = "<header>iradinaGUI version:<reset> " + version
    return msg

  def getConfirmMode(self):
    return self._confirmMode

  def setConfirmMode(self, value):
    self._confirmMode = value

  def getBatchMode(self):
    return self.options.batch

  def getAnswer(self, msg):
    """
    question and user answer (in console)
    if confirm mode and not batch mode.

    | returns 'YES' or 'NO' if confirm mode and not batch mode
    | returns 'YES' if batch mode
    """
    logger = self.getLogger()
    if self.getConfirmMode() and not self.getBatchMode():
      logger.info("\n" + msg)
      rep = eval(input(_("Do you want to continue? [yes/no] ")))
      if rep.upper() == _("YES"):
        return "YES"
      else:
        return "NO"
    else:
      logger.debug(msg) # silent message
      logger.debug("<green>YES<reset> (as automatic answer)")
      return "YES"

  '''useless
  def runGUI_example(self):
    """a main window example"""
    import salomepy.onceQApplication as OQA
    import ribbonpy.ribbonClassFactory as RCF
    import ribbonpy.ribbonQMainWindow as RQM

    logger = self.getLogger()
    config = self.getConfig()
    logger.trace("runGUI config:\n\n%s" % config.writeToStr())
    cfg = config.toCatchAll() # EZ use

    app = OQA.OnceQApplication()
    aJsonValue = RCF.getExampleJsonRibbon()
    fen = RQM.QMainWindowForRibbon(setFromJson=aJsonValue)
    fen.setWindowTitle(cfg.MainWindow.title)
    fen.resize(int(cfg.MainWindow.sizex), int(cfg.MainWindow.sizey))
    msg = "BEGIN run GUI"
    logger.step(msg)
    fen.show()
    res = app.exec_()
    msg = "END run GUI, return '%s'" % res
    logger.step(msg)
    return RCO.ReturnCode("OK", msg, res)
  '''

  def runGUI(self):
    """a main window for iradinaGUI"""
    import salomepy.onceQApplication as OQA
    app = OQA.OnceQApplication()

    import ribbonpy.ribbonClassFactory as RCF
    import iradinapy.mainWindowIra as MWIRA
    import iradinapy.controllerIra as CTIRA
    import iradinapy.treeViewIra as TVIRA
    import iradinapy.configIra as CFGIRA

    logger = self.getLogger()
    config = self.getConfig()
    # cvw TODO logger.trace("runGUI config:\n\n%s" % config.writeToStr())
    logger.trace("runGUI create window")
    cfg = config.toCatchAll() # EZ access

    # settings from user config 'cfg'
    usermode = cfg.General.usermode
    CFGIRA.setCurrentMode(usermode) # global

    aJsonValue = None #TODO for future: a ribbon RCF.getExampleJsonRibbon()
    desktop = MWIRA.QMainWindowIra(setFromJson=aJsonValue, logger=logger)
    desktop.setWindowTitle(cfg.MainWindow.title)
    sx = int(cfg.MainWindow.sizex)
    sy = int(cfg.MainWindow.sizey)
    logger.debug("iradina main window size %sx%s" % (sx, sy))
    desktop.resize(sx, sy)
    controller = CTIRA.ControllerIra(desktop=desktop)

    msg = "BEGIN run GUI"
    logger.step(msg)
    desktop.show()
    res = app.exec_()
    msg = "END run GUI, return '%s'" % res
    logger.step(msg)
    return RCO.ReturnCode("OK", msg, res)

