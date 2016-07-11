# -------------------------------------------------------------------------
# Name:        Messages
# Purpose:
#
# Author:      burekpe
#
# Created:     16/05/2016
# Copyright:   (c) burekpe 2016
# -------------------------------------------------------------------------


import xml.dom.minidom
import datetime
import time as xtime
import os

from pcraster import*
from pcraster.framework import *

from globals import *


class CWATMError(Exception):
    """
    the error handling class
    prints out an error
    """
    def __init__(self, msg):
        header = "\n\n ========================== CWATM ERROR =============================\n"
        try:
           self._msg = header + msg +"\n" +  sys.exc_info()[1].message
        except:
           self._msg = header + msg +"\n"
    def __str__(self):
        return self._msg

class CWATMFileError(CWATMError):
    """
    the error handling class
    prints out an error
    """
    def __init__(self, filename,msg=""):
        path,name = os.path.split(filename)
        if os.path.exists(path):
            text1 = "path: "+ path + " exists\nbut filename: "+name+ " does not\n"
            text1 +="file name extension can be .nc4 or .nc\n"
        else:
            text1 = "searching: "+filename
            text1 += "\npath: "+ path + " does not exists\n"

        header = "\n\n ======================== CWATM FILE ERROR ===========================\n"
        self._msg = header + msg + text1

class CWATMWarning(Warning):
    """
    the error handling class
    prints out an error
    """
    def __init__(self, msg):
        header = "\n\n ========================== CWATM Warning =============================\n"
        self._msg = header + msg
    def __str__(self):
        return self._msg

class CWATMRunInfo(Warning):
    """
    prints out an error
    """
    def __init__(self, mode, outputDir, Steps = 1, ensMembers=1, Cores=1):
        header = "\n\n ========================== CWATM Simulation Information and Setting =============================\n"
        msg = "   The simulation output as specified in the settings file can be found in "+str(outputDir)+"\n"
        self._msg = header + msg
    def __str__(self):
        return self._msg
