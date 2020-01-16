#!/usr/bin/env python

import time, calendar
from netCDF4 import Dataset, date2num, num2date
import numpy as np
import sys
from datetime import datetime, timedelta

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime(2018, 4, 5)
__modified__ = datetime(2018, 12, 2)
__version__ = "1.5"
__status__ = "Development"


class ConfOpenDrift(object):

    def definesvimpath(self, year):
        if year==2015:
            return '/Volumes/DATASETS/SVIM/2015/ocean_avg_2015%02d.nc4'
        if year==2016:
            return '/Volumes/DATASETS/SVIM/2016/ocean_avg_2016%02d.nc4'
        return None

    def definemainforcingpath(self):
        return '/Volumes/DATASETS/NorFjord/' if self.runLocally else '/global/work/trondk/MALANGEN/'

    def definenorkystpath(self):
        return 'http://thredds.met.no/thredds/dodsC/fou-hi/norkyst800m/NorKyst-800m_ZDEPTHS_avg.an.{:02d}{:02d}{:02d}{:02d}.nc'

    def defineshapefilepath(self):
        return '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/Shapefile05112018/kelpExPol_11sept2018.shp' if self.runLocally else '/global/work/trondk/kelp/Shapefile/Shapefile05112018/kelpExPol_11sept2018.shp'
    
    def __init__(self):
        print('\n--------------------------\n')
        print('Started ' + time.ctime(time.time()))

        self.runLocally = True
        self.experiments = [3,4,5]
        self.releaseParticles = 4
        self.lowDepth, self.highDepth = -4, -2  # in negative meters
        self.verticalBehavior = False   

        self.mainForcingDirectory = self.definemainforcingpath()
        self.apattern = 'norfjords_160m_avg.nc4_*'
        self.argument = "%s%s" % (self.mainForcingDirectory, self.apattern)

        self.use_norkyst = False
        self.allnorkyst=[]
        self.use_svim=True

        self.shapefile = self.defineshapefilepath()

        self.interMethod = 'linearNDFast'  # linearND
        self.mapResolution = "i"  # f, c, i, h
        self.norkyst=self.definenorkystpath()

        self.outputFilename=None
        self.animationFilename=None 
        self.plotFilename=None

        if  self.use_svim:
            self.svimfiles2015=[]
            self.svimfiles2016=[]