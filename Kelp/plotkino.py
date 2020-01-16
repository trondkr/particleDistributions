#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.models.oceandrift import OceanDrift
from kelp.kelpClass import PelagicPlanktonDrift

o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information

#######################
# Preparing readers
#######################

for polygonIndex in [5]:
	startdate='10042010'
	enddate='20042010'

	filename='../results/Kelp_polygon_'+str(polygonIndex)+'_kelp_opendrift_'+str(startdate)+'_to_'+str(enddate)+'_novertical.nc'
	plotfilename='../figures/Kelp_polygon_'+str(polygonIndex)+'_kelp_opendrift_'+str(startdate)+'_to_'+str(enddate)+'_novertical.png'
	plotfilenameColor='../figures/Kelp_polygon_'+str(polygonIndex)+'_kelp_opendrift_'+str(startdate)+'_to_'+str(enddate)+'_novertical_color.png'
	plotfilenameAnime='../figures/Kelp_polygon_'+str(polygonIndex)+'_kelp_opendrift_'+str(startdate)+'_to_'+str(enddate)+'_novertical_color.mp4'
	
	print(filename)
	o.io_import_file(filename)

	o.plot_vertical_distribution()
#	o.plot(linecolor='z',lvmin=-35, lvmax=0,filename=plotfilenameColor)
	o.plot(filename=plotfilename)

	o.animation(filename=plotfilenameAnime)
	 
