#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from kelp.kelpClass import PelagicPlanktonDrift
from opendrift.readers import reader_netCDF_CF_generic
import logging
import gdal
import os
from netCDF4 import Dataset, datetime, date2num,num2date
from numpy.random import RandomState
import random
import glob
import matplotlib.pyplot as plt
try:
    import ogr
    import osr
except Exception as e:
    print(e)
    raise ValueError('OGR library is needed to read shapefiles.')

def setupSeed(hoursBetweenTimestepInROMSFiles,startTime,endTime,startReleaseTime,endReleaseTime,releaseParticles):
    ##################################################
    # Create seed variation as function of day
    ##################################################

    # Make datetime array from start to end at 3 hour interval
    #interval = timedelta(hours=hoursBetweenTimestepInROMSFiles)
    difference=endTime-startTime
    hoursOfSimulation=divmod(difference.total_seconds(), 3600)
     
    difference=endReleaseTime-startReleaseTime
    hoursOfRelease=divmod(difference.total_seconds(), 3600)
     
    #startSimulationJD=startTime.timetuple().tm_yday
    #endSimulationJD=endTime.timetuple().tm_yday
    timeStepsSimulation=int(int(hoursOfSimulation[0])/hoursBetweenTimestepInROMSFiles)
	
    #startReleaseJD=startReleaseTime.timetuple().tm_yday
    #endReleaseJD=endReleaseTime.timetuple().tm_yday
    #timeStepsRelease=int(int(hoursOfRelease[0])/hoursBetweenTimestepInROMSFiles)
	
    print("=>SIMULATION: Drift simulation will run for %s simulation hours" %(timeStepsSimulation))
    print("=>Release: Simulated Release will run for %s simulation hours\n initiated on %s and ending on %s"%(timeStepsSimulation,startReleaseTime,endReleaseTime))

    interval = timedelta(hours=6)
    hoursPerRelease=divmod(interval.total_seconds(), 3600) #hours per Release event
    timeStepsRelease=int(int(hoursOfRelease[0])/int(hoursPerRelease[0])) #number of Release timesteps
    ReleaseTimes = [startReleaseTime + interval*n for n in range(timeStepsRelease)] #times of Release

    num=np.random.normal(releaseParticles,int(releaseParticles/2)-1, size=len(ReleaseTimes)).astype(int)
    num=np.sort(num) #sort particles in increasing order 

    #num=np.concatenate((num[len(num)%2::2],num[::-2]),axis=0) #release the highest number of particles at the midpoint of the Release period
	
    print("Release: Simulated Release will release %s kelp particles"%(np.sum(num)))

    return num, ReleaseTimes

def kelpProperties(num):

    kelpProps=np.zeros((len(num), 4))
    mylist=[0,1,2]

    weights=[0.5069324, 0.4485244,0.10138648]
    volumes=[0.000140226, 0.000753867, 2.80452E-05]
    areas=[0.1149934, 0.05, 0.02299868]
    densities=[2000,1100,300]

    for i in range(len(num)):
        ind=random.choice(mylist)
        kelpProps[i,0]=weights[ind]
        kelpProps[i,1]=volumes[ind]
        kelpProps[i,2]=areas[ind]

        # Calculate the density of the kelp
        density=weights[ind]/volumes[ind]
        kelpProps[i,3]=densities[ind]
      
    return kelpProps

def createOutputFilenames(startTime,endTime,polygonIndex,shapefile,verticalBehavior):
    startDate=''
    if startTime.day<10:
        startDate+='0%s'%(startTime.day)
    else:
        startDate+='%s'%(startTime.day)

    if startTime.month<10:
        startDate+='0%s'%(startTime.month)
    else:
        startDate+='%s'%(startTime.month)

    startDate+='%s'%(startTime.year)

    endDate=''
    if endTime.day<10:
        endDate+='0%s'%(endTime.day)
    else:
        endDate+='%s'%(endTime.day)

    if endTime.month<10:
        endDate+='0%s'%(endTime.month)
    else:
        endDate+='%s'%(endTime.month)

    endDate+='%s'%(endTime.year)
 
    # Special file naming for KINO. Each layer has name 'species.shp' and we want teh species name only.
    head,tail=os.path.split(shapefile)
   
    specie="Kelp"
    if verticalBehavior:
		outputFilename='results/%s_polygon_%s_kelp_opendrift_%s_to_%s_vertical.nc'%(specie,polygonIndex,startDate,endDate)
		animationFilename='figures/%s_polygon_%s_kelp_animation_%s_to_%s_vertical.mp4'%(specie,polygonIndex,startDate,endDate)
		plotFilename='figures/%s_polygon_%s_kelp_plot_%s_to_%s_vertical.png'%(specie,polygonIndex,startDate,endDate)
    else:
		outputFilename='results/%s_polygon_%s_kelp_opendrift_%s_to_%s_novertical.nc'%(specie,polygonIndex,startDate,endDate)
		animationFilename='figures/%s_polygon_%s_kelp_animation_%s_to_%s_novertical.mp4'%(specie,polygonIndex,startDate,endDate)
		plotFilename='figures/%s_polygon_%s_kelp_plot_%s_to_%s_novertical.png'%(specie,polygonIndex,startDate,endDate)
    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('results'):
        os.makedirs('results')
    return outputFilename, animationFilename, plotFilename

   
def createAndRunSimulation(lowDepth,highDepth,endTime,layer,polygonIndex,shapefile,outputFilename,animationFilename,plotFilename,releaseParticles,kinoDirectory,pattern_kino,svimDirectory,pattern_svim,verticalBehavior):

    # Setup a new simulation
    o = PelagicPlanktonDrift(loglevel=0)  # Set loglevel to 0 for debug information

    #######################
    # Preparing readers
    #######################
    reader_basemap = reader_basemap_landmask.Reader(
                       llcrnrlon=16, llcrnrlat=68,
                       urcrnrlon=20, urcrnrlat=72,
                       resolution='f', projection='merc')
    o.add_reader([reader_basemap]) #Do not include basemap when stranding is deactivated

    print([s for s in pattern_kino])
    reader_kino = reader_ROMS_native.Reader([s for s in pattern_kino])
    reader_kino.interpolation = 'linearNDFast' #linearND
    reader_svim = reader_ROMS_native.Reader(svimDirectory+pattern_svim)
    reader_svim.interpolation = 'linearNDFast' #linearND

    #reader_arome = reader_netCDF_CF_generic.Reader('http://thredds.met.no/thredds/dodsC/arome25/arome_metcoop_default2_5km_latest.nc')

    o.add_reader([reader_kino,reader_svim])
    
    num, ReleaseTimes = setupSeed(hoursBetweenTimestepInROMSFiles,startTime,endTime,startReleaseTime,endReleaseTime,releaseParticles)
 
    #######################
    #Adjusting configuration
    #######################
    o.set_config('processes:turbulentmixing', True)
    o.set_config('turbulentmixing:diffusivitymodel','environment')
    o.set_config('turbulentmixing:timestep', 30) # seconds
    o.set_config('turbulentmixing:verticalresolution', 1) # default is 1 meter, but since we have longer timestep we justify it
    o.set_config('processes:verticaladvection', True)
    o.set_config('turbulentmixing:TSprofiles', True)
    
   # o.set_config('drift:scheme', 'euler')
    o.set_config('drift:scheme', 'runge-kutta')
    #del o.fallback_values['x_sea_water_velocity']
    #del o.fallback_values['y_sea_water_velocity']
    o.set_config('general:coastline_action', 'stranding') #Prevent stranding, jump back to previous position

 
    #######################
    # Seed kelp particles
    #######################
        
    kelpProps=kelpProperties(num)
   

    for i, nums in enumerate(num):
        
        if nums <= 0:
            continue
        print("Running i=%s num=%s and polygon=%s"%(i,nums,polygonIndex))
  

        o.seed_from_shapefile(shapefile, nums,featurenum=[polygonIndex], 
            z=np.random.randint(low=lowDepth, high=highDepth, size=np.shape(nums)),
            weight=kelpProps[i,0], 
            volume=kelpProps[i,1],
            diameter=kelpProps[i,2],
            density=kelpProps[i,3],
            time=ReleaseTimes[i])

    #reader_basemap.plot() 

    #########################
    # Run the model
    #########################
    o.run(end_time=endTime, time_step=timedelta(hours=2),
          outfile=outputFilename)
          #export_variables=['lon', 'lat', 'z','temp','length','weight','survival'])
    print(o)

#########################
# SETUP FOR KELP PROJECT
#########################

startTime=datetime(2016,4,10,12,0,0)
endTime=datetime(2016,5,26,23,0,0)
startReleaseTime=startTime
endReleaseTime=datetime(2016,4,12,12,0,0)
releaseParticles=4 # Per timestep multiplied by gaussian bell (so maximum is releaseParticles and minimum is close to zero)
lowDepth, highDepth = -7, -2 # in negative meters
verticalBehavior=False

hoursBetweenTimestepInROMSFiles=1

#kinoDirectory='/work/users/trondk/KINO/FORWARD/Run/RESULTS/'+str(startTime.year)+'/'
kinoDirectory='/work/shared/nn9297k/Nordfjord/'
kinoDirectory='/imr/vol1/NorFjords5/Malangen-160m_AUG2015-AUG2016/'
svimDirectory='/work/shared/imr/SVIM/'+str(startTime.year)+'/'

firstkino = int(date2num(startTime,units="days since 1948-01-01 00:00:00",calendar="standard"))
lastkino = int(date2num(endTime,units="days since 1948-01-01 00:00:00",calendar="standard"))


apattern = 'norfjords_160m_his.nc4_%s*'%(startTime.year)
argument="%s%s"%(kinoDirectory,apattern)
        
pattern_kino = glob.glob(argument)
pattern_kino.sort()
print(pattern_kino)

pattern_svim='ocean_avg_*.nc' 

 
shapefile='/work/shared/nn9297k/Kelp/Shapefile/KelpExPol_utenNASAland.shp'
   
print("=> Using shapefile %s"%(shapefile))
s = ogr.Open(shapefile)

for layer in s:
    polygons=[x+1 for x in range(layer.GetFeatureCount()-1)]
    #polygons=[1,2,3,4,7] #N.Trench,Dogger bank C, Dogger bank, German bight, Viking bank
    #polygons=[2] #N.Trench,Dogger bank C, Dogger bank, German bight, Viking bank
		
    for polygonIndex in polygons:
            

        feature = layer.GetFeature(polygonIndex-1)
           
        print("Area",feature.GetGeometryRef().GetArea())
        geom = feature.GetGeometryRef()
        points = geom.GetGeometryCount()
        ring = geom.GetGeometryRef(0)
        print("jj",polygonIndex, points)

        if ring.GetPointCount() > 3:
            outputFilename, animationFilename, plotFilename = createOutputFilenames(startTime,endTime,polygonIndex,shapefile,verticalBehavior)
              
            print("Result files will be stored as:\nnetCDF=> %s\nmp4=> %s"%(outputFilename,animationFilename))

            createAndRunSimulation(lowDepth,highDepth,endTime,
                    layer,polygonIndex,shapefile,
                    outputFilename,animationFilename,plotFilename,releaseParticles,
                    kinoDirectory,pattern_kino,svimDirectory,pattern_svim,verticalBehavior)
