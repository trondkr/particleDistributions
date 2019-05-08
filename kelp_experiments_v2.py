#!/usr/bin/env python

from datetime import datetime, timedelta
import numpy as np
import time

from opendrift.readers import reader_basemap_landmask
from opendrift.readers import reader_ROMS_native
from opendrift.readers import reader_netCDF_CF_generic
from kelp.kelpClass_v2 import KelpDrift
import os
from netCDF4 import Dataset, date2num, num2date
from datetime import datetime
import random
import math
import glob
from random import randint
from calendar import monthrange
import confOpenDrift
from osgeo import gdal, osr, ogr


def commonKelpProperties():
    # S feces    = 0
    # New blade / lamina = 1
    # Stipe/Thallus     = 2
    # Fragment  = 3
    # M/L feces = 4

    sinkspeed = [ 0.008, 0.074, 0.173, 0.036, 0.013]
    sinkspeedsstd = [0.003, 0.038, 0.068, 0.020, 0.006]
    return sinkspeed, sinkspeedsstd #weights, densities, SDdensities, areas, lengths, volumes, diameters


def commonDateProperties(experiment):
    if experiment == 0:
        # Test experiment
        startTime = datetime(2016, 6, 1, 0, 0, 0)
        endTime = datetime(2016, 6, 20, 0, 0, 0) #8
    if experiment == 1:
        startTime = datetime(2016, 5, 1, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0) #8
    if experiment == 2:
        startTime = datetime(2016, 3, 1, 0, 0, 0)
        endTime = datetime(2016, 5, 15, 0, 0, 0) # 5
    if experiment == 3:
        startTime = datetime(2015, 11, 20, 0, 0, 0)
        endTime = datetime(2016, 4, 1, 0, 0, 0) # 4
    if experiment == 4:
        startTime = datetime(2016, 5, 1, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0) #8
    if experiment == 5:
        startTime = datetime(2015, 8, 3, 0, 0, 0)
        endTime = datetime(2016, 8, 1, 0, 0, 0)

    return startTime, endTime


def kelpProperties(num, kelpTypes):
    # Get the options of weights and densities
    sinkspeeds, sinkspeedsstd = commonKelpProperties()

    # Loop over num release dates and randomly select type from provided kelpTypes list
    # Kelptypes list contents indicate which of the indices/options you have selected in commonKelpProperties
    # If onlye old blades: kelpTypes=[0]

    kelpSinkSpeeds = []
    kelpType = []
    for i in range(len(num)):
        ind = random.choice(kelpTypes)
        kelpType.append(ind)
        # Calculate a random sinking speed based on the mean and std values
        kelpSinkSpeeds.append(np.random.normal(sinkspeeds[ind], sinkspeedsstd[ind] / 1.0))
        print(kelpSinkSpeeds[i], sinkspeeds[ind], sinkspeedsstd[ind])
    return kelpSinkSpeeds, kelpType


def createOutputFilenames(confOD, experiment, polygonIndex):
    startTime, endTime = commonDateProperties(experiment)
    startDate = ''
    if startTime.day < 10:
        startDate += '0%s' % (startTime.day)
    else:
        startDate += '%s' % (startTime.day)

    if startTime.month < 10:
        startDate += '0%s' % (startTime.month)
    else:
        startDate += '%s' % (startTime.month)

    startDate += '%s' % (startTime.year)

    endDate = ''
    if endTime.day < 10:
        endDate += '0%s' % (endTime.day)
    else:
        endDate += '%s' % (endTime.day)

    if endTime.month < 10:
        endDate += '0%s' % (endTime.month)
    else:
        endDate += '%s' % (endTime.month)

    endDate += '%s' % (endTime.year)

    # Special file naming for KINO. Each layer has name 'species.shp' and we want teh species name only.
    head, tail = os.path.split(confOD.shapefile)

    specie = "Kelp"
    confOD.outputFilename = 'results/%s_polygon_%s_experiment_%s_%s_to_%s.nc' % (
        specie, polygonIndex, experiment, startDate, endDate)
    confOD.animationFilename = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.mp4' % (
        specie, polygonIndex, experiment, startDate, endDate)
    confOD.plotFilename = 'figures/%s_polygon_%s_experiment_%s_%s_to_%s.png' % (
        specie, polygonIndex, experiment, startDate, endDate)

    if not os.path.exists('figures'):
        os.makedirs('figures')
    if not os.path.exists('results'):
        os.makedirs('results')
    if os.path.exists(confOD.outputFilename):
            os.remove(confOD.outputFilename)
            print(("Removing file ", confOD.outputFilename))
    if os.path.exists('figures/'+confOD.animationFilename):
            os.remove('figures/'+confOD.animationFilename)
    if os.path.exists('figures/'+confOD.plotFilename):
            os.remove('figures/'+confOD.plotFilename)
    

def createAndRunSimulation(confOD, experiment,layer, polygonIndex, allNum, allReleaseTimes, allKelpSpeeds, allKelpTypes):
 
    o = KelpDrift(loglevel=20)  # Set loglevel to 0 for debug information
    startTime, endTime = commonDateProperties(experiment)

    allKelpSpeeds_flat = [item for sublist in allKelpSpeeds for item in sublist]
    allKelpTypes_flat = [item for sublist in allKelpTypes for item in sublist]
    allNum_flat = [item for sublist in allNum for item in sublist]
    allReleaseTimes_flat=[item for sublist in allReleaseTimes for item in sublist]
    print(("=> Simulation will release a total of %s particles\n" % (np.sum(allNum_flat))))

    # Randomly distribute the particles at depths varying between lowDepth and highDepth
    depths = [randint(confOD.lowDepth, confOD.highDepth) for i in range(len(allNum_flat))]

    #######################
    # Preparing readers
    #######################
    #reader_basemap = reader_basemap_landmask.Reader(
    #    llcrnrlon=15, llcrnrlat=68,
    #    urcrnrlon=23, urcrnrlat=74,
    #    resolution=mapResolution, projection='merc')
    #o.add_reader([reader_basemap])  # Do not include basemap when stranding is deactivated

    #wind = reader_netCDF_CF_generic.Reader('http://oos.soest.hawaii.edu/thredds/dodsC/hioos/model/atm/ncep_global/NCEP_Global_Atmospheric_Model_best.ncd')
    #o.add_reader(wind)
    
    o.set_config('general:use_basemap_landmask', False)
    #o.set_config('general:basemap_resolution','h')

    if confOD.use_svim:
        reader_svim2015 = reader_ROMS_native.Reader(svimfiles2015)
        reader_svim2015.interpolation = confOD.interMethod
        reader_svim2016 = reader_ROMS_native.Reader(svimfiles2016)
        reader_svim2016.interpolation = confOD.interMethod
 
    if confOD.use_norkyst:
        reader_norkyst = reader_netCDF_CF_generic.Reader(confOD.allnorkyst)
        reader_norkyst.interpolation = confOD.interMethod

    if not confOD.use_svim and confOD.use_norkyst:
        o.add_reader([reader_roms, reader_norkyst])

    reader_roms = reader_ROMS_native.Reader(confOD.argument)
    reader_roms.interpolation = confOD.interMethod

    if confOD.use_svim and not confOD.use_norkyst:
        o.add_reader([reader_roms, reader_svim2016, reader_svim2015])

    if confOD.use_svim and confOD.use_norkyst:
        o.add_reader([reader_roms, reader_norkyst, reader_svim2016, reader_svim2015])

    #######################
    # Adjusting configuration
    #######################
    o.set_config('processes:turbulentmixing', True)
    o.set_config('turbulentmixing:diffusivitymodel', 'windspeed_Sundby1983')
    o.set_config('turbulentmixing:timestep', 60)  # secondsS
    o.set_config('turbulentmixing:verticalresolution', 1)  # default is 1 meter, but since we have longer timestep we justify it
    o.set_config('processes:verticaladvection', True)
    o.set_config('turbulentmixing:TSprofiles', True)
    #o.set_config('turbulentmixing:max_iterations', 100)

    o.set_config('drift:scheme', 'euler')
    o.set_config('general:coastline_action', 'previous')  # Prevent stranding, jump back to previous position
  

    # depths=[randint(lowDepth,highDepth) for i in xrange(len(allKelpWeights))]
    allNum_flatInt = [int(i) for i in allNum_flat]
    for i, nums in enumerate(allNum_flatInt):

        if nums <= 0:
            continue
        print(("Running i=%s num=%s and polygon=%s" % (i, nums, polygonIndex)))
        print("allKelpSpeeds_flat[i] ", allKelpSpeeds_flat[i])
        o.seed_from_shapefile(confOD.shapefile, allNum_flatInt[i], featurenum=[polygonIndex],
                              z="seafloor+1",  # depths[i],
                              sink_speed=allKelpSpeeds_flat[i],
                              time=allReleaseTimes_flat[i],
                              plantpart=allKelpTypes_flat[i])

    # reader_basemap.plot()

    #########################
    # Run the model
    #########################
   # o.plot()

    o.run(end_time=endTime, time_step=timedelta(hours=1),
          outfile=confOD.outputFilename,
    export_variables=['lon', 'lat', 'z', 'temp','sink_speed','plantpart','sea_floor_depth_below_sea_level'])
    print(o)


def setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime, releaseParticles):
    ##################################################
    # Create seed variation as function of day
    # Called multiple times  from setupSeedsForExperiment
    ##################################################

    difference = endTime - startTime
    hoursOfSimulation = divmod(difference.total_seconds(), 3600)

    difference = endReleaseTime - startReleaseTime
    hoursOfRelease = divmod(difference.total_seconds(), 3600)

    timeStepsSimulation = int(int(hoursOfSimulation[0]) / 3)

    print("=>Release: Simulated Release will run for %s simulation hours\n initiated on %s and ending on %s"%(timeStepsSimulation,startReleaseTime,endReleaseTime))

    interval = timedelta(hours=intervalHours)
    hoursPerRelease = divmod(interval.total_seconds(), 3600)  # hours per Release event
    timeStepsRelease = int(int(hoursOfRelease[0]) / int(hoursPerRelease[0]))  # number of Release timesteps
    ReleaseTimes = [startReleaseTime + interval * n for n in range(timeStepsRelease)]  # times of Release

    # num=np.random.normal(releaseParticles,int(releaseParticles/2)-1, size=len(ReleaseTimes)).astype(int)
    num = [releaseParticles for n in range(timeStepsRelease)]
    # num=np.sort(num) #sort particles in increasing order

    print(("=> Seed episode: %s => Release of %s kelp particles" % (seedCount, np.sum(num))))

    return num, ReleaseTimes


def setupSeedsForExperiment(experiment, releaseParticles):
    print("\nSeed setup started --------")
    seedCount = 1
    allNum = []
    allReleaseTimes = []
    allKelpSpeeds = []
    allKelpTypes = []

    # Batch one : Old lamina (77%) released evenly 6 times a day between 122 and 135. 
    startTime, endTime = commonDateProperties(experiment)
    startReleaseTime = startTime
    endReleaseTime = endTime

    intervalHours = 24
    print(("=> Release: daily: %s to %s" % (startReleaseTime, endReleaseTime)))

    num, ReleaseTimes = setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime,
                                  releaseParticles)
    allNum.append(num)
    allReleaseTimes.append(ReleaseTimes)
    seedCount += 1

    # Release thallus, new, and stipe
    kelpTypes = [0,1,2]
    kelpSpeeds, kelpTypes = kelpProperties(num,kelpTypes)
    allKelpSpeeds.append(kelpSpeeds)
    allKelpTypes.append(kelpTypes)
    # Batch two : feces released evenly 6 times a day once a week 
    # Find the number of weeks between start and stop
    num_of_weeks = int(math.ceil((endTime - startTime).days / 7.0))

    # Release once per week
    for week in range(int(num_of_weeks) - 1):
        dayN = week * 7
        startReleaseTime = startTime + timedelta(days=dayN)
        endReleaseTime = startTime + timedelta(days=dayN + 1)
        print(("=> Release: weekly: %s to %s" % (startReleaseTime, endReleaseTime)))
        intervalHours = 6

        num, ReleaseTimes = setupSeed(seedCount, intervalHours, startTime, endTime, startReleaseTime, endReleaseTime,
                                      releaseParticles)
        allNum.append(num)
        allReleaseTimes.append(ReleaseTimes)

        seedCount += 1
        # Define the properties of kelp to be released in batch 2
        # new blades, stipes, and fragments (in equal counts - check with Eli)
        kelpTypes = [3,4]
        kelpSpeeds, kelpTypes = kelpProperties(num, kelpTypes)
        allKelpSpeeds.append(kelpSpeeds)
        allKelpTypes.append(kelpTypes)

    print("Seed setup done --------\n")
    # Return the total number of particles per release date
    return allNum, allReleaseTimes, allKelpSpeeds, allKelpTypes


#########################
# SETUP FOR KELP PROJECT
#########################

confOD=confOpenDrift.ConfOpenDrift()

for experiment in confOD.experiments:
    startTime, endTime = commonDateProperties(experiment)

    if confOD.use_svim:
    
        svimfiles2015 = [confOD.definesvimpath(2015) % d for d in range(1, 13, 1)]
        svimfiles2016 = [confOD.definesvimpath(2016) % d for d in range(1, 13, 1)]

    if confOD.use_norkyst:
        
        for year in [2016]:
            for month in range(5,7,1):
                daysinmonth=monthrange(year, month)
                for day in range(1,daysinmonth[1],1):
                    hour=12
                    if (month==6 and day==15): hour=11

                    filename=confOD.norkyst.format(year,month,day,hour)
                
                    print("\nYear: {} month {:02d} day {:02d}".format(year,month,day))
                    print(filename)

                    if not (month==7 and day==26) and not (month==12 and day==24) and not (month==12 and day==26):
                        confOD.allnorkyst.append(filename)


    print("=============================================")
    print(("Starting kelp simulations for experiment %s\n" % (experiment)))
    print(("    %s" % (time.strftime("%c"))))
    print("=============================================")
    allNum, allReleaseTimes, allKelpSpeeds, allKelpTypes = setupSeedsForExperiment(experiment, confOD.releaseParticles)

    # Remember to setup/edit kelp properties in function: commonKelpProperties unless
    # using the default values

    # Setup complete --------------------------------------------------------------
    # No edit below

    print(("=> Using shapefile %s" % confOD.shapefile))
    s = ogr.Open(confOD.shapefile)

    # Find all kelp polygons in Shapefile
    for layer in s:
        polygons = [x + 1 for x in range(layer.GetFeatureCount() - 1)] 
   
        print(('Running for layer with %s features)' % (layer.GetFeatureCount())))
        # Loop over all kelp polygons, releasing kelp and tracking their drift, writing results to file
        for polygonIndex in polygons:

            feature = layer.GetFeature(polygonIndex - 1)

            print(("Kelp area %s for polygon %s" % (feature.GetGeometryRef().GetArea(), polygonIndex)))
            geom = feature.GetGeometryRef()
            points = geom.GetGeometryCount()
            ring = geom.GetGeometryRef(0)

            if ring.GetPointCount() > 3:
                createOutputFilenames(confOD, experiment, polygonIndex)
                print(("Result files will be stored as:\nnetCDF=> %s\nmp4=> %s" % (confOD.outputFilename, confOD.animationFilename)))

                print("Starting simulations....")
                createAndRunSimulation(confOD, experiment,layer, polygonIndex,
                                       allNum, allReleaseTimes, allKelpSpeeds, allKelpTypes)
