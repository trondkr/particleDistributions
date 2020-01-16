# coding=utf-8
import matplotlib
matplotlib.use('Agg')
import os, sys
import numpy as np
import glob
import matplotlib
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xarray as xr
from datetime import datetime
from scipy.ndimage.filters import gaussian_filter
import time
import utils
import pandas as pd

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime(2016, 8, 10)
__modified__ = datetime(2016, 8, 10)
__version__  = "1.0"
__status__   = "Production"

   
# Sigma is used in gaussian filter 
# gaussian_filter(weeklyFrequency, sigma) 
# this depends on how noisy your data is, play with it!

# --------
# calculateaverages.py
#
# This script takes the output from opendrift and calculates area averages 
# as a function of time. The total area around North Sea is divided into
# bins of specific resolution and the number of particles within each bin 
# is summed for a specific time period (e.g. 1 month). The total output 
# is a heatmap of where the most particles reside for each time period.
# --------
def distance(lat1, lon1, lat2, lon2):
    p = 0.017453292519943295     #Pi/180
    a = 0.5 - np.cos((lat2 - lat1) * p)/2 + np.cos(lat1 * p) * np.cos(lat2 * p) * (1 - np.cos((lon2 - lon1) * p)) / 2
    return 12742 * np.arcsin(np.sqrt(a))

def createBins(requiredResolution):

    print('func: createBins() => Creating bins for averaging')
    xmin=15.0; xmax=21.0
    ymin=69.0; ymax=72.0
   
    dy = distance(ymin,xmin,ymax,xmin)
    dx = distance(ymin,xmin,ymin,xmax)

    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))
    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))
    
    delta = int(np.round(dx/requiredResolution,0))
    print("delta {}".format(delta))
    xi = np.linspace(np.floor(xmin),np.ceil(xmax),delta)
    deltaX=abs(xi[0]-xi[1])
    print(deltaX)

    # We are only using longitude to calculate the approximate distance in degrees 
    # that is equvalent to requiredResolsution. Then we create latitude grid based on the 
    # same deltaX.
    yi = np.arange(np.floor(ymin),np.ceil(ymax),deltaX)

    print('=> created binned array of domain of grid cell size (%s) with resolution %s'%(delta,requiredResolution))
    
    return xi,yi,deltaX

def calculateAreaAverages(xi,yi,cdf,weeksInYear):

    print('func: calculateAreaAverages() => Calculating averages within bins')
    print('=> binned domain (%2.1f,%2.1f) to (%2.1f,%2.1f)'%(np.min(xi),np.min(yi),np.max(xi),np.max(yi)))

    for week in range(1,weeksInYear):
        weekcdf = cdf.where(cdf.time.dt.week == week,drop = True)
        
        #utils.filter_sedimented()
        #TODO:Filter by only sedimented particles
 
        Xpos = weekcdf.coords['lon'].values.flatten()
        Ypos = weekcdf.coords['lat'].values.flatten()               

        # create average for one week
        H, xedges, yedges = np.histogram2d(Xpos, Ypos, bins=(xi, yi), normed=False)  
        if week == 1: 
            weeklyFrequency=np.zeros((weeksInYear,np.shape(H)[0],np.shape(H)[1]), dtype=np.float32)
           
        weeklyFrequency[week,:,:] = H

    return gaussian_filter(weeklyFrequency, sigma)

def plotDistribution(kelpData,week,baseout,xii,yii,distType,polygons,experiment,plotfile,showgrid=None):
    print("Plotting the distributions for week: %s"%(week))
    plt.clf()
    plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,resolution='h', 
                    projection='merc')

    x, y = mymap(xii, yii)
    
    # = np.arange(10,110,10)  
    z = np.fliplr(np.rot90(kelpData,3))     
    #z = 100*z/np.max(z)  # % of max data 
    #levels = np.arange(0,10,0.1)
    levels=np.linspace(1,np.max(z)/3,20)
    #cblevels = np.arange(0,10,5)    
    CS1 = mymap.contourf(x,y,z,vmin = 1,vmax = np.max(z)/3, levels = 100, cmap=cm.get_cmap('Spectral_r'), extend='max') 
    #levels,,len(levels)-1 Spectral_r #, alpha=1.0levels
    #CS2 = mymap.contour(x,y,z,[1], colors = 'k',linewidths = 0.1)

    plt.colorbar(CS1,orientation='vertical',extend='max', shrink=0.7) #,ticks = cblevels)
    if showgrid:
        print("Creating showgrid")
        mymap.plot(x,y,lineWidth=0.5, color='k',zorder=1)
        mymap.plot(y,x,lineWidth=0.5, color='k',zorder=1)
    
    mymap.drawmapboundary(fill_color='#677a7a')
    mymap.fillcontinents(color='#b8a67d',zorder=2)
    mymap.drawcoastlines()

    plt.title('Kelp week: {} , polygons: {}-{}, kelp Type {}, experiment {}'.format(
                            week,min(polygons),max(polygons),kelpType,experiment))
    
    print("=> Creating plot %s"%(plotfile))             
    plt.savefig(plotfile,dpi=300)
 

def writeToAsciiRaster(dataarray,gisfile,requiredResolution,lons,lats,deltaDegrees):
    """Write data to ascii file"""
    
    if os.path.exists(gisfile): os.remove(gisfile)
    f = open(gisfile,'a')

    z = np.fliplr(np.rot90(dataarray,3)) 

    # Write the results to an ascii file suitable for arcgis input"""
    cellsize         =deltaDegrees
    Xlowerleftcorner =lons.min()
    Ylowerleftcorner =lats.min()
    nodataValue="-9999.0"

    # All the latitude/longitude values given here indicate the mid-point of the
    # location where the value is found. Therefore, to correctly position the cell, we
    # regard the minimum long/lat point as the mid-point and write the lower left corner as
    # half a cellsize down from mid-point.

    x_range=np.shape(z)[0]
    y_range=np.shape(z)[1]

    header1='ncols %i \n'%(x_range)
    header2='nrows %i \n'%(y_range)
    header3='xllcorner %.5f \n'%((Xlowerleftcorner) + (requiredResolution/2.0))
    header4='yllcorner %.5f \n'%((Ylowerleftcorner) + (requiredResolution/2.0))
    header5='cellsize %.5f \n'%(cellsize)
    header6='NODATA_VALUE %s \n'%(nodataValue)

    f.writelines(header1)
    f.writelines(header2)
    f.writelines(header3)
    f.writelines(header4)
    f.writelines(header5)
    f.writelines(header6)

    for y in range(y_range):
        for x in range(x_range):
          #  if y > 1 and y < y_range-1 and x > 1 and x < x_range-1:
          #      print("delta dx {} dy {}".format(lons[x,y]-lons[x,y-1],lats[x,y]-lats[x-1,y]))
            d='%f\t'%(z[x,(y_range-1)-y])
            f.writelines(d)

    f.close()


def writeToAscii(dataarray,xyzfile,requiredResolution,lons,lats,deltaDegrees):
    """Write data to ascii file"""
    
    if os.path.exists(xyzfile): os.remove(xyzfile)
    f = open(xyzfile,'a')
    x_range=np.shape(dataarray)[0]
    y_range=np.shape(dataarray)[1]
    z = np.fliplr(np.rot90(dataarray,3))   

    header1='longitude, latitude, value\n'
    f.writelines(header1)
   
    for x in range(x_range-1):
        for y in range(y_range-1):
            
            d='{},{},{}\n'.format(lons[y,x],lats[y,x],z[y,x])
            f.writelines(d)
        x=0

    f.close()

def getResultNames(distType,kelpType,experiment):
    if distType == "integrated":
        plotfile=baseout+'/Kelp_distribution_all_fullperiod_kelpType{}_experiment_{}.png'.format(kelpType,experiment)
        gisfile=baseout+'/Kelp_distribution_all_fullperiod_kelpType{}_experiment_{}.txt'.format(kelpType,experiment)
        xyzfile=baseout+'/Kelp_distribution_all_fullperiod_kelpType{}_experiment_{}.xyz'.format(kelpType,experiment)
    else:
        plotfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'polygons_{}-{}'.format(min(polygons),max(polygons))+'kelpType_{}.png'.format(kelpType)
        gisfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'polygons_{}-{}'.format(min(polygons),max(polygons))+'kelpType_{}.txt'.format(kelpType)
        xyzfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'polygons_{}-{}'.format(min(polygons),max(polygons))+'kelpType_{}.xyz'.format(kelpType)
    return plotfile, gisfile, xyzfile

def main(shapefile,requiredResolution,experiment,distType,polygons = None,kelpType = None, showGrid = None):
    
    # Which species to calculate for
    # The timespan is part of the filename
    startdate,enddate = utils.ranges[experiment]

    # Create the grid you want to calculate frequency on
    # The resolution of the output grid in km between each binned box
    xi,yi,deltaDegrees = createBins(requiredResolution)
    xii,yii=np.meshgrid(xi[:-1],yi[:-1])

    if polygons == None:
        polygons = utils.get_polygons()
        print(polygons)
    #Created final array for all data of size : (17, 52, 205, 333)
    # 17 polygons in shapefile , 52 dimension for weeks, xi,yi - bins 
    allData=np.zeros((len(polygons),weeksInYear,len(xi)-1,len(yi)-1))
    print("=> Created final array for all data of size :",np.shape(allData))
    infiles = utils.get_paths(polygons,experiment)
   
    plotfile, gisfile, xyzfile = getResultNames(distType,kelpType,experiment)
    print("Kelptype {}".format(kelpType))
    for polygonIndex,path in enumerate(infiles):   
        if os.path.exists(path):
            #print("=> Opening input file: %s"%(path))
            cdf = xr.open_dataset(path)   
            if kelpType != None: # if kelpType is specified, we filter the dataframe by type 
                cdf = cdf.where(cdf.plantpart == kelpType,drop = True) 
        
            allData[polygonIndex,:,:,:] = calculateAreaAverages(xi,yi,cdf,weeksInYear)
            cdf.close()
        else:
            print("==>> Input file %s could not be opened"%(path))

    if distType == "weekly": 
        for week in range(1,weeksInYear,1):
            # Calculate the cumulative distribution for each week
            kelpData=np.squeeze(np.sum(allData[:,week,:,:], axis=0))
            levels=np.arange(np.ma.min(kelpData),np.ma.max(kelpData),0.5)        
            if (len(levels)>2 and np.ma.mean(kelpData)>0):
                print("Kelp data {}".format(np.shape(kelpData)))
                plotDistribution(kelpData,week,baseout,xii,yii,"weekly",polygons,experiment,plotfile,showGrid)

    elif distType == "integrated":     
        week = 'Sum over all time range'  
        # Plot the distribution for all weeks
        # sum data along all polygons and all weeks. 
        kelpData=np.ma.sum(allData[:,:,:,:], axis=(0,1))     
        print("LAST Kelp data {}".format(np.shape(kelpData)))
        # if value is zero for all polygons,weeks the position is masked
        kelpData=np.ma.masked_where(kelpData<=0,kelpData)
        plotDistribution(kelpData,week,baseout,xii,yii,"integrated",polygons,experiment,plotfile,showGrid)
  
    np.ma.set_fill_value(kelpData,-9999.0)
  #  writeToAsciiRaster(kelpData.filled(),gisfile,requiredResolution,xii,yii,deltaDegrees)
    writeToAscii(kelpData.filled(),xyzfile,requiredResolution,xii,yii,deltaDegrees)

if __name__ == "__main__":
    start_time = time.time()

    # Results and storage folders

    global base,baseout,weeksInYear,sigma 
    global kelpType

    base= r'Data'
    baseout='distributionFigures'
    shapefile = 'Shapefile05112018/kelpExPol_11sept2018.shp'
    if not os.path.exists(baseout): os.makedirs(baseout)   

    weeksInYear = 52
    showGrid=True
    ##Parameters: 
    sigma = 0.5  
    kelpType = 1
    alltypes = (0,1,2,4)  
    experiments = (1,2,3,4,5)
    # kelpType is the parameter for chosing, the type of kelp (differs by density?)
    # available types are 0,1,2,4; -32767 is used for masking in the nc file ,
    #  None can be used for all types       
    requiredResolution=1 # km

    polygons = None 
    #polygons = (1,2,3) # #None means all polygons, otherwise create a list of needed polygons 

    distType ="integrated"
    # integrated to plot all positions of particles during recorded time
    # 'weekly' to plot separately positions during all weeks in the recorded time range 

    #requiredResolution = 1   Resolution of the data, used in creating bins 
    for experiment in experiments:
        for kelpType in alltypes:

            main(shapefile, requiredResolution, experiment = experiment, distType = distType,polygons=polygons, kelpType=kelpType, showGrid=showGrid)

    #        kelpType = kelpType
    #        [main(shapefile,experiment = e, distType = distType,polygons = polygons,
    #    requiredResolution = 1,kelpType = kelpType) for e in experiments]'''

    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   
