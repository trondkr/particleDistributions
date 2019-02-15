# coding=utf-8

import os, sys
import numpy as np
import glob
import string
import matplotlib

#matplotlib.use('Agg')
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from pylab import *
from datetime import datetime
from pprint import pprint
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
import ogr
import osr


__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime(2016, 8, 10)
__modified__ = datetime(2016, 8, 10)
__version__  = "1.0"
__status__   = "Production"

# --------
# calculateaverages.py
#
# This script takes the output from opendrift and calculates area averages 
# as a function of time. The total area around North Sea is divided into
# bins of specific resolution and the number of particles within each bin 
# is summed for a specific time period (e.g. 1 month). The total output 
# is a heatmap of where the most particles reside for each time period.
# --------

def createBins(requiredResolution):

    print('func: createBins() => Creating bins for averaging')
    xmin=15.0; xmax=21.0
    ymin=69.0; ymax=72.0

    deg2rad=np.pi/180.
    R = 6371  # radius of the earth in km
    # Distance from minimum to maximim longitude
    x = (xmax*deg2rad - xmin*deg2rad) * cos( 0.5*(ymax*deg2rad+ymax*deg2rad) )
    y =  ymax*deg2rad - ymax*deg2rad
    dx = R * sqrt( x*x + y*y )
    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))

    # Distance from minimum to maximim latitude
    x = (xmax*deg2rad - xmax*deg2rad) * cos( 0.5*(ymax*deg2rad+ymin*deg2rad) )
    y =  ymax*deg2rad - ymin*deg2rad
    dy = R * sqrt( x*x + y*y )

    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))

    ngridx = int(np.round(dx/requiredResolution,0))
    ngridy = int(np.round(dy/requiredResolution,0))
    
    xi = np.linspace(np.floor(xmin),np.ceil(xmax),ngridx)
    yi = np.linspace(np.floor(ymin),np.ceil(ymax),ngridy)

    print('=> created binned array of domain of size (%s,%s) with resolution %s'%(ngridx,ngridy,requiredResolution))

    return xi,yi

def calculateAreaAverages(xi,yi,cdf,weeksInYear):

    print('func: calculateAreaAverages() => Calculating averages within bins')
    print('=> binned domain (%2.1f,%2.1f) to (%2.1f,%2.1f)'%(np.min(xi),np.min(yi),np.max(xi),np.max(yi)))

    timesteps = cdf.variables['time'][:]
    timeunits = cdf.variables["time"].units
    
    print('=> found %s timesteps in input file'%(len(timesteps)))
    newWeek=-9

    for tindex, t in enumerate(timesteps): 
        currentDate = num2date(t, units=timeunits, calendar="gregorian")
        year,weekNumber,DOW = currentDate.isocalendar()
       
        Xpos = cdf.variables['lon'][:,tindex]
        Ypos = cdf.variables['lat'][:,tindex]
        
        H, xedges, yedges = np.histogram2d(Xpos, Ypos, bins=(xi, yi), normed=False)
        
        if (tindex==0):
            weeklyFrequency=np.zeros((weeksInYear,np.shape(H)[0],np.shape(H)[1]), dtype=float32)
            print("Inside t==0")
        if weekNumber != newWeek:
            print("=> Adding data to week: %s (startdate: %s)"%(weekNumber,currentDate))
            weeklyFrequency[weekNumber,:,:]=weeklyFrequency[weekNumber,:,:] + H
            newWeek=weekNumber

    # Create log values and levels for frequencyplot
    #weeklyFrequency=ma.log(weeklyFrequency)
    levels = np.arange(weeklyFrequency.min(),weeklyFrequency.max(),(weeklyFrequency.max()-weeklyFrequency.min())/10)
    print(levels)
    sigma = 0.2 # this depends on how noisy your data is, play with it!

    return gaussian_filter(weeklyFrequency, sigma)

def plotDistribution(kelpData,week,baseout,xii,yii,distType):
    print("Plotting the distributions for week: %s"%(week))
    plt.clf()
    plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,resolution='f', 
                    projection='merc')

    x, y = mymap(xii, yii)
    #levels=np.arange(np.min(kelpData),np.max(kelpData),10)
    levels=np.arange(1,100,5)
                              
    CS1 = mymap.contourf(x,y,np.fliplr(np.rot90(kelpData,3)),
    levels,cmap=cm.get_cmap('Spectral_r',len(levels)-1), 
    extend='max',alpha=1.0)

    plt.colorbar(CS1,orientation='vertical',extend='max', shrink=0.5)

    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey',zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    plt.title('Kelp week: %s'%(week))
    if distType=="integrated":
        plotfile=baseout+'/Kelp_distribution_all_fullperiod.png'
    else:
        plotfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'.png'
    print("=> Creating plot %s"%(plotfile))             
    plt.savefig(plotfile,dpi=300)
 
def main():
    # Which species to calculate for
    # The timespan part of the filename

    shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/Shapefile05112018/kelpExPol_11sept2018.shp'
    experiment = 1

    if experiment == 1:
        startdate = '01052016'
        enddate = '01082016'
    if experiment == 2:
        startdate = '01032016'
        enddate = '15052016'
    if experiment == 3:
        startdate = '01112015'
        enddate = '01042016'
    if experiment == 4:
        startdate = '01052016'
        enddate = '01082016'

    # Results and storage folders
    base='results15012019'
    baseout='distributionFigures'
    
    if not os.path.exists(baseout): os.makedirs(baseout)

    # The resolution of the output grid in kilometers
    requiredResolution = 1.000 # km between each binned box

    # END EDIT ----------------------------------

    # Create the grid you want to calculate frequency on
    xi,yi = createBins(requiredResolution)
    weeksInYear=52
    xii,yii=np.meshgrid(xi[:-1],yi[:-1])

    s = ogr.Open(shapefile)
    for layer in s:
        polygons = [x + 1 for x in range(layer.GetFeatureCount()-1)]
       # polygons=[1,2]
  
    allData=np.zeros((len(polygons),weeksInYear,len(xi)-1,len(yi)-1))
    print("=> Created final array for all data of size :",np.shape(allData))
        
    for polygonIndex, polygon in enumerate(polygons):
        infile = base+'/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, startdate, enddate)

        print("=> Opening input file: %s"%(infile))
        if os.path.exists(infile):
            cdf = Dataset(infile)
            allData[polygonIndex,:,:,:]=calculateAreaAverages(xi,yi,cdf,weeksInYear)
            cdf.close()
        else:
            print("==>> Input file %s could not be opened"%(infile))
    
   # 
        
    for week in range(1,weeksInYear,1):
        # Calculate the cumulative distribution for each week
        kelpData=np.squeeze(np.sum(allData[:,week,:,:], axis=0))
        #print("1", allData[0,week,10,10])
        #print("2", allData[1,week,10,10])
       
      #  print("4", np.ma.sum(allData[:,week,10,10], axis=0))
      #  print(np.ma.min(kelpData),np.ma.max(kelpData))
        levels=np.arange(np.ma.min(kelpData),np.ma.max(kelpData),0.5)
      
        if (len(levels)>2 and np.ma.mean(kelpData)>0):
            print("Kelp data {}".format(np.shape(kelpData)))
        #    plotDistribution(kelpData,week,baseout,xii,yii,"weekly")
            
    # Plot the distribution for all weeks
    kelpData=np.squeeze(np.ma.sum(allData[:,:,:,:], axis=0))
    print("LAST Kelp data {}".format(np.shape(kelpData)))
    kelpData=np.ma.masked_where(kelpData==0,kelpData)
    kelpData=np.squeeze(np.ma.sum(kelpData[:,:,:], axis=0))
    plotDistribution(kelpData,week,baseout,xii,yii,"integrated")

if __name__ == "__main__":
    main()

