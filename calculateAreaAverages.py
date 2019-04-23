# coding=utf-8

import os, sys
import numpy as np
import glob
#import string
import matplotlib
import geopandas
#matplotlib.use('Agg')
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
#import matplotlib.path as mpath
#import matplotlib.patches as mpatches
#from matplotlib.collections import PatchCollection
#from pylab import *
from datetime import datetime
#from pprint import pprint
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
import ogr
#import osr
import time

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
    def get_dist(xx,yy): 
        return R * np.sqrt( xx**2 + yy**2 )

    xmaxrad = xmax*deg2rad
    xminrad = xmin*deg2rad
    ymaxrad = ymax*deg2rad
    yminrad = ymin*deg2rad   

    dx = get_dist((xmaxrad - xminrad) * np.cos(ymaxrad),0)
    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))

    # Distance from minimum to maximim latitude
    dy = get_dist(0,ymaxrad-yminrad)   
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

    timesteps = cdf['time']
    #timeunits = cdf["time"].units    
    print('=> found %s timesteps in input file'%(len(timesteps)))
    newWeek=-9  
    nWeeksWithData = []

    for tindex, t in enumerate(timesteps): 
        currentDate = t.values
        #currentDate = num2date(t, units=timeunits, calendar="gregorian")
        weekNumber = t.dt.week.values #isocalendar()

        # should be faster to add weeknumber as a column
        # then to group by week number and find the sum ?

        Xpos = cdf.coords['lon'][:,tindex].values
        Ypos = cdf['lat'][:,tindex].values               
        H, xedges, yedges = np.histogram2d(Xpos, Ypos, bins=(xi, yi), normed=False)
        
        if (tindex==0):
            weeklyFrequency=np.zeros((weeksInYear,np.shape(H)[0],np.shape(H)[1]), dtype=np.float32)
            print("Inside t==0")

        if weekNumber != newWeek:
            #print("=> Adding data to week: %s (startdate: %s)"%(weekNumber,currentDate))
            weeklyFrequency[weekNumber,:,:] += H
            newWeek=weekNumber
            nWeeksWithData.append(weekNumber)
    
    # Create log values and levels for frequencyplot
    #weeklyFrequency=ma.log(weeklyFrequency)
    ##weekfmin = weeklyFrequency.min()
    ##weekfmax = weeklyFrequency.max()
    ##levels = np.arange(weekfmin,weekfmax,(weekfmax-weekfmin)/10)
    #print(levels)
    return gaussian_filter(weeklyFrequency, sigma),nWeeksWithData

def plotDistribution(kelpData,week,baseout,xii,yii,distType,polygons):
    print("Plotting the distributions for week: %s"%(week))
    plt.clf()
    plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,resolution='h', 
                    projection='merc')

    x, y = mymap(xii, yii)
    #levels=np.arange(np.min(kelpData),np.max(kelpData),10)
    levels = np.arange(1,100,1)
    cblevels = np.arange(10,110,10)  
    z = np.fliplr(np.rot90(kelpData,3))       
    CS1 = mymap.contourf(x,y,z,levels,cmap=cm.get_cmap('Spectral_r',len(levels)-1), extend='max', alpha=1.0)
    CS2 = mymap.contour(x,y,z,cblevels, colors = 'k', linestyle = '--',linewidths = 0.1)

    plt.colorbar(CS1,orientation='vertical',extend='max', shrink=0.5,ticks = cblevels)

    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey',zorder=2)
    #mymap.drawcountries()
    #mymap.drawmapboundary()

    plt.title('Kelp week: {}-{} , polygons: {}-{}'.format(min(week),max(week),min(polygons),max(polygons)))
    if distType == "integrated":
        plotfile=baseout+'/Kelp_distribution_all_fullperiod.png'
    else:
        plotfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'polygons_{}-{}'.format(min(polygons),max(polygons))+'.png'
    print("=> Creating plot %s"%(plotfile))             
    plt.savefig(plotfile,dpi=300)
 
def main(shapefile,experiment,distType = "integrated",polygons = None,requiredResolution = 1.):
    # Which species to calculate for
    # The timespan is part of the filename
    ranges = {1:['01052016','01082016'],2:['01032016','15052016'],
              3:['01112015','01042016'],4:['01052016','01082016']}
    startdate,enddate = ranges[experiment]

    # Create the grid you want to calculate frequency on
    # The resolution of the output grid in km between each binned box
    xi,yi = createBins(requiredResolution)

    xii,yii=np.meshgrid(xi[:-1],yi[:-1])

    if polygons == None:
        s = geopandas.read_file(shapefile)
        polygons = s.index.values[:-2] + 1

    ###s = ogr.Open(shapefile)
    #for layer in s:
    #    polygons = [x + 1 for x in range(layer.GetFeatureCount()-1)]
    #    print ('layer',layer.GetFeatureCount())
    #   # polygons=[1,2]

    allData=np.zeros((len(polygons),weeksInYear,len(xi)-1,len(yi)-1))

    print("=> Created final array for all data of size :",np.shape(allData))
    '''for polygonIndex, polygon in enumerate(polygons): '''   
    import xarray as xr
    infiles = [base+'/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, 
                startdate, enddate) for polygonIndex, polygon in enumerate(polygons)] 

    for index,infile in enumerate(infiles):
    #for polygonIndex, polygon in enumerate(polygons):
    #    infile = base+'/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, startdate, enddate)
    #    
        if os.path.exists(infile):
            print("=> Opening input file: %s"%(infile))
            cdf = xr.open_dataset(infile) 
            #cdf = Dataset(infile)
            allData[index,:,:,:],nWeeksWithData = calculateAreaAverages(xi,yi,cdf,weeksInYear)
            cdf.close()
        else:
            print("==>> Input file %s could not be opened"%(infile))

    if distType == "weekly": 
        for week in range(1,weeksInYear,1):
            # Calculate the cumulative distribution for each week
            kelpData=np.squeeze(np.sum(allData[:,week,:,:], axis=0))

            levels=np.arange(np.ma.min(kelpData),np.ma.max(kelpData),0.5)        
            if (len(levels)>2 and np.ma.mean(kelpData)>0):
                print("Kelp data {}".format(np.shape(kelpData)))
                plotDistribution(kelpData,week,baseout,xii,yii,"weekly",polygons)

    elif distType == "integrated":     
        week = '{}{}'.format(nWeeksWithData[0],nWeeksWithData[-1])       
        # Plot the distribution for all weeks
        kelpData=np.squeeze(np.ma.sum(allData[:,:,:,:], axis=0))
        print("LAST Kelp data {}".format(np.shape(kelpData)))
        kelpData=np.ma.masked_where(kelpData==0,kelpData)
        kelpData=np.squeeze(np.ma.sum(kelpData[:,:,:], axis=0))
        plotDistribution(kelpData,week,baseout,xii,yii,"integrated",polygons)

if __name__ == "__main__":

    # Results and storage folders
    global base,baseout,weeksInYear,sigma 
    base= r'Data'
    baseout='distributionFigures'
    shapefile = 'Shapefile05112018/kelpExPol_11sept2018.shp'
    if not os.path.exists(baseout): os.makedirs(baseout)   

    weeksInYear = 52
    #Sigma is used in gaussian filter 
    # gaussian_filter(weeklyFrequency, sigma) 
    # this depends on how noisy your data is, play with it!
    sigma = 0.2 
 
    start_time = time.time()
    main(shapefile,experiment = 1, distType ="integrated",polygons = None,requiredResolution = 1)    
    #for experiment in range(1,10): 
    #    main(shapefile,experiment = experiment, distType ="integrated")
    ## ---  It took 409.14820981025696 seconds to run the script ---
    # ---  It took 401.2906005382538 seconds to run the script --- Withough prints 
    # ---  It took 341.74021339416504 seconds to run the script --- resolution ,'i',
    #---  It took 468.2105269432068 seconds to run the script --- xarray made it slowlier
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))     