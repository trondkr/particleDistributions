# coding=utf-8

import os, sys
import numpy as np
import glob
import matplotlib
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import xarray as xr
from datetime import datetime
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
#import ogr
#import osr
import time
import utils

__author__   = 'Trond Kristiansen'
__email__    = 'me (at) trondkristiansen.com'
__created__  = datetime(2016, 8, 10)
__modified__ = datetime(2016, 8, 10)
__version__  = "1.0"
__status__   = "Production"

global base,baseout,weeksInYear,sigma 
global kelpType
weeksInYear = 52

##Parameters: 
sigma = 0.2     
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

    return xi,yi,ngridx,ngridy

def calculateAreaAverages(xi,yi,cdf,weeksInYear):

    print('func: calculateAreaAverages() => Calculating averages within bins')
    print('=> binned domain (%2.1f,%2.1f) to (%2.1f,%2.1f)'%(np.min(xi),np.min(yi),np.max(xi),np.max(yi)))

    for week in range(1,weeksInYear):
        weekcdf = cdf.where(cdf.time.dt.week == week,drop = True)
        
        #utils.filter_sedimented()
        #TODO:Filter by only sedimented particles
        #print (weekcdf)
        print ('week',week)
        break
        Xpos = weekcdf.coords['lon'].values.flatten()
        Ypos = weekcdf.coords['lat'].values.flatten()               

        # create average for one week
        H, xedges, yedges = np.histogram2d(Xpos, Ypos, bins=(xi, yi), normed=False)  
        if week == 1: 
            weeklyFrequency=np.zeros((weeksInYear,np.shape(H)[0],np.shape(H)[1]), dtype=np.float32)     
        weeklyFrequency[week,:,:] = H

    return gaussian_filter(weeklyFrequency, sigma)

def plotDistribution(kelpData,week,baseout,xii,yii,distType,polygons,experiment):
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
    print (np.min(kelpData))
    levels=np.linspace(1,np.max(z)/3,20)
    #cblevels = np.arange(0,10,5)    
    CS1 = mymap.contourf(x,y,z,vmin = 1,vmax = np.max(z)/2, levels = 100, cmap=cm.get_cmap('Spectral_r'), extend='max') 
    #levels,,len(levels)-1 Spectral_r #, alpha=1.0levels
    #CS2 = mymap.contour(x,y,z,[1], colors = 'k',linewidths = 0.1)

    plt.colorbar(CS1,orientation='vertical',extend='max', shrink=0.7) #,ticks = cblevels)

    mymap.drawmapboundary(fill_color='#677a7a')
    mymap.fillcontinents(color='#b8a67d',zorder=2)
    mymap.drawcoastlines()
    plt.title('Kelp week: {} , polygons: {}-{}, kelp Type {}, experiment {}'.format(
                            week,min(polygons),max(polygons),kelpType,experiment))
    if distType == "integrated":
        plotfile=baseout+'/Kelp_distribution_all_fullperiod_kelpType{}_experiment_{}.png'.format(kelpType,experiment)
    else:
        plotfile=baseout+'/Kelp_distribution_all_week_'+str(week)+'polygons_{}-{}'.format(min(polygons),max(polygons))+'kelpType_{}.png'.format(kelpType)
    print("=> Creating plot %s"%(plotfile))             
    plt.savefig(plotfile,dpi=300)
 
def main(shapefile,experiment,distType,polygons = None,requiredResolution = 1.,kelpType = None):
    # Which species to calculate for
    # The timespan is part of the filename
    startdate,enddate = utils.ranges[experiment]

    # Create the grid you want to calculate frequency on
    # The resolution of the output grid in km between each binned box
    xi,yi,ngridx,ngridy = createBins(requiredResolution)
    xii,yii=np.meshgrid(xi[:-1],yi[:-1])

    if polygons == None:
        polygons = utils.get_polygons()

    #Created final array for all data of size : (17, 52, 205, 333)
    # 17 polygons in shapefile , 52 dimension for weeks, xi,yi - bins 
    allData=np.zeros((len(polygons),weeksInYear,ngridx-1,ngridy-1))
    print("=> Created final array for all data of size :",np.shape(allData))
    infiles = utils.get_paths(polygons,experiment)
   
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
                plotDistribution(kelpData,week,baseout,xii,yii,"weekly",polygons,experiment)

    elif distType == "integrated":     
        week = 'Sum over all time range'  
        # Plot the distribution for all weeks
        # sum data along all polygons and all weeks. 
        kelpData=np.ma.sum(allData[:,:,:,:], axis=(0,1))     
        print("LAST Kelp data {}".format(np.shape(kelpData)))
        # if value is zero for all polygons,weeks the position is masked
        kelpData=np.ma.masked_where(kelpData==0,kelpData)
        plotDistribution(kelpData,week,baseout,xii,yii,"integrated",polygons,experiment)

if __name__ == "__main__":
    start_time = time.time()

    # Results and storage folders


    base= r'Data'
    baseout='distributionFigures'
    shapefile = 'Shapefile05112018/kelpExPol_11sept2018.shp'
    if not os.path.exists(baseout): os.makedirs(baseout)   


    kelpType = 1
    alltypes = (0,1,2,4)  
    experiments = (1,2,3,4,5)
    # kelpType is the parameter for chosing, the type of kelp (differs by density?)
    # available types are 0,1,2,4; -32767 is used for masking in the nc file ,
    #  None can be used for all types       

    polygons = None 
    #polygons = (1,2,3) # #None means all polygons, otherwise create a list of needed polygons 

    distType ="integrated"
    # integrated to plot all positions of particles during recorded time
    # 'weekly' to plot separately positions during all weeks in the recorded time range 

    #requiredResolution = 1   Resolution of the data, used in creating bins 
    main(shapefile,experiment = 3, distType = distType,polygons = polygons,
        requiredResolution = 1,kelpType = 1)

    '''for kelpType in alltypes:
        kelpType = kelpType
        [main(shapefile,experiment = e, distType = distType,polygons = polygons,
        requiredResolution = 1,kelpType = kelpType) for e in experiments]'''

    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   