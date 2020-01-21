# coding=utf-8
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import utils
import animateScatter
import tracks
import xarray as xr
from mpl_toolkits.basemap import Basemap
#
#
# MAIN program for creating particle track plots and animations of kelp experiments
# Requires tracks.py, amimateScatter.py, laplacefilter.py, and mpl_util.py

def create_map():
    fig = plt.figure(figsize=(10,10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.0,resolution='h', 
                    projection='merc')

    #mymap.drawmapboundary(fill_color='#677a7a')
    mymap.fillcontinents(color='#b8a67d',zorder=2)
    mymap.drawcoastlines()
    return mymap

def to_plot(d,mymap,start,sed):
    # = d.dropna(dim = 'time')
    l = np.ma.masked_greater(d.lon.values[start:sed+1],100)  
    lons = pd.Series(l).dropna()

    l = np.ma.masked_greater(d.lat.values[start:sed+1],100)  
    lats = pd.Series(l).dropna()

    mymap.plot(lons.values,lats.values,linewidth = 0.5,
            latlon = True,alpha = 0.6) 

    mymap.plot(lons.values[-1],lats.values[-1],'ko', markersize = 1,
            latlon = True,alpha = 0.6)        
    mymap.plot(lons.values[0],lats.values[0],'ro', markersize = 1.3,
            latlon = True,alpha = 0.6)                               
    return mymap

def make_map(paths,kelpType,type,experiment,polygons):

    df = xr.open_mfdataset(paths,concat_dim='time')
    #df = df.where(df.lat > 100,drop = True)  
    df = df.where(df.status > -1, drop = True)   
    if kelpType != 'All':
        df = df.where(df.plantpart == kelpType,drop = True)  
    df['z'] = df['z'] * -1.    
    df['dif_depth'] =  df.sea_floor_depth_below_sea_level - df.z  

    mymap = create_map()
    grp = df.groupby('trajectory')
    #counter = 0
    for n,d in grp:
        #if counter > 500:
        #    break
        start,sed,sed_depth = utils.get_start_sed_depth(d)
        if sed != None and start != sed:
            to_plot(d,mymap,start,sed)
            #counter += 1

    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))  
    figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,experiment)
    plt.savefig(figname)
    #plt.show()  

def call_make_map(animation,kelpType,plot_type,experiments,polygons):
    if len(experiments) > 1:
        p = [utils.get_paths(polygons,exp)for exp in experiments]
        paths = [val for sublist in p for val in sublist]
    else:
        paths = utils.get_paths(experiment = experiments[0],polygons = polygons)
    
    make_map(paths,kelpType,plot_type,experiment = experiments,polygons = polygons)

if __name__ == "__main__":
    start_time = time.time()

    #Examples:
    # 1 ) 
    # takes  20 seconds   
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = [1], polygons = [1]) # 'All')

    # 2)
    #Will create map for Kelp 1, polygon 1, all experiment 
    animation=False
    if animation:
        plot_type='animate_particles'
    else:
        plot_type='plot_particletracks'

    call_make_map(animation,kelpType=[1],plot_type=plot_type,experiments = [1], polygons = [1]) #,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17]) # 'All')

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment
    # takes 20 minutes to plot     
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   
    #call_make_map(kelpType = 2,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment    
    #call_make_map(kelpType = 'All',plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')
   

