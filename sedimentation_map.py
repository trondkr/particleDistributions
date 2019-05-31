# coding=utf-8

import os, sys
import numpy as np
import numpy.ma as ma
import glob
import matplotlib
from matplotlib.pyplot import cm 
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import pandas as pd
import xarray as xr
from datetime import datetime
from netCDF4 import Dataset, date2num,num2date
from scipy.ndimage.filters import gaussian_filter
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
#import ogr
#import osr
import time
import utils


def find_depth(data):
    data['z'] = data['z'] * -1.
    # ! find first non nan at first and cut the rest 
    #data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    # find differences between floor depth and particle depth for each trajectory
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data

'''def get_groups(new_df,p_part):
    d = new_df.where(new_df.plantpart == p_part,drop = True)
    # apply method to each trajectory (particle release event)
    return d.groupby(d.trajectory).apply(find_depth)'''

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

def get_sed_pos(paths,kelpType): 
    with xr.open_mfdataset(paths,concat_dim='time') as ds: #
        df = ds.load()
        df = df.where(df.status > -1, drop = True)  
        df['z'] = df['z'] * -1.        

        # Function can plot either all kelp types or one 
        if kelpType != 'All':
            d = df.where(df.plantpart == kelpType,drop = True)     
        d['dif_depth'] =  d.sea_floor_depth_below_sea_level - d.z     
        grp = d.groupby('trajectory')
        loop  = [utils.get_latlon(d) for n,d in grp if utils.get_latlon(d) != None]
        lats = list(map(lambda x : x[0], loop))  
        lons = list(map(lambda x : x[1], loop)) 
    return lats,lons

def createBins(res):
    print('func: createBins() => Creating bins for averaging')
    xmin=15.0; xmax=21.0
    ymin=69.0; ymax=72.0

    deg2rad=np.pi/180.
    R = 6371  # radius of the earth in km
    # Distance from minimum to maximim longitude
    def get_dist(xx,yy): 
        return R * np.sqrt( xx**2 + yy**2 )

    def rad(var):
        return var*deg2rad

    xmaxrad = xmax*deg2rad
    xminrad = xmin*deg2rad
    ymaxrad = ymax*deg2rad
    yminrad = ymin*deg2rad   

    dx = get_dist((xmaxrad - xminrad) * np.cos(ymaxrad),0)
    # number of kilometers .Then if res is 1 we get 1 km resolution 
    print("Distance from minimum to maximim longitude binned area is %s km"%(dx))

    # Distance from minimum to maximim latitude
    dy = get_dist(0,ymaxrad-yminrad)   
    print("Distance from minimum to maximim latitude binned area is %s km"%(dy))

    ngridx = int(dx/res)
    ngridy = int(dy/res)

    db = 1.e-6 # bin padding     
    xi = np.linspace(xmin-db,xmax+db,ngridx)
    yi = np.linspace(ymin-db,ymax+db,ngridy)

    print('=> created binned array of domain of size (%s,%s) with resolution %s'%(ngridx,ngridy,res))

    return xi,yi #,ngridx,ngridy

def get_density(lats, lons,nlevels,cmap):
    # compute appropriate bins to chop up the data: 

    lon_bins,lat_bins = createBins(res = 1)   
    density, _, _ = np.histogram2d(lats, lons, [lat_bins, lon_bins])

    # mask 0 density if needed
    density = ma.masked_where(density == 0, density)

    levels = MaxNLocator(nbins=nlevels).tick_values(0,100)
    #norm .tick_values(density.min(),density.max())
   
    norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # Turn the lon/lat of the bins into 2 dimensional arrays 
    lon_bins_2d, lat_bins_2d = np.meshgrid(lon_bins, lat_bins)
  
    return lon_bins_2d,lat_bins_2d,density,norm

def make_map(paths,kelpType,type,experiment,polygons):
    mymap = create_map()
    lats,lons = get_sed_pos(paths,kelpType)

    if type == 'heatmap':
        nlevels = 50        
        cmap = plt.get_cmap('Spectral_r')        
        lon_bins_2d,lat_bins_2d,density,norm = get_density(lats, lons,nlevels,cmap)
        xs, ys = mymap(lon_bins_2d, lat_bins_2d) # will be plotted using pcolormesh          
        cs = mymap.pcolormesh(xs, ys, density, cmap=cmap, norm=norm)      
        plt.colorbar(cs, shrink=0.7)
        figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,str(experiment)) 
        #TODO:convert to tiff

    elif type == 'scatter':
        x,y = mymap(lons,lats)
        mymap.scatter(x,y,alpha = 0.5,c = 'k',s = 10)
        figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,experiment)         
    plt.savefig(figname,format = 'png',dpi = 300)

def call_make_map(kelpType,plot_type,experiments,polygons):
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
    call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = [1,2,3,4,5], polygons = [1]) # 'All') #,2,3,4,5

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment
    # takes 20 minutes to plot     
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   
    #call_make_map(kelpType = 2,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment    
    #call_make_map(kelpType = 'All',plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')
   
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   
