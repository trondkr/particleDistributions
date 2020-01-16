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
import sys

#
# MAIN program for creating particle track plots and animations of kelp experiments
# Requires tracks.py, amimateScatter.py, laplacefilter.py, and mpl_util.py
def find_depth(data):
    data['z'] = data['z'] * -1.
    # ! find first non nan at first and cut the rest 
    #data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    # find differences between floor depth and particle depth for each trajectory
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data

def get_groups(new_df,p_part):
    d = new_df.where(new_df.plantpart == p_part,drop = True)
    # apply method to each trajectory (particle release event)
    return d.groupby(d.trajectory).apply(find_depth)

def get_pos_function_of_time(paths,kelpType):
 
    df = xr.open_mfdataset(paths,concat_dim='trajectory')
    # Function can plot either all kelp types or one 
    if kelpType != 'All':
        df = df.where(df.plantpart == kelpType,drop = True)     
    df = df.where(df.status > -1, drop = True)
    d = df.groupby(df.trajectory).apply(find_depth)
    
    return d

def make_map(paths,kelpType,type,experiment,polygons):  
    
    first=True
    
    for polygon,path in enumerate(paths):
     
        df = xr.open_dataset(path)
        #df = df.where(df.lat > 100,drop = True)  
        df = df.where(df.status > -1, drop = True)   
    
        if kelpType != 'All':
            df = df.where(df.plantpart == kelpType,drop = True)  
        df['z'] = df['z'] * -1.    
        df['dif_depth'] =  df.sea_floor_depth_below_sea_level - df.z  
        lats=df['lat'][:].values
        lons=df['lon'][:].values
        times=df['time'][:].values
        print(times)
        z=df['z'][:].values

        if type == 'animate_particles':
            figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.mp4'.format(type,kelpType,polygons,experiment)
            anim = animateScatter.AnimatedScatter(lons,lats,z,times,figname)
            anim.saveAnim()

        if type=='plot_particletracks':
            figname = r'{}_for_kelp_type_{}_polygons_{}_experiment_{}.png'.format(type,kelpType,polygons,experiment)
            
            if first:
                trk = tracks.Tracks(lons,lats,z,figname)
                first=False
            
            diff = df.dif_depth.values
            for traj in range(np.shape(lons)[0]):
                message="Plotting trajectory {} polygon number {} ({})".format(traj,polygon,path)
                print(message, end='\r')
                sys.stdout.flush()
                
                start,sed,sed_depth = utils.get_start_sed_depth(diff,z,traj)
              
                if sed != None and start != sed:
                    trk.to_plot(df,start,sed,traj)
        print("")
    print("Saving figure {}".format(trk.figname))
    plt.savefig(trk.figname,format = 'png',dpi = 300)
    first=True


def call_make_map(animation,kelpType,plot_type,experiments,polygons):
   
    paths = utils.get_paths(experiment = experiments,polygons = polygons) 
    make_map(paths,kelpType,plot_type,experiment = experiments,polygons = polygons)
   
if __name__ == "__main__":
    start_time = time.time()

    #Examples:
    # 1 ) 
    # takes  20 seconds   
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = [1], polygons = [1]) # 'All')

    # 2)
    #Will create map for Kelp 1, polygon 1, all experiment 
    animation=True
    if animation:
        plot_type='animate_particles'
    else:
        plot_type='plot_particletracks'
    experiments=[5]
    for exp in experiments:
        call_make_map(animation,kelpType='All',plot_type=plot_type,experiments = exp, polygons = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18]) # 'All')

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment
    # takes 20 minutes to plot     
    #call_make_map(kelpType = 1,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   
    #call_make_map(kelpType = 2,plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')   

    # 3)
    #Will create map for Kelp 1, polygon 1, all experiment    
    #call_make_map(kelpType = 'All',plot_type = 'heatmap',experiments = (1,2,3,4,5), polygons = 'All')
   
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))   
