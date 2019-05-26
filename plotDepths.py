import netCDF4
import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.dates as mdates
import pandas as pd 
import xarray as xr
from pandas.plotting import register_matplotlib_converters
register_matplotlib_converters()
import seaborn as sns
sns.set() 
import utils 
import time
sed_crit = 0.1

'''def find_depth(data):
    # ! find first non nan at first and cut the rest 
    data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    # find differences between floor depth and particle depth for each trajectory
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data['dif_depth']'''

def plt_part(df,p_part,col,axis,norm):
    d = df.where(df.plantpart == p_part,drop = True)
    d['dif_depth'] =  d.sea_floor_depth_below_sea_level - d.z     

    grp = d.groupby('trajectory')
    parts = d.coords['trajectory'].values

    sed_depths = [utils.get_sed_depth2(d) for n,d in grp]
    sed_depths = list(filter(None.__ne__,sed_depths))    
    #starts = [utils.get_start_sed(d.isel(trajectory = n))[0] for n,val in enumerate(parts)]
    #seds =   [utils.get_start_sed(d.isel(trajectory = n))[1] for n,val in enumerate(parts)]


    startdate = np.datetime64('2000-01-01T00:00:00')

    for n,ds in grp:  # loop over trajectories 
        #ds = d.isel(trajectory = n)
        start,sed = utils.get_start_sed(ds) #starts[n],seds[n]
        if start != sed:   
            if norm == True:
                #lifetime = ds.time[stop].values - ds.time[start].values
                dif = ds.time[start].values-startdate
                x = ds.time[start:sed+1] - dif          
                z = ds.z[start:sed+1]
   
            elif norm == False:      
                x = d.time[start:sed+1]
                z = d.z[n][start:sed+1]    

            axis.plot(x,z,'-', color = col,linewidth = 0.3,alpha = 0.5,zorder = 9)             
            axis.plot(x[-1],z[-1].values,'ko', markersize = 0.5,zorder = 10)                                
        
    if norm == True:
        axis.set_title('Distibution of particles (type {}), normalized by time'.format(p_part))    
        frmt = '%M-%d'
    elif norm == False: 
        axis.set_title('Distibution of particles (type {})'.format(p_part))
        frmt = '%b/%d'             
    axis.xaxis.set_major_formatter(mdates.DateFormatter(frmt))   
    axis.set_ylabel('Depth, m')
    axis.set_xlabel('Month,day of the release')
    axis.set_ylim(350,0)    
    return sed_depths,

def plt_part_dif(df,p_part,col,axis):
    d = get_groups(new_df = df, p_part = p_part)

    parts = d.coords['trajectory'].values
    sed_depths = [utils.get_sed_depth2(d,n) for n,val in enumerate(parts)]  
    sed_depths = list(filter(None.__ne__,sed_depths))
    #startdate = np.datetime64('2000-01-01T00:00:00')

    # groupby trajectory ?
    for n,val in enumerate(parts): # loop over trajectories 
        ds = d.isel(trajectory = n)
        start,sed = utils.get_start_sed(n,ds)
        if start != sed:   
            '''if norm == True:
                dif = ds.time[start].values-startdate
                x = ds.time[start:sed+1] - dif
                z = ds.dif_depth[start:sed+1]                
                z = ds.z[start:sed+1]''' 
            x = d.time[start:sed+1]
            z = d.z[n][start:sed+1]    
            axis.plot(x,z,'-', color = col,linewidth = 0.3,alpha = 0.5,zorder = 9)             
            axis.plot(x[0],z[0],'ko', markersize = 0.5,zorder = 10)              
    #if norm == True:
    #    axis.set_title('Distibution of particles (type {}), normalized by time'.format(p_part))    
    #frmt = '%M-%d'
    #elif norm == False: 
    #    axis.set_title('Distibution of particles (type {})'.format(p_part))
    frmt = '%b/%d'             
    axis.xaxis.set_major_formatter(mdates.DateFormatter(frmt))   
    axis.set_ylabel('Depth, m')
    axis.set_xlabel('Month,day of the release')
    #axis.set_ylim(350,0)    
  
def call_make_plot_mf(paths,experiment,normalize):

    fig = plt.figure(figsize=(11.69 , 8.27), dpi=100,
                                            facecolor='white')
    gs = gridspec.GridSpec(3,2, width_ratios=[3, 1])
    gs.update(left=0.08, right=0.98 ,top = 0.96,bottom = 0.08,
            wspace=0.13 ,hspace=0.37) 

    ax1   = fig.add_subplot(gs[0])
    ax1_1 = fig.add_subplot(gs[1]) 
    ax2   = fig.add_subplot(gs[2])
    ax2_1 = fig.add_subplot(gs[3])
    ax3   = fig.add_subplot(gs[4])
    ax3_1 = fig.add_subplot(gs[5])

    with xr.open_mfdataset(paths,concat_dim='time') as ds: #
        df = ds.load()
    df = df.where(df.status > -1, drop = True)
    df['z'] = df['z'] * -1.

    sed_depths1 = plt_part(df,1,'#d65460',ax1,normalize) 
    sed_depths2 = plt_part(df,2,'g',ax2,normalize)
    sed_depths4 = plt_part(df,4,'#006080',ax3,normalize)        

    bins = np.arange(1,200,10)

    ax1_1.hist(sed_depths1,bins = bins,density = True,color = 'k')   
    ax2_1.hist(sed_depths2,bins = bins,density = True,color = 'k')   
    ax3_1.hist(sed_depths4,bins = bins,density = True,color = 'k')

    for axis2 in (ax1_1,ax2_1,ax3_1):  
        axis2.set_title('Sedimentation depths')    
        axis2.set_xlim(0,200)  
    #if normalize == True:
    #    plt.savefig('Figures/Kelp_trajectories_and_sedimentation_norm.png',format = 'png')
    #else:
    #    plt.savefig('Figures/Kelp_trajectories_and_sedimentation.png',format = 'png')            
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))        
    plt.show()

def call_make_plot_dif(paths,experiment):

    fig = plt.figure(figsize=(11.69 , 8.27), dpi=100,
                                            facecolor='white')
    gs = gridspec.GridSpec(3,1)
    gs.update(left=0.08, right=0.98 ,top = 0.96,bottom = 0.08,
            wspace=0.13 ,hspace=0.37) 

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])


    with xr.open_mfdataset(paths,concat_dim='trajectory') as ds:
        df = ds.load()

    df = df.where(df.status > -1, drop = True)
    df['z'] = df['z'] * -1.

    plt_part_dif(df,1,'#d65460',ax1) 
    plt_part_dif(df,2,'g',ax2)
    plt_part_dif(df,4,'#006080',ax3)        
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))  
    if normalize == True:
        plt.savefig('Figures/Kelp_trajectories_and_sedimentation_norm.png',format = 'png')
    else:
        plt.savefig('Figures/Kelp_trajectories_and_sedimentation.png',format = 'png')   
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))            
    #plt.show()

if __name__ == '__main__':
    start_time = time.time()
    polygons = utils.get_polygons()
    experiments = (1,2,3,4,5)
    paths = utils.get_paths(polygons,experiment = 1)[:5]
    # for all experiment: 
    # p = []
    #for exp in experiments:
    #    p.append(utils.get_paths(polygons,exp))
    #paths = np.concatenate(p)

    #call_make_plot_mf(paths,'all',normalize = True)    
    #call_make_plot_mf(paths,experiment = 1,normalize = False) 
    call_make_plot_mf(paths,experiment = 1,normalize = True) 
    #call_make_plot_dif(paths,experiment = 1) 
        
