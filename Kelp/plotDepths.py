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
import glob
def plt_part(df,p_part,col,axis,norm):

    startdate = np.datetime64('2000-01-01T00:00:00')    
    d = df.where(df.plantpart == p_part,drop = True)
    d['dif_depth'] =  d.sea_floor_depth_below_sea_level - d.z     

    grp = d.groupby('trajectory')
    loop  = [[utils.get_start_sed_depth(d),n,d] for n,d in grp if utils.get_start_sed_depth(d) != (None,None,None)]

    s  =            list(map(lambda x : x[0], loop))  
    trajectories  = list(map(lambda x : x[1], loop))   
    ds_all  =       list(map(lambda x : x[2], loop))

    starts =     list(map(lambda x : x[0], s))  
    seds =       list(map(lambda x : x[1], s)) 
    sed_depths = list(map(lambda x : x[2], s))    

    for k,ds in enumerate(ds_all):  # loop over trajectories 
        start,sed = starts[k],seds[k]
        if start != sed:   
            if norm == True:
                #lifetime = ds.time[stop].values - ds.time[start].values
                dif = ds.time[start]-startdate
                x = ds.time[start:sed+1] - dif          
                z = ds.z[start:sed+1]
            elif norm == False:      
                x = d.time[start:sed+1]
                z = d.z[start:sed+1]    
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
    import datetime
    axis.set_xlim(startdate,'2000-02-15T00:00:00') 
    return sed_depths

  
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
    #df = xr.open_mfdataset(paths,concat_dim='time')
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
    #    plt.savefig('Figures/Kelp_trajectories_and_sedimentation_norm_experiment{}.png'.format(experiment),format = 'png')
    #else:
    #    plt.savefig('Figures/Kelp_trajectories_and_sedimentation.png',format = 'png')            
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))        
    plt.show()

if __name__ == '__main__':
    start_time = time.time()
    pol = utils.get_polygons()
    experiments = (1,2,3,4,5)
    #paths = utils.get_paths(polygons,experiment = 1)
    allpaths = (glob.glob("Data/*.nc"))
    
    #for exp in experiments:
    #    call_make_plot_mf(utils.get_paths(pol,experiment = exp),experiment = exp,normalize =True) 
    #call_make_plot_mf(utils.get_paths(pol,experiment = 2),experiment = 2,normalize =True)

    call_make_plot_mf(utils.get_paths(pol,experiment = 5),experiment = 5,normalize =True)   