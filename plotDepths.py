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

def find_depth(data):
    # ! find first non nan at first and cut the rest 
    data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    # find differences between floor depth and particle depth for each trajectory
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data

def get_groups(new_df,p_part):
    d = new_df.where(new_df.plantpart == p_part,drop = True)
    # apply method to each trajectory (particle release event)
    return d.groupby(d.trajectory).apply(find_depth)

def plt_part(df,p_part,col,axis,norm):
    d = get_groups(new_df = df, p_part = p_part)
    not_deposited = 0
    parts =  range(0,len(d.trajectory)-1)
    sed_depths = [utils.get_sed_depth(d,n) for n in parts if utils.is_sedimented(d,n)]  

    for n in parts: # loop over trajectories 
        if utils.is_sedimented(d,n):
            start = utils.get_start(d,n)  
            sed = utils.get_sed(d,n) 
            if start != sed:  
                if norm == True:
                    x = np.arange(0,sed+1-start)
                    z =  d.z[n][start:sed+1]    
                    axis.plot(x,z,'-', color = col,markersize = 1,linewidth = 0.5,alpha = 1,zorder = 9)       
                    axis.plot(sed-start,d.z[n][sed],'ko', markersize = 1,zorder = 10)      
                elif norm == False:      
                    x = d.time[start:sed+1]
                    z = d.z[n][start:sed+1]    
                    axis.plot(x,z,'-', color = col,markersize = 1,linewidth = 0.5,alpha = 1,zorder = 9)       
                    axis.plot(d.time[sed],d.z[n][sed],'ko', markersize = 1,zorder = 10)                     
        else : 
            # if particle is not deposited or all values are nan      
            not_deposited += 1   
    if norm == True:
        axis.set_title('Distibution of particles (type {}), normalized by time'.format(p_part))    
    elif norm == False: 
        axis.set_title('Distibution of particles (type {})'.format(p_part))   
        axis.set_xlabel('Days since the release') 
        frmt = '%b/%d'             
        axis.xaxis.set_major_formatter(mdates.DateFormatter(frmt))   
    axis.set_ylabel('Depth, m')
    axis.set_ylim(350,0)    
    #axis.set_xlim(0,650)       
    return sed_depths

def call_make_plot_mf(paths,experiment,normalize):

    fig = plt.figure(figsize=(11.69 , 8.27), dpi=100,
                                            facecolor='white')
    gs = gridspec.GridSpec(3,2, width_ratios=[3, 1])
    gs.update(left=0.08, right=0.98 ,top = 0.96,bottom = 0.08,
            wspace=0.13 ,hspace=0.37) 

    ax1 = fig.add_subplot(gs[0])
    ax1_1 = fig.add_subplot(gs[1]) 

    ax2 = fig.add_subplot(gs[2])
    ax2_1= fig.add_subplot(gs[3])

    ax3 = fig.add_subplot(gs[4])
    ax3_1 = fig.add_subplot(gs[5])

    sed_depths1,sed_depths2,sed_depths4 = [],[],[]
    for path in paths:
        df = xr.open_dataset(path)
        df['z'] = df['z'] * -1.

        s1 = plt_part(df,1,'#d65460',ax1,normalize) 
        s2 = plt_part(df,2,'g',ax2,normalize)
        s4 = plt_part(df,4,'#006080',ax3,normalize)        

        sed_depths1.extend(s1)
        sed_depths2.extend(s2)
        sed_depths4.extend(s4)
        df.close()  

    bins = np.arange(1,200,10)
    ax1_1.hist(sed_depths1,bins = bins,density = True,color = 'k')   
    ax2_1.hist(sed_depths2,bins = bins,density = True,color = 'k')   
    ax3_1.hist(sed_depths4,bins = bins,density = True,color = 'k')   

    for axis2 in (ax1_1,ax2_1,ax3_1):  
        axis2.set_title('Sedimentation depths')    
        axis2.set_xlim(0,200)  
    if normalize == True:
        plt.savefig('Figures/Kelp_trajectories_and_sedimentation_norm.png',format = 'png')
    else:
        plt.savefig('Figures/Kelp_trajectories_and_sedimentation.png',format = 'png')            
    #plt.show()


if __name__ == '__main__':
    start_time = time.time()
    polygons = utils.get_polygons()
    experiments = (1,2,3,4,5)
    paths = utils.get_paths(polygons,experiment = 1)[:1]
    # for all experiment: 
    # p = []
    #for exp in experiments:
    #    p.append(utils.get_paths(polygons,exp))
    #paths = np.concatenate(p)

    #call_make_plot_mf(paths,'all',normalize = True)    
    call_make_plot_mf(paths,experiment = 1,normalize = False) 
    #call_make_plot_mf(paths,experiment = 1,normalize = True) 

        #for exp in experiments:
    #    p.append(utils.get_paths(polygons,exp))
    #paths = np.concatenate(p)

    '''paths = [utils.get_paths(polygons,exp)[0] for exp in (1,2,3,4,5)]
    p = [utils.get_paths(polygons,exp) for exp in (1,2,3,4,5)]
    paths2 = [val for sublist in p for val in sublist]
    
    [x[0] for x in p]

    print (len(paths),'------', len(p),'\n',len(paths2))'''
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))          