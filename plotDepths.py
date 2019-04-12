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

specie='Kelp'
startdate='01052016'
enddate='01082016'
polygons=[1] #,2,3,4,5]
#densities=[0,1,2,3,4] 
experiment=1

df = xr.open_dataset(r'Kelp_polygon_1_experiment_1_01052016_to_01082016.nc' )
df['z'] = df['z'] * -1.

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

def plt_part(p_part,col,axis):

    d = get_groups(new_df = df, p_part = p_part)
    not_deposited = 0
    for n in range(0,len(d.trajectory)-1): 
        # find indices with nan, after first nan (sed) we assume that particle is sedimented and do not plot anything after 
        try:
            start = int(min(np.argwhere(np.invert(np.isnan(d.dif_depth[n].values)))))    
            sed = int(min(np.argwhere(d.dif_depth[n].values < 0.2)))
            #print ('sed,start',sed,start) 
            if start != sed:   
                axis.plot(d.time[start:sed+1],d.z[n][start:sed+1],'-', color = col,markersize = 1,linewidth = 0.5,alpha = 1,zorder = 9)       
                axis.plot(d.time[sed],d.z[n][sed],'ko', markersize = 3,zorder = 10)              
                axis.plot(d.time[start],d.z[n][start],'bo', markersize = 3,zorder = 10)      
        except : 
            # if particle is not deposited or all values are nan      
            not_deposited += 1                     
    #print ('not_deposited or wrong array',not_deposited)    
    axis.set_ylim(350,0)        
 
def plt_part_normilized(p_part,col,axis,axis2):
    d = get_groups(new_df = df, p_part = p_part)
    not_deposited = 0
    sed_depths = []
    for n in range(0,len(d.trajectory)-1): 
        try:
            # find first nonnan value 
            start = int(min(np.argwhere(np.invert(np.isnan(d.dif_depth[n].values)))))    
            # find sedimentation  time
            sed = int(min(np.argwhere(d.dif_depth[n].values < 0.2)))
            #print ('sed,start',sed,start) 
            sed_depths.append(d.z[n][sed])
            if start != sed:        
                axis.plot(np.arange(0,sed+1-start),d.z[n][start:sed+1],'-', color = col,markersize = 1,linewidth = 0.5,alpha = 1,zorder = 9)       
                axis.plot(sed-start,d.z[n][sed],'ko', markersize = 2,zorder = 10)                   
        except :    
            not_deposited += 1 

    axis2.hist(sed_depths,50,density = True)  
    axis.set_title('Distibution of particles (type {}), normalized by time'.format(p_part))      
    axis2.set_title('Sedimentation depths')
    axis.set_xlabel('Days since the release') 
    axis.set_ylabel('Depth, m')
    axis.set_ylim(350,0)    
    axis.set_xlim(0,650) 
    axis2.set_xlim(0,200)  

def call_make_plot(norm = False):
    fig = plt.figure(figsize=(11.69 , 8.27), dpi=100,
                                            facecolor='white')
    gs = gridspec.GridSpec(3,2, width_ratios=[3, 1])
    gs.update(left=0.08, right=0.98 ,top = 0.96,bottom = 0.08,
            wspace=0.1 ,hspace=0.35) 

    ax1 = fig.add_subplot(gs[0])
    ax1_1 = fig.add_subplot(gs[1]) 

    ax2 = fig.add_subplot(gs[2])
    ax2_1= fig.add_subplot(gs[3])

    ax3 = fig.add_subplot(gs[4])
    ax3_1 = fig.add_subplot(gs[5])

    if norm == True:
        plt_part_normilized(1,'#d65460',ax1,ax1_1)   
        plt_part_normilized(2,'g',ax2,ax2_1)
        plt_part_normilized(4,'#006080',ax3,ax3_1)

    else:            
        plt_part(1,'#d65460',ax1)   
        plt_part(2,'g',ax2)
        plt_part(4,'#006080',ax3)
        #plt_part(5,'k',ax1)

    plt.savefig('sed_depth.png',format = 'png')
    #plt.show()

if __name__ == '__main__':
   call_make_plot(norm = True)
   #call_make_plot(norm = False)
