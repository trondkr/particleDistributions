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
 
def plt_part_normilized(df,p_part,col,axis,axis2):
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
                axis.plot(sed-start,d.z[n][sed],'ko', markersize = 1,zorder = 10)                   
        except :    
            not_deposited += 1 

    axis.set_title('Distibution of particles (type {}), normalized by time'.format(p_part))      
    axis.set_xlabel('Days since the release') 
    axis.set_ylabel('Depth, m')
    axis.set_ylim(350,0)    
    axis.set_xlim(0,650) 

    return sed_depths

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

def call_make_plot_mf(paths,experiment):

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

        s1 = plt_part_normilized(df,1,'#d65460',ax1,ax1_1) 
        s2 = plt_part_normilized(df,2,'g',ax2,ax2_1)
        s4 = plt_part_normilized(df,4,'#006080',ax3,ax3_1)        
        sed_depths1.extend(s1)
        sed_depths2.extend(s2)
        sed_depths4.extend(s4)
        df.close()
    ax1_1.hist(sed_depths1,bins = np.arange(1,200,10),density = True,color = 'k')   
    ax2_1.hist(sed_depths2,50,density = True,color = 'k')   
    ax3_1.hist(sed_depths4,50,density = True,color = 'k')   

    for axis2 in (ax1_1,ax2_1,ax3_1):  
        axis2.set_title('Sedimentation depths')    
        axis2.set_xlim(0,200)  

    #plt.savefig('sed_depth.png',format = 'png')
    plt.show()


if __name__ == '__main__':

    polygons=(1,2) #,2,3,4,5]
    experiment = 1
    paths = utils.get_paths(polygons,experiment)
    call_make_plot_mf(paths,experiment)
   #call_make_plot(norm = False)
