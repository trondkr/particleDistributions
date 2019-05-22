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
np.seterr(divide='ignore', invalid='ignore')

def get_x_z_norm(d,n,start,sed):
    startdate = np.datetime64('2016-01-01T00:00:00.000000000')    
    dif = d.time[start].values - startdate
    x = d.time[start:sed+1] - dif  
    z =  d.z[n][start:sed+1]
    sed_time = d.time[sed] - dif 
    return x,z,sed_time

def get_x_z(d,n,start,sed):
    x = d.time[start:sed+1]
    z = d.z[n][start:sed+1]   
    sed_time = d.time[sed] 
    return x,z,sed_time

def make_plot(axis,d,n,start,sed):
    x,z,sed_time = get_x_z_norm(d,n,start,sed)  
    axis.plot(sed_time,d.z[n][sed],'ko', markersize = 1,zorder = 10)                    
    axis.plot(x,z,'-', color = col,markersize = 1,linewidth = 0.5,alpha = 1,zorder = 9)       

def func(data):
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z   
    data['start'] = np.isnan(data['dif_depth']) #np.amin(np.argwhere())

def plt_part(df,p_part,col,axis,norm):
    d1 = df.where(df.plantpart == p_part,drop = True)
    d = d1.groupby(d1.trajectory).apply(func)
    # def df grouped by trajectories and with dif_depth variable 
    #d = utils.get_groups(new_df = df, p_part = p_part)
    print (d)

    #not_deposited = 0
    #parts =  range(0,len(d.trajectory)-1)
    ####sed_depths = [utils.get_sed_depth(d,n) for n in parts if utils.is_sedimented(d,n)]  
    #s = [(utils.get_start_sed(d,n)) for n,val in enumerate(d.trajectory.values) if utils.is_sedimented(d,n)]
    #for n in s:
    #    print (len(n))
    #    break
    # filter start = stop
    #print (type(s))    
    #s = [n for n in s if n[0]!= n[1]] 

    #for n,val in enumerate(d.trajectory.values): # loop over trajectories 
    #    if utils.is_sedimented(d,n):
    #        start = utils.get_start(d,n)  
    #        sed = utils.get_sed(d,n) 
    #        if start != sed:
    #            make_plot(axis,d,n,start,sed)  
                                    
                        
    #else : 
    #    # if particle is not deposited or all values are nan      
    #    not_deposited += 1 
    # 
    # elif norm == False:      
    #            x,z,sed_time = get_x_z(d,n,start,sed)     
    print("---  It took %s seconds to Loop over trajectories ---" % (time.time() - start_time))   
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
    #print ('not_deposited',not_deposited)  
    #print ('sed_depths',sed_depths)      
    ####return sed_depths



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

    '''dfs = [utils.get_df(path) for path in paths]

    d1s = [utils.get_groups(new_df = df, p_part = 1) for df in dfs]  
    d2s = [utils.get_groups(new_df = df, p_part = 2) for df in dfs]  
    d4s = [utils.get_groups(new_df = df, p_part = 4) for df in dfs] 

    [plt_part(d1, p_part = 1,col = '#d65460',axis = ax1,norm = normalize) for d1 in d1s]
    [plt_part(d2, p_part = 2,col = 'g',      axis = ax2,norm = normalize) for d2 in d2s]
    [plt_part(d4, p_part = 4,col = '#006080',axis = ax3,norm = normalize) for d4 in d4s]  ''' 

    for path in paths:
        df = utils.get_df(path)
        xr.open_mfdataset('my/files/*.nc')
        s1 = plt_part(df,1,'#d65460',ax1,normalize) 
        s2 = plt_part(utils.get_groups(new_df = df, p_part = 2),2,'g',ax2,normalize)
        s4 = plt_part(utils.get_groups(new_df = df, p_part = 4),4,'#006080',ax3,normalize)        

        #sed_depths1.extend(s1)
        #sed_depths2.extend(s2)
        #sed_depths4.extend(s4)
        df.close()  

    #bins = np.arange(1,200,10)
    #ax1_1.hist(sed_depths1,bins = bins,density = True,color = 'k')   
    #ax2_1.hist(sed_depths2,bins = bins,density = True,color = 'k')   
    #ax3_1.hist(sed_depths4,bins = bins,density = True,color = 'k')

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
    #---  It took 97.51012945175171 seconds git to run the script ---
    paths = utils.get_paths(polygons,experiment = 1)[:1]
    
    print("---  It took %s seconds to load paths  ---" % (time.time() - start_time))       
    # for all experiment: 
    # p = []
    #for exp in experiments:
    #    p.append(utils.get_paths(polygons,exp))
    #paths = np.concatenate(p)

    #call_make_plot_mf(paths,'all',normalize = True)    
    #call_make_plot_mf(paths,experiment = 1,normalize = False) 
    #1 hour
    call_make_plot_mf(paths,experiment = 1,normalize = True) 

        #for exp in experiments:
    #    p.append(utils.get_paths(polygons,exp))
    #paths = np.concatenate(p)

    '''paths = [utils.get_paths(polygons,exp)[0] for exp in (1,2,3,4,5)]
    p = [utils.get_paths(polygons,exp) for exp in (1,2,3,4,5)]
    paths2 = [val for sublist in p for val in sublist]
    
    [x[0] for x in p]

    print (len(paths),'------', len(p),'\n',len(paths2))'''
    print("---  It took %s seconds to run the script ---" % (time.time() - start_time))          