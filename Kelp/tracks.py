import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.pyplot import cm 
from mpl_toolkits.basemap import Basemap
import sys
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
import animateScatter
import xarray as xr
import laplacefilter
import mpl_util
import pandas as pd

class Tracks(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self,lons,lats,z,figname):
     
        self.z = z
        self.figname=figname
        self.llcrnrlon=16.0
        self.llcrnrlat=69.4
        self.urcrnrlon=19.0
        self.urcrnrlat=70.2

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.mymap = Basemap(llcrnrlon=self.llcrnrlon, llcrnrlat=self.llcrnrlat,
                    urcrnrlon=self.urcrnrlon, urcrnrlat=self.urcrnrlat,resolution='h', 
                    projection='merc')
        #self.mymap.shadedrelief(zorder=2)
        self.mymap.fillcontinents(color='#b8a67d',zorder=2)
        self.mymap.drawcoastlines(zorder=2)

        etopo = self.addBathymetry()
        polygons = animateScatter.KelpPolygons(self.mymap,self.ax)
       
        self.x,self.y=self.mymap(lons,lats)

    def to_plot(self,d,start,sed,traj):
        
        l = np.ma.masked_greater(d.lon.values[traj,start:sed+1],100)  
        lons = pd.Series(l).dropna()

        l = np.ma.masked_greater(d.lat.values[traj,start:sed+1],100)  
        lats = pd.Series(l).dropna()

        self.mymap.plot(lons.values,lats.values,linewidth = 0.2,
                latlon = True,alpha = 0.4, color='lightgray',zorder=5) 

        self.mymap.plot(lons.values[-1],lats.values[-1],'ko', markersize = 1.2,
                latlon = True,alpha = 0.7, zorder=10)        
        self.mymap.plot(lons.values[0],lats.values[0],'ro', markersize = 1.2,
                latlon = True,alpha = 0.7,zorder=10)                               
      
    def plot_tracks(self):

        for traj_index in range(len(self.x[:,0])):
            print("Plotting trajectory {}".format(traj_index))
            self.mymap.plot(self.x[traj_index,:],self.y[traj_index],c='k',alpha=0.5,linewidth=0.4)
        
        plt.savefig(self.figname,format = 'png',dpi = 300)

    def addBathymetry(self):

        etopo1='/Users/trondkr/Dropbox/NIVA/ETOPO1/ETOPO1_Ice_g_gmt4.grd'
        etopo=xr.open_dataset(etopo1)
        lats=etopo['y'].values
        lons=etopo['x'].values

        res = self.findSubsetIndices(self.llcrnrlat-5,self.urcrnrlat+5,self.llcrnrlon-40,self.urcrnrlon+10,lats,lons)
        res = [ int(x) for x in res ]
     
        bathy = etopo["z"][int(res[2]):int(res[3]),int(res[0]):int(res[1])].values
        bathySmoothed = laplacefilter.laplace_filter(bathy,M=None)
        lon,lat=np.meshgrid(lons[res[0]:res[1]],lats[res[2]:res[3]]) 
        
        levels=[-6000,-5000,-3000, -2000, -1500, -1000,-500, -400, -300, -250, -200, -150, -100, -75, -65, -50, -35, -25, -15, -10, -5, 0]
        x, y = self.mymap(lon,lat)
        etopo_contour = self.mymap.contourf(x,y,bathySmoothed,levels,
                       cmap=mpl_util.LevelColormap(levels,cmap=cm.Blues_r),
                       extend='upper',
                       alpha=1.0,
                       origin='lower',zorder=1)
    
    
    def findSubsetIndices(self,min_lat,max_lat,min_lon,max_lon,lats,lons):

        """Array to store the results returned from the function"""
        res=np.zeros((4),dtype=np.float64)
        minLon=min_lon; maxLon=max_lon

        distances1 = []; distances2 = []
        indices=[]; index=1

        for point in lats:
            s1 = max_lat-point # (vector subtract)
            s2 = min_lat-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        distances1 = []; distances2 = []; index=1

        for point in lons:
            s1 = maxLon-point # (vector subtract)
            s2 = minLon-point # (vector subtract)
            distances1.append((np.dot(s1, s1), point, index))
            distances2.append((np.dot(s2, s2), point, index-1))
            index=index+1

        distances1.sort()
        distances2.sort()
        indices.append(distances1[0])
        indices.append(distances2[0])

        """ Save final product: max_lat_indices,min_lat_indices,max_lon_indices,min_lon_indices"""
        minJ=indices[1][2]
        maxJ=indices[0][2]
        minI=indices[3][2]
        maxI=indices[2][2]

        res[0]=minI; res[1]=maxI; res[2]=minJ; res[3]=maxJ;
        return res
