import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.pyplot import cm
from mpl_toolkits.basemap import Basemap
import sys
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pandas as pd
import geopandas as gpd
import shapely

class KelpPolygons(object):

    def __init__(self,mymap,ax,fill=False):
        self.shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/Shapefile05112018/kelpExPol_11sept2018.shp'
        self.mymap=mymap
        self.ax=ax
        self.fill=fill
        self.add_polygons_to_map()
       
    def getPolygons(self):
        #https://borealperspectives.wordpress.com/2016/03/07/plotting-polygon-shapefiles-on-a-matplotlib-basemap-with-geopandas-shapely-and-descartes/
        lsoas = gpd.read_file(self.shapefile)
        studyarea = shapely.geometry.box(-1., 60, 30., 80.)

        patches = []
        selection = lsoas[lsoas.geometry.intersects(studyarea)]
        for poly in selection.geometry:
            if poly.geom_type == 'Polygon':
                mpoly = shapely.ops.transform(mm, poly)
                patches.append(PolygonPatch(mpoly))
            elif poly.geom_type == 'MultiPolygon':
                for subpoly in poly:
                    mpoly = shapely.ops.transform(mm, poly)
                    patches.append(PolygonPatch(mpoly))
            else:
                print("poly, is neither a polygon nor a multi-polygon. Skipping it")

        return lsoas, patches


    def add_polygons_to_map(self):
        
        print("Adding polygons to plot from file {}".format(self.shapefile))
        lsoas, patches = self.getPolygons()

        print("Plotting {} patches".format(len(patches)))
        p = PatchCollection(patches,alpha=1.0, match_original=True,zorder=3)  # cmap=matplotlib.cm.RdYlBu,alpha=1.0,facecolor='red',lw=2.0,edgecolor='red',zorder=10)

        colors = 100 * np.random.rand(len(patches))
        p.set_array(np.array(colors))
        self.ax.add_collection(p)


class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self,lons,lats,z,times,figname):
     
        self.z = z
        self.times=times
        self.figname=figname
        Writer = animation.writers['ffmpeg']
        self.writer = Writer(fps=15, metadata=dict(artist='Trond.Kristiansen@niva.no'), bitrate=1800)

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.mymap = Basemap(llcrnrlon=16.0, llcrnrlat=69.4,urcrnrlon=19.0, urcrnrlat=70.2,resolution='f',projection='merc')
     #   self.mymap.shadedrelief()
        self.mymap.fillcontinents(color='#b8a67d',zorder=2)
        self.mymap.drawcoastlines()

        polygons = KelpPolygons(self.mymap,self.ax)

        self.x,self.y=self.mymap(lons,lats)
      
        # Then setup FuncAnimation.
        framelength=len(self.z[:,0])-1
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=20,frames=1000, #framelength, 
                                          init_func=self.setup_plot, blit=False)


    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        
        x, y, c = (np.c_[self.x[:,0], self.y[:,0], self.z[:,0]]).T
        self.scat = self.mymap.scatter(x, y, c=c, s=15,
                                    cmap="ocean", edgecolor="k")
        self.ttl = self.ax.text(1.5, 1.05, '', transform = self.ax.transAxes, va='center')

        divider = make_axes_locatable(self.ax)
        cax = divider.append_axes("right", size="5%", pad=0.05)

        plt.colorbar(self.scat, cax=cax)
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat

    def saveAnim(self):
        self.ani.save(self.figname,writer=self.writer)

    def update(self, i):
        """Update the scatter plot."""
        data = (np.c_[np.squeeze(self.x[:,i]), np.squeeze(self.y[:,i])])
        dateobject=self.times[i]
        d = pd.to_datetime(dateobject)
        plt.title('{}.{}.{} {}:{}'.format(d.year,d.month,d.day,d.hour,d.minute))

        print("Updating frame {} for data {}".format(i,np.shape(data)))
        # Set x and y data...
        self.scat.set_offsets(data)

        # Set sizes...
      #  self.scat.set_sizes(20)
        # Set colors..
        self.scat.set_array(self.z[:,i])

        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat


if __name__ == '__main__':
    a = AnimatedScatter("test.mp4")
   # plt.show()