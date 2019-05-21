import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
from matplotlib.pyplot import cm 
from mpl_toolkits.basemap import Basemap
import sys
import ogr
import osr
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
#import geopandas

class KelpPolygons(object):

    def __init__(self,mymap,ax):
        self.shapefile = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/Shapefile05112018/kelpExPol_11sept2018.shp'
        self.mymap=mymap
        self.ax=ax
        self.add_polygons_to_map()
        #df = geopandas.read_file(self.shapefile)
        #ax = df.plot(figsize=(10, 10), alpha=0.5, edgecolor='k')

        #print(df.head())

    def getPathForPolygon(self,ring):
        codes = []
        x = [ring.GetX(j) for j in range(ring.GetPointCount())]
        y = [ring.GetY(j) for j in range(ring.GetPointCount())]

        codes += [mpath.Path.MOVETO] + (len(x) - 1) * [mpath.Path.LINETO]

        pathX, pathY = self.mymap(x, y)
        mymappath = mpath.Path(np.column_stack((pathX, pathY)), codes)

        return mymappath


    def createPathsForPolygons(self):
        # Returns a list of all the patches (polygons) in the Shapefile
        mypatches = []
        s = ogr.Open(self.shapefile)

        for layer in s:

            # get projected spatial reference
            sr = layer.GetSpatialRef()
            # get geographic spatial reference
            geogr_sr = sr.CloneGeogCS()
            # define reprojection
            proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)

            polygons = [x + 1 for x in range(layer.GetFeatureCount() - 1)]
            for polygonIndex, polygon in enumerate(polygons):

                feature = layer.GetFeature(polygonIndex)
                geom = feature.GetGeometryRef()

                ring = geom.GetGeometryRef(0)
                geom.Transform(proj_to_geog)

                if ring.GetPointCount() > 3:
                    print("Looping over polygon index {} with {} points".format(polygonIndex, ring.GetPointCount()))
                    polygonPath = self.getPathForPolygon(ring)
                    path_patch = mpatches.PathPatch(polygonPath, lw=2, edgecolor="black", facecolor='red')
                    mypatches.append(path_patch)
            return mypatches

    def add_polygons_to_map(self):
        mypatches = []
        print("Adding polygons to plot from file {}".format(self.shapefile))
        mypatches = self.createPathsForPolygons()

        print("Plotting {} patches".format(len(mypatches)))
        p = PatchCollection(mypatches,
                            alpha=1.0)  # cmap=matplotlib.cm.RdYlBu,alpha=1.0,facecolor='red',lw=2.0,edgecolor='red',zorder=10)

        colors = 100 * np.random.rand(len(mypatches))
        p = PatchCollection(mypatches, alpha=0.7)
        p.set_array(np.array(colors))
        self.ax.add_collection(p)


class AnimatedScatter(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self,lons,lats,z,figname):
     
        self.z = z
        self.figname=figname
        Writer = animation.writers['ffmpeg']
        self.writer = Writer(fps=15, metadata=dict(artist='Trond.Kristiansen@niva.no'), bitrate=1800)

        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize=(10,10))
        self.mymap = Basemap(llcrnrlon=16.0, llcrnrlat=69.4,
                    urcrnrlon=19.0, urcrnrlat=70.2,resolution='h', 
                    projection='merc')
        self.mymap.shadedrelief()
        self.mymap.fillcontinents(color='#b8a67d',zorder=2)
        self.mymap.drawcoastlines()

        polygons = KelpPolygons(self.mymap,self.ax)

        self.x,self.y=self.mymap(lons,lats)
      
        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, interval=20,frames=len(lons), 
                                          init_func=self.setup_plot, blit=False)


    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        
        x, y, c = (np.c_[self.x[:,0], self.y[:,0], self.z[:,0]]).T
        self.scat = self.mymap.scatter(x, y, c=c, s=10,
                                    cmap="ocean", edgecolor="k")

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
    a = AnimatedScatter()
   # plt.show()