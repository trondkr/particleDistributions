# coding=utf-8

import os, sys
import matplotlib

matplotlib.use('Agg')
import numpy as np
import glob
import string
from matplotlib.pyplot import cm
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap
import matplotlib.path as mpath
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
from pylab import *
import datetime

from netCDF4 import Dataset, datetime, date2num, num2date
from scipy.ndimage.filters import gaussian_filter
import ogr
import osr

__author__ = 'Trond Kristiansen'
__email__ = 'me (at) trondkristiansen.com'
__created__ = datetime(2017, 3, 7)
__modified__ = datetime(2017, 3, 7)
__version__ = "1.0"
__status__ = "Production"


# --------
# plotKelpPolygons.py
#
# This script reads the kelp forest distributions from shapefiles provided by Karen Filbee-Dexter
# --------

def plotDistribution(shapefileList):
    plt.clf()
    plt.figure(figsize=(10, 10), frameon=False)
    ax = plt.subplot(111)

    mymap = Basemap(llcrnrlon=17.45,
                    llcrnrlat=69.5,
                    urcrnrlon=18.3,
                    urcrnrlat=69.7,
                    resolution='f', projection='tmerc', lon_0=15, lat_0=65.)

    mymap.drawcoastlines()
    mymap.fillcontinents(color='grey', zorder=2)
    mymap.drawcountries()
    mymap.drawmapboundary()

    plt.title('Kelp forest 07.03.2017')

    mypatches = []
    for shapefile in shapefileList:
        print "Adding polygons to plot from file %s" % (shapefile)
        mypatches = createPathsForPolygons(shapefile, mymap)

    print "Plotting %s patches" % (len(mypatches))
    p = PatchCollection(mypatches,
                        alpha=1.0)  # cmap=matplotlib.cm.RdYlBu,alpha=1.0,facecolor='red',lw=2.0,edgecolor='red',zorder=10)

    colors = 100 * np.random.rand(len(mypatches))
    p = PatchCollection(mypatches, alpha=0.4)
    p.set_array(np.array(colors))
    ax.add_collection(p)

    plotfile = 'kelp.png'
    print "=> Creating plot %s" % (plotfile)
    plt.savefig(plotfile, dpi=300)


def getPathForPolygon(ring, mymap):
    codes = []
    x = [ring.GetX(j) for j in range(ring.GetPointCount())]
    y = [ring.GetY(j) for j in range(ring.GetPointCount())]

    codes += [mpath.Path.MOVETO] + (len(x) - 1) * [mpath.Path.LINETO]

    pathX, pathY = mymap(x, y)
    mymappath = mpath.Path(np.column_stack((pathX, pathY)), codes)

    return mymappath


def createPathsForPolygons(shapefile, mymap):
    # Returns a list of all the patches (polygons) in the Shapefile
    mypatches = []
    s = ogr.Open(shapefile)

    for layer in s:

        # get projected spatial reference
        sr = layer.GetSpatialRef()
        # get geographic spatial reference
        geogr_sr = sr.CloneGeogCS()
        # define reprojection
        proj_to_geog = osr.CoordinateTransformation(sr, geogr_sr)

        polygons = [x + 1 for x in xrange(layer.GetFeatureCount() - 1)]
        for polygonIndex, polygon in enumerate(polygons):

            feature = layer.GetFeature(polygonIndex)
            geom = feature.GetGeometryRef()

            ring = geom.GetGeometryRef(0)
            geom.Transform(proj_to_geog)

            if ring.GetPointCount() > 3:
                print "Looping over polygon index %s with %s points" % (polygonIndex, ring.GetPointCount())
                polygonPath = getPathForPolygon(ring, mymap)
                path_patch = mpatches.PathPatch(polygonPath, lw=2, edgecolor="purple", facecolor='red')
                mypatches.append(path_patch)
        return mypatches


def main():
    shapefileMerge = '/Users/trondkr/Dropbox/NIVA/KelpFloat/Kelp/Shapefile/ShapefilesHGU/kelpexpol_exp_grazed_combined.shp'
    shapefileList = [shapefileMerge]

    plotDistribution(shapefileList)


if __name__ == "__main__":
    main()
