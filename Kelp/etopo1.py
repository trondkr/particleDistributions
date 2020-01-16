
import os, sys, datetime, string
import numpy as np
from netCDF4 import Dataset
import numpy.ma as ma
from netCDF4 import Dataset, date2num,num2date

__author__ = 'Trond Kristiansen'
__email__ = 'trond.kristiansen@niva.no'
__created__ = datetime.datetime(2017, 5, 16)
__modified__ = datetime.datetime(2017, 5, 16)
__version__ = "1.0"
__status__ = "Development, modified on 16.05.2017"


def findSubsetIndices(min_lat,max_lat,min_lon,max_lon,lats,lons):

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
