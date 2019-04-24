import geopandas
import xarray as xr
import numpy as np

 
ranges = {1:['01052016','01082016'],2:['01032016','15052016'],
          3:['20112015','01042016'],4:['01052016','01082016'],
          5:['03082015','01082016']}

def get_polygons():
    shapefile = r'Shapefile05112018\kelpExPol_11sept2018.shp'
    s = geopandas.read_file(shapefile)
    return(s.index.values[:-2] + 1)

def get_paths(polygons = None,experiment = 1):
    if polygons == None:
        polygons = get_polygons()
    startdate,enddate = ranges[experiment]
    base= r'Data' 
    return [base+'/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, 
                startdate, enddate) for polygonIndex, polygon in enumerate(polygons)] 

if __name__ is '__main__':
    print (get_paths([1,2],experiment = 1))                