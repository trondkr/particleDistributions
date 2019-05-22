import geopandas
import xarray as xr
import numpy as np

sed_criteria = 0.1

ranges = {1:['01052016','01082016'],2:['01032016','15052016'],
          3:['20112015','01042016'],4:['01052016','01082016'],
          5:['03082015','01082016']}

def get_polygons():
    shapefile = r'Shapefile05112018\kelpExPol_11sept2018.shp'
    s = geopandas.read_file(shapefile)
    return(s.index.values[:-2] + 1)

def get_paths(polygons = None,experiment = 1):
    if polygons is None:
        polygons = get_polygons()
    elif polygons is 'All':
        polygons = get_polygons()        
    startdate,enddate = ranges[experiment]
    base= r'Data' 
    return [base+'/Kelp_polygon_%s_experiment_%s_%s_to_%s.nc' % (polygon, experiment, 
                startdate, enddate) for polygonIndex, polygon in enumerate(polygons)] 

def is_sedimented(d,n):
    # check if particle was sedimented
    d = d.fillna(999) 
    sed = np.argwhere(d.dif_depth[n].values < sed_criteria)
    if sed.size == 0:
        return False
    else: 
        return True 

def get_start(d,n):
    # find index of the release event 
    return int(min(np.argwhere(np.invert(np.isnan(d.dif_depth[n].values))))) 

def get_sed(d,n):  
    # find index in array of sedimentations time 
    # first time when difference between seafloor 
    # and particle depth is < 0.1 m (sed_criteria)
    #d = d.fillna(999) 
    return int(min(np.argwhere(d.dif_depth[n].values < sed_criteria)))

def get_start_sed(d,n):
    # find index of the release event 
    try:
        start = np.amin(np.argwhere(np.invert(np.isnan(d.dif_depth[n].values)))) #int(
        d = d.fillna(999) 
        sed = np.amin(np.argwhere(d.dif_depth[n].values < sed_criteria))
        #sed1 = np.argwhere(d.dif_depth[n].values < sed_criteria)[0]
        #print (sed,sed1) 
        return start,sed 
    except:
        pass    

def get_df(path):
    df = xr.open_dataset(path)
    df = df.where(df.status > -1, drop = True)
    df['z'] = df['z'] * -1.
    return df

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

def get_lat(d,n):
    # find lat of the particle at the sedimentation time 
    return d.lat[n][get_sed(d,n)]

def get_lon(d,n):
    # find lat of the particle at the sedimentation time 
    return d.lon[n][get_sed(d,n)]

def get_sed_depth(d,n):
    return d.z[n][get_sed(d,n)].values

def get_latlon(d,n):
    # fin lat of the particle at the sedimentation time 
    s = get_sed(d,n)
    return (d.lat[n][s].values,d.lon[n][s].values)

def get_dif_depth(data):
    # find differences between floor depth and particle depth for each trajectory     
    # ! find first non nan at first and cut the rest 
    data = data.where(data.z != 'nan')
    data = data.where(data.z != np.nan)
    data = data.where(data.sea_floor_depth_below_sea_level != 'nan',drop = True)
    data['dif_depth'] =  data.sea_floor_depth_below_sea_level - data.z 
    return data

def get_df(path):
    df = xr.open_dataset(path)
    df['z'] = df['z'] * -1.
    return df.where(df.status > -1, drop = True)

if __name__ is '__main__':
    print (get_paths([1,2],experiment = 1))                