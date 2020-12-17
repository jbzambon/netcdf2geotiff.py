# Program to convert netCDF output to GeoTIFF
#
# Joseph B. Zambon
#  jbzambon@ncsu.edu
#  19 March 2020

# conda install -c anaconda netcdf4
# conda install -c conda-forge gdal
# conda install -c anaconda scikit-learn

from netCDF4 import Dataset
import numpy as np
import os
from osgeo import osr, gdal
import time

yyyymmdd = input()

domains = np.array([-98.0,    -80.0,   18.0, 31.0])
names = ['GoM']
dimensions = np.array([1660, 1314])
#Input file info
ncfile = '/home/ubuntu/geotiff/jbz_' + yyyymmdd + '.nc'
ncfile = Dataset(ncfile,'r')

lat    = ncfile.variables['lat'][:]
lon    = ncfile.variables['lon'][:]

frames = [0, 1]

for frame in frames:
    sst    = ncfile.variables['sst'][frame,:,:]
    chl    = ncfile.variables['chl'][frame,:,:]
    for domain in range(0,np.size(names)):
        # jbz SST
        lat_index = np.where((lat[:]-domains[domain,2]>=0) & (lat[:]-domains[domain,3]<=0))
        lon_index = np.where((lon[:]-domains[domain,0]>=0) & (lon[:]-domains[domain,1]<=0))
        l_lat = lat_index[0][0]
        u_lat = lat_index[0][-1]
        l_lon = lon_index[0][0]
        u_lon = lon_index[0][-1]
        jbz_sst = sst[l_lat:u_lat,l_lon:u_lon]
        if frame == 0:
            output_file = '/var/www/html/' + folders[domain] + '/jbz-SST_' + names[domain] + '_' + yyyymmdd + '_F16.tif'
        elif frame == 1:
            output_file = '/var/www/html/' + folders[domain] + '/jbz-SST_' + names[domain] + '_' + yyyymmdd + '_F34.tif'
        # Create gtif
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create('sst_out.tif', np.shape(jbz_sst[0,:])[0], np.shape(jbz_sst[:,0])[0], bands=1, eType=gdal.GDT_Byte )
        raster = jbz_sst
        raster[raster>33.8]=0
        raster[raster<0]=0
        raster=raster/0.13286   # converting sst in degree to 0-255 pixels
        raster=abs(raster-255)  # invert color, so lower (higher) temperature has brighter (darker) color
        raster = np.flipud(raster)
            # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform( [ lon[lon_index[0][0]], 0.01, 0, lat[lat_index[0][-1]], 0, -0.01 ] )
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS72")
        srs.SetProjParm("false_easting",0.0),
        srs.SetProjParm("false_northing",0.0),
        dst_ds.SetProjection( srs.ExportToWkt() )
            # write the band
        dst_ds.GetRasterBand(1).WriteArray(raster)
        del(dst_ds)
        os.system(''.join(['gdalwarp -s_srs /home/ubuntu/geotiff/src.prj -t_srs /home/ubuntu/geotiff/roffs_projections/' + \
                            projection_files[domain] + ' sst_out.tif sst_warp.tif']))
        os.system(''.join(['convert -resize ' + np.str(dimensions[domain,0]) + 'x' + np.str(dimensions[domain,1]) + '! sst_warp.tif ' + output_file]))
        os.system('rm -rf sst_out.tif sst_warp.tif')
 
        # jbz CHL-a
        lat_index = np.where((lat[:]-domains[domain,2]>=0) & (lat[:]-domains[domain,3]<=0))
        lon_index = np.where((lon[:]-domains[domain,0]>=0) & (lon[:]-domains[domain,1]<=0))
        l_lat = lat_index[0][0]
        u_lat = lat_index[0][-1]
        l_lon = lon_index[0][0]
        u_lon = lon_index[0][-1]
        jbz_chl = chl[l_lat:u_lat,l_lon:u_lon]
        if frame == 0:
            output_file = '/var/www/html/' + folders[domain] + '/jbz-CHL_' + names[domain] + '_' + yyyymmdd + '_F16.tif'
        elif frame == 1:
            output_file = '/var/www/html/' + folders[domain] + '/jbz-CHL_' + names[domain] + '_' + yyyymmdd + '_F34.tif'
        # Create gtif
        driver = gdal.GetDriverByName("GTiff")
        dst_ds = driver.Create('chl_out.tif', np.shape(jbz_chl[0,:])[0], np.shape(jbz_chl[:,0])[0], bands=1, eType=gdal.GDT_Byte )
        raster = jbz_chl
        raster = 10**jbz_chl
        raster[raster>20]=20
        raster[raster<=0]=0;
        raster= np.log(raster+1)/np.log(10)/0.00519;
        raster = np.flipud(raster)
            # top left x, w-e pixel resolution, rotation, top left y, rotation, n-s pixel resolution
        dst_ds.SetGeoTransform( [ lon[lon_index[0][0]], 0.01, 0, lat[lat_index[0][-1]], 0, -0.01 ] )
        srs = osr.SpatialReference()
        srs.SetWellKnownGeogCS("WGS72")
        srs.SetProjParm("false_easting",0.0),
        srs.SetProjParm("false_northing",0.0),
        dst_ds.SetProjection( srs.ExportToWkt() )
            # write the band
        dst_ds.GetRasterBand(1).WriteArray(raster)
        del(dst_ds)
        os.system(''.join(['gdalwarp -s_srs /home/ubuntu/geotiff/src.prj -t_srs /home/ubuntu/geotiff/roffs_projections/' + \
                            projection_files[domain] + ' chl_out.tif chl_warp.tif']))
        os.system(''.join(['convert -resize ' + np.str(dimensions[domain,0]) + 'x' + np.str(dimensions[domain,1]) + '! chl_warp.tif ' + output_file]))
        os.system('rm -rf chl_out.tif chl_warp.tif')
