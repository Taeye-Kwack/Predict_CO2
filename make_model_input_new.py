import numpy as np
import glob
from netCDF4 import Dataset as D
import matplotlib.pyplot as plt
import sys
from datetime import datetime, timedelta
#import warnings

##################tropomi##############################
dx = 0.25
dy = 0.25
nnx = int(360 / dx)
nny = int(180 / dy) + 1
year = 2020
yyyy = repr(year)
lon_cams = np.arange(-180,180,step=0.25) # n = 1440
lat_cams = np.arange(-90,90.25,step=0.25) # n = 721

#for m in range(1,13):
m=8
fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/t2m/'
fileid = D(fpath + 't2m' + yyyy + str(m).zfill(2) +'.nc','r')
lon_t = fileid.variables['longitude'][:]
lon_t[720:] = lon_t[720:] - 360 #[0,360] -> [0,...,180,-179.75,...,-0.25]
lat_t = fileid.variables['latitude'][:]
t2m = fileid.variables['t2m']
time= fileid.variables['time'][:] # hours since 1900-01-01
delta_hour = (datetime(2010, 1, 1) - datetime(1900, 1, 1)).total_seconds()/3600
time = time - delta_hour #hours since 2010-01-01
time = time.astype(float)

fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/sfcp/'
fileid = D(fpath + 'sfcp' + yyyy + str(m).zfill(2) +'.nc','r')
sfcp = fileid.variables['sp']

fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/pl_wind/'
fileid = D(fpath + 'ERA5_wind_pl_500hPa_' + yyyy + str(m).zfill(2) +'.nc','r')
u = fileid.variables['u']  
v = fileid.variables['v']
w = fileid.variables['w']
         
fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/single_wind/'
fileid = D(fpath + 'ERA5_wind_single_' + yyyy + str(m).zfill(2) +'.nc','r')
u10 = fileid.variables['u10']
v10 = fileid.variables['v10'] 
u100 = fileid.variables['u100']
v100 = fileid.variables['v100']

fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/SST_SKT_Vgtype/'
f = D(fpath + 'sst_skt_vgtype_cover' + yyyy + str(m).zfill(2) + '.nc','r')
sst = f.variables['sst']
skt = f.variables['skt']

for d in range(1,32):
    Allvariables = np.empty((0,16), float)
    fpath = '/home/DATA_ARCHIVE/TROPOMI_HiR/CO/%s/'%(year)
    flist_co = np.sort(glob.glob(fpath + 'S5P_OFFL_L2__CO_____%s%s%s'%(year,str(m).zfill(2),str(d).zfill(2)) + '*'))
    fpath = '/home/DATA_ARCHIVE/TROPOMI_HiR/NO2/%s/'%(year)
    flist_no2 = np.sort(glob.glob(fpath +'S5P_OFFL_L2__NO2____%s%s%s'%(year,str(m).zfill(2),str(d).zfill(2)) + '*'))
    if len(flist_co) == len(flist_co):
        for ii in range(len(flist_co)): 
            x_time = datetime.now()
            gridat = np.zeros((nny,nnx)) ### mean
            gridco = np.zeros((nny,nnx)) ### count
            gridtime = np.zeros((nny,nnx))
            fname_S5P = flist_co[ii]
            tro = D(fname_S5P, 'r');
            var = tro.groups['PRODUCT'].variables
            lonset = var['longitude'][:][0,:,:] #(time,scanline,grid_pixel)
            latset = var['latitude'][:][0,:,:]
            rawdat = var['carbonmonoxide_total_column'][:][0,:,:]
            rawdat = rawdat.filled(np.nan)
            rawdat = np.where(rawdat<0, np.nan, rawdat)
            delta_time = var['delta_time'][0].astype(float)/1000 + var['time'][0] #seconds since 2010-01-01
            qa_value = var['qa_value'][0]

            for i in range(rawdat.shape[0]): #4175
                idx = np.where(~np.isnan(rawdat[i,:]) * (qa_value[i,:] > 0.5));
                lonidx, latidx, rawidx = lonset[i,:][idx], latset[i,:][idx], rawdat[i,:][idx] # 1D  scanline (215,)
                lonidx = np.int16(np.floor((lonidx + 180)/dx))
                idx = np.where(lonidx == nnx)
                lonidx[idx] = 0
                latidx = np.int16(np.floor((latidx + 90)/dx))
                    
                for j in range(rawidx.shape[0]): # if all data is good, rawdat.shape[0] = 215
                    xx = lonidx[j]
                    yy = latidx[j]
                    gridat[yy, xx] = gridat[yy,xx] + rawidx[j]
                    gridco[yy, xx] = gridco[yy,xx] + 1.
                    gridtime[yy,xx] = gridtime[yy,xx] + delta_time[i]

            idx = np.where(gridco >= 1.)
            nidx = np.where(gridco == 0.)
            gridat[idx]  = gridat[idx] / gridco[idx]
            gridtime[idx] = gridtime[idx] / gridco[idx]
            gridat[nidx], gridtime[nidx] = np.nan, np.nan

            gridat_co, gridco_co, gridtime_co = gridat, gridco, gridtime

            fpath = '/home/xodpwkd/oco2/TROPOMI_gridded/CO/'
            saving = D(fpath + 'Grid_' + fname_S5P[39:], 'w', format='NETCDF4')
            saving.createDimension('longitude', 1440)
            saving.createDimension('latitude', 721)

            x = saving.createVariable('longitude', 'f8', ('longitude',))
            y = saving.createVariable('latitude', 'f8', ('latitude',))
            CO =saving.createVariable('co2 concentration', 'f8', ('latitude','longitude'))
            GRIDCO = saving.createVariable('grid_count', 'f8', ('latitude','longitude'))
            GRIDTIME = saving.createVariable('grid_time', 'f8', ('latitude','longitude'))

            x[:] = lon_cams
            y[:] = lat_cams
            CO[:,:] = gridat_co
            GRIDCO[:,:] = gridco_co
            GRIDTIME[:,:] = gridtime_co

        ##################################NO2################################
            gridat = np.zeros((nny,nnx)) ### mean
            gridco = np.zeros((nny,nnx)) ### count
            gridtime = np.zeros((nny,nnx))
            
            fname_S5P = flist_no2[ii]
            tro = D(fname_S5P, 'r');
            var = tro.groups['PRODUCT'].variables
            lonset = var['longitude'][:][0,:,:] #(time,scanline,grid_pixel)
            latset = var['latitude'][:][0,:,:]
            rawdat = var['nitrogendioxide_tropospheric_column'][:][0,:,:] #masked data
            rawdat = rawdat.filled(np.nan) #masked data -> np.nan
            rawdat = np.where(rawdat<0, np.nan, rawdat)
            delta_time = var['delta_time'][0].astype(float)/1000 + var['time'][0]
            qa_value = var['qa_value'][0]

            for i in range(rawdat.shape[0]): #4175
                idx = np.where(~np.isnan(rawdat[i,:]) * (qa_value[i,:] >= 0.75))
                lonidx, latidx, rawidx = lonset[i,:][idx], latset[i,:][idx], rawdat[i,:][idx] # 1D  scanline (215,)
                lonidx = np.int16(np.floor((lonidx + 180)/dx))
                latidx = np.int16(np.floor((latidx + 90)/dx))
                idx = np.where(lonidx == nnx)
                lonidx[idx] = 0

                for j in range(rawidx.shape[0]): # if all data is good, rawdat.shape[0] = 215
                    xx = lonidx[j]
                    yy = latidx[j]
                    gridat[yy, xx] = gridat[yy,xx] + rawidx[j]
                    gridco[yy, xx] = gridco[yy,xx] + 1.
                    gridtime[yy,xx] = gridtime[yy,xx] + delta_time[i]

            idx = np.where(gridco >= 1.)
            nidx = np.where(gridco == 0.)
            gridat[idx]  = gridat[idx] / gridco[idx]
            gridtime[idx] = gridtime[idx] / gridco[idx]
            gridat[nidx], gridtime[nidx] = np.nan, np.nan

            gridat_no2, gridco_no2, gridtime_no2 = gridat, gridco, gridtime

            fpath = '/home/xodpwkd/oco2/TROPOMI_gridded/NO2/'
            saving = D(fpath + 'Grid_' + fname_S5P[40:], 'w', format='NETCDF4')
            saving.createDimension('longitude', 1440)
            saving.createDimension('latitude', 721)

            x = saving.createVariable('longitude', 'f8', ('longitude',))
            y = saving.createVariable('latitude', 'f8', ('latitude',))
            NO2 =saving.createVariable('no2 concentration', 'f8', ('latitude','longitude'))
            GRIDCO = saving.createVariable('grid_count', 'f8', ('latitude','longitude'))
            GRIDTIME = saving.createVariable('grid_time', 'f8', ('latitude','longitude'))

            x[:] = lon_cams
            y[:] = lat_cams
            NO2[:,:] = gridat_no2
            GRIDCO[:,:] = gridco_no2
            GRIDTIME[:,:] = gridtime_no2
        #####################################merge data#################################
            delta_time = (gridtime_co + gridtime_no2)/2
            idx = np.where(~np.isnan(delta_time))
            co = gridat_co[idx]
            no2 = gridat_no2[idx]
            lat = lat_cams[idx[0]]
            lon = lon_cams[idx[1]]
            time_co = delta_time[idx]
            time_co_hour = (delta_time[idx]/3600).astype(int) #hours since 2010-01-01
            time_co_hour = np.where(time_co_hour > np.max(time), time_co_hour-1, time_co_hour)
            time_co_hour = np.where(time_co_hour > np.max(time), time_co_hour-1, time_co_hour)

            for i in range(len(co)):
                x = np.where(time == time_co_hour[i])[0][0]
                y = np.where(lat_t == lat[i])[0][0]
                z = np.where(lon_t == lon[i])[0][0]
                Allvariables = np.append(Allvariables, np.array([[t2m[x,y,z],sfcp[x,y,z],u10[x,y,z],v10[x,y,z],u100[x,y,z],v100[x,y,z],u[x,y,z],v[x,y,z],w[x,y,z],co[i],time_co[i],no2[i],lat[i],lon[i],skt[x,y,z],sst[x,y,z]]]), axis=0)
            print(fname_S5P)
            y_time = datetime.now()
            print('Time spent is ' + str(y_time-x_time))
            print('Allvariables shape is', yyyy, m, d, Allvariables.shape)
            fpath = '/home/xodpwkd/oco2/newgrid/'
            np.save(fpath + 'Model_input_'+ yyyy + str(m).zfill(2) + str(d).zfill(2), Allvariables)

