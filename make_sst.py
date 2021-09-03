import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset as D
#import cartopy.crs as ccrs
#import cartopy
from datetime import datetime, timedelta

#[t2m,sfcp,u10,v10,u100,v100,u,v,w,co,co_time,no2,co2_time,hour,lat,lon,sif,co2]
# 0     1   2   3   4     5  6 7 8  9   10     11    12     13   14  15  16  17

for year in ['2020']:
    for m in range(12,13):
        all = np.load('Allvariables_' + year + str(m).zfill(2) + '.npy')
        idx = np.where(all[:,-2] == all[:,-1]) # my fault
        all[idx,-2] = np.nan
        co_time = all[:,10]
        co2_time = all[:,12]
        mean_time = ((co_time + co2_time)*0.5/3600).astype(int)
        lat = all[:,14]
        lon = all[:,15]

        fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/SST_SKT_Vgtype/'
        f = D(fpath + 'sst_skt_vgtype_cover' + year + str(m).zfill(2) + '.nc','r')
        time = f.variables['time'][:]
        delta_hour = (datetime(2010, 1, 1) - datetime(1900, 1, 1)).total_seconds()/3600
        time = time - delta_hour #hours since 2010-01-01
        time= time.astype(float)
        mean_time = np.where(mean_time > np.max(time), mean_time-1, mean_time)
        lat_t = f.variables['latitude']
        lon_t = f.variables['longitude']
        sst = f.variables['sst']
        skt = f.variables['skt']

        sst_np, skt_np = np.array([]), np.array([])
        for i in range(len(lat)):
            x = np.where(time == mean_time[i])[0][0]
            y = np.where(lat_t == lat[i])[0][0]
            z = np.where(lon_t == lon[i])[0][0]
            sst_np = np.append(sst_np, sst[x,y,z])
            skt_np = np.append(skt_np, skt[x,y,z])
            
        idx = np.where(sst_np ==0)
        sst_np[idx] = np.nan

        ############[t2m,sfcp,u10, ...., lat,lon,sif,skt,sst,co2]##############
        all = np.c_[all,all[:,-1]]
        all = np.c_[all,all[:,-1]]
        all[:,-3] = skt_np
        all[:,-2] = sst_np

        np.save('Allvariables_' + year + str(m).zfill(2), all)
        print('Allvariables_' + year + str(m).zfill(2) + ' is done')

'''
    ax = plt.figure()
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=180))
    ax.coastlines(zorder=1)
    plt.scatter(lon,lat,s=0.1)
    plt.title('202001 ocean data (w/o sif)')
'''
