import numpy as np
from netCDF4 import Dataset as D
from datetime import datetime, timedelta

year = 2018
yyyy=repr(year)

for m in range(8,13):
    fpath = '/home/DATA_ARCHIVE/ERA5/HOURLY/t2m/'
    fileid = D(fpath + 't2m' + yyyy + str(m).zfill(2) +'.nc','r')
    lon_t = fileid.variables['longitude'][:]
    lon_t[720:] = lon_t[720:] - 360 #[0,360] -> [0,...,180,-179.75,...,-0.25]
    lat_t = fileid.variables['latitude'][:]
    t2m = fileid.variables['t2m'][:]
    time= fileid.variables['time'][:] # hours since 1900-01-01
    delta_hour = (datetime(2010, 1, 1) - datetime(1900, 1, 1)).total_seconds()/3600
    time = time - delta_hour #hours since 2010-01-01
    time= time.astype(float)

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

    #print(t2m,sfcp,u,u10)

    fpath = '/home/xodpwkd/oco2/newgrid/'
    overlap = np.load(fpath + 'overlap_' + yyyy + str(m).zfill(2) + '.npy', allow_pickle=True)
    co = overlap[0] 
    co_time = overlap[1]
    no2 = overlap[2]
    co2 = overlap[3]
    co2_time = overlap[4]
    lat = overlap[5] - 0.125 # [-90, -89.75, ..., 89.75]
    lon = overlap[6] - 0.125 # [-180, -179.75, ..., 179.75]

    mean_time = ((co_time + co2_time)*0.5/3600).astype(int) #meantime (hour) since 2010-01-01
    #mean_time = np.where(mean_time >= 85440, 85439, mean_time)
    mean_time = np.where(mean_time > np.max(time), mean_time-1, mean_time)

    Allvariables = np.empty((0,18), float)
    for i in range(len(co)):
        x = np.where(time == mean_time[i])[0][0]
        y = np.where(lat_t == lat[i])[0][0]
        z = np.where(lon_t == lon[i])[0][0]
        Allvariables = np.append(Allvariables, np.array([[t2m[x,y,z],sfcp[x,y,z],u10[x,y,z],v10[x,y,z],u100[x,y,z],v100[x,y,z],u[x,y,z],v[x,y,z],w[x,y,z],co[i],co_time[i],no2[i],co2_time[i],lat[i],lon[i],skt[x,y,z],sst[x,y,z],co2[i]]]),axis=0)
    print('Allvariables shape is', yyyy, m, Allvariables.shape)
    fpath = '/home/xodpwkd/oco2/newgrid/'
    np.save(fpath + 'Allvariables_'+ yyyy + str(m).zfill(2), Allvariables)
