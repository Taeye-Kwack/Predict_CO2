import numpy as np
import joblib
import glob
from sklearn.preprocessing import StandardScaler

fpath = '/home/xodpwkd/oco2/Model_input/'
flist = np.sort(glob.glob(fpath+'Model_input_2020*'))
model = joblib.load('/home/xodpwkd/oco2/RF_with_time.pkl')

for i in flist:
    data = np.load(i)
    data = np.delete(data,-1,axis=1) # remove sst
    data = pd.DataFrame(data)
    data.columns = ['t2m','sfcp','u10','v10','u100','v100','u','v','w','co','co_time','no2','lat','lon','skt']
    #data['sin(lat)'] = np.sin(data['lat']*np.pi/180)
    #data['sin(lon)'] = np.sin(data['lon']*np.pi/180)
    #data = data[['t2m', 'sfcp', 'u10', 'v10', 'u100', 'v100', 'u', 'v', 'w', 'co', 'co_time', 'no2', 'lat', 'lon','sin(lat)', 'sin(lon)','skt']]

    #data = np.delete(data,10,axis=1) # remove co_time : ['t2m','sfcp','u10','v10','u100','v100','u','v','w','co','no2','lat','lon','skt']
    #scaler = StandardScaler()
    #scaler.fit(data)
    #data_scaled = scaler.transform(data)
    y_pred = model.predict(data)
    y_pred = y_pred.reshape(-1,1)
    data = np.append(data, y_pred, axis=1)
    np.save('/home/xodpwkd/oco2/Model_output/RF_with_time/Model_output_' + i[-12:-4], data)
    print('/home/xodpwkd/oco2/Model_output/RF_with_time/Model_output_' + i[-12:-4])

