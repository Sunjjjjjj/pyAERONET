# -*- coding: utf-8 -*-
"""
reading AERONET data

@author: sunj
"""


######### load python packages
import sys, os
import time
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import glob
import pandas as pd
import datetime
from scipy import interpolate
from scipy.signal import argrelextrema
from numpy import trapz
import re
import pickle
#from otherFunctions import *

dataOutputDir = '/nobackup/users/sunj/'
dataInputDir = '/nobackup_1/users/sunj/'


pwd = os.getcwd()
dataOutputDir = pwd + '/'
dataInputDir = pwd + '/'

# =============================================================================
# inversion product
# =============================================================================
def AERONETinversion(caseName, startdate, enddate, parameter = 'all'): 
    """
    Function to read AERONET Inversion version 3 product.
    
    caseName: the name of the folder containing AEROPNET sites to be used. 
    starttime/endtime: select data within a peroid (closed interval), format 
    in "YYYY-MM-DD".
    parameter: if 'all' (default), then retrieve all parameters in the original 
    AERONET data file; else, choose parameters of interest. The following list 
    contains the most used parameters. If users are interested in other 
    parameters, please find the full parameters list returned by the function. 
    most used parameters:
        ['AOD_Extinction-Total', 'AOD_Extinction-Fine','AOD_Extinction-Coarse', 
         'Single_Scattering_Albedo', 'Absorption_AOD',
         'Refractive_Index-Real_Part', 'Refractive_Index-Imaginary_Part', 
         'Asymmetry_Factor-Total', 'Asymmetry_Factor-Fine', 
         'Asymmetry_Factor-Coarse', 'Size_Distribution_Function'
         ...]
    
    Return: 
    output: a dataframe containing all sites during the selected period. 
    The columns are parameters.  
    support_info: contains the full parameter list and 
         
    
    
    @author: Sunji
    Last updated date: 2019-10-24
    """
# =============================================================================
#   Initialization
# =============================================================================
    caseDir = dataInputDir + 'AERONET/%s/' % caseName
    filelist = glob.glob(caseDir + '*.all')
    output = pd.DataFrame()

    for isite, ff in enumerate(sorted(filelist)[:]): 
        sys.stdout.write('\r Reading AERONET inversion version 3 # %i/%i sites' % (isite + 1, len(filelist)))
        data = pd.read_csv(ff, sep = ",", header = 6)
# =============================================================================
#     Date and time
# =============================================================================
        aerDate = np.array(data['Date(dd:mm:yyyy)'])
        aerTime = np.array(data['Time(hh:mm:ss)'])
        lag = round(data['Longitude(Degrees)'] / 15)

        dateTime = []
        dataTimeLocal = []
        for i in range(len(data)):
            temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
            dateTime.append(temp)
            dataTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
        data['dateTime'] = dateTime
        data['dateTimeLocal'] = dataTimeLocal
        data['timeStamp'] = data['dateTime'].values.astype(np.int64) // 10 ** 9
        data['timeStampLocal'] = data['dateTimeLocal'].values.astype(np.int64) // 10 ** 9
        del data['Date(dd:mm:yyyy)'], data['Time(hh:mm:ss)']
# =============================================================================
#     Select data in time period using local date and time
# =============================================================================
        mask = (data['dateTimeLocal'] >= startdate) & (data['dateTimeLocal'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))
        data = data[mask] 
# =============================================================================
#     Select parameters of interest
# =============================================================================
        # retrieve all parameters in the original AERONET data file
        parameterList = list(data.keys())
        # initialization
        POI = []
        # retrieve all parameters
        if parameter == 'all':
            POI = parameterList.copy()
        # retrieve parameters of interest (POI)
        else:
            # site information
            paraSite = ['Site', 'Latitude(Degrees)', 'Longitude(Degrees)', 'Elevation(m)', \
                        'dateTime', 'dateTimeLocal', 'timeStamp', 'timeStampLocal']
            # check whether each parameter is in POI list
            for ipara in parameter:
                # retrieve particle size if required by users
                if ipara == 'Size_Distribution_Function':
                    for iparaList in parameterList: 
                        try:
                            temp = float(iparaList)
                            POI.append(iparaList)
                        except:
                            pass
                # retrieve other parameters other than particle size
                else:
                    for iparaList in parameterList: 
                        if ipara in iparaList: 
                            POI.append(iparaList)
            POI += paraSite
# =============================================================================
#    retrieve wavelength and particle size 
# =============================================================================
        # initialization
        wvl = []
        particleSize = []
         
        for iparaList in parameterList: 
            # retrieve wavelength
            if 'Single_Scattering_Albedo' in iparaList:
                try: 
                    wvl.append(float(re.findall("\d+", iparaList)[0]))
                except:
                    pass 
             # retrieve particle size 
            try:
                particleSize.append(float(iparaList))
            except:
                pass
        
        output = output.append(data[POI])
# =============================================================================
#     Output
# =============================================================================
    # remove outliers
    output[output == -999] = np.nan
    support_info = {'parameter': parameterList, 'wavelength': wvl, 'particleSize': np.array(sorted(particleSize))}
    return output.reset_index(drop = True), support_info
    



# =============================================================================
# direct sun products
# =============================================================================
def AERONETdirectSun(caseName, startdate, enddate, parameter = 'all'):
    """
    Function to read AERONET direct sun version 3 product.
    
    caseName: the name of the folder containing AEROPNET sites to be used. 
    starttime/endtime: select data within a peroid (closed interval), format 
    in "YYYY-MM-DD".
    parameter: if 'all'(default), then retrieve all parameters in the original 
    AERONET data file; else, choose parameters of interest. The following list 
    contains the most used parameters. If users are interested in other 
    parameters, please find the full parameters list returned by the function. 
    most used parameters:
        ['AOD', 'Angstrom_Exponent', 
         ...]
    
    Return: 
    output: a dataframe containing all sites during the selected period. 
    The columns are parameters.  
    support_info: contains the full parameter list and 
         
    
    @author: Sunji
    Last updated date: 2019-10-24
    """

# =============================================================================
#     Initialization
# =============================================================================
    aerCaseDir = dataInputDir + 'AERONET/%s/' % caseName    
    filelist = glob.glob(aerCaseDir + '*lev*')
    output = pd.DataFrame()
    
    for isite, ff in enumerate(sorted(filelist)[:]):
        sys.stdout.write('\r Reading AERONET direct sun version 3 # %i/%i sites' % (isite + 1, len(filelist)))
        data = pd.read_csv(ff, sep = ",", header = 6)    

        data = data.rename(columns = {'AERONET_Site_Name': 'Site',
                             'Site_Latitude(Degrees)': 'Latitude(Degrees)',
                             'Site_Longitude(Degrees)': 'Longitude(Degrees)', 
                             'Site_Elevation(m)': 'Elevation(m)'})
# =============================================================================
#    Date and time  
# =============================================================================
        aerDate = np.array(data['Date(dd:mm:yyyy)'])
        aerTime = np.array(data['Time(hh:mm:ss)'])
        lag = round(data['Longitude(Degrees)'] / 15)

        dateTime = []
        dataTimeLocal = []
        for i in range(len(data)):
            temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
            dateTime.append(temp)
            dataTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
        data['dateTime'] = dateTime
        data['dateTimeLocal'] = dataTimeLocal
        data['timeStamp'] = data['dateTime'].values.astype(np.int64) // 10 ** 9
        data['timeStampLocal'] = data['dateTimeLocal'].values.astype(np.int64) // 10 ** 9
        del data['Date(dd:mm:yyyy)'], data['Time(hh:mm:ss)']
# =============================================================================
#    Select data in time period using local date and time
# =============================================================================
        mask = (data['dateTimeLocal'] >= startdate) & (data['dateTimeLocal'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))
        data = data[mask] 
# =============================================================================
#    Select parameter of interest
# =============================================================================
        # retrieve all parameters in the original AERONET data file
        parameterList = list(data.keys())
        # initialization
        POI = []
        wvl = []
        # retrieve all parameters
        if parameter == 'all':
            POI = parameterList.copy()
            # retrieve particle size if required by users
            for iparaList in parameterList: 
                if (iparaList[:4] == 'AOD_') & (iparaList[-2:] == 'nm'):
                    try:
                        wvl.append(float(re.findall("\d+", iparaList)[0]))
                    except:
                        pass
        # retrieve parameters of interest (POI)
        else:
            # site information
            paraSite = ['Site', 'Latitude(Degrees)', 'Longitude(Degrees)', 'Elevation(m)', \
                        'dateTime', 'dateTimeLocal', 'timeStamp', 'timeStampLocal', 'Solar_Zenith_Angle(Degrees)']
            # check whether each parameter is in POI list
            for ipara in parameter:
                # retrieve particle size if required by users
                if ipara == 'AOD':
                    for iparaList in parameterList: 
                        if (iparaList[:4] == 'AOD_') & (iparaList[-2:] == 'nm'):
                            try:
                                wvl.append(float(re.findall("\d+", iparaList)[0]))
                                POI.append(iparaList)
                            except:
                                pass
                    # retrieve other parameters other than particle size
                else:
                    for iparaList in parameterList: 
                        if ipara in iparaList: 
                            POI.append(iparaList)
            POI += paraSite

        output = output.append(data[POI])
# =============================================================================
#         Output
# =============================================================================
    output[output == -999] = np.nan
    support_info = {'parameter': parameterList, 'wavelength': np.array(sorted(wvl))}
    return output.reset_index(drop = True), support_info



# =============================================================================
# AERONET time process
# =============================================================================
def AERONETtimeProcess(data, freq = 'day', window = False, **kwargs):
    """
    Function to process AERONET.
    
    data: outputs of AERONETinversion or AERONETdirectSun.
    freq: time processing frequency, chosse from 'original', 'month', 'day' and
    'hour'.
    window: whether use a time window in case of, e.g, satellite 
    overpass period. If False, then use all records.
    span: if period is True, then specify the time window in '**kwarg', e.g. if 
    the time window is 13 p.m -14 p.m., then specify as ['13:00:00', '14:00:00'].
    Note the time is LOCAL time, NOT the UTC.
    
    Return:
    Temporal mean and std of input AERONET input data.
    
    @author: Sunji
    Last updated date: 2019-10-25
    """

# =============================================================================
#     Initialization
# =============================================================================
    timedf = pd.DataFrame()
    # use local date and time
    timedf['year'] = data.dateTimeLocal.dt.year
    timedf['YYMM'] = data.dateTimeLocal.dt.to_period('M')
    timedf['date'] = data.dateTimeLocal.dt.to_period('D')
    timedf['time'] = data.dateTimeLocal.dt.time
    timedf['hour'] = data.dateTimeLocal.dt.to_period('H')
# =============================================================================
#    Select records within period of interest
# =============================================================================
    # no time window
    if ~window:
        pass
    # specified time winder
    if window:
        span = kwargs['span']
        starttime = pd.to_datetime(span[0], format = '%H:%M:%S')
        endtime = pd.to_datetime(span[1], format = '%H:%M:%S')
        mask = (timedf['time'] >= starttime.time()) & (timedf['time'] <= endtime.time())
        data = data[mask] 
# =============================================================================
#    Temporal processed
# =============================================================================
    # no process
    if freq == 'original':
        data_mean = data.copy()
        data_std = None
    # process based on frequency of hour, day or month
    else:
        if freq == 'hour':
            data_mean = data.groupby([timedf['hour'], 'Site']).mean()
            data_std = data.groupby([timedf['hour'], 'Site']).std()
            count = data.groupby([timedf['hour'], 'Site']).count()
        
        if freq == 'day':
            data_mean = data.groupby([timedf['date'], 'Site']).mean()
            data_std = data.groupby([timedf['date'], 'Site']).std()
            count = data.groupby([timedf['date'], 'Site']).count()
    
        if freq == 'month':
            data_mean = data.groupby([timedf['YYMM'], 'Site']).mean()
            data_std = data.groupby([timedf['YYMM'], 'Site']).std()
            count = data.groupby([timedf['YYMM'], 'Site']).count()

        # reset index, dateTime and time Stamp
        data_mean.reset_index(inplace = True)
        data_std.reset_index(inplace = True)
        data_mean['num'] = count.values[:, 0]
        data_std['num'] = count.values[:, 0]
        
        data_mean['dateTime'] = pd.to_datetime(data_mean['timeStamp'] * 10 ** 9)
        data_std['dateTime'] = pd.to_datetime(data_mean['timeStamp'] * 10 ** 9)
        data_mean['dateTimeLocal'] = pd.to_datetime(data_mean['timeStampLocal'] * 10 ** 9)
        data_std['dateTimeLocal'] = pd.to_datetime(data_mean['timeStampLocal'] * 10 ** 9)
# =============================================================================
#    Output
# =============================================================================
    return data_mean, data_std

# =============================================================================
# AERONET wavelength interpolation/extrapolation
# =============================================================================
def AERONETwvlProcess(data, support_info, parameter, wavelength, method = 'linear'):
    """
    Function to process AERONET product.
    data: outputs of AERONETinversion or AERONETdirectSun.
    parameter: parameters of interest. Select from the following parameters: 
    [AOD_Extinction-Total, AOD_Extinction-Coarse, AOD_Extinction-Fine, 
    symmetry_Factor-Total, Asymmetry_Factor-Coarse, Asymmetry_Factor-Fine, 
    Absorption_AOD, Single_Scattering_Albedo, 
    Refractive_Index-Imaginary_Part, Refractive_Index-Real_Part, 
    Scattering_Angle_Bin_3.2_to_<6_degrees, 
    Scattering_Angle_Bin_30_to_<80_degrees,
    Scattering_Angle_Bin_6_to_<30_degrees,
    
    Depolarization_Ratio, Lidar_Ratio,
    AOD_Coincident_Input,]
    wavelength: wavelengths of interest to be interpolated/extrapolated.
    method: interpolation/extrapolation methods,
    choose among 'linear' (default), 'quadratic', 'cubic.
    If the parameter is AOT or AOTabsp, the Angstorm Exponent will be used instead of given method.
    
    
    
    Return: 
    a dictionary contains all sites with processed data. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2019-10-25
    """

# =============================================================================
#     Initialization
# =============================================================================
    parameterList, wvl = support_info['parameter'], support_info['wavelength'] 
    
    parameter_ = []
    for ipara in parameterList:
        if parameter in ipara:
            parameter_.append(ipara)
            
    data_subset = data[parameter_]        
    data_subset.columns = wvl
    
    WOI = sorted(list(set(wavelength) - set(wvl)) + wvl)
    data_subset = data_subset.reindex(columns = WOI)
    
    data_subset = data_subset.T
    
#            paradata = paradata.dropna(axis = 1, how = 'all')
# =============================================================================
#   interpolation
# =============================================================================
#    data_subset.interpolate(method = method, limit_direction='forward')
#            """
#            Interpolation/extrapolation
#            """
#            x = wvl
#            y = paradata
#            if len(y) > 0: 
#                f = interpolate.interp1d(x, y, kind = wvlProcessMethod, fill_value = 'extrapolate', bounds_error = False)
#                for iwvl in WOI:
#                    paraname = columns[0]
#                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
#                    processed['%s' % (paraname)] = f(iwvl)
#                
#                """
#                Angstorm exponent:
#                If the parameter is AOT or AOTabsp, ignore the given method, use Angstorm exponent (AE) instead. 
#                Use nearby AOT/AOTabsp to calculate AE and predict the target wavelength.
#                """
#                wvl = sorted(wvl)
#                if ipara.find('AOT') >= 0: 
#                    for iwvl in WOI:
#                        for i in range(0, len(wvl) - 1):
#                            if int(iwvl) in range(int(wvl[i]), int(wvl[i + 1])): 
#                                if ipara == 'AOT':  
#                                    AE = Angstorm(wvl[i], paradata['AOT_%i' % (wvl[i])], wvl[i + 1], paradata['AOT_%i' % (wvl[i + 1])])
#                                    paraname = columns[0]
#                                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
#                                    processed['%s' % (paraname)] = wvldepAOD(wvl[i], paradata['AOT_%i' % (wvl[i])], iwvl, AE)
#                                else:
#                                    AE = Angstorm(wvl[i], paradata['AOTAbsp%i-T' % (wvl[i])], wvl[i + 1], paradata['AOTAbsp%i-T' % (wvl[i + 1])])
#                                    paraname = columns[0]
#                                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
#                                    processed['%s' % (paraname)] = wvldepAOD(wvl[i], paradata['AOTAbsp%i-T' % (wvl[i])], iwvl, AE)
#                                
#
# =============================================================================
#         output
# =============================================================================
#        output[isite]['data'] = processed                
    return pd.DataFrame.from_dict(output)


# =============================================================================
# Test code
# =============================================================================
#def main(): 
#caseName = 'CA201712'
#ROI = {'S':25, 'N': 50, 'W': -130, 'E': -110}
#
#startdate = '%4i-%02i-%02i' % (2017, 12, 1)
#enddate   = '%4i-%02i-%02i' % (2017, 12, 31) 
#ctime = datetime.datetime(2017,12,12,19,50,1)
#dtime = 1800
#
## read data
#INV = AERONETinversion(caseName, startdate, enddate)
##DS = AERONETdirectSun(caseName, startdate, enddate)
### time process 
##inv_mm = AERONETtimeProcess(INV, 'monthly')
##INV_pp = AERONETtimeProcess(INV, 'period', ctime, dtime)
##INV_dd = AERONETtimeProcess(INV, 'daily')
## wavelength interpolation
#INV_int = AERONETwvlProcess(INV, ['SSA', 'AOTAbsp'], [388, 550], wvlProcessMethod = 'linear')
##DS_int = AERONETwvlProcess(DS, ['AOT'], [388, 500, 550], 'linear')
#
#
#if __name__ == '__main__':
#    main()
    

ROI = {'S':-90, 'N': 90, 'W': -180, 'E': 180}

t1 = time.time()
caseName = 'CA2017-18'

startdate = '%4i-%02i-%02i' % (2019, 5, 1)
enddate   = '%4i-%02i-%02i' % (2019, 5, 31) 

parameter = ['AOD_Extinction-Total', 'AOD_Extinction-Fine','AOD_Extinction-Coarse', 
         'Single_Scattering_Albedo', 'Absorption_AOD'] 

parameter = ['AOD_Extinction-Total', 
         'Single_Scattering_Albedo', 'Absorption_AOD'] 

INV, support_info = AERONETinversion(caseName, startdate, enddate, parameter = 'all')
INV_mean, INV_std = AERONETtimeProcess(INV, freq = 'day', window = True, span = ['12:00:00', '15:00:00'])


#parameter = ['AOD', 'Angstrom_Exponent']
#DS, support_info = AERONETdirectSun(caseName, startdate, enddate, parameter = 'all')
#DS_mean, DS_std = AERONETtimeProcess(DS, freq = 'day', window = True, span = ['12:00:00', '15:00:00'])
#


#INV_int = AERONETwvlProcess(INV, ['SSA', 'AOTAbsp'], [388, 440, 500, 532, 550], 'linear')
#DS_int = AERONETwvlProcess(DS, ['AOT'], [388, 440, 500, 532, 550], 'linear')

#with open(casedir + 'INV_%4i.pickle' % (iyear), 'wb') as handle:
#    pickle.dump(INV_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open(casedir + 'DS_%4i.pickle' % (iyear), 'wb') as handle:
#    pickle.dump(DS_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
t2 = time.time()
print('Time: %1.2f s' % (t2 - t1))        
