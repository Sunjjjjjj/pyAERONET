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
import matplotlib.pyplot as pltport 
import matplotlib.dates as mdates
import glob
import pandas as pd
import datetime
from scipy import interpolate
from scipy.signal import argrelextrema
from numpy import trapz
import re
import pickle
from otherFunctions import *

dataOutputDir = '/nobackup/users/sunj/'
dataInputDir = '/nobackup_1/users/sunj/'


#pwd = os.getcwd()
#dataOutputDir = pwd + '/'
#dataInputDir = pwd + '/'

# =============================================================================
# inversion product
# =============================================================================
def AERONETinversion(caseName, dates, parameter = 'all'): 
    """
    Function to read AERONET Inversion version 3 product.
    
    caseName: the name of the folder containing AEROPNET sites to be used. 
    
    starttime/endtime: select data within a peroid (closed interval), format 
    in "YYYY-MM-DD".
    
    parameter: a list of parameters of interest. If 'all' (default), then 
    retrieve all parameters in the original AERONET data file; else, choose 
    parameters of interest. The following list contains the most used 
    parameters. If users are interested in other parameters, please find the 
    full parameters list returned by the function. 
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
    support_info = {}
    
    try:
        startdate, enddate = dates[0], dates[1]
    except:
        startdate = enddate = dates[0]
        
    
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
        dateTimeLocal = []
        for i in range(len(data)):
            temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
            dateTime.append(temp)
            dateTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
        data['dateTime'] = dateTime
        data['dateTimeLocal'] = dateTimeLocal
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
    try:
        output[output == -999] = np.nan
        support_info = {'parameter': parameterList, 'wavelength': sorted(wvl), 'particleSize': np.array(sorted(particleSize))}
    except:
        print('Error: no AERONET site available!')
    return output.reset_index(drop = True), support_info
    



# =============================================================================
# direct sun products
# =============================================================================
def AERONETdirectSun(caseName, dates, parameter = 'all'):
    """
    Function to read AERONET direct sun version 3 product.
    
    caseName: the name of the folder containing AEROPNET sites to be used.
    
    starttime/endtime: select data within a peroid (closed interval), format 
    in "YYYY-MM-DD".
    
    parameter: a list of parameters of interest. If 'all'(default), then 
    retrieve all parameters in the original AERONET data file; else, choose 
    parameters of interest. The following list contains the most used 
    parameters. If users are interested in other parameters, please find the 
    full parameters list returned by the function. 
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
    support_info = {}
    try:
        startdate, enddate = dates[0], dates[1]
    except:
        startdate = enddate = dates[0]
    
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
        dateTimeLocal = []
        for i in range(len(data)):
            temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
            dateTime.append(temp)
            dateTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
        data['dateTime'] = dateTime
        data['dateTimeLocal'] = dateTimeLocal
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
    support_info = {'parameter': parameterList, 'wavelength': sorted(wvl)}
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
    Dataframes contain temporal mean and std of each AERONET site.
    
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
def AERONETwvlProcess(data, support_info, parameter, wvl_int, method = 'linear'):
    """
    Function to process AERONET product.
    data and support_info: outputs of AERONETinversion or AERONETdirectSun.
    
    parameter: a parameter of interest. The parameter should be 
    wavelength-dependent. The following list contains the most used parameters. 
    If users are interested in other parameters, please find the full 
    parameters list returned by AERONETinversion or AERONETdirectSun. 
    most used parameters:
    [AOD_Extinction-Total, AOD_Extinction-Coarse, AOD_Extinction-Fine, 
    symmetry_Factor-Total, Asymmetry_Factor-Coarse, Asymmetry_Factor-Fine, 
    Absorption_AOD, Single_Scattering_Albedo, 
    Refractive_Index-Imaginary_Part, Refractive_Index-Real_Part, 
    ...]
    
    wvl_int: a list of wavelengths of interest to be interpolated/extrapolated. It can
    be a single band or a list of wavelengths.
    
    
    method: interpolation/extrapolation methods. Choose among 
    'linear' (default), 'nearest', 'quadratic', 'cubic', and 'AngstromExponent' 
    (only applicable to AOD or Absorption_AOD).
    
    
    Return: 
    a dictionary contains all sites with processed data. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2019-11-02
    """
# =============================================================================
#    interpolation functinos
# =============================================================================
    def interpAOD(tau0, lambda0, lambda1, alpha):
#        tau0, lambda1, lambda0, alpha
        return tau0 * (lambda1 / lambda0) ** (-alpha)
    
    def interpOthers(y):
        x = wvl
        f = interpolate.interp1d(x, y, kind = method, fill_value = 'extrapolate', bounds_error = False)
        return f(wvl_int)
# =============================================================================
#     Initialization
# =============================================================================
    # retrieve information
    parameterList, wvl = support_info['parameter'], support_info['wavelength'] 
    # retrieve the columns of the wavelength-dependent parameter
    parameter_ = []
    for ipara in parameterList:
        # retrieve the parameter of interest at different wavelength
        if ipara.startswith(parameter) & ('nm' in ipara):
            parameter_.append(ipara)
#        # retrieve Angstrom Exponent
#        if ipara.endswith('Angstrom_Exponent') | ipara.startswith('Absorption_Angstrom_Exponent'):
#            if (~data[ipara].isnull()).any():
#                AE.append(ipara)
    data_subset = data[parameter_]        
    data_subset.columns = wvl
    
    # a list of wavelength at which parameter is to be interpolated
    
    wvl_int = list(map(float, wvl_int))
    
    # combine wavelength and wvl
    wvl_all = sorted(list(set(wvl_int) - set(wvl)) + list(wvl))
    data_subset = data_subset.reindex(columns = wvl_all)
# =============================================================================
#   Interpolation / extrapolation
# =============================================================================
    """
    Notice: 
    For using AngstromExponent to interpolate AOD or Absorption_AOD: 
        AERONET measuring bands may differ from one site to another, the 
        Angstrom Exponent may also available in different wavelength windows. 
        Thus, the wavelength processing for each site should be processed 
        individually. For each wavelengths of interest (iwvl), the step is as 
        following:
            1) Check whether iwvl is in the range of measuring band.
            2) Check whether iwvl is in the range of Angstrom Exponent 
            wavelength windows. 
        It is also common that a iwvl is within over one wavelength window. The 
        rule of selecting AE is: choose the minimum the total range of iwvl to the 
        window boundary
    
    For using other interpolation method:
        Each record is interpolated independently.
    
    """
    if method == 'AngstromExponent': 
        # AngstromExponent is only applicable for AOD and Absorption_AOD
        if 'AOD' in parameter:
            # interpolate a site at each time
            for isite in list(set(data.Site)):
                mask = (data.Site == isite)
                # retrieve columns containing Angstrom Exponent
                AE = []
                for ipara in parameterList:
                    # Angstrom_Exponent is for direct sun and Absorption_Angstrom_Exponent is for inversion
                    if ipara.endswith('Angstrom_Exponent') | ipara.startswith('Absorption_Angstrom_Exponent'):
                        if (~data[mask][ipara].isnull()).any():
                            AE.append(ipara)
                # convert AE into a dataframe containing the boundary of wavelength windows 
                temp = []
                for ipara in AE: 
                    temp.append([float(re.findall("\d+", ipara)[0]), float(re.findall("\d+", ipara)[1])])
                AE = pd.DataFrame(temp, index = AE, columns = ['bin1', 'bin2'])
                # wavelengths containing valid data (non-NAN records)
                wvl_obs = np.array(sorted(data_subset[mask].columns[(~data_subset[mask].isnull()).any()]))
                # For each wavelength in wvl_int, check whether it is within range of Angstrom Exponent window
                for iwvl in wvl_int:
                    if (iwvl < AE.bin1.min()):
                        ipara = AE[AE.bin1 == AE.bin1.min()].index
                        if len(ipara) > 1:
                            window_width = (AE.loc[ipara] - iwvl).sum(axis = 1)
                            ipara = AE.loc[ipara][window_width == window_width.min()].index[0]
                        else:
                            ipara = ipara[0]
                        print('Warning: %.1f is outside the Angstrom Exponent window of %s from %.1f to %.1f nm, use %s.' \
                              % (iwvl, isite, AE.bin1.min(), AE.bin2.max(), ipara))
                        alpha = data[mask][ipara]
                        idx = np.argmin(abs(wvl_obs - iwvl))
                        if iwvl < wvl_obs.min():
                            print('Warning: %.1f is outside the measuring band of %s from %.1f to %.1f nm, interpolated by %.1f.' \
                                  % (iwvl, isite, wvl_obs.min(), wvl_obs.max(), wvl_obs[idx]))

                        tau0 = data_subset[mask][float(wvl_obs[idx])]
                        lambda0 = wvl_obs[idx]
                        lambda1 = iwvl
                        data_subset.loc[mask, float(iwvl)] = interpAOD(tau0, lambda0, lambda1, alpha)
                        
                    if (iwvl > AE.bin2.max()):
                        ipara = AE[AE.bin2 == AE.bin2.max()].index
                        if len(ipara) > 1:
                            window_width = (AE.loc[ipara] - iwvl).sum(axis = 1)
                            ipara = AE.loc[ipara][window_width == window_width.min()].index[0]
                        else:
                            ipara = ipara[0]
                        print('Warning: %.1f is outside the Angstrom Exponent window of %s from %.1f to %.1f nm, use %s.' \
                              % (iwvl, isite, AE.bin1.min(), AE.bin2.max(), ipara))
                        alpha = data[mask][ipara]
                        idx = np.argmin(abs(wvl_obs - iwvl))
                        if iwvl > wvl_obs.max():
                            print('Warning: %.1f is outside the measuring band of %s from %.1f to %.1f nm, interpolated by %.1f nm.' \
                                  % (iwvl, isite, wvl_obs.min(), wvl_obs.max(), wvl_obs[idx]))
                        tau0 = data_subset[mask][float(wvl_obs[idx])]
                        lambda0 = wvl_obs[idx]
                        lambda1 = iwvl
                        data_subset.loc[mask, float(iwvl)] = interpAOD(tau0, lambda0, lambda1, alpha)
                        
                    if (iwvl >= AE.bin1.min()) & (iwvl <= AE.bin2.max()):
                        ipara = AE[(iwvl >= AE.bin1) & (iwvl <= AE.bin2)].index
                        if len(ipara) > 1:
                            window_width = (AE.loc[ipara] - iwvl).sum(axis = 1)
                            ipara = AE.loc[ipara][window_width == window_width.min()].index[0]
                        else:
                            ipara = ipara[0]
                        alpha = data[mask][ipara]

                        idx = np.argmin(abs(wvl_obs - iwvl))
#                        print('%.1f nm is interpolated by %s and %.1f nm.' \
#                              % (iwvl, ipara, wvl_obs[idx]))

                        tau0 = data_subset[mask][float(wvl_obs[idx])]
                        lambda0 = wvl_obs[idx]
                        lambda1 = iwvl
                        data_subset.loc[mask, float(iwvl)] = interpAOD(tau0, lambda0, lambda1, alpha)
        else:
            raise Exception('AngstromExponent method is only applicable to AOD and Absorption_AOD!')
    else:
        data_subset[wvl_int] = list(map(interpOthers, data_subset[wvl].values))
#
# =============================================================================
#         output
# =============================================================================
#        output[isite]['data'] = processed                
    output = pd.DataFrame.from_dict(data_subset)
    output.reset_index(drop = True)
    return output


# =============================================================================
# 
# =============================================================================
def AERONETcollocation(data1, data2, data3, dates, timeWindow, Range):
    """
    Collocate AERONET with satellite / model data.
    
    @author: Sunji
    Last updated date: 2019-11-13
    """
    try:
        startdate, enddate = dates[0], dates[1]
    except:
        startdate = enddate = dates[0]

    data1 = data1.rename(columns={"dateTime": "dateTime(INV)", "timeStamp": "timeStamp(INV)"})
    DSparameter = list(set(data2) - set(['dateTime', 'dateTimeLocal', 'timeStamp', 'timeStampLocal', 'Site', 'Latitude(Degrees)', 'Longitude(Degrees)', 'Elevation(m)']))

    COL = pd.DataFrame()
    dates = pd.date_range(startdate, enddate)
    for idate in dates:
        mask = (data1['dateTime(INV)'].dt.date == idate.date())
        data1_day = data1[mask].reset_index(drop = True)
        mask = (data2['dateTime'].dt.date == idate.date())
        data2_day = data2[mask].reset_index(drop = True)
        mask = (data3['dateTime'].dt.date == idate.date())
        data3_day = data3[mask].reset_index(drop = True)
        
        if len(data1_day) * len(data2_day) * len(data3_day) > 0:
            for i in range(len(data1_day)):
                lat1, lon1 = data1_day.iloc[i: i + 1]['Latitude(Degrees)'].values, data1_day.iloc[i: i + 1]['Longitude(Degrees)'].values
                lat2, lon2 =  data3_day.lat.values, data3_day.lon.values
                distance = geoDistance(lat1, lon1, lat2, lon2) # unit: km
                dTime = abs(data1_day.iloc[i: i + 1]['timeStamp(INV)'].values - data3_day['timeStamp'].values) # unit: s
                data3_day['dTime'] = dTime
                mask1 = (distance <= Range) & (dTime <= timeWindow) # radius: 50 km, time window: +/-3 hr
                
                if any(mask1):
                    DS_subset = data2_day[data2_day.Site == data1_day.iloc[i].Site].reset_index(drop = True)
                    dTime =  abs(data1_day.iloc[i: i + 1]['timeStamp(INV)'].values - DS_subset.timeStamp.values) # unit: s 
                    
                    mask2 = (dTime <= 0.5 * 3600)  # time winder: +/-30 min
                    if any(mask2):
                        DS_subset = DS_subset[mask2].mean().to_frame().T[DSparameter]
                    
                        COL = COL.append(pd.concat([data3_day[mask1].mean().to_frame().T,\
                                                    data1_day.iloc[i: i + 1].reset_index(drop = True), \
                                                    DS_subset], axis = 1))
    COL = COL.reset_index(drop = True)
    return COL
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
    


##ROI = {'S':-90, 'N': 90, 'W': -180, 'E': 180}
##
#t1 = time.time()
##caseName = 'CA2019-01'
#caseName = 'Global_2005-2018_v3'
#
#startdate = '%4i-%02i-%02i' % (2017, 1, 1)
#enddate   = '%4i-%02i-%02i' % (2019, 12, 31) 
#
#parameter = ['Single_Scattering_Albedo', 'Absorption_AOD', 'Asymmetry_Factor-Total', 'Angstrom_Exponent',
#             'Refractive_Index-Real_Part', 'Refractive_Index-Imaginary_Part'] 
#INV, support_info1 = AERONETinversion(caseName, [startdate, enddate], parameter = parameter)
#INV_mean, INV_std = AERONETtimeProcess(INV, freq = 'day', window = False, span = ['12:00:00', '15:00:00'])
##para = AERONETwvlProcess(INV, support_info1, 'Single_Scattering_Albedo', [380, 550.], method = 'linear')
#INV['Single_Scattering_Albedo[550nm]'] = AERONETwvlProcess(INV, support_info1, 'Single_Scattering_Albedo', [550.], method = 'linear')[float(550)]
#INV['Absorption_AOD[550nm]'] = AERONETwvlProcess(INV, support_info1, 'Absorption_AOD', [550.], method = 'AngstromExponent')[float(550)]
#INV['Asymmetry_Factor-Total[550nm]'] = AERONETwvlProcess(INV, support_info1, 'Asymmetry_Factor-Total', [550.], method = 'linear')[float(550)]
#INV['Refractive_Index-Real_Part[550nm]'] = AERONETwvlProcess(INV, support_info1, 'Refractive_Index-Real_Part', [550.], method = 'linear')[float(550)]
#INV['Refractive_Index-Imaginary_Part[550nm]'] = AERONETwvlProcess(INV, support_info1, 'Refractive_Index-Imaginary_Part', [550.], method = 'linear')[float(550)]
#INV.to_pickle(dataInputDir + 'AERONET/INV_2017-onwards.pickle')
#
#
#parameter = ['AOD', 'Angstrom_Exponent', 'Solar_Zenith_Angle(Degrees) ']
#DS, support_info2 = AERONETdirectSun(caseName, [startdate, enddate], parameter = parameter)
#DS_mean, DS_std = AERONETtimeProcess(DS, freq = 'day', window = True, span = ['12:00:00', '15:00:00'])
##para = AERONETwvlProcess(DS, support_info2, 'AOD', [340, 380, 500, 550], method = 'AngstromExponent')
#DS['AOD_550nm'] = AERONETwvlProcess(DS, support_info2, 'AOD', [550.], method = 'AngstromExponent')[float(550)]
#DS.to_pickle(dataInputDir + 'AERONET/DS_2017-onwards.pickle')
#
#
#
#t2 = time.time()
#print('Time: %1.2f s' % (t2 - t1))        
