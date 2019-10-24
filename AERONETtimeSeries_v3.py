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
from otherFunctions import *

dataOutputDir = '/nobackup/users/sunj/'
dataInputDir = '/nobackup_1/users/sunj/'


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

        dateTime = []
#        timeStamp = []
#        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(data)):
#            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            dateTime.append(pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S'))
        data['dateTime'] = dateTime
        del data['Date(dd:mm:yyyy)'], data['Time(hh:mm:ss)']
#        data['timeStamp'] = pd.timeStamp(dateTime)
# =============================================================================
#     Select data in time period 
# =============================================================================
        mask = (data['dateTime'] >= startdate) & (data['dateTime'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))
        data = data[mask] 
# =============================================================================
#     Select parameters of interest
# =============================================================================
        # parameter names are same for all AERONET site, thus only retrieve  parameters of the first site
#        if isite == 0: 
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
            paraSite = ['Site', 'Latitude(Degrees)', 'Longitude(Degrees)', 'Elevation(m)', 'dateTime']
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
        
        # pass this procedue for other site
#        else:
#            pass
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

# =============================================================================
#    Date and time  
# =============================================================================
        aerDate = np.array(data['Date(dd:mm:yyyy)'])
        aerTime = np.array(data['Time(hh:mm:ss)'])

        dateTime = []
#        timeStamp = []
#        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(data)):
#            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            dateTime.append(pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S'))
        data['dateTime'] = dateTime
        del data['Date(dd:mm:yyyy)'], data['Time(hh:mm:ss)']
#        data['timeStamp'] = pd.timeStamp(dateTime)
# =============================================================================
#    Select data in time period 
# =============================================================================
        mask = (data['dateTime'] >= startdate) & (data['dateTime'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))
        data = data[mask] 
# =============================================================================
#    Select parameter of interest
# =============================================================================
        # parameter names are same for all AERONET site, thus only retrieve parameters of the first site
#        if isite == 0: 
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
            paraSite = ['AERONET_Site_Name', 'Site_Latitude(Degrees)', 'Site_Longitude(Degrees)', 'Site_Elevation(m)', 'dateTime', \
                        'Solar_Zenith_Angle(Degrees)']
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
        # pass this procedue for other site
#        else:
#            pass
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
def AERONETtimeProcess(Data, timeProcessMethod = 'daily', *arg):
    """
    Function to process AERONET product.
    
    -Data: outputs of AERONETinversion or AERONETdirectSun.
        
    -timeProcessMethod: time processing methods.
    choose among 'period' (average over a certain period), 'daily' (daily mean, default), 'monthly' (monhly mean).
    
    -*arg: if processing method is 'period', time period is need to be specified by ctime and dtime. 
    period = ctime (format in datetime or timestamp) +/- dtime (in second) 
    
    Return: 
    a dictionary contains all sites with processed data. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2018-10-18
    """

    """
    Initialization
    """
    output = {}
    for i, isite in enumerate(Data.keys()):
        sys.stdout.write('\r Time processing # %i/%i sites' % (i + 1, len(Data.keys())))
        output[isite] = {'lat': Data[isite]['lat'], 'lon': Data[isite]['lon'], 'elev': Data[isite]['elev']}
        processed = pd.DataFrame()
        data = Data[isite]['data'].copy()
        
        if len(data) > 0: 
            """
            Monthly mean
            """
            if timeProcessMethod == 'monthly':
                data['YY-MM'] = data.index.to_period('M')
                processed = data.groupby('YY-MM').mean()
            else:
                data['YY-MM-DD'] = data.index.to_period('D')
                """
                Daily mean
                """            
                if timeProcessMethod == 'daily':
                    try:
                        processed = data.groupby('YY-MM-DD').mean()
                    except:
                        pass
                """
                Co-time-period mean
                """
                if timeProcessMethod == 'period':
                    ctime, dtime = arg[0], arg[1]
                    try:
                        timestart = ctime - dtime
                        timeend = ctime + dtime
                        tsrange = (data.timeStamp >= timestart) & (data.timeStamp < timeend)
                        ctime = datetime.datetime.fromtimestamp(ctime)
                    except:
                        timestart = ctime - datetime.timedelta(seconds=dtime)
                        timeend = ctime + datetime.timedelta(seconds=dtime)
                        tsrange = (data.index >= timestart) & (data.index < timeend)
                    if tsrange.sum() > 0:
                        processed = data[tsrange].mean()
                        processed = processed.to_frame().T
                        processed.index = pd.Series(ctime)
        """
        Output
        """
        output[isite]['data'] = processed
                    
    return pd.DataFrame.from_dict(output)

# =============================================================================
# AERONET wavelength interpolation/extrapolation
# =============================================================================
def AERONETwvlProcess(Data, POI, WOI, wvlProcessMethod = 'linear'):
    """
    Function to process AERONET product.
    
    -Data: outputs of AERONETinversion or AERONETdirectSun.
    
    -POI: parameters of interest.
    
    -WOI: wavelengths of interest to be interpolated/extrapolated.
    
    -wvlProcessMethod: interpolation/extrapolation methods,
    choose among 'linear' (default), 'quadratic', 'cubic.
    If the parameter is AOT or AOTabsp, the Angstorm Exponent will be used instead of given method.
    
    
    
    Return: 
    a dictionary contains all sites with processed data. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2018-10-18
    """

    """
    Initialization
    """
    output = {}
    for i, isite in enumerate(Data.keys()):
        sys.stdout.write('\r Wavelength processing # %i/%i sites' % (i + 1, len(Data.keys())))
        output[isite] = {'lat': Data[isite]['lat'], 'lon': Data[isite]['lon'], 'elev': Data[isite]['elev']}
        data = Data[isite]['data'].copy()
        processed = data.copy()
        for ipara in POI:
            columns = []
            for ikey in data.keys():
                if ikey.find(ipara) == 0: 
                    columns.append(ikey)
            paradata = data[columns]
            """
            Remove NaN columns (no measurement at this wavelength)
            """
            paradata = paradata.dropna(axis = 1, how = 'all')
            """
            Update wavelengths
            """         
            wvl = []
            for ikey in paradata:
                wvl.append(float(re.findall(r'\d+', ikey)[0]))
            """
            Interpolation/extrapolation
            """
            x = wvl
            y = paradata
            if len(y) > 0: 
                f = interpolate.interp1d(x, y, kind = wvlProcessMethod, fill_value = 'extrapolate', bounds_error = False)
                for iwvl in WOI:
                    paraname = columns[0]
                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
                    processed['%s' % (paraname)] = f(iwvl)
                
                """
                Angstorm exponent:
                If the parameter is AOT or AOTabsp, ignore the given method, use Angstorm exponent (AE) instead. 
                Use nearby AOT/AOTabsp to calculate AE and predict the target wavelength.
                """
                wvl = sorted(wvl)
                if ipara.find('AOT') >= 0: 
                    for iwvl in WOI:
                        for i in range(0, len(wvl) - 1):
                            if int(iwvl) in range(int(wvl[i]), int(wvl[i + 1])): 
                                if ipara == 'AOT':  
                                    AE = Angstorm(wvl[i], paradata['AOT_%i' % (wvl[i])], wvl[i + 1], paradata['AOT_%i' % (wvl[i + 1])])
                                    paraname = columns[0]
                                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
                                    processed['%s' % (paraname)] = wvldepAOD(wvl[i], paradata['AOT_%i' % (wvl[i])], iwvl, AE)
                                else:
                                    AE = Angstorm(wvl[i], paradata['AOTAbsp%i-T' % (wvl[i])], wvl[i + 1], paradata['AOTAbsp%i-T' % (wvl[i + 1])])
                                    paraname = columns[0]
                                    paraname = paraname.replace(re.findall(r'\d+', columns[0])[0], str(iwvl))
                                    processed['%s' % (paraname)] = wvldepAOD(wvl[i], paradata['AOTAbsp%i-T' % (wvl[i])], iwvl, AE)
                                

        """
        output
        """
        output[isite]['data'] = processed                
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
casedir = dataInputDir + 'AERONET/Global_2005-2018_v3/'

t1 = time.time()
caseName = 'CA2017-18'

startdate = '%4i-%02i-%02i' % (2019, 1, 1)
enddate   = '%4i-%02i-%02i' % (2019, 12, 31) 

parameter = ['AOD_Extinction-Total', 'AOD_Extinction-Fine','AOD_Extinction-Coarse', 
         'Single_Scattering_Albedo', 'Absorption_AOD'] 


INV, support_info = AERONETinversion(caseName, startdate, enddate, parameter = parameter)
parameter = ['AOD', 'Angstrom_Exponent']

DS, support_info = AERONETdirectSun(caseName, startdate, enddate, parameter = parameter)
#
#INV_int = AERONETwvlProcess(INV, ['SSA', 'AOTAbsp'], [388, 440, 500, 532, 550], 'linear')
#DS_int = AERONETwvlProcess(DS, ['AOT'], [388, 440, 500, 532, 550], 'linear')

#with open(casedir + 'INV_%4i.pickle' % (iyear), 'wb') as handle:
#    pickle.dump(INV_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
#with open(casedir + 'DS_%4i.pickle' % (iyear), 'wb') as handle:
#    pickle.dump(DS_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
t2 = time.time()
print('Time: %1.2f s' % (t2 - t1))        
