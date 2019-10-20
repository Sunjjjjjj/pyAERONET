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


# =============================================================================
# inversion product
# =============================================================================
def AERONETinversion(caseName, startdate, enddate): 
    """
    Function to read AERONET Inversion product.
    
    -caseName: name of the folder contain AEROPNET sites for a specify case. 
    
    -starttime/endtime: select data from start date to end date, format in "YYYY-MM-DD".
    
    Return: 
    a dictionary contains all sites. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2018-10-18
    """
    
    
    """
    Initialization
    """
    aerCaseDir = '/nobackup/users/sunj/AERONET/%s/' % caseName     
    filelist = glob.glob(aerCaseDir + '*.all')
    Data = pd.DataFrame()

    for i, ff in enumerate(sorted(filelist)[:]): 
        sys.stdout.write('\r Reading AERONET inversion # %i/%i sites' % (i + 1, len(filelist)))
        """
        Site information: name, lat, lon, elevation    
        """
        aerData = pd.read_csv(ff, sep=",", header = 6)    
    
        """
        Date and time  
        """
        aerDate = np.array(aerData['Date(dd:mm:yyyy)'])
        aerTime = np.array(aerData['Time(hh:mm:ss)'])

        dateTime = []
        timeStamp = []
        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(aerDatetime)):
            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            year, month, day, hour, minute, second = aerDatetime[i][6:10], aerDatetime[i][3:5], aerDatetime[i][:2], aerDatetime[i][11:13], aerDatetime[i][14:16], aerDatetime[i][17:19]
            dateTime.append(pd.to_datetime('%s-%s-%s %s:%s:%s' % (year, month, day, hour, minute, second)))
        aerData['timeStamp'] = timeStamp
        aerData['dateTime'] = dateTime
        keywords = list(aerData.keys())

        """
        Select data in time period from start-date to end-date        
        """       
        mask = (aerData['dateTime'] >= startdate + ' 00:00:00') & (aerData['dateTime'] <= enddate + ' 23:59:59')
        aerData = aerData[mask] 
        
        """
        keys of interest
        """
        keywords = keywords[0:1] + keywords[5:75] + keywords[-8:]

        """
        Output  
        """
        Data = Data.append(aerData[keywords])

#        """
#        Wavelength and complex refractive index, SSA, AAOD, asymmetry factor
#        """
#        refr = []
#        refi = []
#        ssa = []
#        aaod = []
#        asy = []
#        Angstorm = []
#        size = []
#        wvl = []
#        for ikey in keywords:
#            if ikey.find('REFR') == 0:
#                refr.append(ikey)
#                wvl.append(float(re.findall(r'\d+', ikey)[0]))
#            if ikey.find('REFI') == 0:
#                refi.append(ikey)
#            if ikey.find('SSA') == 0:
#                ssa.append(ikey)
#            if (ikey.find('AOTAbsp') == 0):
#                aerData[ikey][aerData[ikey]<0] = np.nan
#                aaod.append(ikey)
#            if ikey.find('ASYM') == 0:
#                asy.append(ikey)
#            if ikey.find('Angstrom') >= 0:
#                Angstorm.append(ikey)
#            try:
#                size.append(float(ikey))
#            except ValueError:
#                pass
#
#        parameters = refr + refi + ssa + aaod + asy + Angstorm
#        
#        """
#        Size distribution function
#        """
#        def float2str(data):
#            return '%1.6f' % (data)
#        sizekey = list(map(float2str, size))
#
#        prob = aerData[sizekey].values
#        size = np.array(size)
#        # seperate fine and coarse mode         
#        fine = np.ones(len(prob)) * np.nan
#        coarse = np.ones(len(prob)) * np.nan
#        Cv_f = np.ones(len(prob)) * np.nan
#        Cv_c = np.ones(len(prob)) * np.nan
#        for k in range(len(prob)):
#            peaks = argrelextrema(prob[k,:], np.greater)[0]
#            if len(peaks) > 0: 
#                fine[k] = size[peaks[0]]
#                coarse[k] = size[peaks[-1]]
#                fcidx = argrelextrema(prob[k,:],np.less)[0]
#                if len(fcidx) == 1:
#                    if fine[k] < 1: 
##                        print('Only fine mode', k)
#                        Cv_c[k] = np.nan
#                        Cv_f[k] = 1 - Cv_c[k]
#                        coarse[k] = np.nan
#                    else:
##                        print('Only coarse mode', k)
#                        Cv_f[k] = np.nan
#                        Cv_c[k] = 1
#                        fine[k] = np.nan
#                if len(fcidx) > 1:    
#                    fcidx = argrelextrema(prob[k,:],np.less)[0][0]
#                    Cv_f[k] = trapz(prob[k, 0:fcidx+1],dx = 0.01)
#                    Cv_c[k] = trapz(prob[k, fcidx:-1], dx = 0.01) 
#            else:
#                # no data
#                fine[k] = np.nan
#                coarse[k] = np.nan
#                Cv_f[k] = np.nan
#                Cv_c[k] = np.nan
#        # convert volume to number density
#        sigma_f = 1.5
#        sigma_c = 2.0 
#        r_f = fine / np.exp(3*(np.log(sigma_f)**2))
#        r_c = coarse / np.exp(3*(np.log(sigma_c)**2))
#        
#        try: 
#            rfc_Cn = (Cv_f / Cv_c) * (r_c / r_f)**3 * np.exp(-4.5 * (sigma_f**2 - sigma_c**2))  
#            rfc_Cv = Cv_f / Cv_c
#            wfc_n = rfc_Cn / (rfc_Cn + 1.) 
#            wfc_v = rfc_Cv / (rfc_Cv + 1.)
#        except RuntimeWarning:
#            print('Divided by zero!')
#        
#        temp = np.c_[r_f, r_c, wfc_n, wfc_v]
#        sizefunc = pd.DataFrame(temp, columns = ['rf', 'rc', 'wnum', 'wvol'])

#        data = pd.concat([timeseries, aerData[parameters], sizefunc], axis = 1).iloc[idxDate]
#        
#        
#        data.index = data['dateTime']
#        del data['dateTime']
#        
#        data = data.dropna(axis = 1, how = 'all')
#        para = {'lat': aerlat, 'lon': aerlon, 'elev': aerelv, 'data': data}
#        Data[aername] = para
    Data[Data == -999] = np.nan
    return Data
    



# =============================================================================
# direct sun products
# =============================================================================
def AERONETdirectSun(caseName, startdate, enddate):
    """
    Function to read AERONET direct sun product.
    
    -caseName: name of the folder contain AEROPNET sites for a specify case. 
    
    -starttime/endtime: select data from start date to end date, format in "YYYY-MM-DD".
    
    Return: 
    a dictionary contains all sites. 
    Each site is a sub-dictionary contains elements of name, lat, lon, elevation, wavelengths and data.
    Data is in format of dataframe.
    
    
    @author: Sunji
    Last updated date: 2018-10-18
    """

    """
    Initialization
    """
    aerCaseDir = '/nobackup/users/sunj/AERONET/%s/' % caseName    
    filelist = glob.glob(aerCaseDir + '*.lev15*')
    Data = pd.DataFrame()
    
    for i, ff in enumerate(sorted(filelist)[:]):
        sys.stdout.write('\r Reading AERONET direct sun # %i/%i sites' % (i + 1, len(filelist)))
        """
        Site information: name, lat, lon, elevation
        """
        aerData = pd.read_csv(ff, sep=",", header = 6)    

        """
        Date and time  
        """
        aerDate = np.array(aerData['Date(dd:mm:yyyy)'])
        aerTime = np.array(aerData['Time(hh:mm:ss)'])
        
        dateTime = []
        timeStamp = []
        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(aerDatetime)):
            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            year, month, day, hour, minute, second = aerDatetime[i][6:10], aerDatetime[i][3:5], aerDatetime[i][:2], aerDatetime[i][11:13], aerDatetime[i][14:16], aerDatetime[i][17:19]
            dateTime.append(pd.to_datetime('%s-%s-%s %s:%s:%s' % (year, month, day, hour, minute, second)))
        aerData['timeStamp'] = timeStamp
        aerData['dateTime'] = dateTime
        keywords = list(aerData.keys())

        """
        Select data in time period from start-date to end-date        
        """       
        mask = (aerData['dateTime'] >= startdate + ' 00:00:00') & (aerData['dateTime'] <= enddate + ' 23:59:59')
        aerData = aerData[mask] 
        
        """
        keys of interest
        """
        for ikey in sorted(keywords):
            if (ikey.find('Triplet_Variability') >= 0) | (ikey.find('Exact_Wavelengths_') >= 0):
                keywords.remove(ikey)
        """
        Output
        """
        Data = Data.append(aerData[keywords])
    Data[Data == -999] = np.nan
    return Data



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
    

#ROI = {'S':-90, 'N': 90, 'W': -180, 'E': 180}
#casedir = '/nobackup/users/sunj/AERONET/AERONET_global_yearly/'
#
#t1 = time.time()
#for iyear in np.arange(2005, 2017):
#    caseName = 'Global_2005-2017'
#    
#    startdate = '%4i-%02i-%02i' % (iyear, 1, 1)
#    enddate   = '%4i-%02i-%02i' % (iyear, 12, 31) 
#    
#    INV = AERONETinversion(caseName, startdate, enddate)
#    DS = AERONETdirectSun(caseName, startdate, enddate)
#    
#    INV_int = AERONETwvlProcess(INV, ['SSA', 'AOTAbsp'], [388, 440, 500, 532, 550], 'linear')
#    DS_int = AERONETwvlProcess(DS, ['AOT'], [388, 440, 500, 532, 550], 'linear')
#    
#    with open(casedir + 'INV_%4i.pickle' % (iyear), 'wb') as handle:
#        pickle.dump(INV_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
#    with open(casedir + 'DS_%4i.pickle' % (iyear), 'wb') as handle:
#        pickle.dump(DS_int, handle, protocol=pickle.HIGHEST_PROTOCOL)
#        
#t2 = time.time()
#print('Time: %1.2f s' % (t2 - t1))        
