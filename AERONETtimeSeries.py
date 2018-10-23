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
    filelist = glob.glob(aerCaseDir + '*.dubovik')
    Data = {}

    for i, ff in enumerate(sorted(filelist)[:]): 
        sys.stdout.write('\r Reading AERONET inversion # %i/%i sites' % (i + 1, len(filelist)))
        """
        Site information: name, lat, lon, elevation    
        """
        aerData = pd.read_csv(ff, sep=",", header = 3)    
        keywords = aerData.keys() 
    
        aerHeader = open(ff, 'r')
        lines = aerHeader.readlines()
        siteInfo = lines[0]
    
        idx1 = siteInfo.find('Locations')
        idx2 = siteInfo.find(',Nmeas')
        siteInfo = siteInfo[idx1:idx2]
        idx0 = siteInfo.find('Locations=')
        idx1 = siteInfo.find('long=')
        idx2 = siteInfo.find('lat=')
        idx3 = siteInfo.find('elev=')
        
        aername = siteInfo[idx0 + 10 : idx1 - 1]
        aerlon = float(siteInfo[idx1 + 5 : idx2 - 1])
        aerlat = float(siteInfo[idx2 + 4 : idx3 - 1])
        aerelv = float(siteInfo[idx3 + 5 :])
        """
        Date and time  
        """
        aerDate = np.array(aerData['Date(dd-mm-yyyy)'])
        aerTime = np.array(aerData['Time(hh:mm:ss)'])
        jDay = np.array(aerData['Julian_Day'])
        
        dateTime = []
        timeStamp = []
        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(aerDatetime)):
            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            dateTime.append(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S'))

        jDay = pd.DataFrame(jDay, columns = ['julian'])
        dateTime = pd.DataFrame(dateTime, columns = ['dateTime'])
        timeStamp = pd.DataFrame(timeStamp, columns = ['timeStamp'])
        timeseries = pd.concat([dateTime, jDay, timeStamp], axis = 1)
        """
        Select data in time period from start-date to end-date        
        """       
        dates = pd.date_range(startdate, enddate).to_datetime()

        idxDate = [] 
        for idate in dates:
            idx = np.where(aerDate == '%02i:%02i:%4i' % (idate.day, idate.month, idate.year))[0]
            if not idx.any():
                pass
#                print('\n %s: there is no measurements on %02i:%02i:%4i!' % (aername, idate.day, idate.month, idate.year))
            else:
                idxDate += list(idx)

        """
        Wavelength and complex refractive index, SSA, AAOD, asymmetry factor
        """
        refr = []
        refi = []
        ssa = []
        aaod = []
        asy = []
        Angstorm = []
        size = []
        wvl = []
        for ikey in keywords:
            if ikey.find('REFR') == 0:
                refr.append(ikey)
                wvl.append(float(re.findall(r'\d+', ikey)[0]))
            if ikey.find('REFI') == 0:
                refi.append(ikey)
            if ikey.find('SSA') == 0:
                ssa.append(ikey)
            if (ikey.find('AOTAbsp') == 0):
                aaod.append(ikey)
            if ikey.find('ASYM') == 0:
                asy.append(ikey)
            if ikey.find('Angstrom') >= 0:
                Angstorm.append(ikey)
            try:
                size.append(float(ikey))
            except ValueError:
                pass

        parameters = refr + refi + ssa + aaod + asy + Angstorm
        
        """
        Size distribution function
        """
        def float2str(data):
            return '%1.6f' % (data)
        sizekey = list(map(float2str, size))

        prob = aerData[sizekey].values
        size = np.array(size)
        # seperate fine and coarse mode         
        fine = np.ones(len(prob)) * np.nan
        coarse = np.ones(len(prob)) * np.nan
        Cv_f = np.ones(len(prob)) * np.nan
        Cv_c = np.ones(len(prob)) * np.nan
        for k in range(len(prob)):
            peaks = argrelextrema(prob[k,:], np.greater)[0]
            if len(peaks) > 0: 
                fine[k] = size[peaks[0]]
                coarse[k] = size[peaks[-1]]
                fcidx = argrelextrema(prob[k,:],np.less)[0]
                if len(fcidx) == 1:
                    if fine[k] < 1: 
#                        print('Only fine mode', k)
                        Cv_c[k] = np.nan
                        Cv_f[k] = 1 - Cv_c[k]
                        coarse[k] = np.nan
                    else:
#                        print('Only coarse mode', k)
                        Cv_f[k] = np.nan
                        Cv_c[k] = 1
                        fine[k] = np.nan
                if len(fcidx) > 1:    
                    fcidx = argrelextrema(prob[k,:],np.less)[0][0]
                    Cv_f[k] = trapz(prob[k, 0:fcidx+1],dx = 0.01)
                    Cv_c[k] = trapz(prob[k, fcidx:-1], dx = 0.01) 
            else:
                # no data
                fine[k] = np.nan
                coarse[k] = np.nan
                Cv_f[k] = np.nan
                Cv_c[k] = np.nan
        # convert volume to number density
        sigma_f = 1.5
        sigma_c = 2.0 
        r_f = fine / np.exp(3*(np.log(sigma_f)**2))
        r_c = coarse / np.exp(3*(np.log(sigma_c)**2))
        
        try: 
            rfc_Cn = (Cv_f / Cv_c) * (r_c / r_f)**3 * np.exp(-4.5 * (sigma_f**2 - sigma_c**2))  
            rfc_Cv = Cv_f / Cv_c
            wfc_n = rfc_Cn / (rfc_Cn + 1.) 
            wfc_v = rfc_Cv / (rfc_Cv + 1.)
        except RuntimeWarning:
            print('Divided by zero!')
        
        temp = np.c_[r_f, r_c, wfc_n, wfc_v]
        sizefunc = pd.DataFrame(temp, columns = ['rf', 'rc', 'wnum', 'wvol'])


        """
        Output  
        """
        data = pd.concat([timeseries, aerData[parameters], sizefunc], axis = 1).iloc[idxDate]
        data.index = data['dateTime']
        del data['dateTime']
        para = {'lat': aerlat, 'lon': aerlon, 'elev': aerelv, 'wvl': np.array(wvl), 'data': data}
        Data[aername] = para
        
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
    filelist = glob.glob(aerCaseDir + '*.lev*')
    Data = {}

    for i, ff in enumerate(sorted(filelist)[:]):
        sys.stdout.write('\r Reading AERONET direct sun # %i/%i sites' % (i + 1, len(filelist)))
        """
        Site information: name, lat, lon, elevation
        """
        aerData = pd.read_csv(ff, sep=",", header = 4)    
        keywords = aerData.keys() 
    
        aerHeader = open(ff, 'r')
        lines = aerHeader.readlines()
        siteInfo = lines[2]
        idx1 = siteInfo.find('Location')
        idx2 = siteInfo.find(',Nmeas')
        siteInfo = siteInfo[idx1:idx2]
        idx0 = siteInfo.find('Location=')
        idx1 = siteInfo.find('long=')
        idx2 = siteInfo.find('lat=')
        idx3 = siteInfo.find('elev=')
        
        aername = siteInfo[idx0 + 9 : idx1 - 1]
        aerlon = float(siteInfo[idx1 + 5 : idx2 - 1])
        aerlat = float(siteInfo[idx2 + 4 : idx3 - 1])
        aerelv = float(siteInfo[idx3 + 5 :])
        """
        Date and time  
        """
        aerDate = np.array(aerData['Date(dd-mm-yy)'])
        aerTime = np.array(aerData['Time(hh:mm:ss)'])
        jDay = np.array(aerData['Julian_Day'])
        
        dateTime = []
        timeStamp = []
        aerDatetime = aerDate + ' ' + aerTime
        for i in range(len(aerDatetime)):
            timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
            dateTime.append(pd.to_datetime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S')))
            
        jDay = pd.DataFrame(jDay, columns = ['julian'])
        dateTime = pd.DataFrame(dateTime, columns = ['dateTime'])
        timeStamp = pd.DataFrame(timeStamp, columns = ['timeStamp'])
        timeseries = pd.concat([dateTime, jDay, timeStamp], axis = 1)
        """
        Select data in time period from start-date to end-date        
        """       
        dates = pd.date_range(startdate, enddate).to_datetime()
        idxDate = [] 
        for idate in dates:
            idx = np.where(aerDate == '%02i:%02i:%4i' % (idate.day, idate.month, idate.year))[0]
            if not idx.any():
                pass
#                print('\n %s: there is no measurements on %02i:%02i:%4i!' % (aername, idate.day, idate.month, idate.year))
            else:
                idxDate += list(idx)
        """
        Wavelengths
        """
        Angstrom = []
        aod = []
        wvl = []
        for ikey in keywords:
            if ikey.find('AOT_') == 0:
                aod.append(ikey)
                wvl.append(float(re.findall(r'\d+', ikey)[0]))
            if ikey.find('Angstrom') >= 0:
                Angstrom.append(ikey)
        
        parameters = aod + Angstrom
        """
        Output
        """
        data = pd.concat([timeseries, aerData[parameters]], axis = 1).iloc[idxDate]
        data.index = data['dateTime']
        del data['dateTime']
        para = {'lat': aerlat, 'lon': aerlon, 'elev': aerelv, 'wvl': np.array(wvl), 'data': data}
        Data[aername] = para
    
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
        sys.stdout.write('\r Time processing # %i/%i sites' % (i + 1, len(Data)))
        output[isite] = {'lat': Data[isite]['lat'], 'lon': Data[isite]['lon'], 'elev': Data[isite]['elev'], 'wvl': Data[isite]['wvl']}
        processed = pd.DataFrame()
        data = Data[isite]['data'].copy()
        
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
                processed = data.groupby('YY-MM-DD').mean()
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
        sys.stdout.write('\r Wavelength processing # %i/%i sites' % (i + 1, len(Data)))
        output[isite] = {'lat': Data[isite]['lat'], 'lon': Data[isite]['lon'], 'elev': Data[isite]['elev'], 'wvl': Data[isite]['wvl']}
        data = Data[isite]['data'].copy()
        processed = data.copy()
        for ipara in POI:
            columns = []
            for ikey in data.keys():
                if ikey.find(ipara) == 0: 
                    columns.append(ikey)
            paradata = data[columns]
            """
            Interpolation/extrapolation
            """
            x = Data[isite]['wvl']
            y = paradata
            try: 
                f = interpolate.interp1d(x, y, kind = wvlProcessMethod)
                for iwvl in WOI:
                    processed['%s%i' % (ipara, iwvl)] = f(iwvl)
            except: 
                pass
    """
    output
    """
    output[isite]['data'] = processed                
    return pd.DataFrame.from_dict(output)


# =============================================================================
# Test code
# =============================================================================
def main(): 
    caseName = 'CA201712'
    ROI = {'S':25, 'N': 50, 'W': -130, 'E': -110}
    
    startdate = '%4i-%02i-%02i' % (2017, 11, 1)
    enddate   = '%4i-%02i-%02i' % (2017, 12, 31) 
    ctime = datetime.datetime(2017,12,12,19,50,1)
    dtime = 1800
    
    # read data
    inv = AERONETinversion(caseName, startdate, enddate)
    ds = AERONETdirectSun(caseName, startdate, enddate)
    # time process 
    inv_mm = AERONETtimeProcess(inv, 'monthly')
    inv_pp = AERONETtimeProcess(inv, 'period', ctime, dtime)
    ds_dd = AERONETtimeProcess(ds, 'daily')
    # wavelength interpolation
    inv_dd_int = AERONETwvlProcess(inv_pp, ['SSA', 'AOTAbsp'], [550], wvlProcessMethod = 'linear')
    ds_int = AERONETwvlProcess(ds_dd, ['AOT'], [550], wvlProcessMethod = 'linear')


if __name__ == '__main__':
    main()
