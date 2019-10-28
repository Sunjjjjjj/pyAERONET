# -*- coding: utf-8 -*-
"""
reading AERONET data

@author: sunj
"""


######### load python packages
import os, sys
import shutil
import numpy as np
import glob
import pandas as pd
import subprocess
import time
import datetime




# =============================================================================
# Case information 
# =============================================================================
#caseName = 'Global_2005-2018_v3'
#caseName = 'CA201712'
#caseName = 'AU201903'
#caseName = 'CA201811'
#caseName = 'SA201901'
caseName = 'CAN201905'

caseDir = '/nobackup/users/sunj/AERONET/%s/' % (caseName)
if not os.path.isdir(caseDir):
    os.makedirs(caseDir)   

#ROI = {'S':-90, 'N': 90, 'W': -180, 'E': 180}
#ROI = {'S':30, 'N': 42.5, 'W': -130, 'E': -117.5}
#ROI = {'S': -50, 'N': -20, 'W': 130, 'E': 160}
#ROI = {'S': 40, 'N': 60, 'W': -125, 'E': -80}
#ROI = {'S': 20, 'N': 50, 'W': -140, 'E': -115}
#ROI = {'S': -15, 'N': 15, 'W': -20, 'E': 15}
ROI = {'S': 40, 'N': 60, 'W': -120, 'E': -80}

#timeStart = '%4i-%02i-%02i %02i:%02i:%02i' % (2019, 5, 29, 0, 0, 0)
#timeEnd   = '%4i-%02i-%02i %02i:%02i:%02i' % (2019, 5, 30, 23, 59, 59)
startdate = '%4i-%02i-%02i' % (2019, 5, 29)
enddate   = '%4i-%02i-%02i' % (2019, 5, 30)

# =============================================================================
# AERONET site information
# =============================================================================
start = time.time() 
dataOutputDir = '/nobackup/users/sunj/'
dataInputDir = '/nobackup_1/users/sunj/'

invDir = dataInputDir + 'AERONET/INV_Level15_All_Points_V3/INV/LEV15/ALL/ALL_POINTS/'
dsDir = dataInputDir + 'AERONET/AOD_Level15_All_Points_V3/AOD/AOD15/ALL_POINTS/'


filelistINV = glob.glob(invDir + '*.all')
filelistDS = glob.glob(dsDir + '*lev*')



sys.stdout.write('\r Searching site from AERONET inversion product...' )
for ifile in filelistINV[:]:
    data = pd.read_csv(ifile, sep=",", header = 6)    

    site = data['Site'][0]
    lon = data['Longitude(Degrees)'][0]
    lat = data['Latitude(Degrees)'][0]
    aerDate = data['Date(dd:mm:yyyy)']
    aerTime = data['Time(hh:mm:ss)']
    lag = round(data['Longitude(Degrees)'] / 15)
    
    
    dateTime = []
    dateTimeLocal = []
    for i in range(len(data)):
        temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
        dateTime.append(temp)
        dateTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
    data['dateTime'] = dateTime
    data['dateTimeLocal'] = dateTimeLocal
# =============================================================================
# Copy the sites within the ROI to case folder
# =============================================================================
    if (lon >= ROI['W']) & (lon <= ROI['E']) & (lat >= ROI['S']) & (lat <= ROI['N']):
        if any((data['dateTimeLocal'] >= startdate) & (data['dateTimeLocal'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))):
#            shutil.copy2(ifile, caseDir)
            print(ifile)

## =============================================================================
## Delete the sites outside TOI from case folder
## =============================================================================
#filelist = glob.glob(caseDir + '*.all')
#count = 0
#for i, ff in enumerate(sorted(filelist)[:]):
#    
#    aerData = pd.read_csv(ff, sep=",", header = 6)    
#
#    aerDate = aerData['Date(dd:mm:yyyy)']
#    aerTime = aerData['Time(hh:mm:ss)']
#    
#    timeStamp = []
#    aerDatetime = aerDate + ' ' + aerTime
#    for i in range(len(aerDatetime)):
#        timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
#    
#    tsstart = np.array(time.mktime(datetime.datetime.strptime(timestart, '%Y-%m-%d %H:%M:%S').timetuple()))
#    tsend  =  np.array(time.mktime(datetime.datetime.strptime(timeend, '%Y-%m-%d %H:%M:%S').timetuple()))
#    tsrange = (timeStamp >= tsstart) & (timeStamp <= tsend)
#    
#    if len(aerData[tsrange]) == 0: 
#        os.remove(ff)
#        count += 1 
#
#num = len(filelist) - count
#sys.stdout.write('\r Total %i inversion sites selected' %(num) )
#
## =============================================================================
## copy the AOT data of the same site
## =============================================================================
#filelist = glob.glob(caseDir + '*.all')
#dsName = []
#print('Searching direct sun product...')
#for i in sorted(filelist):
#    idx = i.find('19930101') 
#    dsName.append(i[idx + 18: -4])
#    
#count = 0
##for i, ff in enumerate(sorted(filelist)[:]):
##    aerData = pd.read_csv(ff, sep=",", header = 6)    
##    siteName = aerData['Site'][0]
##    print(siteName)
##
##    if set([siteName]) & set(dsName) == set([siteName]): 
#temp = []
#for siteName in dsName:
#    idx = dslist[0].find('19930101')
#    dsf = dslist[0][:idx+18] + siteName + '.lev15'
#    idx = filelist[0].find('19930101')
#    invf = filelist[0][:idx+18] + siteName + '.all'
#    try:
#        shutil.copy2(dsf, caseDir)
#    except:
#        os.remove(invf)
## =============================================================================
## check rudundent Inversion data
## =============================================================================
#dslist = glob.glob(caseDir + '*lev*')
#invlist = glob.glob(caseDir + '*.all')
#
#num += count            
#print('Total %i sites selected' %(num))
#print('Number of direct sun data %i' % (len(dslist)))  
#print('Number of inversion data %i' % (len(invlist)))
#end = time.time() 
#print('Time of selecting:',end - start,'s')
#
#    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    