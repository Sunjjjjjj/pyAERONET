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





dataOutputDir = '/nobackup/users/sunj/'
dataInputDir = '/nobackup_1/users/sunj/'


def selectAERONET(caseName, startdate, enddate, ROI):
# =============================================================================
#   Initialization
# =============================================================================
    invDir = dataInputDir + 'AERONET/INV_Level15_All_Points_V3/INV/LEV15/ALL/ALL_POINTS/'
    dsDir = dataInputDir + 'AERONET/AOD_Level15_All_Points_V3/AOD/AOD15/ALL_POINTS/'
    filelistINV = glob.glob(invDir + '*.all')
    filelistDS = glob.glob(dsDir + '*lev*')
    # make case directory to restore selected AERONET site
    caseDir = dataInputDir + 'AERONET/%s/' % (caseName)
    if not os.path.isdir(caseDir):
        os.makedirs(caseDir)   
# =============================================================================
#   Select inversion
# =============================================================================
    t1 = time.time()
    sys.stdout.write('\r Searching site from AERONET inversion product...')
    sites = []
    for ifile in filelistINV:
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
#   Copy the sites within the ROI and TOI into the caseDir
# =============================================================================
        if (lon >= ROI['W']) & (lon <= ROI['E']) & (lat >= ROI['S']) & (lat <= ROI['N']):
            if any((data['dateTimeLocal'] >= startdate) & (data['dateTimeLocal'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))):
                shutil.copy2(ifile, caseDir)
                sites.append(site)
    
    sys.stdout.write('\r Total %i sites of inversion product selected' % (len(sites)))
# =============================================================================
#   Select direct sun. 
# =============================================================================
    sys.stdout.write('\r Searching site from AERONET direct sun product...')
    for ifile in filelistDS:
        data = pd.read_csv(ifile, sep=",", header = 6)
        site = data['AERONET_Site_Name'][0]
        
        if site in sites: 
            aerDate = data['Date(dd:mm:yyyy)']
            aerTime = data['Time(hh:mm:ss)']
            lag = round(data['Site_Longitude(Degrees)'] / 15)
            
            dateTime = []
            dateTimeLocal = []
            for i in range(len(data)):
                temp = pd.to_datetime('%s %s' %(aerDate[i], aerTime[i]), format ='%d:%m:%Y %H:%M:%S')
                dateTime.append(temp)
                dateTimeLocal.append(temp + np.sign(lag[i]) * pd.Timedelta(abs(lag[i]), 'h'))
            data['dateTime'] = dateTime
            data['dateTimeLocal'] = dateTimeLocal
            
            if any((data['dateTimeLocal'] >= startdate) & (data['dateTimeLocal'] <= pd.to_datetime(enddate) + pd.Timedelta(days = 1))):
                shutil.copy2(ifile, caseDir)
            #If direct sun is not available for the same site, then the corresponding inversion is removed.
            else:
                subprocess.call('rm ' + caseDir + '*%s*' % (site), shell = True)
                pass
        #If direct sun is not available for the same site, then the corresponding inversion is removed.
#        else:
#            subprocess.call('rm ' + caseDir + '*%s*' % (site), shell = True)
# =============================================================================
# check rudundent Inversion data
# =============================================================================
    dslist = glob.glob(caseDir + '*lev*')
    invlist = glob.glob(caseDir + '*.all')
    
    print('Total %i sites selected' %(len(dslist) + len(invlist)))
    print('Number of direct sun data %i' % (len(dslist)))  
    print('Number of inversion data %i' % (len(invlist)))
    t2 = time.time() 
    print('Time of selecting:',t2 - t1,'s')
    

# =============================================================================
# Case information 
# =============================================================================
#caseName = 'Global_2005-2018_v3'
#caseName = 'CA201712'
#caseName = 'AU201903'
#caseName = 'CA201811'
#caseName = 'SA201901'
#cases = ['ECN', 'NAF', 'CPC', 'ATA', 'SAF']


#ROI = {'S':-90, 'N': 90, 'W': -180, 'E': 180}
#ROI = {'S':30, 'N': 42.5, 'W': -130, 'E': -117.5}
#ROI = {'S': -50, 'N': -20, 'W': 130, 'E': 160}
#ROI = {'S': 40, 'N': 60, 'W': -125, 'E': -80}
#ROI = {'S': 20, 'N': 50, 'W': -140, 'E': -115}
#ROI = {'S': -15, 'N': 15, 'W': -20, 'E': 15}
ROI = {'S': 40, 'N': 60, 'W': -120, 'E': -80}

ROIs = np.load(dataOutputDir + 'P3_output/ROIs.npy').reshape(-1)[0]




startdate = '%4i-%02i-%02i' % (2006, 1, 1)
enddate   = '%4i-%02i-%02i' % (2016, 12, 31)

# =============================================================================
# AERONET site information
# =============================================================================
for caseName in ['CPC']:
    ROI = ROIs[caseName]
    selectAERONET(caseName, startdate, enddate, ROI)

