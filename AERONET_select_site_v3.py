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
timestart = '%4i-%02i-%02i' % (2019, 5, 29)
timeend   = '%4i-%02i-%02i' % (2019, 5, 30)

# =============================================================================
# AERONET site information
# =============================================================================
start = time.time() 
invDir = '/nobackup/users/sunj/AERONET/INV_Level15_All_Points_V3/INV/LEV15/ALL/ALL_POINTS/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340
dsDir = '/nobackup/users/sunj/AERONET/AOD_Level15_All_Points_V3/AOD/AOD15/ALL_POINTS/'     #AERONET Santiago_Beaichef (S33.46 W7066 560m) @340


invlist = glob.glob(invDir + '*.all')
dslist = glob.glob(dsDir + '*lev*')

#subprocess.call('rm ' + caseDir + '*.lev', shell = True)
#subprocess.call('rm ' + caseDir + '*.lev*', shell = True)
sys.stdout.write('\r Searching site from AERONET inversion product...' )
for invf in invlist[:]:
    aerData = pd.read_csv(invf, sep=",", header = 6)    

    siteInfo = aerData['Site'][0]
    aerlon = aerData['Longitude(Degrees)'][0]
    aerlat = aerData['Latitude(Degrees)'][0]
    aerelv = aerData['Elevation(m)'][0]

# =============================================================================
# Copy the sites within the ROI to case folder
# =============================================================================
    if (aerlon > ROI['W']) and (aerlon < ROI['E']) and (aerlat > ROI['S']) and (aerlat < ROI['N']) :
        shutil.copy2(invf, caseDir)                         

# =============================================================================
# Delete the sites outside TOI from case folder
# =============================================================================
filelist = glob.glob(caseDir + '*.all')
count = 0
for i, ff in enumerate(sorted(filelist)[:]):
    
    aerData = pd.read_csv(ff, sep=",", header = 6)    

    aerDate = aerData['Date(dd:mm:yyyy)']
    aerTime = aerData['Time(hh:mm:ss)']
    
    timeStamp = []
    aerDatetime = aerDate + ' ' + aerTime
    for i in range(len(aerDatetime)):
        timeStamp.append(time.mktime(datetime.datetime.strptime(aerDatetime[i], '%d:%m:%Y %H:%M:%S').timetuple()))
    
    tsstart = np.array(time.mktime(datetime.datetime.strptime(timestart, '%Y-%m-%d %H:%M:%S').timetuple()))
    tsend  =  np.array(time.mktime(datetime.datetime.strptime(timeend, '%Y-%m-%d %H:%M:%S').timetuple()))
    tsrange = (timeStamp >= tsstart) & (timeStamp <= tsend)
    
    if len(aerData[tsrange]) == 0: 
        os.remove(ff)
        count += 1 

num = len(filelist) - count
print('Total %i inversion sites selected' %(num))
# =============================================================================
# copy the AOT data of the same site
# =============================================================================
filelist = glob.glob(caseDir + '*.all')
dsName = []
print('Searching direct sun product...')
for i in sorted(filelist):
    idx = i.find('19930101') 
    dsName.append(i[idx + 18: -4])
    
count = 0
#for i, ff in enumerate(sorted(filelist)[:]):
#    aerData = pd.read_csv(ff, sep=",", header = 6)    
#    siteName = aerData['Site'][0]
#    print(siteName)
#
#    if set([siteName]) & set(dsName) == set([siteName]): 
temp = []
for siteName in dsName:
    idx = dslist[0].find('19930101')
    dsf = dslist[0][:idx+18] + siteName + '.lev15'
    idx = filelist[0].find('19930101')
    invf = filelist[0][:idx+18] + siteName + '.all'
    try:
        shutil.copy2(dsf, caseDir)
    except:
        os.remove(invf)
# =============================================================================
# check rudundent Inversion data
# =============================================================================
dslist = glob.glob(caseDir + '*lev*')
invlist = glob.glob(caseDir + '*.all')

num += count            
print('Total %i sites selected' %(num))
print('Number of direct sun data %i' % (len(dslist)))  
print('Number of inversion data %i' % (len(invlist)))
end = time.time() 
print('Time of selecting:',end - start,'s')

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    