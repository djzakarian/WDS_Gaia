
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  9 11:56:21 2022

Read in the WDS using pandas, turn into an astropy table, find secondary coords, 
account for proper motion etc.

Final table will be used to query Gaia -- most importantly, we have primary 
and secondary stars' coordinates in deg

@author: djzak
"""



import pandas as pd
from astropy.table import Table
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord, FK5
from astropy.io import ascii
import math
import numpy as np



""" Read in Table: fixed width file to pandas to astropy """

# manually list WDS txt file column widths
columns = ((0,13),(13,17),(17,23),(23,28),(28,32),(32,37),
           (37,42),(42,45),(45,51),(51,57),(57,64),(64,68),
           (69,80),(80,84),(84,88),(88,93),(93,98),
           (98,101),(101,106),(106,112),(112,121),(121,131))

# read fixed width file into pandas 
    # easier to work with this fwf in pandas than astropy
wdspd = pd.read_fwf(r"C:\Users\djzak\OneDrive\Documents\NOFS\wds_full.txt",
                    colspecs=columns, skiprows=2)

# pandas table -> astropy table
wdstab = Table.from_pandas(wdspd)


# RENAME COLUMNS:
    # make lists of the old and new names
    # use for loop to go through and replace old names
oldnames=wdstab.colnames 
newnames= ['WDS Identifier', 'Discovr', 'Comp', 'Epoch-1', 'Epoch-2', '#', 
           'Theta-1', 'Theta-2', 'Rho-1', 'Rho-2', 'Mag-pri', 'Mag-sec',
           'SpecType','PMpri-RA', 'PMpri-DEC', 'PMsec-RA', 'PMsec-DEC', 'DM', 
           'Desig', 'Note', 'WDS_RA"', 'WDS_DEC"']
counter = 0
for oldname in oldnames:
    wdstab[oldname].name = newnames[counter]
    counter+=1
    
    
    
""" REMOVE UNNECESSARY COLUMNS """
wdstab.remove_columns(['Discovr', 'Comp','#','SpecType','DM', 
'Desig', 'Note',])



###########################################################
# !!!  Don't run this after they have been removed !!!

# """ REMOVE THE ROWS WITHOUT COORDINATES """
# c = []
# r=0
# for row in wdstab:
#     if row['WDS_RA"']== '.':
#         c.append(r)
#     r+=1
        
# wdstab.remove_rows(c)


    


""" ASSIGN UNITS TO THE COLUMNS """

unitlist = [u.yr, u.yr, u.deg, u.deg, u.arcsec, u.arcsec, u.mag, u.mag, 
            u.mas/u.yr, u.mas/u.yr, u.mas/u.yr, u.mas/u.yr,'' ,'' ]

counter = 0
for column in wdstab.colnames:
    wdstab[column].unit = unitlist[counter]
    counter +=1


# Establish format for RA and DEC

# make new column to store updated RA and DEC -- is there a better way to
    # make column dtype longer so I don't have to manually type 15 spaces???
wdstab['RAprideg']=0.0
wdstab['DECprideg']=0.0
wdstab['RAhms']='               '
wdstab['DECdms']='               '
 

# assign units to the new columns for coordinates in degrees
wdstab['RAprideg'].unit=u.deg
wdstab['DECprideg'].unit=u.deg


# loop through the ra and dec of each row and make new columns for coords in 2 formats:
    # hms/dms format and degrees
for row in wdstab:
    try:
        ra1 = row['WDS_RA"']
        dec1 = row['WDS_DEC"']
        rastr = ra1[:2] + "h" + ra1[2:4]+"m"+ra1[4:]+"s"
        decstr = dec1[:3] + "d" + dec1[3:5]+"m"+dec1[5:]+'s'
        
        
        # put the strings of ra and dec into the table
        row['RAhms']= rastr
        row['DECdms']= decstr
        
        # put the coordinates into degree form - should work with GAIA
        coo = ICRS(rastr,decstr)
        # put the new coordinates into a column with a float dtype
        row['RAprideg']= coo.ra.deg
        row['DECprideg']=coo.dec.deg
    # skip rows that don't have coordinates
    except ValueError:
        pass




""" FIND COORDINATE OF SECONDARY  """

# make new columns for secondary coords
wdstab['RAsecdeg']=0.0
wdstab['DECsecdeg']=0.0
wdstab['RAsecdeg'].unit=u.deg
wdstab['DECsecdeg'].unit=u.deg


## OLD WAY
# # loop through rows and calculate approximate coordinates for the secondary star
# for row in wdstab:
#     try:
#         rho, theta = row['Rho-2']*u.arcsec, row['Theta-2']*u.deg
#         ra,dec = row['RAprideg'], row['DECprideg']
#         rhodeg = rho.to(u.deg).value
#         thetarad = theta.to(u.rad).value
#         row['RAsecdeg']= ra + rhodeg*math.sin(thetarad)
#         row['DECsecdeg'] = dec +rhodeg*math.cos(thetarad)
#     except ValueError:
#         pass
    

c=0
for row in wdstab:
    prira, pridec = row['RAprideg'], row['DECprideg']
    angle, sep = row['Theta-2']*u.deg, row['Rho-2']*u.arcsec
    pricoord = SkyCoord(prira*u.deg, pridec*u.deg, frame='icrs')
    seccoord= pricoord.directional_offset_by(angle, sep)
    row['RAsecdeg']= seccoord.ra.deg
    row['DECsecdeg'] = seccoord.dec.deg
    c+=1
    if c % 1000==0:
        print(c)





""" ACCOUNT FOR PROPER MOTION """
# Change units of pm to deg/year
wdstab['PMpri-RAdeg'] = wdstab['PMpri-RA'].to(u.deg/u.year)
wdstab['PMpri-DECdeg'] = wdstab['PMpri-DEC'].to(u.deg/u.year)
wdstab['PMsec-RAdeg'] = wdstab['PMsec-RA'].to(u.deg/u.year)
wdstab['PMsec-DECdeg'] = wdstab['PMsec-DEC'].to(u.deg/u.year)

# Make new columns for the final, Gaia-ready coordinates
wdstab['RApri-prepped']=0.0
wdstab['DECpri-prepped']=0.0
wdstab['RAsec-prepped']=0.0
wdstab['DECsec-prepped']=0.0

# assign unit (degrees) to the new columns
wdstab['RApri-prepped'].unit=u.deg
wdstab['DECpri-prepped'].unit=u.deg
wdstab['RAsec-prepped'].unit=u.deg
wdstab['DECsec-prepped'].unit=u.deg








#### OLD METHOD, EUCLIDEAN GEO ####
# # carry out primary and secondary coords in different loops
# # to avoid skipping primary blocks when secondary coords dont exist 
    

# # PRIMARY     
# for row in wdstab:
#     # bring in all proper motions [deg] and coordinates [deg]
#     pmpriRA, pmpriDEC = row['PMpri-RAdeg'], row['PMpri-DECdeg'],
    
#     # if there isn't a pm, just use J200 coords
#     # this is done for ra and dec separately because a few coords only have pm 
#     # in one 'direction'
#     if str(pmpriRA) == 'nan':
#         row['RApri-prepped'] = row['RAprideg'] 
#         if str(pmpriDEC) == 'nan':
#             row['DECpri-prepped'] = row['DECprideg']
#     elif str(pmpriDEC) == 'nan':
#         row['DECpri-prepped'] = row['DECprideg']

#     # otherwise, update for pm
#     else:
#     # updated coordinates will be of the foorm coord-prepped = coord +16*propermotion
        
#         row['RApri-prepped'] = row['RAprideg'] + 16 * row['PMpri-RAdeg']
#         row['DECpri-prepped'] = row['DECprideg'] + 16 * row['PMpri-DECdeg']


        

# # SECONDARY      
# for row in wdstab:
#     # bring in all proper motions [deg] and coordinates [deg]
#     pmsecRA, pmsecDEC = row['PMsec-RAdeg'], row['PMsec-DECdeg']
    
#     # if there isn't a pm, just use 2016 coords
#     if str(pmsecRA) == 'nan':
#         row['RAsec-prepped'] = row['RAsecdeg'] 
#         if str(pmsecDEC) == 'nan':
#             row['DECsec-prepped'] = row['DECsecdeg']
#     elif str(pmsecDEC) == 'nan':
#         row['DECsec-prepped'] = row['DECsecdeg']
        
#     # otherwise, update for pm
#     else:
#     # updated coordinates will be of the foorm coord-prepped = coord +16*propermotion
        
#         row['RAsec-prepped'] = row['RAsecdeg'] + 16 * row['PMsec-RAdeg']
#         row['DECsec-prepped'] = row['DECsecdeg'] + 16 * row['PMsec-DECdeg']








## pm_ra_cosdec ????? maybe right?

# my_coord = SkyCoord(0.027817777777777775, 75.48330000000001, unit='deg', pm_ra_cosdec=34*u.mas/u.year, pm_dec=5*u.mas/u.year)
# newcoord = my_coord.apply_space_motion(dt=16*u.yr)
# newcoord.ra.deg



# use the proper motions to precess coordinates to J2016 if pm info is available in Gaia
# PRIMARY

c=0   
for row in wdstab:
    # bring in all proper motions [deg] and coordinates [deg]
    pmpriRA, pmpriDEC = row['PMpri-RAdeg'], row['PMpri-DECdeg'],

    c+=1
    if c%1000 == 0:
        print(c)
            
    # if there isn't a pm, just use J200 coords
    # this is done for ra and dec separately because a few coords only have pm 
    # in one 'direction'
    if str(pmpriRA) == 'nan':
        row['RApri-prepped'] = row['RAprideg'] 
        if str(pmpriDEC) == 'nan':
            row['DECpri-prepped'] = row['DECprideg']
    elif str(pmpriDEC) == 'nan':
        row['DECpri-prepped'] = row['DECprideg']

    # otherwise, update for pm using SkyCoord
    else:
        coord = SkyCoord(row['RAprideg'], row['DECprideg'], unit='deg', 
                         pm_ra_cosdec = pmpriRA *u.deg/u.year, pm_dec = pmpriDEC *u.deg/u.year)
        newcoord = coord.apply_space_motion(dt=16*u.yr)
        
        row['RApri-prepped'], row['DECpri-prepped'] = newcoord.ra.deg, newcoord.dec.deg



# use the proper motions to precess coordinates to J2016 if pm info is available in Gaia
# SECONDARY       
c=0   
     
for row in wdstab:
    # bring in all proper motions [deg] and coordinates [deg]
    pmsecRA, pmsecDEC = row['PMsec-RAdeg'], row['PMsec-DECdeg']
    
    c+=1
    if c%1000 == 0:
        print(c)
    
    # if there isn't a pm, just use 2016 coords
    if str(pmsecRA) == 'nan':
        row['RAsec-prepped'] = row['RAsecdeg'] 
        if str(pmsecDEC) == 'nan':
            row['DECsec-prepped'] = row['DECsecdeg']
    elif str(pmsecDEC) == 'nan':
        row['DECsec-prepped'] = row['DECsecdeg']
        
    # otherwise, update for pm
    else:
        coord = SkyCoord(row['RAsecdeg'], row['DECsecdeg'], unit='deg', 
                         pm_ra_cosdec = pmsecRA *u.deg/u.year, pm_dec = pmsecDEC *u.deg/u.year)
        newcoord = coord.apply_space_motion(dt=16*u.yr)
        
        row['RAsec-prepped'], row['DECsec-prepped'] = newcoord.ra.deg, newcoord.dec.deg



# DR3 vs EDR3

# skycoord matching

### Space motion !!! 










""" SAVE TABLE AS ECSV and CSV """

ascii.write(wdstab, r'C:\Users\djzak\OneDrive\Documents\NOFS\wdstab6-27.ecsv', format='ecsv')
ascii.write(wdstab, r'C:\Users\djzak\OneDrive\Documents\NOFS\wdstab6-27.csv', format='csv')






""" OPEN AND READ ECSV """

from astropy.table import Table

# # My Computer 
# path = r"C:\Users\djzak\OneDrive\Documents\NOFS\wdstab6-27.ecsv"

# # NOFS Laptop
path = '/home/daphne/NOFS/wdstab6-27.ecsv'
wdstab = Table.read(path, header_start=0, data_start=1) 









