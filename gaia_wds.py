# -*- coding: utf-8 -*-
"""gaia_wds

__author__="Daphne Zakarian, Bob Zavala" 
__credits__=["Stephen Williams, Rachel Matson"]
__copyright__="A work of the United States Government"
__version__="0.1.5"
__maintainer_="Bob Zavala"
__email__="robert.t.zavala.civ@us.navy.mil"
__status__="Prototype"

Module contains code used for a project to match the WDS catalog against 
the Gaia DR3 catalog. The goal is to use Gaia data to provide the astrometric 
and photometric data to determine which entries in the WDS are true physical 
systems and which systems or components in those systems are optical components. 

Only WDS entries with precise (arcsec) coordinates will be matched against Gaia.  

Coding began in the summer of 2022 with REU student Daphne Zakarian (DZ, then an 
undergraduate at Truman State University) working with Stephen Williams and Bob
Zavala (BZ) at NOFS. 

BZ undertook subsequent development to improve and develop upon DZ's work. 

Module consists of:

    1) gaiaTable: read in the WDS catalog and after some database preparation
    writes out files in formats suitable for the project especially a VOTable
    to query against Gaia DR3. 

"""

import pandas as pd
from astropy.table import Table, Column, MaskedColumn
from astropy import units as u
from astropy.coordinates import ICRS, SkyCoord, FK5
from astropy.io import ascii
import math
import numpy as np
import os

# to time task completion only
import time

def gaiaTable():

    """
    Original Created on Thu 09 Jun 11:56:21 MDT 2022
    @author: Daphne Zakarian
    @author: Bob Zavala

    Read in the WDS using pandas, turn into an astropy table, find secondary coords, 
    account for proper motion, and write out the table as CSV, ECSV and a
    VOTable. The intent is that the VOTable is used for the query against Gaia
    and the CSV and ECSV are provided for ease of access. 

    Final table will be used to query Gaia -- most importantly, we have primary 
    and secondary stars' coordinates in degrees as the WDS provides precision 
    coordinates in sexagesimal format. 

    Changes began on Thu 30 May 09:08:00 MDT 2024 
    by Bob Zavala to include Discoverer and Components in order to 
    pass unique identifier information along. Coordinates and WDS ID alone
    do not uniquely identify the WDS entries. 

    Requirements: 
    In the current working directory the code looks for an ASCII
    copy of the WDS called: wdsweb_summ2.txt

    Output:
    Three (3) files in total saved to the current working directory.

    wdstab_new.ecsv: A CSV file with metadata 
    wdstab_new.csv:  A CSV file, no metadata
    wdstab_new.vot   A VOTable for the query against Gaia DR3 suitable for 
                     use for example on the ESA Gaia data server. 

    """



    this_path = os.getcwd()

    # Time the code
    startTime    = time.time()
    startCpuTime = time.process_time()

    print('\n')
    print('Welcome to the gaiaTable module! This takes about 5 minutes (+/-).')
    print('Starting to build the WDS tables for the Gaia queries.\n')
    print('**********')
    print('  Any ErfaWarnings from function pmsafe are normal and appear because the WDS lacks distances.')
    print('**********')

    # Read in Table: fixed width file to pandas to astropy 

    # manually list WDS txt file column widths to read in 
    # the formatted data. A change made by Daphane is to 
    # separate the DM column into two parts and the WDS 
    # 2000 arcsecond coordinates into RA and DEC portions. 
    ## IMPORTANT: colspecs are HALF-OPEN intervals: [,)
    ## SO THE RIGHT SIDE IS NOT INCLUSIVE!!!
    columns = ((0,10),(10,17),(17,23),(23,28),(28,32),(33,37),
            (37,42),(42,45),(45,51),(51,57),(57,64),(64,68),
            (69,80),(80,84),(84,88),(88,93),(93,98),(98,101),
            (101,106),(107,111),(112,121),(121,130))

    oldnames= ['WDS_Identifier', 'Discovr', 'Comp', 'Epoch-1', 'Epoch-2', '#', 
            'Theta-1', 'Theta-2', 'Rho-1', 'Rho-2', 'Mag-pri', 'Mag-sec',
            'SpecType','PMpri-RA', 'PMpri-DEC', 'PMsec-RA', 'PMsec-DEC', 'DM', 
            'Desig', 'Note', 'WDS_RA', 'WDS_DEC']

    # read fixed width file into pandas 
        # easier to work with this fixed-width-formatted (fwf) file 
        # in pandas than astropy
    wdspd = pd.read_fwf(this_path+"/wdsweb_summ2.txt",
                        colspecs=columns,header=None,names=oldnames,skiprows=3)

    # pandas table -> astropy table
    wdstab = Table.from_pandas(wdspd)

    print('\n')
    print('An astropy table has been created from the WDS data.')
    print('Converting precise WDS coordinates, and some other tasks. This takes a bit.')

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

    # Remove the WDS entries without precise (arcsec-level) 
    # coordinates

    c = []
    r=0
    for row in wdstab:
        if row['WDS_RA']== '.':
            c.append(r)
        r+=1

    wdstab.remove_rows(c)

    print('\n')
    print('I found ',len(c),' WDS entries without precise (arcsec) coordinates.')
    print('These imprecise entries will be excluded from the output.\n')


    # loop through the ra and dec of each row and make new columns for coords in 2 formats:
        # hms/dms format and degrees
    # setup counter for no coordinates
    noPrecise = 0
    for row in wdstab:
        try:
            ra1 = row['WDS_RA']
            dec1 = row['WDS_DEC']
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
            noPrecise+=1
            pass
    
    print('Coordinate conversion of primaries completed. \n')
        

    # FIND COORDINATE OF SECONDARY

    # make new columns for secondary coords
    wdstab['RAsecdeg']=0.0
    wdstab['DECsecdeg']=0.0
    wdstab['RAsecdeg'].unit=u.deg
    wdstab['DECsecdeg'].unit=u.deg

    # The for loop that the vector algorithm below replaced took 6.5 minutes 
    # on medea to run. 

    pricoord=SkyCoord(wdstab[0:]['RAprideg'],wdstab[0:]['DECprideg'], frame='icrs')
    angle = wdstab[0:]['Theta-2']*u.deg
    sep = wdstab[0:]['Rho-2']*u.arcsec
    seccoord = pricoord.directional_offset_by(angle, sep)

    wdstab['RAsecdeg'] = seccoord.ra.deg
    wdstab['DECsecdeg'] = seccoord.dec.deg

    print('Coordinates (RA/DEC) of secondaries found via offsets from primaries.\n')
    
    # Insert some missing units into the table
    # The for loop seems necessary and takes little time
    col_need_units = ['Epoch-1','Epoch-2','Theta-1','Theta-2','Rho-1','Rho-2','Mag-pri','Mag-sec',
                    'PMpri-RA','PMpri-DEC','PMsec-RA','PMsec-DEC']
    units_for_cols = [u.yr,u.yr,u.deg,u.deg,u.arcsec,u.arcsec,u.mag,u.mag,u.mas/u.yr,u.mas/u.yr,
                    u.mas/u.yr,u.mas/u.yr]
    counter = 3
    for a_col in col_need_units:
        # print(a_col)
        wdstab[a_col].unit = units_for_cols[counter-3]
        counter +=1

    # Change the '--' to 0.0 in the proper motion columns 
    # in order to use zero proper motions in PM calculations
    wdstab['PMpri-RA'].fill_value = 0.0
    wdstab['PMpri-RA'] = wdstab['PMpri-RA'].filled(0.0)
    wdstab['PMpri-DEC'].fill_value = 0.0
    wdstab['PMpri-DEC'] = wdstab['PMpri-DEC'].filled(0.0)
    wdstab['PMsec-RA'].fill_value=0.0
    wdstab['PMsec-RA'] = wdstab['PMsec-RA'].filled(0.0)
    wdstab['PMsec-DEC'].fill_value=0.0
    wdstab['PMsec-DEC'] = wdstab['PMsec-DEC'].filled(0.0)

    # ACCOUNT FOR PROPER MOTION
    # Need to add units for PMpri-RA,PMpri-DEC,PMsec-RA,PMsec-DEC as 
    # these did not make it into the table.

    # Change units of pm to deg/year
    wdstab['PMpri-RAdeg'] = wdstab['PMpri-RA'].to(u.deg/u.year)
    wdstab['PMpri-DECdeg'] = wdstab['PMpri-DEC'].to(u.deg/u.year)
    wdstab['PMsec-RAdeg'] = wdstab['PMsec-RA'].to(u.deg/u.year)
    wdstab['PMsec-DECdeg'] = wdstab['PMsec-DEC'].to(u.deg/u.year)

    # Apply proper motions to primaries -> 2016.0
    priRA    = wdstab['RAprideg']
    priDEC   = wdstab['DECprideg']
    pmpriRA  = wdstab['PMpri-RAdeg']
    pmpriDEC = wdstab['PMpri-DECdeg']

    coord = SkyCoord(priRA, priDEC, unit='deg', 
                    pm_ra_cosdec = pmpriRA, pm_dec = pmpriDEC)
    newcoord = coord.apply_space_motion(dt=16*u.yr)

    # While newcoord above IS a SkyCoord (astropy) object the conversion
    # to degrees from the hms default transforms newcoord -> numpy array.
    # I need to convert this numpy array back to an astropy column 
    # AND THEN insert the proper-motion adjusted coordinates into the 
    # table. Arrgh. 
    # The other thing I could not do was create the RApri-prepped 
    # and DECpri-prepped columns, prefill with zeroes and insert the new 
    # newcoord columns. 
    new_prira_col = Column(newcoord[0:].ra.deg)
    new_pridec_col = Column(newcoord[0:].dec.deg)
    # Here is how I finally entered the PM applied coordinates. 
    # Except I still need to fix the 'nan' proper motion prepped coordinates
    wdstab['RApri-prepped'] = new_prira_col
    wdstab['DECpri-prepped'] = new_pridec_col

    # Apply proper motions to secondaries -> 2016.0
    secRA    = wdstab[0:]['RAsecdeg']
    secDEC   = wdstab[0:]['DECsecdeg']
    pmsecRA  = wdstab[0:]['PMsec-RAdeg']
    pmsecDEC = wdstab[0:]['PMsec-DECdeg']

    coord = SkyCoord(secRA, secDEC, unit='deg', 
                    pm_ra_cosdec = pmsecRA, pm_dec = pmsecDEC)
    newcoord = coord.apply_space_motion(dt=16*u.yr)

    new_secra_col = Column(newcoord[0:].ra.deg)
    new_secdec_col = Column(newcoord[0:].dec.deg)

    wdstab['RAsec-prepped']  = new_secra_col
    wdstab['DECsec-prepped'] = new_secdec_col

    # Insert more missing units into the table
    col_need_units = ['RApri-prepped','DECpri-prepped','RAsec-prepped','DECsec-prepped','RAsecdeg','DECsecdeg']

    the_unit = u.deg

    for a_col in col_need_units:
        #print(a_col)
        #print(wdstab[a_col].unit)
        wdstab[a_col].unit = the_unit
        #print(wdstab[a_col].unit)

    print('\n')
    print('Proper motions applied to 2016.0 for querying against Gaia DR3.\n')

    # Save this for the END as a P code indicates a Proper Motion note and 
    # other such details we may need later.
    wdstab.remove_columns(['#','SpecType','DM', 'Desig', 'Note',])

    print('Writing the three (3) output tables. Standby please.\n')

    # SAVE TABLE AS ECSV and CSV and as a VOTable
    # Two different ways of coding the ascii-write are 
    # shown as an example for how to write these and one is a 
    # little more concise.
    ascii.write(wdstab, this_path+'/wdstab_new.ecsv', format='ecsv', overwrite=True)
    wdstab.write(this_path+'/wdstab_new.csv', overwrite=True)
    wdstab.write(this_path+'/wdstab_new.vot', format = 'votable', overwrite=True)

    print('\n')
    print('Tables built and written as CSV, ECSV and VOTables.\n')

    endTime    = time.time()
    endCpuTime = time.process_time()

    totalElapTime = endTime - startTime
    totalCpuTime  = endCpuTime - startCpuTime

    print('Elapsed time for building and writing tables: {:.3f} seconds.'.format(totalElapTime))
    print('Elapsed CPU for building and writing tables: {:.3f} seconds.\n'.format(totalCpuTime))
    print('Finished, thank you for using gaiaTable!')

    print('\n')









