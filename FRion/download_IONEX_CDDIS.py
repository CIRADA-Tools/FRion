#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This file contains the code necessary to download IONEX files from CDDIS.
It is intended to replace the getIONEX functionality of RMextract, since 
that code does not currently support CDDIS (which requires authentication).

The majority of the code is taken directly from RMextract. This has been 
modified by Art Davydov and Anna Ordog to be compatible with CDDIS, and
further modified by me to streamline the process.

To fit with the existing RMextract integration, the main predict module
will attempt to pre-download the IONEX files using this code before invoking
RMextract. That way, RMextract will find the already-downloaded IONEX files
and will not try to download them itself, thus dodging all the authentication
problems.

Created on Tue Aug  1 14:09:25 2023
@author: cvaneck
"""

import datetime
import logging
import os
# AO, AD added these:
#import requests

logging.basicConfig(level=logging.ERROR)


def get_CDDIS_IONEXfile(time="2023/03/23/02:20:10.01",
                    prefix="jplg",
                    outpath='./',
                    overwrite=False):
    """Get IONEX file with prefix from server for a given day

    Downloads files with given prefix from CDDIS, unzips and stores
    the data.
    
    CDDIS requires authentication to download files. Thus must be done as a
    .netrc file with the credentials in them.
    (https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html)

    Args:
        time (string or list) : date of the observation
        prefix (string) : prefix of the IONEX files (case insensitive)
        outpath (string) : path where the data is stored
        overwrite (bool) : Do (not) overwrite existing data
    """
    
    server="https://cddis.nasa.gov"
    
    
    # AO, AD changed this to lower from upper to comply with cddis:
    prefix=prefix.lower()
    if outpath[-1] != "/":
        outpath += "/"
    if not os.path.isdir(outpath):
        try:
            os.makedirs(outpath)
        except:
            print("cannot create output dir for IONEXdata: %s",
                          outpath)

    try:
        yy = int(time[2:4])
        year = int(time[:4])
        month = int(time[5:7])
        day = int(time[8:10])
    except:
        year = time[0]
        yy = year - 2000
        month = time[1]
        day = time[2]
    mydate = datetime.date(year, month, day)
    dayofyear = mydate.timetuple().tm_yday
    if not overwrite and os.path.isfile("%s%s%03d0.%02dI"%(outpath,prefix.upper(),dayofyear,yy)):
        logging.info("FILE exists: %s%s%03d0.%02dI",outpath,prefix,dayofyear,yy)
        return "%s%s%03d0.%02dI"%(outpath,prefix,dayofyear,yy)
   
    #If proxy url is given, enable proxy using pysocks
    try:
        from urllib import request
    except ImportError:
        import urllib2 as request


    #Naming conventions changed, at different times for different data sets.
    #Depending on the data source/prefix, the filename has to be completely 
    #different for different date ranges. This needs to be manually coded,
    #which means supporting a limited range of data sources. I'll try to support
    #the main ones stored at CDDIS.
    
    if prefix == 'jplg' and mydate > datetime.date(2023,8,7):
        filename=f"{prefix[0:3].upper()}0OPSFIN_{year}{dayofyear}0000_01D_02H_GIM.INX.gz"
    elif prefix == 'codg' and mydate > datetime.date(2022,11,26):
        filename=f"{prefix[0:3].upper()}0OPSFIN_{year}{dayofyear}0000_01D_01H_GIM.INX.gz"
    elif prefix == 'igsg' and mydate > datetime.date(2022,11,26):
        filename=f"{prefix[0:3].upper()}0OPSFIN_{year}{dayofyear}0000_01D_02H_GIM.INX.gz"
    elif prefix == 'casg' and mydate > datetime.date(2022,12,31):
        filename=f"{prefix[0:3].upper()}0OPSFIN_{year}{dayofyear}0000_01D_30M_GIM.INX.gz"
    elif prefix == 'esag':
        if mydate <= datetime.date(2023,2,4):
            raise Exception("ESA did not publish GIM.INX.gz files prior to 2023 Feb 04.")
        else:
            filename=f"{prefix[0:3].upper()}0OPSFIN_{year}{dayofyear}0000_01D_02H_GIM.INX.gz"
    else:
        filename=f"{prefix}{dayofyear:03d}0.{yy:02d}i.Z"


    url = "https://cddis.nasa.gov/archive/gnss/products/ionex/%4d/%03d/%s"%(year,dayofyear,filename)

    # Download IONEX file, make sure output format is old-format name and uppercase.
    fname = outpath+'/'+(f"{prefix}{dayofyear:03d}0.{yy:02d}i.Z").upper()
    print("Downloading ",url)


    os.system(f'wget --auth-no-challenge -O {fname} "{url}"')

    
    ###### gunzip files
    if fname[-2:].upper()==".Z":
        command = "gunzip -dc %s > %s" % (fname, fname[:-2])
        retcode = os.system(command)
        if retcode:
            raise RuntimeError("Could not run '%s'" % command)
        else:
            os.remove(fname)
        fname=fname[:-2]
    #returns filename of uncompressed file

    return fname




