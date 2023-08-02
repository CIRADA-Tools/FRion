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
import requests

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

    # Url of the primary server has the syntax "ftp://ftp.aiub.unibe.ch/CODE/YYYY/CODGDOY0.YYI.Z" where DOY is the day of the year, padded with leading zero if <100, and YY is the last two digits of year.
    # Url of the backup server has the syntax "ftp://cddis.gsfc.nasa.gov/gnss/products/ionex/YYYY/DOY/codgDOY.YYi.Z where DOY is the day of the year, padded with leading zero if <100, and YY is the last two digits of year.
    #try primary url

    try:
        primary = request.urlopen(server,timeout=30)
    except:
        logging.error('Server not responding')


    url = "https://cddis.nasa.gov/archive/gnss/products/ionex/%4d/%03d/%s%03d0.%02di.Z"%(year,dayofyear,prefix,dayofyear,yy)

    # Download IONEX file, make sure it is always uppercase
    fname = outpath+'/'+(url.split('/')[-1]).upper()
    print("Downloading ",url)

    # AO, AD put in this block:
    try:
        r = requests.get(url)
    except:
        logging.info("No files found on %s for %s",server,fname)
        return -1

    if (r.text[0:15] == '<!DOCTYPE html>') or (r.text[0:15] == 'HTTP Basic: Acc'):
        logging.info('CDDIS is requesting authentication; download failed.')
        return -1

    with open(fname, 'wb') as fd:
        for chunk in r.iter_content(chunk_size=1000):
            fd.write(chunk)

    
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





