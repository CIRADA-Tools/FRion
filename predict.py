#!/usr/bin/env python3

"""
Tools for predicting time-integrated ionospheric Faraday rotation corrections.
Uses RMextract package to get time-dependent corrections for a given 
observation, then performs the time-integration to work out the effective
change in polarization angle, and the effective depolarization.


"""



import RMextract.getRM as RME
from astropy.time import Time
import numpy as np
from astropy.coordinates import SkyCoord, EarthLocation
import astropy.units as u


def get_telescope_coordinates(telescope_name):
    """Return the astropy.coordinates EarthLocation object associated
    with the position of the telescope named. Each telescope must be 
    hardcoded in.
    """
    if telescope_name == 'ASKAP':
        lat = -1*26+42/60+15/3600 #degree
        long = +1*116+39/60+32/3600 # degree
        height = 381.0 #
        return EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=height*u.m)
    else:
        raise Exception('Telescope position not known.')


def integrate_single_timestep(time1,time2,RM1,RM2):
    pass
    return ""



def write_correction(correction_dictionary,filename):
    pass



def calculate_correction(start_time, end_time, freq_array, telescope_location,
                         ra,dec, timestep=600.,ionexPath='./IONEXdata/'):
    """Calculate the ionospheric FR correction, as a function of frequency,
    for a given observation (time, location, target direction).
    Inputs:
        start_time (string readable by astropy.time.Time): starting
                    time of observation
        end_time (string readable by astropy.time.Time): ending
                    time of observation
            example time string: '2010-01-02T00:00:00'
        freq_array (array-like): vector of channel frequencies (in Hz)
        telescope_location (astropy.coordinates EarthLocation or string):
            location of telescope, or name of telescope known to 
            get_telescope_coordinates() function.
        ra (float): right ascension of observation center (in deg, J2000)
        dec (float): declination of observation center (in deg, J2000)
        timestep (float): time in seconds between ionosphere FR estimates
    Returns: dictionary containing information about correction, containing:
        ???
    """

    #If necessary, convert telescope name into telescope location object:    
    if type(telescope_location) == str:
        telescope_location=get_telescope_coordinates(telescope_location)

    #RMextract wants time ranges in MJD seconds:
    timerange=[Time(start_time).mjd*86400.0,
               Time(end_time).mjd*86400.0] 

    #Extract telescope coordinates into expected format (geodetic x,y,z):
    telescope_coordinates=[telescope_location.x.value,
                           telescope_location.y.value,
                           telescope_location.z.value]
    
    #Get RMExract to generate it's RM predictions
    predictions=RME.getRM(ionexPath=ionexPath, radec=np.deg2rad([ra,dec]), 
                          timestep=timestep, 
                          timerange = timerange, 
                          stat_positions=[telescope_coordinates,])
    #predictions dictionary contains STEC, Bpar, BField, AirMass, elev, azimuth
    #  RM, times, timestep, station_names, stat_pos, flags, reference_time
    times=predictions['times']
    RMs=np.squeeze(predictions['RM']['st1'])


def main():
    """When invoked from the command line, generate ionospheric RM predictions
    and save the predictions, as a function of frequency, to a text file.
    """
    import argparse
    
    



if __name__ == "__main__":
    main()









