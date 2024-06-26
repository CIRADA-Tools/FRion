#!/usr/bin/env python3

"""
Functions for predicting ionospheric Faraday rotation effects, both
time-dependent and time-averaged.
Uses RMextract package to get time-dependent ionospheric RMs for a given
observation, then performs the time-integration to work out the effective
change in polarization angle, and the effective depolarization (together,
called the ionospheric *modulation*\ , which is called Theta in the derivation).

The ionospheric prediction is currently derived from RMExtract (https://github.com/lofar-astron/RMextract/).
Other ionosphere RM codes are available (ionFR, ALBUS) are available, but
RMextract was selected for its comparative ease of install and use.

RMextract relies on external maps of Total Electron Content (TEC). As of 
version 1.1, the default is to get the TEC data from CDDIS
(https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html).
This requires an account and a .netrc file with the credentials
(https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html).
Advanced users can change the TEC source using the kwargs for RMextract's
getRM function, with the caveat that this is not supported for the command line
tools.
The default TEC maps are now the JPL global maps, based on the results of
Porayko et al. 2019 who found they performed the best.


Two command line functions are available and documented below:
    - predict(): time-averaged Faraday rotation.
    - timeseries(): time-dependent Faraday rotation.


"""


try:
    import RMextract.getRM as RME
    from RMextract import getIONEX as ionex
except:
    #This is the easiest solution to the documentation problem:
    #This code needs to be importable when RMextract isn't installed, for
    #ReadTheDocs to work. Getting RMextract to install properly in RTD is too
    #much work, so this is my workaround.
    print('Cannot import RMextract. Continuing import, but will fail if called.')
from astropy.time import Time,TimeDelta
import numpy as np
from scipy.integrate import simpson as simps
from astropy.coordinates import EarthLocation,SkyCoord,Angle, UnknownSiteException
import astropy.units as u
from FRion.correct import find_freq_axis
import astropy
from FRion.download_IONEX_CDDIS import get_CDDIS_IONEXfile

C = 2.99792458e8 # Speed of light [m/s]


def predict():
    """Wrapper for command line interface, for time-averaging mode.
    Gets command line arguments, calculates ionospheric effects, and saves
    modulation and/or figure if specified.

    If the package is installed, it can be called as frion_predict,
    otherwise as python predict.py

    Call with the -h flag to see command line options.
    """
    import argparse
    descStr = """
    Calculate ionospheric Faraday rotation and predict time-integrated effect
    as a function of frequency.
    Can determine the frequency, location, direction, and observation time
    parameters from a supplied FITS cube or PSRFITS file, if it has the correct
    keywords, otherwise from those parameters must be supplied on the command
    line.
    
    By default, gets ionosphere TEC data from CDDIS (https://cddis.nasa.gov/Data_and_Derived_Products/GNSS/atmospheric_products.html)
    This requires an account and a .netrc file with the credentials to run.
    (https://cddis.nasa.gov/Data_and_Derived_Products/CreateNetrcFile.html)
    """

    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-F",dest="fits",default=None,metavar='FILENAME',
                        help="FITS cube relevant information in header.")
    parser.add_argument("-d", dest="times", nargs=2,type=str,default=None,
                        metavar=('START','END'),
                        help="start and end time strings.")
    parser.add_argument("-c", dest=("freq_parms"),nargs=3,default=None,type=float,
                        metavar=('MINFREQ','MAXFREQ','CHANNELWIDTH'),
                        help=("Generate channel frequencies (in Hz): \n"
                              "    minfreq, maxfreq, channel_width"))
    parser.add_argument("-t",dest='telescope_name',type=str,default=None,
                        help="Telescope name")
    parser.add_argument("-T",dest='telescope_coords',nargs=3,type=float,default=None,
                        metavar=('LONG','LAT','ALT'),
                        help="Telescope coordinates:\n    lat[deg],long[deg], altitude[m].")
    parser.add_argument('-p',dest='pointing',nargs=2,type=float,default=None,
                        metavar=('RA','DEC'),
                        help="Pointing center: RA[deg], Dec[deg]")
    parser.add_argument("-s", dest='savefile',type=str,default=None,metavar='POLFILE',
                        help="Filename to save ionosphere data to.")
    parser.add_argument("-S", dest='savefig',type=str,default=None,metavar='FIGFILE',
                        help="Filename to save the plots to. Entering 'screen' plots to the screen.")
    parser.add_argument("--timestep",dest='timestep',default=600.,type=float,
                        help="Timestep for ionospheric prediction, in seconds. Default = 600")
    parser.add_argument("-P", dest="tec_prefix",type=str,default='uqrg',metavar='TEC_SOURCE')
    args = parser.parse_args()


    start_time=None
    end_time=None
    freq_arr=None
    telescope=None
    ra=None
    dec=None


    #If a FITS file is present, try to fill in any missing keywords.
    #But since FITS headers can be very different, it may be that not all
    #keywords can be found.
    if args.fits is not None:
        start_time,end_time,freq_arr,telescope,ra,dec=get_parms_from_FITS(args.fits)

    #Any parametrs taken from FITS header can be overridden by manual inputs:
    if args.times is not None:
        start_time=Time(args.times[0])
        end_time=Time(args.times[1])

    if args.freq_parms is not None:
        freq_arr=np.arange(args.freq_parms[0],args.freq_parms[1],args.freq_parms[2])

    if args.telescope_name is not None:
        telescope = get_telescope_coordinates(args.telescope_name)
    if args.telescope_coords is not None:
        telescope = get_telescope_coordinates(args.telescope_coords)

    if args.pointing is not None:
        ra=args.pointing[0]
        dec=args.pointing[1]

    #Check that all parameters are set:
    missing_parms=[]
    if (start_time is None): missing_parms.append('Start time')
    if (end_time is None): missing_parms.append('End time')
    if (freq_arr is None): missing_parms.append('Frequency array')
    if (telescope is None): missing_parms.append('Telescope')
    if (ra is None): missing_parms.append('Pointing center')

    if len(missing_parms) > 0:
        print("\n\nMissing parameters:",missing_parms,'\n')
        raise Exception("Cannot continue without parameters listed above.")

    times,RMs,theta=calculate_modulation(start_time, end_time, freq_arr, telescope,
                         ra,dec, timestep=args.timestep,
                         ionexPath='./IONEXdata/', prefix=args.tec_prefix)

    if args.savefile is not None:
        write_modulation(freq_arr,theta,args.savefile)

    if args.savefig is not None:
        generate_plots(times,RMs,theta,freq_arr,position=[ra,dec],savename=args.savefig)


def calculate_modulation(start_time, end_time, freq_array, telescope_location,
                         ra,dec, timestep=600.,ionexPath='./IONEXdata/',
                         prefix='uqrg', **kwargs):
    """Calculate the ionospheric FR modulation (time-averaged effect),
    as a function of frequency, for a given observation (time, location, target direction).



    Args:
        start_time (astropy.time Time, or string readable by astropy.time Time):
            Starting time of observation.
            Example time string: '2010-01-02T00:00:00'
        end_time (astropy.time.Time, or string readable by astropy.time.Time):
            ending time of observation.
        freq_array (array-like, either floats or Astropy Quantity):
            vector of channel frequencies (in Hz)
        telescope_location (astropy.coordinates EarthLocation or string):
            location of telescope, or name of telescope known to
            get_telescope_coordinates() function.
        ra (astropy.coordinates Angle, astropy.units Quantity, or float):
            right ascension of observation center (if float: in deg, J2000)
        dec (astropy.coordinates Angle, astropy.units Quantity, or float):
            declination of observation center (if float: in deg, J2000)
        timestep (astropy.time TimeDelta, astropy.units Quantity, or float):
            time between ionosphere FR estimates. If float, time must be in seconds.
        ionexPath (str, default='./IONEXdata/'): path to download IONEX files to
            for ionosphere calculations.
        **kwargs: additional keyword arguments to pass to RMextract.getRM()
        
    Returns:
        tuple containing

        -times (astropy Time array): vector of times of each ionospheric RM calculation

        -RMs (array): vector of RM values computed for each time step

        -theta (array) vector containing the (complex) ionospheric polarization
        for each frequency channel.

    """

    #Convert frequencies to have units if needed.
    if type(freq_array) == astropy.units.quantity.Quantity:
        frequencies=freq_array
    else:
        frequencies=freq_array*u.Hz

    #Calculation of the time-dependent RMs.
    times,RMs=get_RM(start_time, end_time, telescope_location,
                         ra,dec, timestep=timestep,ionexPath=ionexPath,
                         prefix=prefix,**kwargs)


    #Compute the time-integrated change in polarization.
    theta=numeric_integration(times.mjd*86400.,RMs,frequencies.to(u.Hz).value)

    #Verify that we are not in a regime where numerical instabilities might occur.
    check_numeric_problems(RMs, frequencies.to(u.Hz).value,theta)


    return times,RMs, theta


def get_RM(
    start_time, 
    end_time, 
    telescope_location,
    ra,
    dec, 
    timestep=600.,
    ionexPath='./IONEXdata/',
    pre_download=True,
    prefix='uqrg',
    server='http://cddis.gsfc.nasa.gov',
    **kwargs
):
    """
    Calculate the ionospheric Faraday rotation as a function of time,
    for a given observation (time, location, target direction).

    Args:
        start_time (astropy.time Time, or string readable by astropy.time Time):
            Starting time of observation.
            Example time string: '2010-01-02T00:00:00'
        end_time (astropy.time.Time, or string readable by astropy.time.Time):
            ending time of observation.
        telescope_location (astropy.coordinates EarthLocation or string):
            location of telescope, or name of telescope known to
            get_telescope_coordinates() function.
        ra (astropy.coordinates Angle, astropy.units Quantity, or float):
            right ascension of observation center (if float: in deg, J2000)
        dec (astropy.coordinates Angle, astropy.units Quantity, or float):
            declination of observation center (if float: in deg, J2000)
        timestep (astropy.time TimeDelta, astropy.units Quantity, or float):
            time between ionosphere FR estimates. If float, time must be in seconds.
        ionexPath (str, default='./IONEXdata/'): path to download IONEX files to
            for ionosphere calculations.
        pre_download (bool, default=True): if True, will pre-download the IONEX
            files from CDDIS before running RMextract.
        prefix (str, default='uqrg'): prefix for IONEX files to download.
        server (str, default='http://cddis.gsfc.nasa.gov'): server to download
            IONEX files from.
        **kwargs: additional keyword arguments to pass to RMextract.getRM()
        
    Returns:
        tuple containing

        -times (astropy Time array): vector of times of each ionospheric RM calculation

        -RMs (array): vector of RM values computed for each time step

    """
    #If necessary, convert telescope name into telescope location object:
    if type(telescope_location) != EarthLocation:
        telescope_location=get_telescope_coordinates(telescope_location)

    #RMextract wants time ranges in MJD seconds:
    timerange=[Time(start_time).mjd*86400.0,
               Time(end_time).mjd*86400.0]

    #Extract telescope coordinates into expected format (geodetic x,y,z):
    telescope_coordinates=[telescope_location.x.value,
                           telescope_location.y.value,
                           telescope_location.z.value]

    #Handle all forms of angle input:
    if type(ra) == float or type(ra) == int:
        ra_angle=Angle(ra,'deg')
    elif type(ra) == astropy.units.quantity.Quantity:
        ra_angle=Angle(ra)
    elif type(ra) == astropy.coordinates.angles.Angle:
        ra_angle=ra
    else:
        raise Exception("""RA input object type not recognized.
                        Only astropy.coordinates Angle, astropy.units Quantity, or float (in deg) allowed.""")

    if type(dec) == float or type(dec) == int:
        dec_angle=Angle(dec,'deg')
    elif type(dec) == astropy.units.quantity.Quantity:
        dec_angle=Angle(dec)
    elif type(dec) == astropy.coordinates.angles.Angle:
        dec_angle=dec
    else:
        raise Exception("""Dec input object type not recognized.
                        Only astropy.coordinates Angle, astropy.units Quantity, or float (in deg) allowed.""")

    #Handle all forms of time input:
    if type(timestep) == float or type(timestep) == int:
        timestep_Delta=TimeDelta(timestep*u.second)
    elif type(timestep) == astropy.units.quantity.Quantity:
        timestep_Delta=TimeDelta(timestep)
    elif type(timestep) == astropy.time.core.TimeDelta:
        timestep_Delta=timestep
    else:
        raise Exception("""Timestep input object type not recognized.
                        Only astropy.time TimeDelta, astropy.units Quantity, or float (in seconds) allowed.""")


    #Pre-download the IONEX data from CDDIS, to work around the RMextract lack
    #of support for CDDIS downloads.    
    if pre_download:
        _predownload_CDDIS(start_time, end_time,prefix,outpath=ionexPath)

    #Get RMExtract to generate its RM predictions
    predictions=RME.getRM(ionexPath=ionexPath,
                          radec=[ra_angle.rad,dec_angle.rad],
                          timestep=timestep_Delta.sec,
                          timerange = timerange,
                          stat_positions=[telescope_coordinates,],
                          prefix=prefix,
                          server=server,
                          **kwargs,
                        )
    
    #predictions dictionary contains STEC, Bpar, BField, AirMass, elev, azimuth
    #  RM, times, timestep, station_names, stat_pos, flags, reference_time
    times=Time(predictions['times']/86400.,format='mjd')
    RMs=np.squeeze(predictions['RM']['st1'])


    return times, RMs


def _predownload_CDDIS(start_time,end_time,prefix='uqrg',outpath='./IONEXdata/'):
    """Downloads the IONEX maps for the required dates from CDDIS.
    This must occur before RMextract is invoked, otherwise RMextract will try
    to use its own download tool which is broken.
    
    Will try to download all the days between the start and end dates.
    Is slightly greedy (downloading more days than) may be necessary) as a safety
    against edge cases.
    
    """
    from math import ceil, floor
    
    
    #Work out how many days need to be downloaded:
    start_date=Time(start_time)
    end_date=Time(end_time)
    
    #Download each day one by one:
    for day_mjd in range(floor(start_date.mjd)-1,ceil(end_date.mjd)+1):
        day = Time(day_mjd,format='mjd')
        fname=get_CDDIS_IONEXfile(time=day.to_value('isot'),
                            prefix=prefix,
                            outpath=outpath,
                            overwrite=False)

    if fname == -1:
        raise Exception('Something has gone wrong in downloading IONEX data.')



def numeric_integration(times,RMs,freq_array):
    """Numerical integration of the time-varying ionospheric polariation
    modulation. Testing has shown that numerical integration is accurate
    to better than 1% accuracy except where depolarization is extreme (>99%).

    Args:
        times (array): ionosphere sampling times (in MJD seconds, as used in RMextract)
        RMs (array): ionospheric RMs at each time (in rad/m^2)
        freq_array (array): channel frequencies (in Hz)

    Returns: array: time-integrated ionospheric modulation per channel.
    """

    l2_arr=(C/freq_array)**2
    z=np.exp(2.j*np.outer(l2_arr,RMs))

    #Scipy's numerical integrators can't handle complex numbers, so the
    # integral needs to be broken into real and complex components.
    real=simps(z.real,times,axis=1)
    imag=simps(z.imag,times,axis=1)
    theta=(real+1.j*imag)/(times[-1]-times[0])
    return theta



def check_numeric_problems(RMs, freq_array,theta):
    """Checks for conditions that might cause numeric instability in the
    time-integration, and warns the user if there might be concerns.

    Specifically, checks for extreme jumps in polarization angle between
    timesteps (will cause integration errorrs),
    and for extreme depolarization (high liklihood of large errors).

    Args:
        RMs (array): ionospheric RMs per time step
        freq_array (array): channel frequencies (in Hz)
        theta (array): ionospheric modulation per channel.
    """
    import warnings

    #Check for large jumps in RM/polarization angle between steps.
    #These can cause the numeric integrator to not catch angle wraps.
    longest_l2=(C/np.min(freq_array))**2
    max_deltaRM=np.max(np.diff(RMs))
    max_delta_polangle=longest_l2*max_deltaRM  #in radians
    if max_delta_polangle > 0.5:
        warnings.warn(("\nLarge variations in RM between points, which may "
                       "introduce numerical errors.\n"
                       "Consider trying a smaller timestep."))

    #Warn about very low values of theta (very strong depolarization)
    # as these can probably not be corrected reliably.
    if np.min(np.abs(theta)) < 0.02:
        warnings.warn(("\nExtreme depolarization predicted (>98%). "
                       "Corrected polarization will almost certainly not "
                       "be trustworthy in affected channels."))
    elif np.min(np.abs(theta)) < 0.1:
        warnings.warn(("\nSignificant depolarization predicted (>90%). "
                       "Errors in corrected polarization are likely to be "
                       "very large in some channels."))



def write_modulation(freq_array,theta,filename):
    """Saves predicted ionospheric modulation to a text file.
    File has two columns, whitespace-delimited.

    Args:
        freq_array (array): channel frequencies (in Hz)
        theta (array): ionospheric (complex) modulation at each frequency
        filename (str): file path to save data to.
"""
    np.savetxt(filename, list(zip(freq_array,theta.real,theta.imag)))




def get_telescope_coordinates(telescope):
    """Return the astropy.coordinates EarthLocation object associated
    with the position of the telescope.
    The input must be either a string with the name of a telescope listed in
    Astropy's sites list
    (https://github.com/astropy/astropy-data/blob/gh-pages/coordinates/sites.json),
    an EarthLocation object (which is passed through unchanged),
    or a 3-component tuple with the latitude [deg], longitude [deg],
    and height [m].

    Since ASKAP is not yet listed in the Astropy sites list, it's position is
    manually coded in as a temporary measure.

    """
    if type(telescope) == EarthLocation: #Pass EarthLocations through without processing
        return telescope
    elif type(telescope) == str: #Hardcoded coordinates for some telescopes.
        if telescope == 'ASKAP':
            lat = -1*26+42/60+15/3600 #degree
            long = +1*116+39/60+32/3600 # degree
            height = 381.0 #
        else:
            return EarthLocation.of_site(telescope)
        return EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=height*u.m)
    elif (type(telescope) == tuple) or (type(telescope) == list):
        return EarthLocation(lat=telescope[0]*u.deg, lon=telescope[1]*u.deg,
                             height=telescope[2]*u.m)






def generate_plots(times,RMs,theta,freq_array,position=None,savename=None):
    """Makes a figure with two plots: the RM variation over time,
    and the (modulus of the) modulation as a function of frequency.
    If savename contains a string it will save the plots to that filename,
    otherwise the plots are not saved.
    If position ([ra,dec]) is given, it will be printed above the plots.

    Args:
        times (astropy Time array): ionosphere sampling times
        RMs (array): ionospheric RMs at each time (in rad/m^2)
        theta (array): ionospheric (complex) modulation at each frequency
        freq_array (array): channel frequencies (in Hz)
        position (list): [ra, dec] in degrees, left blank if not supplied.
        savename (str): File path to save plot to; if 'screen' will send to display.

    """
    from matplotlib import pyplot as plt
    from matplotlib import dates as mdates
    plot_times=times.plot_date


    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(8,8))
    ax1.plot_date(plot_times,RMs,fmt='k.')
    locator=mdates.AutoDateLocator(minticks=3,maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax1.xaxis.set_major_locator(locator)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.set_ylabel(r'$\phi_\mathrm{ion}$ [rad m$^{-2}$]')

    ax2.plot(freq_array,np.abs(theta),'k.')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel(r'|$\Theta(\lambda^2)$|')
    if position is not None:
        ax1.set_title("RA: {:.2f}°, Dec: {:.2f}°".format(position[0],position[1]))

    if savename is not None:
        if savename == "screen":
            plt.show()
        else:
            plt.savefig(savename,bbox_inches='tight')





def get_parms_from_FITS(filename):
    """Extract relevant parameters (start time, end time, telescope, ra, dec,
    and frequencies) from a FITS file, possibly a PSRFITS file.
    If the missing keywords are not found in the FITS file, they are left blank
    and the user

    Args:
        filename (str): path to the FITS file.

    Returns:
        start_time
        end_time
        freq_arr
        telescope
        ra
        dec
    """
    import astropy.io.fits as pf
    hdulist=pf.open(filename)
    header=hdulist[0].header


    start_time=None
    end_time=None
    freq_arr=None
    telescope=None
    ra=None
    dec=None

    #PSRFITS files have a different format, so split prpcesing depending on
    #PSRFITS vs FITS image:
    if 'FITSTYPE' in header.keys() and header['FITSTYPE']=='PSRFITS':
        #PRSFITS code adapted from example by David Kaplan.
        #I'm trying to avoid using a PSRFITS reader module, to minimize dependencies,
        #at risk of not handling all variations of PSRFITS properly.

        if 'STT_IMJD' in header.keys():
            start_time=Time(float(header['STT_IMJD']) +
                        (float(header['STT_SMJD'])+
                         float(header['STT_OFFS']))/float(86400),
                        format='mjd')

            try:
                end_time = start_time + np.sum(hdulist['SUBINT'].data['TSUBINT']) * u.s
            except:
                pass

        if 'TELESCOP' in header.keys():
            try:
                telescope = EarthLocation.of_site(header['TELESCOP'])
            except UnknownSiteException:
                telescope = EarthLocation.from_geocentric(header['ANT_X'],
                                                     header['ANT_Y'],
                                                     header['ANT_Z'],
                                                     unit=u.m)

        if 'RA' in header.keys():
            ra=Angle(hdulist[0].header['RA'],unit='hour')
        if 'DEC' in header.keys():
            dec=Angle(hdulist[0].header['DEC'],unit='deg')

        if 'OBSFREQ' in header.keys() and 'OBSNCHAN' in header.keys() and 'OBSBW' in header.keys():
            chan_bw=header['OBSBW']/header['OBSNCHAN']

            startchan= header['OBSFREQ'] - (header['OBSNCHAN']/2 - 1)*chan_bw
            freq_arr=np.arange(header['OBSNCHAN'])*chan_bw+startchan


    #FITS images/cubes:
    else:
        if 'DATE-OBS' in header.keys():
            start_time=Time(header['DATE-OBS'])
        if 'DATE-OBS' in header.keys() and 'DURATION' in header.keys():
            end_time=Time(header['DATE-OBS'])+TimeDelta(header['DURATION'],format="sec")

        #For coordinates, the code will always use the middle pixel to derive
        #the RA and Dec. Assumes position coordinates are in first 2 axes.
        if 'RA' in header['CTYPE1']:
            ra=header['CRVAL1']+header['CDELT1']*(header['NAXIS1']/2-header['CRPIX1'])
        if 'DEC' in header['CTYPE2']:
            dec=header['CRVAL2']+header['CDELT2']*(header['NAXIS2']/2-header['CRPIX2'])
        if 'GLON' in header['CTYPE1'] and 'GLAT' in header['CTYPE2']:
        #Support galactic coordinates, just in case:
            gl=header['CRVAL1']+header['CDELT1']*(header['NAXIS1']/2-header['CRPIX1'])
            gb=header['CRVAL2']+header['CDELT2']*(header['NAXIS2']/2-header['CRPIX2'])
            position=SkyCoord(gl,gb,frame='galactic',unit='deg')
            ra=position.fk5.ra.deg
            dec=position.fk5.dec.deg

        if 'TELESCOP' in header.keys():
            telescope=header['TELESCOP']

        freq_axis=str(find_freq_axis(header))
        if freq_axis != '0':
            chan0=header['CRVAL'+freq_axis]-header['CDELT'+freq_axis]*(header['CRPIX'+freq_axis]-1)
            chan_final=chan0+(header['NAXIS'+freq_axis]-1)*header['CDELT'+freq_axis]
            freq_arr=np.linspace(chan0,chan_final,header['NAXIS'+freq_axis])

    return start_time,end_time,freq_arr,telescope,ra,dec







def timeseries():
    """
    Wrapper for command line interface, for time-series mode.
    Gets command line arguments, feeds them into RMextract, and saves the
    output to the file if specified or prints to terminal.

    If the package is installed, it can be called as frion_timeseries,
    otherwise as python -c 'import FRion.predict; predict.timeseries()'

    Call with the -h flag to see command line options.
    """
    import argparse
    descStr = """
    Calculate ionospheric Faraday rotation as a function of time.
    Can determine the observation time, direction, and location parameters
    from a supplied FITS cube or PSRFITS file, if it has the correct keywords,
    otherwise from those parameters must be supplied on the command line.
    """

    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("-F",dest="fits",default=None,metavar='FILENAME',
                        help="FITS cube relevant information in header.")
    parser.add_argument("-d", dest="times", nargs=2,type=str,default=None,
                        metavar=('START','END'),
                        help="start and end time strings.")
    parser.add_argument("-t",dest='telescope_name',type=str,default=None,
                        help="Telescope name")
    parser.add_argument("-T",dest='telescope_coords',nargs=3,type=float,default=None,
                        metavar=('LONG','LAT','ALT'),
                        help="Telescope coordinates:\n    lat[deg],long[deg], altitude[m].")
    parser.add_argument('-p',dest='pointing',nargs=2,type=float,default=None,
                        metavar=('RA','DEC'),
                        help="Pointing center: RA[deg], Dec[deg]")
    parser.add_argument("-s", dest='savefile',type=str,default=None,metavar='POLFILE',
                        help="Filename to save ionosphere data to.")
    parser.add_argument("-f", dest='timeformat',type=str,default='mjd',
                        metavar='TIMEFORMAT',
                        help="Format for times, must be one from list at https://docs.astropy.org/en/stable/time/index.html#time-format \n Default is mjd, for human-readable try fits.")
    parser.add_argument("-S", dest='savefig',type=str,default=None,metavar='FIGFILE',
                        help="Filename to save the plots to. Entering 'screen' plots to the screen.")
    parser.add_argument("--timestep",dest='timestep',default=600.,type=float,
                        help="Timestep for ionospheric prediction, in seconds. Default = 600")
    args = parser.parse_args()


    start_time=None
    end_time=None
    freq_arr=None
    telescope=None
    ra=None
    dec=None


    #If a FITS file is present, try to fill in any missing keywords.
    #But since FITS headers can be very different, it may be that not all
    #keywords can be found.
    if args.fits is not None:
        start_time,end_time,freq_arr,telescope,ra,dec=get_parms_from_FITS(args.fits)

    #Any parametrs taken from FITS header can be overridden by manual inputs:
    if args.times is not None:
        start_time=Time(args.times[0])
        end_time=Time(args.times[1])


    if args.telescope_name is not None:
        telescope = get_telescope_coordinates(args.telescope_name)
    if args.telescope_coords is not None:
        telescope = get_telescope_coordinates(args.telescope_coords)

    if args.pointing is not None:
        ra=args.pointing[0]
        dec=args.pointing[1]

    #Check that all parameters are set:
    missing_parms=[]
    if (start_time is None): missing_parms.append('Start time')
    if (end_time is None): missing_parms.append('End time')
    if (telescope is None): missing_parms.append('Telescope')
    if (ra is None): missing_parms.append('Pointing center')

    if len(missing_parms) > 0:
        print("\n\nMissing parameters:",missing_parms,'\n')
        raise Exception("Cannot continue without parameters listed above.")

    times,RMs=get_RM(start_time, end_time, telescope,
                         ra,dec, timestep=args.timestep,ionexPath='./IONEXdata/')

    if args.savefile is not None:
        write_timeseries(times,RMs,args.savefile,timeformat=args.timeformat)
    else:
        print('Times:      RM:')
        for tm, rm in zip(times.to_value(args.timeformat),RMs):
            print(tm,rm)


    if args.savefig is not None:
        generate_plots(times,RMs,position=[ra,dec],savename=args.savefig)






def write_timeseries(times,RMs,filename,timeformat='mjd'):
    """Saves the predicted ionospheric RMs as a function of time to a text file.
    File has two columns, whitespace-delimited (Be aware that some time formats
    will add whitespace to the time column.)

    Args:
        times (astropy Time array): ionosphere sampling times
        RMs (array): ionospheric RMs at each time (in rad/m^2)
        filename (str): file path to save data to.
        timeformat (str): preferred format for writing times.
            Possible values can be found at
            https://docs.astropy.org/en/stable/time/index.html#time-format
            default is 'mjd' (e.g., ‘51544.0’).

    """
    fout = open(filename, "w")
    for tm, rm in zip(times.to_value(timeformat),RMs):
            print(tm,rm, file=fout)





def plot_timeseries(times,RMs,position=None,savename=None):
    """Makes a figure plotting the RM variation over time,
    If savename contains a string it will save the plots to that filename
    (unless the name is 'screen', in which case it will plot on screen),
    otherwise the plots are not saved.
    If position ([ra,dec]) is given, it will be printed above the plots.

    Args:
        times (astropy Time array): ionosphere sampling times
        RMs (array): ionospheric RMs at each time (in rad/m^2)
        position (list): [ra, dec] in degrees, left blank if not supplied.
        savename (str): File path to save plot to; if 'screen' will send to display.

    """
    from matplotlib import pyplot as plt
    from matplotlib import dates as mdates
    plot_times=times.plot_date


    fig,ax1=plt.subplots(1,1,figsize=(8,8))
    ax1.plot_date(plot_times,RMs,fmt='k.')
    locator=mdates.AutoDateLocator(minticks=3,maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax1.xaxis.set_major_locator(locator)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.set_ylabel(r'$\phi_\mathrm{ion}$ [rad m$^{-2}$]')

    if position is not None:
        ax1.set_title("RA: {:.2f}°, Dec: {:.2f}°".format(position[0],position[1]))

    if savename is not None:
        if savename == "screen":
            plt.show()
        else:
            plt.savefig(savename,bbox_inches='tight')



if __name__ == "__main__":
    predict()









