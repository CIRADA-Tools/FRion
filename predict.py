#!/usr/bin/env python3

"""
Tools for predicting time-integrated ionospheric Faraday rotation corrections.
Uses RMextract package to get time-dependent corrections for a given 
observation, then performs the time-integration to work out the effective
change in polarization angle, and the effective depolarization.


"""



import RMextract.getRM as RME
from astropy.time import Time,TimeDelta
import numpy as np
from astropy.coordinates import EarthLocation,SkyCoord
import astropy.units as u
from apply import find_freq_axis

C = 2.99792458e8 # Speed of light [m/s]



def get_telescope_coordinates(telescope):
    """Return the astropy.coordinates EarthLocation object associated
    with the position of the telescope.
    The input must be either a string with the name of a pre-programmed telescope,
     a 3-component tuple with the latitude [deg], longitude [deg], and height [m],
     or an EarthLocation object.
    """
    if type(telescope) == EarthLocation: #Pass EarthLocations through without processing
        return telescope
    elif type(telescope) == str: #Hardcoded coordinates for some telescopes.
        if telescope == 'ASKAP':
            lat = -1*26+42/60+15/3600 #degree
            long = +1*116+39/60+32/3600 # degree
            height = 381.0 #
        else:
            raise Exception("Telescope name not recognized.")
        return EarthLocation(lat=lat*u.deg, lon=long*u.deg, height=height*u.m)
    elif (type(telescope) == tuple) or (type(telescope) == list):
        return EarthLocation(lat=telescope[0]*u.deg, lon=telescope[1]*u.deg, 
                             height=telescope[2]*u.m)






def write_correction(freq_array,corr,filename):
    np.savetxt(filename, list(zip(freq_array,corr.real,corr.imag)))


def numeric_integration(times,RMs,freq_array):
    """Numerical integration of the time-varying ionospheric polariation 
    modulation. Testing has shown that numerical integration is accurate 
    to better than 1% accuracy except where depolarization is extreme (>99%).
    """

    from scipy.integrate import simps
    l2_arr=(C/freq_array)**2
    z=np.exp(2.j*np.outer(l2_arr,RMs))
    
    #Scipy's numerical integrators can't handle complex numbers, so the
    # integral needs to be broken into real and complex components.
    real=simps(z.real,times,axis=1)
    imag=simps(z.imag,times,axis=1)
    corr=(real+1.j*imag)/(times[-1]-times[0])
    return corr
    

def check_numeric_problems(RMs, freq_array,corr):
    """Checks for conditions that might cause numeric instability in the
    correction, and warns the user if there might be concerns.
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
    
    #Warn about very low values of the correction (very strong depolarization)
    # as these can probably not be corrected reliably.
    if np.min(np.abs(corr)) < 0.02:
        warnings.warn(("\nExtreme depolarization predicted (>98%). "
                       "Corrected polarization will almost certainly not "
                       "be trustworthy in affected channels."))
    elif np.min(np.abs(corr)) < 0.1:
        warnings.warn(("\nSignificant depolarization predicted (>90%). "
                       "Errors in corrected polarization are likely to be "
                       "very large in some channels."))
    
    


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
    Returns: times (array): vector of times (in MJD seconds) of each 
                            ionospheric RM calculation
             RMs (array): vector of RM values computed for each time step
             correction (array): vector containing the (complex) ionospheric
                polarization (Theta) for each frequency channel.
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
    
    #Get RMExract to generate it's RM predictions
    predictions=RME.getRM(ionexPath=ionexPath, radec=np.deg2rad([ra,dec]), 
                          timestep=timestep, 
                          timerange = timerange, 
                          stat_positions=[telescope_coordinates,])
    #predictions dictionary contains STEC, Bpar, BField, AirMass, elev, azimuth
    #  RM, times, timestep, station_names, stat_pos, flags, reference_time
    times=predictions['times']
    RMs=np.squeeze(predictions['RM']['st1'])
    
    #Compute the time-integrated correction. 
    correction=numeric_integration(times,RMs,freq_array)
    
    return times,RMs, correction
    
    

def generate_plots(times,RMs,corr,freq_array,position=None,savename=None):
    """Makes a figure with two plots: the RM variation over time, 
    and the (modulus of the) correction as a function of frequency.
    If savename contains a string it will save the plots to that filename,
    otherwise the plots are not saved.
    If position ([ra,dec]) is given, it will be printed above the plots.
    """
    from matplotlib import pyplot as plt
    from matplotlib import dates as mdates
    plot_times=Time(times/86400.0,format='mjd').plot_date
    
    
    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(8,8))
    ax1.plot_date(plot_times,RMs,fmt='k.')
    locator=mdates.AutoDateLocator(minticks=3,maxticks=7)
    formatter = mdates.ConciseDateFormatter(locator)
    ax1.xaxis.set_major_locator(locator)
    ax1.xaxis.set_major_formatter(formatter)
    ax1.set_ylabel('$\phi_\mathrm{ion}$ [rad m$^{-2}$]')
    
    ax2.plot(freq_array,np.abs(corr),'k.')
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('|$\Theta(\lambda^2)$|')
    if position is not None:
        ax1.set_title("RA: {:.2f}°, Dec: {:.2f}°".format(position[0],position[1]))
    
    if savename is not None:
        if savename == "screen":
            plt.show()    
        else:
            plt.savefig(savename,bbox_inches='tight')



def main():
    """When invoked from the command line, generate ionospheric RM predictions
    and save the predictions, as a function of frequency, to a text file.
    """
    import argparse
    descStr = """
    Calculate ionospheric Faraday rotation and predict time-integrated effect
    as a function of frequency.
    Can determine the frequency and observation time parameters from a 
    supplied FITS cube, if it has the correct keywords, 
    otherwise from those parameters must be supplied on the command line.
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
    parser.add_argument("-s", dest='savefile',type=str,default=None,metavar='CORRFILE',
                        help="Filename to save correction data to.")
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
        import astropy.io.fits as pf
        header=pf.getheader(args.fits)
        if 'DATE-OBS' in header.keys():
            start_time=header['DATE-OBS']
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
        print("Missing parameters:",missing_parms)
        raise Exception("Cannot continue without those parameters.")
    
    times,RMs,correction=calculate_correction(start_time, end_time, freq_arr, telescope,
                         ra,dec, timestep=args.timestep,ionexPath='./IONEXdata/')
    
    if args.savefile is not None:
        write_correction(freq_arr,correction,args.savefile)

    if args.savefig is not None:
        generate_plots(times,RMs,correction,freq_arr,position=[ra,dec],savename=args.savefig)
    
    
    

if __name__ == "__main__":
    main()









