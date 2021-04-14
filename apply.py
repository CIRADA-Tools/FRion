#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Tools for applying ionospheric Faraday rotation corrections to images and
image cubes. This file contains the functions needed to apply the correction
to data cubes (either in memory, or in FITS files.)
"""

import numpy as np
from astropy.io import fits as pf






def read_correction(filename):
    (frequencies,real,imag)=np.genfromtxt(filename,unpack=True)
    correction=real+1.j*imag
    return frequencies, correction


def find_freq_axis(header):
    """Finds the frequency axis in a FITS header.
    Input: header: a Pyfits header object.
    Returns the axis number (as recorded in the FITS file, **NOT** in numpy ordering.)
    Returns 0 if the frequency axis cannot be found.
    """
    freq_axis=0 #Default for 'frequency axis not identified'
    #Check for frequency axes. Because I don't know what different formatting
    #I might get ('FREQ' vs 'OBSFREQ' vs 'Freq' vs 'Frequency'), convert to
    #all caps and check for 'FREQ' anywhere in the axis name.
    for i in range(1,header['NAXIS']+1):  #Check each axis in turn.
        try:
            if 'FREQ' in header['CTYPE'+str(i)].upper():
                freq_axis=i
        except:
            pass #The try statement is needed for if the FITS header does not
                 # have CTYPE keywords.
    return freq_axis


def readData(Qfilename,Ufilename):
    """
    Open the Stokes Q and U input cubes (from the supplied 
    file names) and return data-access variables and the header. 
    Axes are re-ordered so that frequency is first, beyond that the number
    and ordering of axes doesn't matter.
    Uses the memmap functionality so that data isn't read into data; variables
    are just handles to access the data on disk.
    Returns the header from the Q file, the U file's header is ignored.
    """    
    
    hdulistQ=pf.open(Qfilename,memmap=True)
    header=hdulistQ[0].header
    Qdata=hdulistQ[0].data
    hdulistU=pf.open(Ufilename,memmap=True)
    Udata=hdulistU[0].data
    
    
    N_dim=header['NAXIS'] #Get number of axes

    freq_axis=find_freq_axis(header) 
    #If the frequency axis isn't the last one, rotate the array until it is.
    #Recall that pyfits reverses the axis ordering, so we want frequency on
    #axis 0 of the numpy array.
    if freq_axis != 0 and freq_axis != N_dim:
        Qdata=np.moveaxis(Qdata,N_dim-freq_axis,0)
        Udata=np.moveaxis(Udata,N_dim-freq_axis,0)

    
    return Qdata, Udata, header


def write_corrected_cubes(Qoutputname,Uoutputname,Qcorr,Ucorr,header,overwrite=False):
    """
    Write the corrected Q and U data to FITS files. Copies the supplied header,
    adding a note to the history saying that the correction was applied.
    """
    output_header=header.copy()
    output_header.add_history('Corrected for ionospheric Faraday rotation using ionosphereFR_correct.')

    #Get data back to original axis order, if necessary.
    N_dim=output_header['NAXIS'] #Get number of axes
    freq_axis=find_freq_axis(output_header)
    if freq_axis != 0:
        Qcorr=np.moveaxis(Qcorr,0,N_dim-freq_axis)
        Ucorr=np.moveaxis(Ucorr,0,N_dim-freq_axis)


    pf.writeto(Qoutputname,Qcorr,output_header,overwrite=overwrite)
    pf.writeto(Qoutputname,Qcorr,output_header,overwrite=overwrite)

    


def correct_cubes(Qdata,Udata,correction):
    """
    Applies the ionospheric Faraday rotation correction to the Stokes Q/U
    data, derotating the polarization angle and renormalizing to remove
    depolarization. Note that this will amplify the noise present in the data,
    particularly if the depolarization is large (|corr| is small).
    Note that the 'correction' variable is the effect of the ionosphere 
    (Theta, in the derivation), so the correction process is to divide by this 
    variable to produce the output cubes.
    """
    
    Pdata=Qdata+1.j*Udata #Input complex polarization
    Pcorr=Pdata/correction 
    Qcorr=Pcorr.real
    Ucorr=Pcorr.imag
    
    return Qcorr,Ucorr


def apply_correction_to_files(Qfile,Ufile,correctionfile,Qoutfile,Uoutfile,overwrite=False):
    """
    This function combines all the individual steps needed to apply a 
    correction to a set of Q and U FITS cubes and save the results.
    The user should supply the paths to all the files as specified.
    The output files can probably(?) be the same as the input files, if you
    want to overwrite the existing files.
    """
    
    #Get all data:
    correction=read_correction(correctionfile)
    Qdata,Udata,header=readData(Qfile,Ufile)
    
    #Checks for data consistency.
    if (Qdata.shape != Udata.shape):
        raise Exception("Q and U files don't have same dimensions.")
    if Qdata.shape[0] != correction.size:
        raise Exception("Correction file does not have same number of channels as FITS cube.")

    #Apply correction
    Qcorr,Ucorr=correct_cubes(Qdata,Udata,correction)
    
    #Save results
    write_corrected_cubes(Qoutfile,Uoutfile,Qcorr,Ucorr,header,overwrite=overwrite)




def main():
    """
    When invoked from the command line, parse the input options to get the
    filenames and other parameters, then invoke apply_correction_to_files
    to run all the steps and save the output cubes.
    """
    
    import argparse
    import os
    #TODO
    descStr = """
    Apply correction for ionospheric Faraday rotation to Stokes Q and U FITS
    cubes. Requires the file names for the input cubes, output cubes, and the
    correction file (which contains the ionospheric effect per channel 
    in the cubes).
    """

    parser = argparse.ArgumentParser(description=descStr,
                                 formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("fitsQ",dest="fitsQ",
                        help="FITS cube containing (uncorrected) Stokes Q data.")
    parser.add_argument("fitsU",dest="fitsU",
                        help="FITS cube containing (uncorrected) Stokes U data.")
    parser.add_argument("correctionFile",dest="correctionFile",
                        help="File containing ionospheric correction to be applied.")
    parser.add_argument("Qcorrected",dest="outQ",
                        help="Output filename for corrected Stokes Q cube.")
    parser.add_argument("Ucorrected",dest="outU",
                        help="Output filename for corrected Stokes U cube.")

    parser.add_argument("-o",dest="overwrite",action="store_true",
                        help="Overwrite exising output files? [False]")

    args = parser.parse_args()

    #Check for file existence.
    if not os.path.isfile(args.fitsQ):
        raise Exception("Stokes Q file not found.")
    if not os.path.isfile(args.fitsU):
        raise Exception("Stokes U file not found.")
    
    #Pass file names into do-everything function.
    apply_correction_to_files(args.fitsQ,args.fitsU,args.correctionFile,
                              args.outQ,args.outU,overwrite=args.overwrite)








if __name__ == "__main__":
    main()


