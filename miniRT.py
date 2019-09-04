#! /usr/bin/env python

import sys, os
import argparse
sys.path.append(os.path.dirname(__file__) + "/lib/")
import calcflux as cf
import plotter  as p

"""
miniRT is controlled via command-line arguments.

Example execution from the miniRT directory:
miniRT.py --atm inputs/broadening.atm --par inputs/oneline.par --fname broadening --wnlo 4366 --wnhi 4370 --wnsamp 0.005 --osamp 1080 --outdir output/broadening/

You can run this from a file. See run.sh in the 'examples' subdirectory.

"""


if __name__ == '__main__':
    # Command-line argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--atm",    dest="atm",    type=str,   default=None, 
                        help="Atmospheric model file")
    parser.add_argument("--par",    dest="par",    type=str,   default=None, 
                        help="HITRAN line list files", action='append')
    parser.add_argument("--wnlo",   dest="wnlo",   type=float, default=None, 
                        help="Minimum wavenumber for the spectrum")
    parser.add_argument("--wnhi",   dest="wnhi",   type=float, default=None, 
                        help="Maximum wavenumber for the spectrum")
    parser.add_argument("--wnsamp", dest="wnsamp", type=float, default=1.0,  
                        help="Wavenumber sampling interval")
    parser.add_argument("--osamp",  dest="osamp",  type=int,   default=2160, 
                        help="Wavenumber oversampling amount for high-res " + \
                             "opacity spectra")
    parser.add_argument("--fname",  dest="fname",  type=str,   default=None, 
                        help="Filename prefix for output")
    parser.add_argument("--saveflux", dest="saveflux", action='store_true', 
                        default=True, help="Boolean to save flux spectrum")
    parser.add_argument("--saveopa",  dest="saveopa",  action='store_true', 
                        default=True, help="Boolean to save opacity")
    parser.add_argument("--outdir",   dest="outdir", type=str, 
                        default='./output/', help="Directory to store outputs")
    parser.add_argument("--makeplots",  dest="makeplots",  action='store_true', 
                        default=True, help="Boolean to plot spectrum")
    # Parse the inputs
    known, unknown = parser.parse_known_args(sys.argv[1:])

    # Set them into variables
    for var in dir(known):
        exec("{:s} = known.{:s}".format(var, var))

    # Make sure `outdir` is properly set
    if outdir[-1] != '/':
        outdir = outdir + '/'
    if outdir[0] != '/' and outdir[0] != '.':
        outdir = './' + outdir
    # Create `outdir` if it does not already exist
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    # Set the names for the output, if requested
    if saveflux:
        saveflux = fname
    if saveopa:
        saveopa  = fname

    # Compute the spectrum
    out = cf.calcflux((wnlo, wnhi), atm, par, wnsamp=wnsamp, osamp=osamp, 
                      saveflux=saveflux, saveopa=saveopa, outdir=outdir)

    if makeplots:
        p.plotspec(outdir+saveflux+'.dat', outdir+saveflux+'.png')





