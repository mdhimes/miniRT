import sys, os
import numpy as np
import scipy.interpolate as si
import matplotlib.pyplot as plt

plt.ion()

"""
Functions to plot the output of miniRT.

plotspec(): Plots a single spectrum.

compspecs(): Plots multiple spectra, and computes the difference between 
             the first spectrum and each of the others.
"""


def plotspec(spec, outname, wl=True):
    """
    Plots a caclulated spectrum.

    Inputs
    ------
    spec   : string. Path to text file containing the spectrum.
                     spec[:,0] provides the wavelengths/wavenumbers.
                     spec[:,1] provides the spectrum in erg/s/cm.
    outname: string. Savename for the plot.
    wl     : bool.   Determines whether `spec` has wavelengths (True) 
                     or waveumbers (False)
    """
    spec = np.loadtxt(spec)
    plt.plot(spec[:,0], spec[:,1])
    if wl:
        plt.xlabel(u'Wavelength (\u00b5m)')
    else:
        plt.xlabel('Wavenumber (cm-1)')
    plt.ylabel('Flux (erg/s/cm)')
    plt.savefig(outname)
    plt.close()


def compspecs(specs, labels, outname, wl=True):
    """
    Plots multiple spectra, and difference(s) w.r.t. the first spectrum.
    If the X axis values do not match, spectra are linearly interpolated to 
    the first spectrum's sampling for difference calculation(s).

    Inputs
    ------
    specs  : list, strings. Text files w/ spectrum data.
    labels : list, strings. Names for each spectrum.
    outname: string. Savename for the plot.
    wl     : bool.   Determines whether `specs` has wavelengths (True) 
                     or waveumbers (False)
    """
    # Load the data
    xax = []
    yax = []
    for foo in specs:
        dat = np.loadtxt(foo)
        xax.append(dat[:,0])
        yax.append(dat[:,1])

    fig0 = plt.figure(0, (6,5))
    plt.clf()
    # Spectra
    frame1 = fig0.add_axes((.14, .3, .8, .65))
    for i in range(len(specs)):
        plt.plot(xax[i], yax[i], label=labels[i])
    plt.legend(loc='best')
    plt.ylabel("Flux (erg s$^{-1}$ cm$^{-1}$)")
    frame1.set_xticklabels([])
    # Residuals
    frame2 = fig0.add_axes((.14, .1, .8, .2))
    resid = []
    for i in range(1, len(specs)):
        if xax[0] != xax[i]:
            f   = si.interp1d(xax[i], yax[i])
            yre = f(xax[0])
        else:
            yre = yax[i]
        plt.plot(xax[0], (yax[0] - yre) / yax[0] * 100, ls=":", label=labels[i])
    plt.ylabel("Difference (%)")
    if wl:
        plt.xlabel(u"Wavelength  (\u00b5m)")
    else:
        plt.xlabel('Wavenumber (cm-1)')
    plt.legend(loc='best')
    plt.savefig(outname)
    plt.close()


