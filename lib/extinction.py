import numpy as np
import scipy.constants   as const
import scipy.interpolate as si
import sys
import voigt

"""
This file contains functions related to atmospheric extinction.

extinction(): Function to calculate the extinction coefficients.

resamp(): Function to resample an array at different values.

downsamp(): Function to downsample an array to a lower resolution.

"""


def extinction(nu, P, T, molmass, gf, Elow, dia, dias, Z, atm, osamp, 
               molind, ratio):
    """
    This function calculates the extinction coefficients.

    Inputs
    ------
    nu     : float. Central wavenumber of the line (in cm-1).
    P      : array of floats. Pressure of layer in the atmosphere in cgs units.
    T      : array of floats. Temperature of layer in the atmosphere.
    molmass: float. Mass of the line-producing species.
    gf     : float. Weighted oscillator strength of the line.
    Z      : array of floats. Partition function.
    Elow   : float. Energy of the lower state.
    dia    : float. Diameter  of the  molecule.
    dias   : array. Diameters of each molecule in the atmosphere.
    atm    : array. Contains information of layer in atmospheric file. Format 
                    is [[abundance1, molarmass1], [abundance2, molarmass2], etc]
    osamp  : int.   Oversampling factor for the Voigt profile.
    molind : int.   Index of line-producing molecule in `atm` array.
    ratio  : float. Ratio of isotope of line-producing molecule.

    Outputs
    -------
    K    : array of floats. Broadened opacity spectrum.
    nurng: array of floats. Wavenumbers associated with `K`.
    """
    # Constants in CGS
    c    = const.c   *  100.
    e    = const.e   *   10. * const.c #Coulomb -> statC conv
    me   = const.m_e * 1000.
    R    = const.R   *    1e7
    h    = const.h   *    1e7
    k    = const.k   *    1e7
    amu  = const.physical_constants['atomic mass constant'][0] * 1000.
    Nava = const.N_A

    # Molecule information:
    mass     = molmass * amu        # mass in g
    nd       = P / (k * T)   # number density in cm-3

    # Gaussian HWHM: Goody (1995), Cubillos et al. BART paper eqn 21
    alpha = nu / c * (2. * np.log(2.) * k * T / mass)**0.5

    # Lorentzian HWHM: Goody (1995), Cubillos et al. BART paper eqn 22
    gamma = np.sum((dia/2. + dias/2.)**2 *                         \
                   (2. * k * T / np.pi)**0.5 / c *              \
                   nd * atm[:,0] *                              \
                   (1./mass + 1./(atm[:,1] * amu))**0.5)

    # Generate range of values to calculate the Voigt profile at
    if   (alpha >= gamma):
        nurng = np.linspace(nu - 20.*alpha, nu + 20.*alpha, num=osamp)
    elif (gamma >  alpha):
        nurng = np.linspace(nu - 20.*gamma, nu + 20.*gamma, num=osamp)

    # Calculate the broadened opacity, Cubillos et al. BART paper eqns 19, 20
    K = const.pi * e**2 / c**2 / me            *     \
        nd * ratio                             *     \
        gf / Z * np.exp(-h*c*Elow/k/T)         *     \
        (1 - np.exp(-h*c*nu/k/T))              *     \
        voigt.V(nurng, alpha, gamma, shift=nu)

    # Divide out the density
    K /= (P * molmass / k / T / Nava)

    return K, nurng


def resamp(K, nurng, nuspec, shift=False):
    """
    This function takes the computed opacity from extinction() and resamples it 
    to the desired resolution.

    Note that this function shifts the peak of the opacity spectrum to coincide 
    with the nearest value in `nuspec`

    Inputs
    ------
    K     : array of floats. Broadened opacity spectrum.
    nurng : array of floats. Wavenumbers associated with `K`.
    nuspec: array of floats. Wavenumbers of the desired output spectrum.
    shift : bool.            Determines whether line peaks should be shifted 
                             to a `nuspec` value.

    Outputs
    -------
    Kspec: array of floats. `K` resampled to correspond to `nuspec`.
    """
    # Make cubic spline w/ 0 outside `nurng`, evaluate at `nuspec`
    opa    = si.interp1d(nurng-shift, K, bounds_error=False, fill_value=0)
    Kspec  = opa(nuspec)

    return Kspec


def downsamp(K, nurng, nuspec):
    """
    This function takes the computed opacity from extinction() and downsamples
    it to the desired resolution.

    Inputs
    ------
    K     : array of floats. Broadened opacity spectrum.
    nurng : array of floats. Wavenumbers associated with `K`.
    nuspec: array of floats. Wavenumbers of the desired output spectrum.

    Outputs
    -------
    Kspec: array of floats. `K` resampled to correspond to `nuspec`.
    """
    # Find indices where nurng fits into nuspec
    Kspec = np.zeros(nuspec.shape)
    inds  = np.digitize(nurng, nuspec)
    m     = 0
    for wn in inds:
        if wn > 0 and wn < len(nuspec):
            Kspec[wn] += K[m]
        m += 1

    return Kspec


