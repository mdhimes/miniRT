import sys, os
import numpy as np
import matplotlib.pyplot as plt
#plt.ion()
import scipy.constants as const
import scipy.integrate as integrate

import reader as R
import extinction as E
sys.path.append(os.path.dirname(__file__) + "/../../BART/modules/transit/pylineread/src/pytips")
import pytips

"""
This file contains functions related to calculating the flux emitted by planets.

planck(): Function to compute the Planck function at a given wavenumber and 
          temperature.

calcflux(): Function to compute the flux emitted by a planet.

"""


# CONSTANTS
c    = const.c * 100 #cm s-1
k    = const.k * 1e7 #erg K-1
h    = const.h * 1e7 #erg s
Nava = const.N_A

def planck(wavenum, temp):
    """
    This function calculates the Planck function for an array of wavenumbers 
    at some temperature(s).

    Inputs
    ------
    wavneum: array, floats. Wavenumbers in cm-1 to calculate Planck function.
                            Can also be a single float.
    temp:    array, floata. Temperature of blackbody in Kelvin. 
                            Can also be a single float.

    Outputs
    -------
    array of Planck function values at `wavenum` and `temp`.
    Shape is according to the inputs. If `wavenum` is a float, 
    shape is 1D array of Planck functions for `wavenum` at various temperatures.
    If 'temp' is a float, shape is 1D array of Planck functions for `temp` at 
    various wavenumbers.
    If both `wavenum` and `temp` is a float, then the shape is (wavenum, temp).
    """
    if type(wavenum)==np.ndarray and type(temp)==np.ndarray:
        wavenum = wavenum.reshape(-1, 1)
        temp    = temp   .reshape( 1,-1)

    return 2. * h * c**2 * wavenum**3 / (np.exp(h*wavenum*c/k/temp) - 1)


def calcflux(wnrng, atm, parfile, 
             molfile=os.path.dirname(__file__) + '/../inputs/molecules.dat', 
             hitfile=os.path.dirname(__file__) + '/../inputs/litran.dat',  
             wnsamp=1., osamp=2160, angles=np.array([0, 20, 40, 60, 80]), 
             toomuch=1., wl=True, saveflux=False, saveopa=False, 
             outdir='./', verb=0):
    """
    This function calculates the flux by summing the intensity for an array 
    of angles.
    Equation used is equation 18 in Patricio Cubillos's 2017 BART paper.

    Inputs
    ------
    wnrng  : tuple.  (wavenumber_min, wavenumber_max). Range of wavenumbers for 
                     which to calculate the flux.
    atm    : string. Path/to/file for the atmospheric file
    parfile: list of strings. List of paths/to/files for the line list.
                              Note: MUST be HITRAN format! MUST be a list!
    molfile: string. Path/to/file for the file containing molecule info.
    hitfile: string. Path/to/file for the file containing HITRAN info.
    wnsamp : float.  Sampling interval for wavenumbers.
    osamp  : int.    Oversampling factor when calculating the Voigt profile.
    angles : array.  Angles to calculate the intensity at.
    wl     : bool.   If True, outputs flux with respect to wavelength. 
                     If False, outputs with respect to wavenumber.
    saveflux: bool/str. If False, does not save file w/ wavelengths and flux. If a string, 
                        saves wavelengths and flux to file with name `saveflux`. Do not 
                        include file extension.
    saveopa:  bool/str. Same as saveflux, but for the opacity array.
    verb   : int.       Verbosity parameter to determine debug output messages (0 -- 5)

    Outputs
    -------
    flux: array. flux values for each of the input wavenumbers.
    """   
    # Check number of supplied angles
    if len(angles) <= 2:
        # Force user to specify at least 3 angles
        print("Results are not worthy of your trust!\n")
        print("Specify more angles.\n")
        sys.exit()
    # Convert angles to radians
    else:
        radang = angles * np.pi / 180.

    # Load the atmospheric file
    try:
        mols, atminfo, ur, up = R.readatm(atm, verb)
        # Convert the radii and pressures to CGS using `ur` and `up`
        atminfo[:, 0] *= ur
        atminfo[:, 1] *= up
    except Exception as e:
        print(e)
        print("Unable to read the atmospheric file.\n")
        sys.exit()

    # Load the line database(s)
    try:
        molisoID, wavenum, eincoA, elow, lowstat = R.readpar(parfile, verb)
    except Exception as e:
        print(e)
        print("Unable to read HITRAN .par file(s).\n")
        sys.exit()

    # Load the molecule database
    try:
        moldict = R.readmol(molfile, verb)
    except Exception as e:
        print(e)
        print("Unable to read the molecule database file.\n")
        sys.exit()

    # Load HITRAN database
    try:
        hitdict = R.readhit(hitfile, verb)
    except Exception as e:
        print(e)
        print("Unable to read the HITRAN molecule database file.\n")
        sys.exit()

    # Array to hold wavelengths and flux data
    flux = np.zeros((2, int(round((wnrng[1] - wnrng[0] + wnsamp)/wnsamp))), dtype=float)

    # Store the wavelengths or wavenumbers; endpoints inclusive
    wnums = np.arange(wnrng[0], wnrng[1]+wnsamp/2., wnsamp, dtype=float)

    if wl:
        flux[0] = 10000. / wnums
    else:
        flux[0] = wnums

    # Calculate weighted oscillator strength for each line
    C1 = 4. * const.epsilon_0 * const.m_e * const.c**2 / const.e**2 * 0.01
    gf = eincoA * lowstat * C1 /                                               \
         (8.0 * np.pi * const.c * 100.0) / wavenum**2.0

    # Array to hold opa (opacity)
    opa = np.zeros((len(molisoID), len(atminfo[:,2]),      \
                    len(wnums)), dtype=float)

    # Set up `atmlayers` with molecular masses (abundances filled in later)
    atmlayers = np.zeros((atminfo.shape[0], len(mols), 2), dtype=float)
    for i in range(len(mols)):
        atmlayers[:, i, 1] = moldict[mols[i]][0]

    # Calculate the extinction coefficient for all lines
    print("Calculating extinction coefficients...")
    for i in range(molisoID.shape[0]):
        # Calculate opacity for each layer of the atmosphere
        for lay in range(atminfo.shape[0]):
            # Calculate the partition function for this molecule at this temp
            Z = pytips.tips(int(str(molisoID[i])[0]),     # molecule HITRAN ID
                            hitdict[str(molisoID[i])][1], # isotope ID
                            atminfo[lay, 2])              # temperature
            
            # Atmosphere layer `lay` abundance info
            atmlayers[lay, :, 0] = atminfo[lay][-len(mols):]

            # Line molecule
            thismol     = hitdict[str(molisoID[i])][0] #isotope name
            molisoratio = hitdict[str(molisoID[i])][2] #isotope ratio
            molind      = mols.index(thismol)          #index of this molecule

            # Array of diameters for collisional broadening
            dias = np.zeros(len(mols))
            for j in range(len(mols)):
                dias[j] = moldict[mols[j]][1]
            
            # Calculate the high-res opacity spectrum for this layer
            opaHR, rng = E.extinction(
                           wavenum[i],        #wavenumber->freq
                           atminfo[lay, 1], atminfo[lay, 2],     #pressure, temp
                           hitdict[str(molisoID[i])][3],         #isotope mass
                           gf[i],                  #weighted oscillator strength
                           elow[i],                #energy of lower state
                           moldict[hitdict[str(molisoID[i])][0]][1],   #diameter
                           dias,                          #diameters of all mols
                           Z,                             #partition function
                           atmlayers[lay],                #atmosphere layer info
                           osamp,                         #oversampling factor
                           molind, molisoratio)           #index of molecule in 
                                                          #atm file, iso ratio          
            # Resample onto `wnums` grid
            opa[i, lay] += E.resamp(opaHR, rng, wnums)

    # Save the opacity array, if requested
    if saveopa:
        np.save(outdir + 'opacity_' + saveopa + '.npy', opa)

    # Calculate the optical depth until `toomuch`
    bot    = int(len(atminfo[:,2])-1)
    layopa = np.zeros((opa.shape[2], atminfo[:,0].shape[0]))
    tau    = np.zeros((opa.shape[2], atminfo[:,0].shape[0]))
    # assume it will reach bottom of atm
    itau   = np.zeros(opa.shape[2], dtype=int) + bot
    print("Calculating optical depth until toomuch...")
    for wn in range(opa.shape[2]):
        for lay in range(6, len(atminfo[:,2])):
            P      = atminfo[lay, 1]
            T      = atminfo[lay, 2]
            # opacity for all lines at this wn
            for i in range(len(molisoID)):
                if opa[i,lay,wn] == 0:
                    continue
                else:
                    thismol = hitdict[str(molisoID[i])][0] #isotope name
                    molind  = mols.index(thismol)          #molecule index
                    abun    = atmlayers[lay, molind, 0]
                    molmass = atmlayers[lay, molind, 1]
                    # Add the opacity * mass density
                    layopa[wn,lay] += opa[i, lay, wn] * \
                                      (P * abun * molmass / k / T / Nava)
            # Top of atmosphere, nothing to integrate
            if lay==0:
                continue
            # Integrate optical depth for this layer
            # Number of layers traversed, and radius array
            nrad = lay + 1
            # Need at least 3 entries -- take mean of the two layers
            if nrad==2:
                rad  = np.array([atminfo[1,0], (atminfo[1,0]+atminfo[0,0])/2, 
                                 atminfo[0,0]])
                nrad = 3
            else:
                rad = atminfo[:,0][:nrad][::-1]
            # Integrate tau
            laytau = integrate.simps(layopa[wn,:lay+1][::-1], rad, even='last')
            tau[wn,lay] += laytau
            # Stop at `toomuch`
            if tau[wn,lay] >= toomuch:
                itau[wn] = lay # update the bottom-most layer to calc until
                continue

    # Calculate Planck func and dtau for each layer
    B    = planck(wnums, atminfo[:,2])
    dtau = np.exp(-tau.reshape(len(wnums),len(atminfo[:,2]), 1) / 
           np.cos(radang.reshape(1,1,-1)))

    # Calculate intensity
    print("Calculating intensity at various points on the planet...")
    intens = np.zeros((tau.shape[0], len(radang)))
    for wn in range(len(wnums)):
        for i in range(len(radang)):
            intens[wn, i] = B[wn,itau[wn]] * dtau[wn,itau[wn],i] - \
                            np.trapz(B[wn,:itau[wn]], dtau[wn,:itau[wn], i])

    # Grid of angles
    area_grid = np.zeros(len(radang) + 1)
    area_grid[ 0] = 0
    area_grid[-1] = np.pi/2
    for i in range(1, len(area_grid)-1):
        area_grid[i] = (radang[i-1] + radang[i]) / 2
    # Grid of areas
    area = np.zeros(len(radang))
    for i in range(len(area)):
        area[i] = np.sin(area_grid[i+1])**2 - np.sin(area_grid[i])**2

    # Calculate the fulx
    print("Calculating the flux spectrum...")
    for wn in range(opa.shape[2]):
        for i in range(len(area)):
            flux[1, wn] += np.sum(np.pi * intens[wn,i] * area[i])

    # Save the wavenumber vs flux
    if saveflux:
        if wl:
            np.savetxt(outdir + saveflux + '.dat', flux.T, fmt='%.11e', 
                       header='Wavelength (um) Flux (erg/s/cm)')
        else:
            np.savetxt(outdir + saveflux + '.dat', flux.T, fmt='%.13e', 
                       header='Wavenumber (cm-1) Flux (erg/s/cm)')

    return flux


