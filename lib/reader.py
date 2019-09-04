import sys
import numpy as np

"""
This file contains functions related to reading inputs for miniRT.

readatm()
  Reads an atmospheric model file into a Numpy array. 
  Format follows Jasmina Blecic's TEA.

readpar()
  Reads HITRAN .par file and extracts necessary info into a Numpy array.

readmol()
  Reads Transit's molecules.dat file, creates a dictionary with molecules as 
  keywords that link to the mass and diameter of the molecule.

readhit()
  Reads the HITRAN file that links molecules and isotopologues to numbers.
"""

def readatm(atmfile, verb=0):
    """
    This function reads an atmospheric file into a Numpy array.

    Inputs
    ------
    atmfile: string. Path/to/file for the atmospheric file used.
    verb   : int. Verbosity parameter to determine debug output messages (0 -- 5)

    Returns
    -------
    mols   : list of strings. Molecules present in the atmosphere.
    atminfo: array. Shape is (nlayers, nmolecules + 3). For a given layer, 
                    the data stored is radius, pressure, temperature, 
                    mol_1, mol_2, ...
                    Layers are ordered top to bottom. That is, the 0th index 
                    is the top-most layer.

    Example
    -------
    This example assumes that the directory containing Transit is parallel to 
    the current.
    >>> atm = "../transit/transit/examples/demo/HD209458b_demo.atm"
    >>> mols, atminfo, ur, up = atmreader.readatm(atm)

    Revisions
    ---------
    2017-10-02 mhimes@knights.ucf.edu   Initial implementation.
    2018-10-09 mhimes                   Added verbosity parameter.
    """
    # Open the atm file and read it
    if verb > 3:
        print("Reading .atm file...")
    atmfoo = open(atmfile, 'r')
    lines  = atmfoo.readlines()

    if verb > 0:
        print("Successfully read .atm file")

    # Store conversion factors
    ur   = float(lines[4].split()[1]) # radius   conversion factor
    up   = float(lines[5].split()[1]) # pressure conversion factor
    # Store the molecule info
    mols = lines[9].split()

    # Find where the layer data starts
    for i in range(len(lines)):
        line = lines[i].split()
        try:
            if line[0] == '#Radius':
                break
        except:
            continue

    if verb > 3:
        print(".atm layer data begins on line ", i)

    # Trim atm file
    lines = lines[i + 1:]

    # Array to hold atm file info
    atminfo = np.zeros((len(lines), len(mols) + 3), dtype=float)
    # Note that +3 is to hold radius, pressure, and temperature

    # Read in data
    if verb > 3:
        print("Reading .atm data into arrays...")

    for i in range(len(lines)):
        # Split info for layer `i`
        line = lines[i].split()
        for j in range(3 + len(mols)): # rad, press, temp, mol1, mol2, ...
            # Read in each column
            atminfo[i, j] = float(line[j])

    if verb > 0:
        print(".atm data read into array.")

    # Check if atmosphere is ordered from bottom to top. If so, reverse it
    if atminfo[0, 0] < atminfo[1, 0]:
        if verb > 0:
            print("Atmosphere is ordered bottom to top. Reversing it to be top to bottom.")
        for i in range(atminfo.shape[1]):
            atminfo[:, i] = atminfo[:, i][::-1]

    return mols, atminfo, ur, up


def readpar(parfile, verb=0):
    """
    This function reads HITRAN .par file(s) into a Numpy array.

    Inputs
    ------
    parfile: list of strings. Paths/to/files for HITRAN line lists.
    verb   : int. Verbosity parameter to determine debug output messages (0 -- 5)

    Returns
    -------
    molisoID: array. Molecule and isotope ID number as given in HITRAN .par file
    wavenum : array. Wavenumber of the transition.
    eincoA  : array. Einstein A coefficient.
    elow    : array. Energy of the lower state.
    lowstat : array. Statistical weight of the lower state.

    Example
    -------
    

    Revisions
    ---------
    2017-10-02 mhimes@knights.ucf.edu   Initial implementation.
    2017-10-03  mhimes                  Replaced 2D data structure w/ 1D arrays
                                        for clarity.
    2018-10-09 mhimes                   Added verbosity parameter.
    """
    # Read the file(s)
    if verb > 3:
        print("Reading .par file...")

    lines = []
    for par in parfile:
        parfoo = open(par, 'r')
        lines.append(parfoo.readlines())
        if verb > 0:
            print(".par file " + par + " read successfully.")

    # Extract the relevant information
    # Note that lines is a list of lists, each with potentially varying 
    # numbers of lines
    size     = sum(len(lines[i])      for i in range(len(lines)))
    size    -= sum(lines[i][-1]=='\n' for i in range(len(lines)))
    molisoID = np.zeros(size, dtype=int)
    wavenum  = np.zeros(size, dtype=float)
    eincoA   = np.zeros(size, dtype=float)
    elow     = np.zeros(size, dtype=float)
    lowstat  = np.zeros(size, dtype=float)

    if verb > 3:
        print("Arrays for .par info allocated successfully")

    # Read the values into the array
    n = 0
    for i in range(len(lines)):
        for j in range(len(lines[i])):
            line        = lines[i][j]
            if line=='\n':
                continue
            molisoID[n] = int(  line[  0:  3].strip(' ')) #molecule and iso number
            wavenum[n]  = float(line[  3: 15].strip(' ')) #transition wavenumber
            eincoA[n]   = float(line[ 25: 35].strip(' ')) #Einstein A coefficient
            elow[n]     = float(line[ 45: 55].strip(' ')) #energy of lower-state
            lowstat[n]  = float(line[153:160].strip(' ')) #lower stat weight
            n          += 1

    if verb > 0:
        print("Arrays for .par info filled successfully")

    return molisoID, wavenum, eincoA, elow, lowstat


def readmol(molfile, verb=0):
    """
    This function reads Transit's molecules.dat file and creates a dictionary 
    with molecules as keywords and links each to the mass and diameter of the 
    molecule.

    Inputs
    ------
    molfile: list of strings. Paths/to/file for molecule data file.
    verb   : int. Verbosity parameter to determine debug output messages (0 -- 5)

    Returns
    -------
    moldict: dictionary. Dictionary keywords are molecules (i.e. H2O). Each is 
             associated with a tuple of molar mass [g/mol] and diameter [A].

    Example
    -------
    

    Revisions
    ---------
    2017-10-02 mhimes@knights.ucf.edu   Initial implementation.
    2018-10-09 mhimes                   Added verbosity parameter.
    """
    # Read the file
    if verb > 3:
        print("Reading molecule file...")
    molfoo = open(molfile, 'r')
    lines  = molfoo.readlines()

    if verb > 0:
        print("Molecule file successfully read.")

    # Find where to cut off headers
    for i in range(len(lines)):
        if (lines[i][0] == '#') or (lines[i] == '\n'):
            continue
        else:
            break

    if verb > 3:
        print("Molecule file headers end on line ", i)

    # Trim off the headers
    lines = lines[i:]

    # Array to hold info to be converted to dictionary
    molinfo = np.zeros((len(lines), 3), dtype=object)

    # Read in relevant info
    if verb > 3:
        print("Reading molecule file data into arrays...")
    for i in range(len(lines)):
        line = lines[i].split()
        molinfo[i, 0] = line[1] # molecule name
        molinfo[i, 1] = line[2] # molecule mass
        molinfo[i, 2] = line[3] # diameter

    if verb > 0:
        print("Molecule file data successfully read into arrays.")

    # Create a dictionary from this
    # Diameter is multiplied by 1e-8 to convert from Angstroms to cm
    if verb > 3:
        print("Creating dictionary of molecule file data...")
    moldict = dict(zip(molinfo[:, 0], zip(molinfo[:, 1].astype(float), \
                                          molinfo[:, 2].astype(float)*1e-8)))
    if verb > 0:
        print("Created dictionary of molecule file data.")

    return moldict


def readhit(hitfile, verb=0):
    """
    This function reads a file with HITRAN information relating molecule number 
    and isotopologue number to a molecule.

    Inputs
    ------
    hitfile: list of strings. Paths/to/file for HITRAN info file.
    verb   : int. Verbosity parameter to determine debug output messages (0 -- 5)

    Returns
    -------
    hitdict: dictionary. Dictionary keywords are molecule IDs (i.e. 11 for the 
             first isotope of H2O). Each is associated with a tuple of molecule 
             name, isotope's name, isotope's ratio, and the isotope's mass.

    Example
    -------
    

    Revisions
    ---------
    2017-10-02 mhimes@knights.ucf.edu   Initial implementation.
    2018-10-09 mhimes                   Added verbosity parameter.
    """
    # Read the file
    if verb > 3:
        print("Reading HITRAN file...")
    hitfoo = open(hitfile, 'r')
    lines  = hitfoo.readlines()

    if verb > 0:
        print("Successfully read HITRAN file.")

    # Find where to cut off headers
    for i in range(len(lines)):
        if (lines[i][0] == '#') or (lines[i] == '\n'):
            continue
        else:
            break

    if verb > 3:
        print("HITRAN file headers end on line ", i)

    # Trim off the headers
    lines = lines[i:]

    # Temp array to hold info for dictionary
    hitinfo = np.zeros((len(lines), 5), dtype=object)

    # Read in relevant info
    if verb > 3:
        print("Reading HITRAN file data into arrays...")
    for i in range(len(lines)):
        line          = lines[i].split()
        hitinfo[i, 0] = line[0] # molecule ID
        hitinfo[i, 1] = line[1] # molecule name
        hitinfo[i, 2] = line[2] # iso ID
        hitinfo[i, 3] = line[4] # iso ratio
        hitinfo[i, 4] = line[5] # iso mass

    if verb > 0:
        print("HITRAN file data successfully read into arrays.")

    # Modify molecule ID to match HITRAN format
    if verb > 3:
        print("Modifying molecule IDs to match HITRAN format...")
    hitinfo[0, 0] += '1' # set the first one
    val = 1
    for i in range(1, len(lines)):
        if hitinfo[i, 1] == hitinfo[i - 1, 1]:
            val += 1
        else:
            val = 1
        hitinfo[i, 0] += str(val)

    if verb > 0:
        print("Successfully modified molecule IDs to match HITRAN format.")

    # Build the dictionary
    if verb > 3:
        print("Creating dictionary of HITRAN file data...")
    hitdict = dict(zip(hitinfo[:, 0], zip(hitinfo[:, 1],               \
                                          hitinfo[:, 2].astype(int),   \
                                          hitinfo[:, 3].astype(float), \
                                          hitinfo[:, 4].astype(float))))

    if verb > 0:
        print("Created dictionary of HITRAN file data.")

    return hitdict


