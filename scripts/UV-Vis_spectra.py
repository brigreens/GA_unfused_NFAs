# make UV-Vis spectra

import sys
import imp
import numpy as np
import matplotlib.pyplot as plt
import cclib
import os
import glob
import csv

# script given by Geoff
def spectra(etens, etoscs, low = 0.35, high = 10.0, resolution = 0.01, smear = 0.04):
    """Return arrays of the energies and intensities of a Lorentzian-blurred spectrum"""

    maxSlices = int((high - low) / resolution) + 1
    peaks = len(etens)

    spectraEV = []
    spectraNM = []
    spectraIntensity = []

    # eV = wavenumbers * 1.23981e-4
    # nm = 1.0e7 / wavenumbers

    for i in range(0, maxSlices):
        # in eV
        energy = float(i * resolution + low)
        wavenumber = energy / 1.23981e-4
        intensity = 0.0
        for trans in range(0, len(etens)):
            this_smear = smear / 0.2 * (-0.046 * etoscs[trans] + 0.20)
            #            print this_smear
            deltaE = etens[trans] * 1.23981e-4 - energy
            intensity = intensity + etoscs[trans] * this_smear**2 / (deltaE**2 + this_smear**2)

        spectraEV.append(energy)
        spectraNM.append(float(1.0e7 / wavenumber))
        spectraIntensity.append(intensity)
        
    return spectraEV, spectraNM, spectraIntensity

# code found at http://pychem.rcnelson.com/UV/uv_final.py
# some lines modified by Brianna

def uv_vis_plot(spectraNM, spectraIntensity, filename, max_nm):
    plt.clf()
    plt.rc('font', size=16)
    plt.rc('axes', linewidth=1.5)
    plt.rc('xtick.major', size=8, width=1.5)
    plt.rc('ytick.major', size=8, width=1.5)

    spectraNM = np.array(spectraNM, dtype=float)
    spectraIntensity = np.array(spectraIntensity, dtype=float)

    plt.plot(spectraNM, spectraIntensity, color = '#CA4327', lw=2)
    plt.xlim(250,max_nm + 250)
    plt.ylim(-0.005, 2)
    #plt.xlim(250,600) #for zoomed in on fullerene spectra
    #plt.ylim(-0.005, 0.25)

    plt.title('UV-Vis Spectrum of ' + filename)
    plt.xlabel('Wavelength (nm)')
    plt.ylabel('Oscillator Strength')
    
    plt.gcf().subplots_adjust(bottom=0.25)
    plt.tight_layout
    plt.savefig('/ihome/ghutchison/blp62/GA/running_GA/v5/absorption_spectra/' + filename + '.png')
    
    #plt.show()
def make_csv(filename, spectra_NM, spectra_Intensity):
    
    path = '/ihome/ghutchison/blp62/GA/running_GA/opt_GA/top10/' + filename + '.csv'

    with open(path, "w") as csvoutput:
        writer = csv.writer(csvoutput, lineterminator = '\n')
        writer.writerow(["Wavelength (nm)", "Oscillator strength"])
        
        for i in range(len(spectra_NM)):
            data = []
            data.append(round(spectra_NM[i], 2))
            data.append(spectra_Intensity[i])
            #print(data)
            writer.writerow(data)


def parse_sTDDFT(filename):
    '''
    Parses through sTD-DFT output files

    Parameters
    -----------
    filename: str
        path to output file
    donor_or_acc: str
        define if molecule is acceptor ('acc') or donor ('don')

    Returns
    -------
    acceptors_sTD[filename] or donors_sTD[filename]: list
        adds list of sTDDFT descriptors to a dictionary for the molecule
    '''
    outputs = []
    with open(filename, 'r', encoding = 'utf-8') as file:
        line = file.readline()
        oscs = []
        wavelength = []
        energyEV = []
        wavenumber = []
        while line:

            if 'ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS' in line:
                for x in range(5):
                    line = file.readline()
                opt_bg = float(line[6:15])                    
                while line != '\n':
                    oscs.append(float(line[25:37]))
                    wavelength.append(float(line[17:23]))
                    energyEV.append(float(line[6:15]) * 1.2398e-4)
                    wavenumber.append(float(line[6:15]))
                    line = file.readline()

                break


            line = file.readline()  
        line = file.readline()
        return [wavenumber, oscs]


# for directory of molecules to plot

for file in glob.iglob('/ihome/ghutchison/blp62/GA/running_GA/v5/top10/*.out'):
    props = parse_sTDDFT(file)
    wavenumbers = props[0]
    oscs = props[1]

    min_wavenumber = min(wavenumbers)
    max_nm = 10000000 / min_wavenumber

    filename = file.split('/')[-1].split('.')[0]

    (spectraEV, spectraNM, spectraIntensity) = spectra(wavenumbers, oscs)

    uv_vis_plot(spectraNM, spectraIntensity, filename, max_nm)


# for individual files

'''filename = '/ihome/ghutchison/blp62/OPEP/Github/OPEP/output_files/sTDDFT/sTDDFT_donors/PBDB-T-SF_olig.out'
#file = cclib.io.ccopen(filename)

props = parse_sTDDFT(filename)
wavenumbers = props[0]
oscs = props[1]

max_nm = 1000


filename = 'PBDB-T-SF'

(spectraEV, spectraNM, spectraIntensity) = spectra(wavenumbers, oscs)

uv_vis_plot(spectraNM, spectraIntensity, filename, max_nm)'''


