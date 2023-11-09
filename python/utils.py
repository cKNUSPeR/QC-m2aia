
# namespace and import declaration
import m2aia as m2
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import random as rnd
import statistics as stat
import scipy.stats as SST
import scipy.signal as SSI
import skimage.measure as skim
import pandas as pd
import warnings
from scipy.signal import argrelextrema
# hide the warnings, nice for FutureWarnings
warnings.filterwarnings('ignore')
import matplotlib.backends.backend_pdf
#from reportlab.lib.pagesizes import letter
#from reportlab.pdfgen import canvas



# M2aia-independent tools
def evaluate_formats(file_format  # metadata_string
                      ):
    """Evaluates a file format string and return a dict of flags."""

    # instance flags dict
    flags = dict()

    # check for different flags
    if "Profile" in file_format:
        flags["profile"] = True
        flags["centroid"] = False
    elif "Centroid" in file_format:
        flags["profile"] = False
        flags["centroid"] = True
    else:
        raise ValueError(
            "The loaded file has an undefined spectrum type.\n Please check for accessions 'MS:1000128' or 'MS:1000127'")

    if "Processed" in file_format:
        flags["processed"] = True
        flags["continuous"] = False
    elif "Continuous" in file_format:
        flags["processed"] = False
        flags["continuous"] = True
    else:
        raise ValueError(
            "The loaded file has an undefined alignment type.\n Please check for accessions 'IMS:1000030' or 'IMS:1000031'")

    return flags


def evaluate_image_corners(ndarray):
    """Givel the values of the corners of the pixel-filled ndarray """
    pix_pos = np.argwhere(ndarray)

    # get the corners of data entry, is useful to set limits of plotting
    x_min = pix_pos[np.argmin(pix_pos[:, 1])][1]
    x_max = pix_pos[np.argmax(pix_pos[:, 1])][1]
    y_min = pix_pos[np.argmin(pix_pos[:, 0])][0]
    y_max = pix_pos[np.argmax(pix_pos[:, 0])][0]

    return (x_min, x_max), (y_min, y_max)


def make_subsample(samplenumber: int, percent_sample:float) -> list:
    """Makes a subsample out of the samplenumber with the given percentage (as float)"""
    # Determine the size of a batch
    batch_size = int(samplenumber * percent_sample)
    # get random numbers according to the batch size
    return rnd.sample(range(0, samplenumber), batch_size)


def mask_bad_image(key_list,  # an iterable of valid pixel indices,
                   val_list,  # an iterable of projected Intensities, matched to key_list
                   image  # An array-like object with the given distribution of key variables
                   ):
    """make a mask approach to plot any feature based on mapping onto existing image array with a translation apporach.
    It transfers pixels from 0 (in binary image input ) to NaN, which allows them to be set to bad"""
    # set up a translational dictionary
    trans_dict = dict(zip(key_list, val_list))
    # Important zero-index conversion, otherwise rounding gives error
    trans_dict[0] = np.nan

    # defines the callable function (juhu, we love functional programming
    translate = np.vectorize(lambda ele: trans_dict.get(ele, ele))

    return translate(image)


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]



def read_calibrants(filepath: str) -> tuple(list,list):
    """Reads calibrant files and gives a list of name and thr. mz values
    INVARIANTS: Needs a header column with 'Name' and a col with 'Theoretical m/z'."""
    cal = pd.read_csv(filepath, sep=';', header=0)
    cal_masses = cal["Theoretical m/z"].tolist()
    cal_names = cal["Name"].tolist()
    return (cal_names,cal_masses)

def extract_calibrant_spectra(Image, cal_masses, subsample, mz_bin):
    """Read the list of cal masses and collects the spectral data in the given mz bin."""
    accu_list = []
    for i in range(len(cal_masses)):
        accu_list.append(np.array([[],[]]))

    # and looping over calibrants
    for i, calmass in enumerate(cal_masses):
        # looping over sample
        for ind in subsample:
            mass, intensity = Image.GetSpectrum(ind)
            try:
                #mindex needs to be inclusive
                mindex = min(np.where(mass > (calmass-mz_bin))[0])
                mindex_flag = True
            except:
                mindex_flag = False
            try:
                #maxdex needs to be exclusive
                maxdex = min(np.where(mass > (calmass+mz_bin))[0])
                maxdex_flag = True
            except:
                maxdex_flag = False

            # pixels are only written if there is data present in the specified part
            if maxdex_flag and mindex_flag:
                # collecting of masses and intensities
                adder = np.array((mass[mindex:maxdex],intensity[mindex:maxdex]))
                accu_list[i] = np.concatenate((accu_list[i], adder), axis=1)
            elif mindex_flag and not maxdex_flag:
                adder = np.array((mass[mindex:], intensity[mindex:]))
                accu_list[i] = np.concatenate((accu_list[i], adder), axis=1)
            elif not mindex_flag and maxdex_flag:
                adder = np.array((mass[:maxdex], intensity[:maxdex]))
                accu_list[i] = np.concatenate((accu_list[i], adder), axis=1)

    return accu_list


def collect_calibrant_stats(cal_spectra, cal_masses ):
    """collects bulk statistics of the calibrants. return a tuple of 4 lists:
    0) binary mask if an entry for the calibrant was found
    1) the mostabundant peak in the specified bin
    2) the intensity-weighted mz value avergae in the bin
    3) the nearest local maxima to the calibrant mass."""
    cal_mass_mask = []
    accu_mab_list = []
    accu_center_list = []
    accu_locmax_list = []

    # extraction of most Abundant PEaks and peak centers and theeir validity in mask for cal_masses
    for i in range(len(cal_masses)):
        if len(cal_spectra[i][1]) > 0:
            cal_mass_mask.append(True)

            # peak with hightest intensity
            most_abundant_peak = cal_spectra[i][0][np.where(cal_spectra[i][1] == max(cal_spectra[i][1]))][0]
            # weighted average of mz values weighted by their intensity
            center, _ = np.average(cal_spectra[i][0], weight=cal_spectra[i][1])
            # nearest local maxima of intensities
            loc_max_ind = argrelextrema(cal_spectra[i][1], np.greater)
            loc_max = find_nearest(cal_spectra[i][0][loc_max_ind][0])

            accu_mab_list.append(most_abundant_peak)
            accu_center_list.append(center)
            accu_locmax_list.append(loc_max)
        else:
            cal_mass_mask.append(False)
            accu_mab_list.append(0)  # mab is set to 0 for not found peaks to retain list shape
            accu_center_list.append(0)  # center is set to 0 for not found peaks
            accu_locmax_list.append(0)



# M2aia-dependant tools

#2DO: maybe implement a coverage here aswell
def collect_image_stats(Image # m2aia-ImzML Reader object
                        ):
    """ Expensive function to call. iterates over full spectrum and returns the following metrics:
    0) the index of each pixel, useful for plotting
    1) the number of data points present in each pixel
    2) the TIC of each pixel
    3) the median intensity of each pixel
    4) the maximum intensity that was recorded in the pixel
    5) the minimum intesity that was recorded in the pixel
    6) the highest mz value present
    7) the smallest mz value present
    8) the most abundant mz value in each pixel (Base Peak)"""
    index_nr = tuple()  # container for pixel index, corrected for 0-index
    peak_nr = tuple()  # container for number of features loaded in a single pixel
    tic_nr = tuple()  # container of TIC for each single pixel
    median_nr = tuple()  # Container for the median intnsity in each pixel
    max_int_nr = tuple()  # Container for the max intnsity in each pixel
    min_int_nr = tuple()  # Container for the min intnsity in each pixel
    max_mz_nr = tuple()  # Container for the maximal mz value in each pixel
    min_mz_nr = tuple()  # Container for the minimal mz  value in each pixel
    max_abun_nr = tuple()  # Container for the most abundant mz in each pixel

    for ind, mass, inten in Image.SpectrumIterator():  # loop to run over full imzML dataset
        # values for a single pixel are recorded in their specific tuples
        index_nr = index_nr + (ind + 1,)  # pixel order is 0 for not-recorded pixels
        peak_nr = peak_nr + (len(inten),)
        tic_nr = tic_nr + (sum(inten),)
        median_nr += (stat.median(inten),)
        max_int_nr += (max(inten),)
        min_int_nr += (min(inten),)
        max_mz_nr += (max(mass),)
        min_mz_nr += (min(mass),)
        max_abun_nr += (mass[np.where(inten == max(inten))[0][0]],)

    return (index_nr, peak_nr, tic_nr, median_nr, max_int_nr, min_int_nr, max_mz_nr, min_mz_nr, max_abun_nr)



def generate_table_data(Image, x_limits, y_limits, im_stats):
        table = [
        ["Spectral Type:", str(Image.GetSpectrumType())],
        ["Numeric shape of Image (x, y, z):", str(Image.GetShape())],
        ["Number of recorded Pixels:", str(Image.GetNumberOfSpectra())],
        ["Number of unrecorded Pixels:", str(np.abs(np.product(Image.GetShape()) - Image.GetNumberOfSpectra()))],
        ["Recorded x-range:", str(x_limits)],
        ["Recorded y-range:", str(y_limits)],
        ["Number of individual mz features:", str(np.sum(im_stats[1]))],
        ["Mean TIC ± sd:", str(f"{int(stat.mean(im_stats[2]))} ± {int(stat.stdev(im_stats[2]))}")],
        ["Median TIC ± MAD:", str(f"{int(stat.median(im_stats[2]))} ± {int(SST.median_abs_deviation(im_stats[2]))}")],
        ["Mean number of mz features per spectrum ± sd:", str(f"{int(stat.mean(im_stats[1]))} ± {int(stat.stdev(im_stats[1]))}")],
        ["Median number of mz features per spectrum ± MAD:", str(f"{int(stat.median(im_stats[1]))} ± {int(SST.median_abs_deviation(im_stats[1]))}")],
        ["Range of median intensities per pixel:", str((min(im_stats[3]), max(im_stats[3])))],
        ["Range of Maximal Intensity per pixel:", str((min(im_stats[4]), max(im_stats[4])))],
        ["Range of most abundant mz per pixel:", str((min(im_stats[8]), max(im_stats[8])))],
        ["mz range:", str((min(im_stats[7]), max(im_stats[6])))],
        ["Spacing:", str(Image.GetSpacing())],
        ["m/z Bins:", str(Image.GetXAxisDepth())],
        ["Intensity range:", str((min(im_stats[5]), max(im_stats[4])))],
        ]
        return table