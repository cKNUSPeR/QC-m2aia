
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


