
# namespace and import declaration
import matplotlib.backends as mpb
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








def calc_accuraciues(found_mz, theo_mz, mask):
    """Calculate the ppm accuracy of a list of found peaks comparative to a set of theroetical masses and a binary mask"""
    ppm = []
    for o, t, m in zip(found_mz, theo_mz, mask):
        if m:
            ppm.append(np.abs((o - t) / t * 10 ** 6))
    return ppm


# M2aia-dependant tools


def label_connected_region(Image):
    labeled_image = skim.label(Image.GetMaskArray()[0], connectivity=1)

    # shape of image array:
    rows, cols = labeled_image.shape
    # get a meshed grid to make x-y accesible
    x_coords, y_coords = np.meshgrid(range(cols), range(rows))

    # make a dataframe
    df = pd.DataFrame({'x': x_coords.flatten(), 'y': y_coords.flatten(), 'annotation_value': labeled_image.flatten()})
    # remove 0-entries (they represent empty pixels)
    df = df.loc[df["annotation_value"] > 0]

    # number of regions found
    max_regions = df['annotation_value'].max()

    return df, labeled_image, max_regions


def parse_regionfile(file, annotation_group, x_lims, y_lims):
    """parses the regions tsv file. Assumes that the regions are annotated
    within the coordinates of the imaging file (the point of origin is at (0,0))."""

    # mark invariants of "x", "y" , "annotation"
    df = pd.read_csv(file, sep="\t", header=0)

    # translate the annotation to values starting from 1
    df['annotation_values'] = pd.factorize(df[annotation_group])[0] + 1

    # the maximum number of regions (counting starts at 1)
    max_regions = df['annotation_values'].max()

    # set up empty array with xlims any ylims ranges
    labeled_image = np.zeros((1 + y_lims[1] - y_lims[0], 1 + x_lims[1] - x_lims[0]))  # 1-indexed for inclusiveness

    # fill the labeled image with the annotations:
    for index, row in df.iterrows():
        x, y = row['x'], row['y']
        labeled_image[y, x] = row['annotation_values']  # Set pixel to annotation group

    return df, labeled_image, max_regions


def write_region_tsv(df, path):
    """writes a region pd.df to a tsv for reimport later."""
    file_name = path + "annotated_regions.tsv"
    df.to_csv(file_name, sep="\t", columns=["x", "y", "annotation_value"], index=False)


def group_region_stat(labeled_image, index_image, label_nr, image_stats, keyword):
    """Groups the statistics of a region into a list of lists.
    Input:
    - labeled_image: An image-like array containing the regions labeled with a non-zero int
    - index_image: An image-like array containing the pixel index
    - image_stats: a dict of image statistics per pixel as tuple, accessible by keywords
    -keyword: the keyword for which data is collected (eg. tic_nr)
    - label_nr: the number of labeled regions

    Output:
    - tuple(list[int], list[list]: Tuple with names and statistics:
    ordered in list order of keywords (currently hardcoded)

    """

    # linear reshaping of pixel index image and segmented image
    lab_ar = np.reshape(labeled_image, -1)
    ind_ar = np.reshape(index_image, -1)

    # arrayization of pixel index counting and TICperPixel counting
    # via np.asarray(image_stats["index_nr"]), best to do inplace

    # collectors for plotable boxplot data
    # collectors for plotable boxplot data
    stat_coll_boxplot = []
    name_coll_boxplot = []

    # loop over all segments
    for seg in range(1, label_nr + 1):
        stat_nr_arr = np.asarray(image_stats[keyword])
        ind_nr_array = np.asarray(image_stats["index_nr"])

        pindex = ind_ar[np.where(lab_ar == seg)]  # extracion of pixel indices per segment

        col = stat_nr_arr[np.isin(ind_nr_array, pindex)]  # extraction of tics from pixel index
        col = np.log2(col)

        stat_coll_boxplot.append(col)
        name_coll_boxplot.append(seg)

    return name_coll_boxplot, stat_coll_boxplot


def average_cont_spectra(Image, pixels):
    """
    Input:
    - sequence of pixel indices
    :return:
    array of mz, array of intensities
    """
    # get lngth
    n = len(pixels)

    # get first element
    mz, ints = Image.GetSpectrum(pixels[0])
    ints = ints / n

    # iterate over remaining elements
    for idx in pixels[1:]:
        _, intensity = Image.GetSpectrum(idx)
        intensity = intensity / n
        ints = np.add(ints, intensity)

    return mz, ints


def average_processed_spectra(Image, pixels):
    """Averages processed spectra by their pixels.
    Calculates a mz window with start, end and stepsize and bins the data into this.
    Masses are returned to the head of their respective bin.
    Input:
        - Image: an m2aia imzML reader
        - pixels: a sequence of pixel indices
    returns:
        mzs, ints
        two arrays with mz values and intensity values of the respective average
    """
    # get the normalization numbers
    n = len(pixels)

    # get the mz value range
    mz_start, mz_end = min(Image.GetXAxis()) - 0.0001, max(Image.GetXAxis())  # minimum is reduced by a bit.

    # get the number of bins (estimation by getting the largest number of spectra in the file)
    bin_nr = 1  # the left index of each mass bin
    for ids in pixels:
        cbin_nr = Image.GetSpectrumDepth(ids)
        if cbin_nr > bin_nr:
            bin_nr = cbin_nr

    # binning is enhanced by factor of 10, to counteract spectral pixelation
    bin_nr = bin_nr * 10

    # set up the collection df for the binned ranges:
    bins = np.linspace(mz_start, mz_end, bin_nr)

    full_df = pd.DataFrame({'mz_bins': bins})
    full_df["collect"] = np.nan

    # make a loop:
    for idx in pixels:
        # mzs and ints of the pixel
        mz, ints = Image.GetSpectrum(idx)

        # make andbin the dat in a dataframe
        little_df = pd.DataFrame({'mz': mz, "intensity": ints})
        little_df['binned'] = pd.cut(mz, bins)

        # means of grouped bins and normalized to pixels
        full_df["single"] = little_df.groupby(['binned'])['intensity'].mean() / n
        full_df["collect"] = full_df[['collect', 'single']].sum(axis=1, min_count=1)

    # filter out NaN values:
    full_df = full_df.dropna(subset=['collect'])

    # return as arrays
    return full_df['mz_bins'].to_numpy(), full_df["collect"].to_numpy()


def collect_image_stats(Image, statistic_keywords):
    """ Expensive function to call. iterates over the full spectrum and returns the specified metrics.
    Input:
        -Image: a m2aia ImzML reader
        -statistic_keywords: a controlled set of keyword strings
            The following keywords are supported:
            ['index_nr', 'peak_nr', 'tic_nr', 'median_nr', 'max_int_nr', 'min_int_nr', 'max_mz_nr', 'min_mz_nr', 'max_abun_nr']

    Output:
        -dict[keyword] -> tuple: a dict with tuples of the required statistic per pixel.
    """
    # Create a dictionary to store the results for each statistic
    statistics_result = {keyword: () for keyword in statistic_keywords}

    for ind, mass, inten in Image.SpectrumIterator():  # loop to run over the full imzML dataset
        for keyword in statistic_keywords:
            if keyword == 'index_nr':
                statistics_result[keyword] += (ind + 1,)
            elif keyword == 'peak_nr':
                statistics_result[keyword] += (len(inten),)
            elif keyword == 'tic_nr':
                statistics_result[keyword] += (sum(inten),)
            elif keyword == 'median_nr':
                statistics_result[keyword] += (stat.median(inten),)
            elif keyword == 'max_int_nr':
                statistics_result[keyword] += (max(inten),)
            elif keyword == 'min_int_nr':
                statistics_result[keyword] += (min(inten),)
            elif keyword == 'max_mz_nr':
                statistics_result[keyword] += (max(mass),)
            elif keyword == 'min_mz_nr':
                statistics_result[keyword] += (min(mass),)
            elif keyword == 'max_abun_nr':
                max_abun_index = np.where(inten == max(inten))[0][0]
                statistics_result[keyword] += (mass[max_abun_index],)

    return statistics_result


#decrepateted
def old_collect_image_stats(Image # m2aia-ImzML Reader object
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
        ["Number of individual mz features:", str(np.sum(im_stats['peak_nr']))],
        ["Mean TIC ± sd:", str(f"{int(stat.mean(im_stats['tic_nr']))} ± {int(stat.stdev(im_stats['tic_nr']))}")],
        ["Median TIC ± MAD:", str(f"{int(stat.median(im_stats['tic_nr']))} ± {int(SST.median_abs_deviation(im_stats['tic_nr']))}")],
        ["Mean number of mz features per spectrum ± sd:", str(f"{int(stat.mean(im_stats['peak_nr']))} ± {int(stat.stdev(im_stats['peak_nr']))}")],
        ["Median number of mz features per spectrum ± MAD:", str(f"{int(stat.median(im_stats['peak_nr']))} ± {int(SST.median_abs_deviation(im_stats['peak_nr']))}")],
        ["Range of median intensities per pixel:", str((min(im_stats['median_nr']), max(im_stats['median_nr'])))],
        ["Range of Maximal Intensity per pixel:", str((min(im_stats['max_int_nr']), max(im_stats['max_int_nr'])))],
        ["Range of most abundant mz per pixel:", str((min(im_stats['max_abun_nr']), max(im_stats['max_abun_nr'])))],
        ["mz range:", str((min(im_stats['min_mz_nr']), max(im_stats['max_mz_nr'])))],
        ["Spacing:", str(Image.GetSpacing())],
        ["m/z Bins:", str(Image.GetXAxisDepth())],
        ["Intensity range:", str((min(im_stats['min_int_nr']), max(im_stats['max_int_nr'])))],
        ]
        return table