from utils import *

def read_calibrants(filepath: str):
    """Reads calibrant files and gives a list of name and thr. mz values
    INVARIANTS: Needs a header column with 'name' and a col with 'mz'."""
    cal = pd.read_csv(filepath, sep=';', header=0)
    cal["found"] = np.NaN
    cal["value_wavg"] = np.NaN
    cal["distance_wavg"] = np.NaN
    cal["value_map"] = np.NaN
    cal["distance_map"] = np.NaN
    cal["coverage"] = np.NaN
    return cal


def make_subsample(samplenumber: int, percent_sample:float) -> list:
    """Makes a subsample out of the samplenumber with the given percentage (as float)"""
    # Determine the size of a batch
    batch_size = int(samplenumber * percent_sample)
    # get random numbers according to the batch size
    return rnd.sample(range(0, samplenumber), batch_size)

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def find_nearest_loc_max(mzs, intensites, value, distance=0):
    """Finds the nearest local maxima of a profile line to a certain mass.
    It disqualifies local maxima that are below 1 % of the highest point found inside the specified range.
    If no distance is specified, the search is applied to the whole spectrum.

    """

    if distance == False:
        lindex = min(np.where(mzs > (value - distance))[0])
        hindex = min(np.where(mzs > (value + distance))[0])
    else:
        lindex = 0
        hindex = len(mzs)

    # ensure arrayztion for mapping and slice array for shorter calc time
    intensites = np.asarray(intensites)[lindex:hindex]
    mzs = np.asarray(mzs)[lindex:hindex]

    # get local maxima indices
    max_index = SSI.argrelextrema(intensites, np.greater)
    # and detuple this
    max_index = max_index[0]

    # get all the values from max_index
    locmax_ints = intensites[max_index]
    # and the maximaum intensity found
    max_intensity = max(locmax_ints)
    # get all the intensity indices where the locmax_intensity surpasses 1% of max intensity
    max_index = max_index[np.where(locmax_ints >= max_intensity*0.01)]

    # get values of those local maxima
    mzs = mzs[max_index]
    idx = (np.abs(mzs - value)).argmin()
    return mzs[idx]

def extract_calibrant_spectra(Image, cal_mass, subsample, mz_bin):
    """Read the full image. Collects the spectral data for a given mass in the given mz bin."""
    accu_list = np.array([[],[]])

    # looping over sample
    for ind in subsample:
        mass, intensity = Image.GetSpectrum(ind)
        try:
            #mindex needs to be inclusive
            mindex = min(np.where(mass > (cal_mass-mz_bin))[0])
            mindex_flag = True
        except:
            mindex_flag = False
        try:
            #maxdex needs to be exclusive
            maxdex = min(np.where(mass > (cal_mass+mz_bin))[0])
            maxdex_flag = True
        except:
            maxdex_flag = False

        # pixels are only written if there is data present in the specified part
        if maxdex_flag and mindex_flag:
            # collecting of masses and intensities
            adder = np.array((mass[mindex:maxdex],intensity[mindex:maxdex]))
            accu_list = np.concatenate((accu_list, adder), axis=1)
        elif mindex_flag and not maxdex_flag:
            adder = np.array((mass[mindex:], intensity[mindex:]))
            accu_list = np.concatenate((accu_list, adder), axis=1)
        elif not mindex_flag and maxdex_flag:
            adder = np.array((mass[:maxdex], intensity[:maxdex]))
            accu_list = np.concatenate((accu_list, adder), axis=1)

    return accu_list



def collect_calibrant_stats(cal_spectra, calibrant_df, index):
    """collects bulk statistics of the calibrants. Adds to the provided df the following infos:
    0) df["found"]: Whether spectral data was found for the mass
    1) cal["value_wavg"]: the value of the weighted average
    2) cal["distance_wavg"]: the distance in ppm of weight. avg to the theo. mz
    3) cal["value_map"]: the value of the most abundant peak in interval
    4) cal["distance_map"]: the distance in ppm of m.a.p. to the theo. mz
    # defunc.) the nearest local maxima to the calibrant mass.
    """

    # deep-copy the df (it gets mutated over function call)
    calibrant_df = calibrant_df.copy(deep=True)

    # extraction of most Abundant Peaks and peak centers and their validity
    if len(cal_spectra[1]) > 0:

        # peak with hightest intensity
        most_abundant_peak = cal_spectra[0][np.where(cal_spectra[1] == max(cal_spectra[1]))][0]

        # weighted average of mz values weighted by their intensity
        wavg = np.average(cal_spectra[0], weights=cal_spectra[1])


        # update the dataframe
        calibrant_df.loc[index, "found"] = True
        calibrant_df.loc[index, "value_wavg"] = wavg
        calibrant_df.loc[index, "value_map"] = most_abundant_peak

        # calculate distane ppm
        calibrant_df.loc[index, "distance_wavg"] = calculate_ppm(wavg, calibrant_df.loc[index, "mz"])
        calibrant_df.loc[index, "distance_map"] = calculate_ppm(most_abundant_peak, calibrant_df.loc[index, "mz"])


    else:
        calibrant_df.loc[index, "found"] = False
        # values are not updated, NaN signifies non-found peaks

    return calibrant_df

def calculate_ppm(exp_mass,
                  theo_mass: float):
    """Calulates ppm of na experimental mass againt a theoretical mass.
    Input:
        - exp_mass: observed mz as float
        - theo_mass: theoretical mz value as flaot
    :return
        - ppm value: float
    """
    return ((exp_mass - theo_mass) / theo_mass)*1e6


def collect_accuracy_stats(Image, calibrants_df, dist, format_dict):
    """ Finds and collects the nearest signals around all provided calibrant masses.
    Input:
        - Image: ImzMLReader object
        - calibrants_df: dataframe of calibrant information
        - format_dict: dict of imzML formats to handle signal evaluation

    :returns  accuracies_ar, index_nr
    accuracies_ar: the array for the data of accurasies per pixels, linearized images
    index_nr: tuple of pixel indix sorted to match the shape of accuracies_ar

    """
    # make a matrix for each pixel and
    accuracies_ar = np.zeros((Image.GetNumberOfSpectra(), len(calibrants_df["name"])))
    index_nr = tuple()  # container for pixel index, corrected for 0-index

    if format_dict["centroid"]:
         for ind, mass, inten in Image.SpectrumIterator():  # loop to run over full imzML dataset
             # get nearest elements
             accuracies_ar[ind] =[find_nearest(mass, calmass) for calmass in calibrants_df["mz"]]
             # collect image index in order of iteration
             index_nr = index_nr + (ind + 1,)  # pixel order is 0 for not-recorded pixels


    elif format_dict["profile"]:
        for ind, mass, inten in Image.SpectrumIterator():  # loop to run over full imzML dataset
            # get nearest loc, max
            accuracies_ar[ind] = [find_nearest_loc_max(mass,inten, calmass, dist) for calmass in calibrants_df["mz"]]
            # collect image index in order of iteration
            index_nr = index_nr + (ind + 1,)  # pixel order is 0 for not-recorded pixels

    # transpose to match  ppm calcs form
    accuracies_ar = accuracies_ar.T

    # convert the mass into ppm ranges
    for i, mass in enumerate(calibrants_df["mz"]):
        accuracies_ar[i] = calculate_ppm(accuracies_ar[i], mass)

    return accuracies_ar, index_nr

def collect_calibrant_converage(accuracy_images, calibrants_df, accuracy_cutoff):
    """Evalualtes how many pixels within a accuracy image fall outsde of the defined accuracy cutoff.
    saves these into the calibrant_df as converage, normed on the amount of pixels."""
    # deep-copy the df
    calibrant_df = calibrants_df.copy(deep=True)


    # loop over the calibrants
    for i, mass in enumerate(calibrant_df["mz"]):
        #count how many values are smaller than specified cutoff:
        low_dist = np.sum(accuracy_images[i] < accuracy_cutoff)
        high_dist = np.sum(accuracy_images[i] > accuracy_cutoff)

        # add to calibrants_df
        calibrant_df.loc[i, "coverage"] = (low_dist+high_dist)

    # divide over number of pixels
    pixel_nr = len(accuracy_images[0])
    # normalize over pixel number
    calibrant_df["coverage"] = calibrant_df["coverage"]/pixel_nr

    return calibrant_df
