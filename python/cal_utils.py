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

def find_nearest_loc_max(array, value, lindex, hindex):
    # ensure arrayztion for mapping
    array = np.asarray(array)
    # slice array for shorter calc time
    array = array[lindex:hindex]
    # get local maxima indices
    max_index = SSI.argrelmax(array, np.greater)
    # get values of those local maxima
    array = array[max_index]

    idx = (np.abs(array - value)).argmin()
    return array[idx]

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

def calculate_ppm(exp_mass: float,
                  theo_mass: float) -> float:
    """Calulates ppm of na experimental mass againt a theoretical mass.
    Input:
        - exp_mass: observed mz as float
        - theo_mass: theoretical mz value as flaot
    :return
        - ppm value: float
    """
    return ((exp_mass - theo_mass) / theo_mass)*1e6


def collect_accuracy_stats(Image, calibrants_df, dist, format_dict):

     # make a matrix for each pixel and
     accuracies_ar = np.zeros((Image.GetNumberOfSpectra(), len(calibrants_df["name"])))
     if format_dict["centroid"]:
         for ind, mass, inten in Image.SpectrumIterator():  # loop to run over full imzML dataset
             # get nearest elements TODO: change to nearest local maxima
             tm = [find_nearest(mass, calmass) for calmass in calibrants_df["mz"]]
             accuracies_ar[ind] = tm
     elif format_dict["profile"]:
         find_nearest_loc_max()