
from utils import *

def read_calibrants(filepath: str):
    """Reads calibrant files and gives a list of name and thr. mz values
    INVARIANTS: Needs a header column with 'Name' and a col with 'Theoretical m/z'."""
    cal = pd.read_csv(filepath, sep=';', header=0)
    cal_masses = cal["Theoretical m/z"].tolist()
    cal_names = cal["Name"].tolist()
    return cal_names, cal_masses


def make_subsample(samplenumber: int, percent_sample:float) -> list:
    """Makes a subsample out of the samplenumber with the given percentage (as float)"""
    # Determine the size of a batch
    batch_size = int(samplenumber * percent_sample)
    # get random numbers according to the batch size
    return rnd.sample(range(0, samplenumber), batch_size)



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
    decrep.) the nearest local maxima to the calibrant mass."""
    cal_mass_mask = []
    accu_mab_list = []
    accu_center_list = []
    #accu_locmax_list = []

    # extraction of most Abundant PEaks and peak centers and theeir validity in mask for cal_masses
    for i in range(len(cal_masses)):
        if len(cal_spectra[i][1]) > 0:
            cal_mass_mask.append(True)

            # peak with hightest intensity
            most_abundant_peak = cal_spectra[i][0][np.where(cal_spectra[i][1] == max(cal_spectra[i][1]))][0]
            # weighted average of mz values weighted by their intensity
            print(cal_spectra[i][0])
            print(type(cal_spectra[i][1]))
            center, _ = np.average(cal_spectra[i][0], weights=cal_spectra[i][1])
            # nearest local maxima of intensities
            # = argrelextrema(cal_spectra[i][1], np.greater)
            #loc_max = find_nearest(cal_spectra[i][0][loc_max_ind][0])

            accu_mab_list.append(most_abundant_peak)
            accu_center_list.append(center)
            #accu_locmax_list.append(loc_max)
        else:
            cal_mass_mask.append(False)
            accu_mab_list.append(0)  # mab is set to 0 for not found peaks to retain list shape
            accu_center_list.append(0)  # center is set to 0 for not found peaks
            #accu_locmax_list.append(0)
    return cal_mass_mask, accu_mab_list, accu_center_list