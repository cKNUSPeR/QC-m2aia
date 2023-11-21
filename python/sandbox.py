import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import nectar_msi as nc

import m2aia as m2
def bgin_of_bin():
    n = 2

    # bins must start lower than actual values
    bins = list(range(0,12,2))


    full_df = pd.DataFrame({'mz_bins': bins})
    full_df["collect"] = np.nan



    # loops
    mz, ints = np.array([1, 7, 5, 4, 6, 3]), np.array([2, 4, 6, 8, 10, 4])
    little_df = pd.DataFrame({'mz': mz, "intensity": ints})
    little_df['binned'] = pd.cut(mz, bins)

    # means of grouped bins and normalized to pixels
    full_df["single"] = little_df.groupby(['binned'])['intensity'].mean()/n
    full_df["collect"] = full_df[['collect', 'single']].sum(axis=1, min_count=1)

    # loop finished
    full_df = full_df.dropna(subset=['collect'])


    print(full_df['mz_bins'].to_numpy(), full_df["collect"].to_numpy())

#from region_utils import average_processed_spectra, average_cont_spectra

from cal_utils import find_nearest_loc_max



file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\Example_Continuous_imzML\Example_Continuous.imzML"
I = m2.ImzMLReader(file_name)

mzs, ints = I.GetSpectrum(6)


fig = plt.figure(figsize=[10, 6])
ax = plt.subplot(111)

ax.set_title('Averaged Profile Mass Spectrum')
ax.set_xlabel('m/z')
ax.set_ylabel('Intensity')


ax.plot(mzs, ints, linewidth=0.8)

plt.show()

signal, intensity = mzs, ints

# Find indices of the top 10 maximum values in the "signal" array
top_indices = np.argpartition(intensity, -10)[-10:]

# Select corresponding values from the "intensity" array
selected_intensity = intensity[top_indices]

# Display the results
print("Top 10 signal values:", signal[top_indices])
print("Corresponding intensity values:", selected_intensity)