import pandas as pd

import numpy as np

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

from region_utils import average_processed_spectra, average_cont_spectra

file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\Example_Continuous_imzML\Example_Continuous.imzML"
I = m2.ImzMLReader(file_name)

mzs, ints = average_processed_spectra(I, [2,3,4])

print(mzs)
print(ints)