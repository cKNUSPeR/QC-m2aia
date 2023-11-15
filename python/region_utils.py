import matplotlib.pyplot as plt
import skimage.measure as skim
import pandas as pd
import numpy as np


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


def parse_regionfile(file, annotation_group,  x_lims, y_lims):
    """parses the regions tsv file. Assumes that the regions are annotated
    within the coordinates of the imaging file (the point of origin is at (0,0)."""

    # mark invariants of "x", "y" , "annotation"
    df = pd.read_csv(file, sep="\t", header=0)

    # translate the annotation to values starting from 1
    df['annotation_values'] = pd.factorize(df[annotation_group])[0] + 1

    # the maximum number of regions (counting starts at 1)
    max_regions = df['annotation_values'].max()

    # set up empty array with xlims any ylims ranges
    labeled_image = np.zeros((1+y_lims[1] - y_lims[0], 1+x_lims[1]-x_lims[0])) # 1-indexed for inclusiveness

    # fill the labeled image with the annotations:
    for index, row in df.iterrows():
        x, y = row['x'], row['y']
        labeled_image[y, x] = row['annotation_values']  # Set pixel to annotation group

    return df, labeled_image, max_regions


def write_region_tsv(df, path):
    """writes a region pd.df to a tsv for reimport later."""
    file_name = path + "annotated_regions.tsv"
    df.to_csv(file_name, sep="\t", columns=["x","y","annotation_value"], index=False)


def group_region_stat(labeled_image, index_image, label_nr, image_stats, keyword):
    """Groups the statistics of a region into a list of lists.
    Input:
    - labeled_image: An image-like array containing the regions labeled with an non-zero int
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
    ind_ar = np.reshape(index_image,-1)

    # arrayization of pixel index counting and TICperPixel counting
    # via np.asarray(image_stats["index_nr"]), best to do inplace


    # collectors for plotable boxplot data
    # collectors for plotable boxplot data
    stat_coll_boxplot = []
    name_coll_boxplot = []


    # loop over all segments
    for seg in range(1, label_nr+1):

        stat_nr_arr = np.asarray(image_stats[keyword])
        ind_nr_array = np.asarray(image_stats["index_nr"])

        pindex = ind_ar[np.where(lab_ar==seg)] # extracion of pixel indices per segment

        col = stat_nr_arr[np.isin(ind_nr_array, pindex)] # extraction of tics from pixel index
        col = np.log2(col)

        # debug blog
        print(len(stat_nr_arr))
        print(ind_ar)
        print(pindex)
        print(ind_nr_array[0:10])
        print(len(col))

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
    ints = ints/n

    # iterate over remaining elements
    for idx in pixels[1:]:
        _, intensity = Image.GetSpectrum(idx)
        intensity = intensity / n
        ints = np.add(ints,intensity)

    return mz, ints

def average_processed_spectra(Image, pixels):
    """"""
    # get the bin numbers (estimation by getting the largest number of spectra in the file)
    n = len(pixels)
    # get the mz value range
    bins = [] # the left index of each mass bin

    # collecting the binned ranges: as NaN
    collector = pd.df()

    # make a loop:
    for idx in pixels:
        # get the values from one pixel
        mz, intensity = Image.GetSpectrum(idx)
        # divide them over n
        intensity = intensity / n

        # add them as col to df
        df["raw"] = pd.cut(df.data, bins)

        # add the numbers to the collected df, make sure NaN + NaN is Nan
        df["collected"] = df["collected"] + df


    # filter out NaN values:

    # retun df