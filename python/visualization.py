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

# custom colormaps with white backgrounds (via out-of-lower-bound)
my_vir = cm.get_cmap('viridis').copy()
my_vir.set_under('white')  # Color for values less than vmin

my_rbw = cm.get_cmap('gist_rainbow').copy()
my_rbw.set_under('white')  # Color for values less than vmin

my_coolwarm = cm.get_cmap('coolwarm').copy()
my_coolwarm.set_under('white')  # Color for values less than vmin
my_coolwarm.set_over('darkorange')

my_cw = cm.get_cmap('coolwarm').copy()
my_cw.set_under('purple')  # Color for values less than vmin
my_cw.set_over('darkorange')
my_cw.set_bad(color='white', alpha=1.0)


def image_full_binary(Image, pdf):
    """plots a binary image of the imaging run with the origin coordinates
        Saves this plot to a pdf."""

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_title('Full view of binary image from origin')
    ax.imshow(Image,
              cmap=my_vir, vmin=0.1,
              interpolation='none')  # attention, large images tend to get a smoohing under the hood by plt
    pdf.savefig(fig)
    plt.close()

def image_cropped_binary(Image, pdf, x_limits, y_limits):
    """generates a plot of binary image cropped to size.
        Saves the plot to a pdf"""

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_xlim(x_limits[0],x_limits[1])
    ax.set_ylim(y_limits[0],y_limits[1])
    ax.set_title('Cropped view of binary image within pixel limits')
    ax.imshow(Image,
              cmap=my_vir, vmin=0.1,
              interpolation='none')
    pdf.savefig(fig)
    plt.close()


def image_pixel_index(Image, pdf, x_limits, y_limits):
    """generates a plot of the index of each pixel. Image cropped to size.
        Saves the plot to a pdf"""

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_title('Pixel Index')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_xlim(x_limits[0], x_limits[1])
    ax.set_ylim(y_limits[0], y_limits[1])

    im = ax.imshow(Image,
                   cmap=my_rbw, vmin=0.1, interpolation='none')
    fig.colorbar(im, extend='min')

    pdf.savefig(fig)
    plt.close()

def plot_basic_scatter(x,y,
                       title, x_lab, y_lab,
                       pdf):
    """makes a simple scatterplot, functional template"""
    fig = plt.figure(figsize=[7,5])
    ax = plt.subplot(111)

    ax.set_title(title)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.grid(visible=True, c='lightgray', ls="--")

    ax.scatter(x, y, color='k', marker=".", zorder=-1)
    ax.set_rasterization_zorder(0)

    pdf.savefig(fig)
    plt.close()

def image_basic_heatmap(Image,
                        title, x_lab, y_lab,
                        pdf, x_limits, y_limits):
    """makes basic heatmap, intended as functional template"""
    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_title(title)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_xlim(x_limits[0], x_limits[1])
    ax.set_ylim(y_limits[0], y_limits[1])

    im = ax.imshow(Image, cmap=my_vir, vmin=0.1)
    fig.colorbar(im, ax=ax, extend='min')

    pdf.savefig(fig)
    plt.close()


def plot_feature_number(indices, values, pdf):
    """plot a scatterplot for the number of feautes per pixel"""
    plot_basic_scatter(indices, values,
                       "Number of Peaks per spectrum",
                       "Index of Spectrum",
                       "Number of Peaks",
                       pdf)

def image_feature_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the number of features. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'Number of Peak Projection',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_tic_number(indices, values, pdf):
    """plot a scatterplot for the Total Ion Count per pixel"""
    plot_basic_scatter(indices, values,
                       "TIC per spectrum",
                       "Index of Spectrum",
                       "Intensity",
                       pdf)

def image_tic_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the TIC. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'TIC per pixel projection',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_max_abun_number(indices, values, pdf):
    """plot a scatterplot for the Highest abundance mz value per pixel"""
    plot_basic_scatter(indices, values,
                       "Highest abundance mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_max_abun_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the Highest abundance mz  value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'Highest abundance mz value per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)

def plot_median_number(indices, values, pdf):
    """plot a scatterplot for the median intensity per pixel"""
    plot_basic_scatter(indices, values,
                       "median intensity per spectrum",
                       "Index of Spectrum",
                       "Intensity",
                       pdf)


def image_median_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the median intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'Median Intensity per Spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_max_int_number(indices, values, pdf):
    """plot a scatterplot for the maximal intensity per pixel"""
    plot_basic_scatter(indices, values,
                       "maximum intensity per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)

def image_max_int_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the maximal intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'maximum intensity per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)

def plot_min_int_number(indices, values, pdf):
    """plot a scatterplot for the minimal intensity per pixel"""
    plot_basic_scatter(indices, values,
                       "minimal intensity per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_min_int_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the minimal intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'minimal intensity per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)#


def plot_max_mz_number(indices, values, pdf):
    """plot a scatterplot for the largest mz value per pixel"""
    plot_basic_scatter(indices, values,
                       "largest mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_max_mz_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the largest mz value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'largest mz value per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_min_mz_number(indices, values, pdf):
    """plot a scatterplot for the smallest mz value per pixel"""
    plot_basic_scatter(indices, values,
                       "smallest mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_min_mz_number(Image, pdf, x_limits, y_limits):
    """Images a heatmap of the smallest mz value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(Image,
                        'smallest mz value per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)

def plot_centroid_spectrum(mz_axis, spectrum_data, pdf):
    fig = plt.figure(figsize=[10, 6])
    ax = plt.subplot(111)

    ax.set_title('Averaged Centroid Mass Spectrum')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_xlim(min(mz_axis).round(0), max(mz_axis).round(0))

    ax.vlines(mz_axis, 0, spectrum_data, linewidth=0.8)

    pdf.savefig(fig)
    plt.close()


def plot_profile_spectrum(mz_axis, spectrum_data, pdf):
    fig = plt.figure(figsize=[10, 6])
    ax = plt.subplot(111)

    ax.set_title('Averaged Centroid Mass Spectrum')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_xlim(min(mz_axis).round(0), max(mz_axis).round(0))

    ax.plot(mz_axis, spectrum_data, linewidth=0.8)

    pdf.savefig(fig)
    plt.close()


def write_summary_table(table, pdf):
    # Create a figure and add the table
    fig = plt.figure(figsize=[10, 10])
    ax = plt.subplot(111)
    ax.axis("off")  # Turn off axis
    table = ax.table(cellText= table,
                     colLabels=["Property", "Values"],
                     loc="center", cellLoc="left")

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    table.scale(1.2, 1.2)  # Adjust table scale for better layout
    # weird error, where some text is not getting passed

    pdf.savefig(fig, bbox_inches="tight")
    plt.close()