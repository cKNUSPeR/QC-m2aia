
from dependencies import *
from utils import mask_bad_image, average_cont_spectra, average_processed_spectra

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

# dictionary to keep formatting consistency



def make_pdf_backend(report_path, title):
    pdf_file_path = report_path + title + ".pdf"
    pdf_pages = mpb.backend_pdf.PdfPages(pdf_file_path)
    return pdf_pages


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
              interpolation='none',  # attention, large images tend to get a smoohing under the hood by plt
              origin='lower')
    pdf.savefig(fig)
    plt.close()


def image_cropped_binary(Image, pdf, x_limits, y_limits):
    """generates a plot of binary image cropped to size.
        Saves the plot to a pdf"""

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')
    ax.set_xlim(x_limits[0], x_limits[1])
    ax.set_ylim(y_limits[0], y_limits[1])
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


def image_regions(Image, regionarray, pdf, x_limits, y_limits):
    """Images the annotated regions image as colorful blops-"""
    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_title('Connected Objects Analysis')
    ax.set_xlabel('x axis')
    ax.set_ylabel('y axis')

    im = ax.imshow(regionarray, cmap=my_rbw, vmin=0.1, interpolation='none', origin='lower')
    # extent=[x_limits[0], x_limits[1], y_limits[0], y_limits[1]])

    fig.colorbar(im, extend='min', format=lambda x, _: f"{int(x)}")

    pdf.savefig(fig)
    plt.close()



def plot_basic_scatter(x, y,
                       title, x_lab, y_lab,
                       pdf):
    """makes a simple scatterplot, functional template"""
    fig = plt.figure(figsize=[7, 5])
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


def plot_feature_number(image_stats, pdf):
    """plot a scatterplot for the number of feautes per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["peak_nr"],
                       "Number of Peaks per spectrum",
                       "Index of Spectrum",
                       "Number of Peaks",
                       pdf)


def image_feature_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the number of features. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["peak_nr"], index_image),
                        'Number of Peak Projection',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_tic_number(image_stats, pdf):
    """plot a scatterplot for the Total Ion Count per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["tic_nr"],
                       "TIC per spectrum",
                       "Index of Spectrum",
                       "Intensity",
                       pdf)


def image_tic_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the TIC. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["tic_nr"], index_image),
                        'TIC per pixel projection',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_max_abun_number(image_stats, pdf):
    """plot a scatterplot for the Highest abundance mz value per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["max_abun_nr"],
                       "Highest abundance mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_max_abun_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the Highest abundance mz  value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["max_abun_nr"], index_image),
                        'Highest abundance mz value per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_median_number(image_stats, pdf):
    """plot a scatterplot for the median intensity per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["median_nr"],
                       "median intensity per spectrum",
                       "Index of Spectrum",
                       "Intensity",
                       pdf)


def image_median_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the median intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["median_nr"], index_image),
                        'Median Intensity per Spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_max_int_number(image_stats, pdf):
    """plot a scatterplot for the maximal intensity per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["max_int_nr"],
                       "maximum intensity per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_max_int_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the maximal intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["max_int_nr"], index_image),
                        'maximum intensity per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_min_int_number(image_stats, pdf):
    """plot a scatterplot for the minimal intensity per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["min_int_nr"],
                       "minimal intensity per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_min_int_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the minimal intensity. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["min_int_nr"], index_image),
                        'minimal intensity per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)  #


def plot_max_mz_number(image_stats, pdf):
    """plot a scatterplot for the largest mz value per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["max_mz_nr"],
                       "largest mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_max_mz_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the largest mz value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["max_mz_nr"], index_image),
                        'largest mz value per spectrum',
                        "x axis",
                        "y axis",
                        pdf, x_limits, y_limits)


def plot_min_mz_number(image_stats, pdf):
    """plot a scatterplot for the smallest mz value per pixel"""
    plot_basic_scatter(image_stats["index_nr"], image_stats["min_mz_nr"],
                       "smallest mz value per spectrum",
                       "Index of spectrum",
                       "Intensity",
                       pdf)


def image_min_mz_number(image_stats, index_image, pdf, x_limits, y_limits):
    """Images a heatmap of the smallest mz value. Image cropped to size.
        Saves the plot to a pdf"""
    image_basic_heatmap(mask_bad_image(image_stats["index_nr"], image_stats["min_mz_nr"], index_image),
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
    ax.set_ylim(bottom=0)

    pdf.savefig(fig)
    plt.close()


def plot_profile_spectrum(mz_axis, spectrum_data, pdf):
    fig = plt.figure(figsize=[10, 6])
    ax = plt.subplot(111)

    ax.set_title('Averaged Profile Mass Spectrum')
    ax.set_xlabel('m/z')
    ax.set_ylabel('Intensity')
    ax.set_xlim(min(mz_axis).round(0), max(mz_axis).round(0))

    ax.plot(mz_axis, spectrum_data, linewidth=0.8)
    ax.set_ylim(bottom=0)

    pdf.savefig(fig)
    plt.close()


def write_summary_table(table, pdf):
    # Create a figure and add the table
    fig = plt.figure(figsize=[10, 10])
    ax = plt.subplot(111)
    ax.axis("off")  # Turn off axis
    table = ax.table(cellText=table,
                     colLabels=["Property", "Values"],
                     loc="center", cellLoc="left")

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    table.scale(1.2, 1.2)  # Adjust table scale for better layout
    # weird error, where some text is not getting passed

    pdf.savefig(fig, bbox_inches="tight")
    plt.close()


def write_calibrant_summary_table(data_frame, pdf):
    # Create a figure and add the table
    fig = plt.figure(figsize=[10, 10])
    ax = plt.subplot(111)
    ax.axis("off")  # Turn off axis
    table = ax.table(cellText=data_frame.to_numpy(),
                     colLabels=data_frame.columns,
                     loc="center", cellLoc="center")

    # Style the table
    table.auto_set_font_size(False)
    table.set_fontsize(14)
    table.scale(1.2, 1.2)  # Adjust table scale for better layout
    # weird error, where some text is not getting passed

    pdf.savefig(fig, bbox_inches="tight")
    plt.close()

def plot_boxplots(name_boxplot, stat_boxplot, pdf):
    # 2DO: scaling adjusted to 20, also parametrized with titles, and mabe make a subfunction for plotting
    len_b20 = len(name_boxplot) // 20
    if (len(name_boxplot) % 20) > 0:
        len_b20 = len_b20 + 1

    # plotting functions based on single-line or multi-line plotting:
    if len_b20 > 1:
        fig, ax = plt.subplots(len_b20, figsize=(10, len_b20 * 4))
        fig.suptitle('Boxplots of Pixelwise TIC per Segment')

        for j in range(1, len_b20 + 1):  # change to 1-base index
            ax[j - 1].boxplot(stat_boxplot[(j - 1) * 20:20 * j],
                              labels=name_boxplot[(j - 1) * 20:20 * j])
            ax[j - 1].set_xlabel('Segmented Group')
            ax[j - 1].set_ylabel('log10 of Pixel TIC')

    else:
        fig = plt.figure(figsize=[10, len_b20 * 4])
        ax = plt.subplot(111)
        ax.set_title('Boxplots of Pixelwise TIC per Segment')

        ax.boxplot(stat_boxplot[:],
                   labels=name_boxplot)
        ax.set_xlabel('Segmented Group')
        ax.set_ylabel('log10 of Pixel TIC')

    plt.tight_layout()
    pdf.savefig(fig)
    plt.close()


def plot_regions_average(Image, format_dict, regions_image, region_number, pdf):
    """plot the average spectrum of each region of the regioned image as a spectrum plot.
    
    Input: 
        - image
        - format_flag
        - ragions_image
        
    Output:
        plot of mean in region, adapted to format flag. 
        Also, additional plotting of full mean spectrum in background (for later)
    """

    lab_ar = np.reshape(regions_image, -1)
    ind_ar = np.reshape(Image.GetIndexArray()[0], -1)

    for index in range(1, region_number + 1):
        # get the index per segment
        pindex = ind_ar[np.where(lab_ar == index)]  # extracion of pixel indices per segment

        # make averages
        if format_dict["continuous"]:
            avg_mz, avg_ints = average_cont_spectra(Image, pindex)

            if format_dict["centroid"]:
                plot_centroid_spectrum(avg_mz, avg_ints, pdf)
            elif format_dict["profile"]:
                plot_profile_spectrum(avg_mz, avg_ints, pdf)

        elif format_dict["processed"]:
            avg_mz, avg_ints = average_processed_spectra(Image, pindex)

            if format_dict["centroid"]:
                plot_centroid_spectrum(avg_mz, avg_ints, pdf)
            elif format_dict["profile"]:
                plot_profile_spectrum(avg_mz, avg_ints, pdf)


# plot functions for calibrant QC


def plot_calibrant_spectra(cal_spectra, calibrant_df, index, format_dict, dist, pdf):
    # differentiante the plotting :
    # 1) with profile or centriod  map&wavg
    # 2) only data points + map&wavg
    # 2.2) zoom of 150% around both metrics with only data points
    # 3) zoom on minimal and maximal data points ()

    if calibrant_df.loc[index, "found"]:
        if format_dict["centroid"]:
            # plot centr_ calibrant
            plot_calibrant_centroid_spectra(cal_spectra, calibrant_df, index, dist, pdf)

        elif format_dict["profile"]:
            plot_calibrant_profile_spectra(cal_spectra, calibrant_df, index, dist, pdf)
    else:
        # plot an empty box
        plot_empty_peak(calibrant_df.loc[index, "mz"], calibrant_df.loc[index, "name"], pdf)


def plot_calibrant_centroid_spectra(cal_spectra,
                                    calibrants_df, index,
                                    dist, pdf):
    """ Cal spectrum is the sliced variable of cal_spectra[i]
        # differentiante the plotting :
    # 1) with profile or centriod  map&wavg
    # 2) only data points + map&wavg
    # 2.2) zoom of 150% around both metrics with only data points
    # 3) zoom on minimal and maximal data points ()"""

    name = calibrants_df.loc[index, "name"]
    mass = calibrants_df.loc[index, "mz"]
    mapeak = calibrants_df.loc[index, "value_map"]
    wavg = calibrants_df.loc[index, "value_wavg"]

    fig, axs = plt.subplots(2, 2, figsize=[10, 10])

    # plot of full data as centroid spectrum --------------------------------------------------------------
    axs[0, 0].set_title(f'centroid spectrum of {name}\n({mass})')
    axs[0, 0].set_xlabel('m/z')
    axs[0, 0].set_ylabel('Intensity')
    # set axis limits and style
    axs[0, 0].set_xlim(mass - dist, mass + dist)
    axs[0, 0].ticklabel_format(useOffset=False, )
    axs[0, 0].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[0, 0])

    # plot centroid spectrum
    axs[0, 0].vlines(cal_spectra[0], 0, cal_spectra[1], color='b', linewidth=0.8, zorder=-1)
    axs[0, 0].scatter(cal_spectra[0], cal_spectra[1], s=4, color='b', marker=".", zorder=-1)
    # adjust y limits
    axs[0, 0].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[0, 0].set_rasterization_zorder(0)

    # plot full spectra with only data points--------------------------------------------------------------
    axs[0, 1].set_title(f'spectrum of {name}\n({mass}), only data points')
    axs[0, 1].set_xlabel('m/z')
    axs[0, 1].set_ylabel('Intensity')
    # set the axis range and styles
    axs[0, 1].set_xlim(mass - dist, mass + dist)
    axs[0, 1].ticklabel_format(useOffset=False, )
    axs[0, 1].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[0, 1])

    # scatter centroid spectrum
    axs[0, 1].scatter(cal_spectra[0], cal_spectra[1], color='k', marker="x", zorder=-1)
    # adjust y limits
    axs[0, 1].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[0, 1].set_rasterization_zorder(0)

    # plot the zoom to minimal and maximal data points --------------------------------------------------------------
    axs[1, 0].set_title(f'centroid spectrum of {name}\n({mass}), zoomed to values')
    axs[1, 0].set_xlabel('m/z')
    axs[1, 0].set_ylabel('Intensity')
    # set the axis range and styles
    axs[1, 0].set_xlim(min(cal_spectra[0]), max(cal_spectra[0]))
    axs[1, 0].ticklabel_format(useOffset=False, )
    axs[1, 0].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[1, 0])

    # plot centroid spectrum
    axs[1, 0].vlines(cal_spectra[0], 0, cal_spectra[1], linewidth=0.8, zorder=-1)
    axs[1, 0].scatter(cal_spectra[0], cal_spectra[1], s=4, color='b', marker=".", zorder=-1)
    # adjust yaxis bottom
    axs[1, 0].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[1, 0].set_rasterization_zorder(0)

    # plot zoom with all metrics  -------------------------------------------------------------------------
    # get closest metric
    metrics = [calibrants_df.loc[index, "distance_map"], calibrants_df.loc[index, "distance_wavg"]]
    # get the farthest bulk metric
    closest = max(metrics, key=abs)

    # get the interval width (overscaled to 150%)
    interval = abs((mass * (closest * 1e-6 + 1) - mass)*1.5)

    axs[1, 1].set_title(f'spectrum of {name}\n({mass}), zoomed to metrics')
    axs[1, 1].set_xlabel('m/z')
    axs[1, 1].set_ylabel('Intensity')
    # set the axis range and styles
    axs[1, 1].set_xlim(mass - interval, mass + interval)
    axs[1, 1].ticklabel_format(useOffset=False, )
    axs[1, 1].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[1, 1])

    # scatter centroid spectrum
    axs[1, 1].vlines(cal_spectra[0], 0, cal_spectra[1], linewidth=0.8, zorder=-1)
    axs[1, 1].scatter(cal_spectra[0], cal_spectra[1], s=4, color='b', marker=".", zorder=-1)
    # adjust y limits
    axs[1, 1].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[1, 1].set_rasterization_zorder(0)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close()


def plot_calibrant_profile_spectra(cal_spectra,
                                   calibrants_df, index,
                                   dist, pdf):
    """ Cal spectrum is the sliced variable of cal_spectra[i]
        # differentiante the plotting :
    # 1) with profile or centriod  map&wavg
    # 2) only data points + map&wavg
    # 2.2) zoom of 150% around both metrics with only data points
    # 3) zoom on minimal and maximal data points ()"""

    name = calibrants_df.loc[index, "name"]
    mass = calibrants_df.loc[index, "mz"]
    mapeak = calibrants_df.loc[index, "value_map"]
    wavg = calibrants_df.loc[index, "value_wavg"]

    fig, axs = plt.subplots(2, 2, figsize=[10, 10])

    # plot of full data as centroid spectrum --------------------------------------------------------------
    axs[0, 0].set_title(f'centroid spectrum of {name}\n({mass})')
    axs[0, 0].set_xlabel('m/z')
    axs[0, 0].set_ylabel('Intensity')
    # set axis limits and style
    axs[0, 0].set_xlim(mass - dist, mass + dist)
    axs[0, 0].ticklabel_format(useOffset=False, )
    axs[0, 0].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[0, 0])

    # plot profile spectrum
    axs[0, 0].plot(cal_spectra[0], cal_spectra[1], linewidth=0.5, zorder=-1)
    # adjust y limits
    axs[0, 0].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[0, 0].set_rasterization_zorder(0)


    # plot full spectra with only data points --------------------------------------------------------------
    axs[0, 1].set_title(f'spectrum of {name}\n({mass}), only data points')
    axs[0, 1].set_xlabel('m/z')
    axs[0, 1].set_ylabel('Intensity')
    # set the axis range and styles
    axs[0, 1].set_xlim(mass - dist, mass + dist)
    axs[0, 1].ticklabel_format(useOffset=False, )
    axs[0, 1].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[0, 1])

    # scatter centroid spectrum
    axs[0, 1].plot(cal_spectra[0], cal_spectra[1], linewidth=0.5, zorder=-1)
    # adjust yaxis bottom
    axs[0, 1].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[0, 1].set_rasterization_zorder(0)

    # plot the zoom to minimal and maximal data points --------------------------------------------------------------
    axs[1, 0].set_title(f'centroid spectrum of {name}\n({mass}), zoomed to values')
    axs[1, 0].set_xlabel('m/z')
    axs[1, 0].set_ylabel('Intensity')
    # set the axis range and styles
    axs[1, 0].set_xlim(min(cal_spectra[0]), max(cal_spectra[0]))
    axs[1, 0].ticklabel_format(useOffset=False, )
    axs[1, 0].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[1, 0])

    # plot centroid spectrum
    axs[1, 0].plot(cal_spectra[0], cal_spectra[1], linewidth=0.5, zorder=-1)
    # adjust yaxis bottom
    axs[1, 0].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[1, 0].set_rasterization_zorder(0)

    # plot zoom with all metrics  ------------------------------------------------------------------------------
    # get closest metric
    metrics = [calibrants_df.loc[index, "distance_map"], calibrants_df.loc[index, "distance_wavg"]]
    # get the nearest bulk metric
    closest = max(metrics, key=abs)

    # get the interval width (overscaled to 150%)
    interval = abs((mass * (closest * 1e-6 + 1) - mass) * 1.5)

    axs[1, 1].set_title(f'spectrum of {name}\n({mass}), zoomed to metrics')
    axs[1, 1].set_xlabel('m/z')
    axs[1, 1].set_ylabel('Intensity')
    # set the axis range and styles
    axs[1, 1].set_xlim(mass - interval, mass + interval)
    axs[1, 1].ticklabel_format(useOffset=False, )
    axs[1, 1].ticklabel_format(axis="y", style='sci', scilimits=(0, 0))

    # draw metrics and masses
    draw_vertical_lines(mass, mapeak, wavg, axs[1, 1])

    # scatter centroid spectrum
    axs[1, 1].plot(cal_spectra[0], cal_spectra[1], linewidth=0.5, zorder=-1)
    # adjust yaxis bottom
    axs[1, 1].set_ylim(bottom=0)

    # rasterisazion for better user exerience
    axs[1, 1].set_rasterization_zorder(0)

    fig.tight_layout()
    pdf.savefig(fig)
    plt.close()


def plot_empty_peak(cal_mass, cal_name, pdf):
    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)
    # offset for text annotations
    ax.set_xlim(0,2)
    ax.set_ylim(0,2)
    ax.set_xticks([])
    ax.set_yticks([])

    ax.set_title(f'Spectrum of {cal_mass} ({cal_name})')
    ax.text(1, 1, f'Peak for {cal_mass} m/z \n not found',
            ha='center', fontsize=12)
    pdf.savefig(fig)
    plt.close()


def draw_vertical_lines(mass, mapeak, wavg, axes):
    # make a line of theoretical mass
    axes.axvline(mass, c='r', ls=(0, (1, 3)))
    # make a line for most abundant peak
    axes.axvline(mapeak, color='green', ls="--")
    # make a line for weighted average
    axes.axvline(wavg, c='purple', ls="-.")


def plot_accuracy_barplots(calibrant_df, pdf):
    """plots a barplot for different metrics:
        currently supported: - map,
                            -wavg
                            """

    kw_list = ["distance_map", "distance_wavg"]
    color_list = ['green', "purple"]
    title_list = ["most abundant peak", "weigthed average"]

    for i,key in enumerate(kw_list):
        # drop invalid rows
        df = calibrant_df.copy(deep=True)
        df.dropna(subset=[key])

        #plot the accuracy plots
        plot_accu_barplot(df["name"],df[key],
                          title_list[i], color_list[i],
                          pdf)


def plot_accu_barplot(names, values, metric_name, color, pdf):
    """Makes a bar plot of a given accuracy metric.
    only plots existing values."""

    y_pos = np.arange(len(names))

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_title(f'mass accuracy of calibrants ({metric_name} vs theoretical)')
    ax.set_xlabel('Calibrant')
    ax.set_ylabel('Mass accuracy in ppm')
    ax.set_xticks(y_pos)
    ax.set_xticklabels(names, rotation=45, fontsize=8)

    bars = ax.bar(y_pos, values, color=color)

    # making the bar chart on the data
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}', xy=(bar.get_x() + bar.get_width() / 2, height), xytext=(0, 3),
                    textcoords="offset points", ha='center', va='bottom')

    pdf.savefig(fig)
    plt.close()


def barplot_addlabels(pos,value, axes):
    """writes the rounded value in a barplot above their respecitve position"""
    for i in range(len(pos)):
        axes.text(i, value[i]+0.5, np.round(value[i],5),
                ha='center', fontsize=8)


def plot_accuracy_images(Image, accuracy_images, calibrants_df, index_nr, accuracy_cutoff, x_limits, y_limits, pdf):
    """Makes accuracy heatmaps per pixel ofthe found calibrant accuracy."""
    # loop over the calibrants
    for i, mass in enumerate(calibrants_df["mz"]):
        img = mask_bad_image(index_nr, accuracy_images[i] , Image.GetIndexArray()[0])

        # plot each image
        fig = plt.figure(figsize=[7, 5])
        ax = plt.subplot(111)
        ax.set_title(f'Mass Accuracy of {calibrants_df.loc[i, "name"]}, {calibrants_df.loc[i, "mz"]}')
        ax.set_xlabel('x')
        ax.set_ylabel('y')

        ax.set_xlim(x_limits[0], x_limits[1])
        ax.set_ylim(y_limits[0], y_limits[1])
        im = ax.imshow(img, cmap=my_cw, vmin=-accuracy_cutoff, vmax=accuracy_cutoff)
        fig.colorbar(im, extend='both', label="ppm")

        pdf.savefig(fig)
        plt.close()



def plot_coverage_barplot(names, data, pdf):
    """Makes a bar plot of a given spectral coverage.
    """

    y_pos = np.arange(len(names))

    fig = plt.figure(figsize=[7, 5])
    ax = plt.subplot(111)

    ax.set_title(f'spectral covoverage of mean spectrum)')
    ax.set_xlabel('mz bin')
    ax.set_ylabel('contribution to Total Ion Signal')
    ax.set_xticks(y_pos)
    ax.set_xticklabels(names, rotation=45, fontsize=8)

    bars = ax.bar(y_pos, data, width=0.95, color="blue")

    # making the bar chart on the data
    for bar in bars:
        height = bar.get_height()
        ax.annotate(f'{height:.4f}', xy=(bar.get_x() + bar.get_width() / 2, height), xytext=(0, 3),
                    textcoords="offset points", ha='center', va='bottom')

    pdf.savefig(fig)
    plt.close()