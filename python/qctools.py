from dependencies import *
from visualization import *
from utils import *


def report_agnostic_qc(I,  # m2.imzMLReader (passing by ref allows faster computation)
                       outfile_path: str,  # path for output file
                       ):
    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "_agnostic_QC")

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    image_full_binary(I.GetMaskArray()[0],
                      pdf_pages)

    image_cropped_binary(I.GetMaskArray()[0],
                         pdf_pages, x_lims, y_lims)

    image_pixel_index(I.GetIndexArray()[0],
                      pdf_pages, x_lims, y_lims)

    image_stats = collect_image_stats(I,
                                      ['index_nr', 'peak_nr', 'tic_nr', 'median_nr', 'max_int_nr', 'min_int_nr',
                                       'max_mz_nr', 'min_mz_nr', 'max_abun_nr'])

    # visualize the feature numbers
    plot_feature_number(image_stats, pdf_pages)
    image_feature_number(image_stats, I.GetIndexArray()[0],
                         pdf_pages, x_lims, y_lims)

    # vis the tic metrics
    plot_tic_number(image_stats, pdf_pages)
    image_tic_number(image_stats, I.GetIndexArray()[0],
                     pdf_pages, x_lims, y_lims)

    # vis the mab metrics
    plot_max_abun_number(image_stats, pdf_pages)
    image_max_abun_number(image_stats, I.GetIndexArray()[0],
                          pdf_pages, x_lims, y_lims)

    # vis the median metrics
    plot_median_number(image_stats, pdf_pages)
    image_median_number(image_stats, I.GetIndexArray()[0],
                        pdf_pages, x_lims, y_lims)

    # vis the max intensitsy metrics
    plot_max_int_number(image_stats, pdf_pages)
    image_max_int_number(image_stats, I.GetIndexArray()[0],
                         pdf_pages, x_lims, y_lims)

    # vis the  min intensitsy metrics
    plot_min_int_number(image_stats, pdf_pages)
    image_min_int_number(image_stats, I.GetIndexArray()[0],
                         pdf_pages, x_lims, y_lims)

    # vis the max intensitsy metrics
    plot_max_mz_number(image_stats, pdf_pages)
    image_max_mz_number(image_stats, I.GetIndexArray()[0],
                        pdf_pages, x_lims, y_lims)

    # vis the  min intensitsy metrics
    plot_min_mz_number(image_stats, pdf_pages)
    image_min_mz_number(image_stats, I.GetIndexArray()[0],
                        pdf_pages, x_lims, y_lims)

    # visualize the mean spectra
    if format_flags["centroid"]:
        plot_centroid_spectrum(I.GetXAxis(), I.GetMeanSpectrum(), pdf_pages)
    elif format_flags["profile"]:
        plot_profile_spectrum(I.GetXAxis(), I.GetMeanSpectrum(), pdf_pages)

    # get spectral coverage data:
    mean_bin, mean_coverage = calculate_spectral_coverage(I.GetXAxis(), I.GetMeanSpectrum())



    # plot spectral coverage data
    plot_coverage_barplot(mean_bin, mean_coverage, pdf_pages)


    write_summary_table(generate_table_data(I, x_lims, y_lims, image_stats),
                        pdf_pages)

    pdf_pages.close()
    print("QC sussefully generated at: ", outfile_path+"_agnostic_QC.pdf")

def report_calibrant_qc(I, # m2.imzMLReader (passing by ref allows faster computation)
                        outfile_path: str,  # path for output file
                        calfile_path: str,  # path to tsv file for calibrants
                        dist: float, # allowed distance to check for bulk metrics around theo. masses
                        ppm: float, # +- ppm cutoff for accuracy determination
                        sample_size: float = 1 # coverage of sample to be used for bulk calc, between 0 and 1
                        ):

    #  read in the calibrants
    calibrants = read_calibrants(calfile_path)

    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "_calibrant_QC")

    #Make a subsample to test accuracies on
    randomlist = make_subsample(I.GetNumberOfSpectra(), sample_size)

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    # per calibrant, bulk data is calculated inside the randomlist subsample
    for i in range(len(calibrants)):
        # adressing the field in df: calibrants.loc[i, "name"]

        # Create the data points for calibrant bulk accuracy cals
        cal_spectra = extract_calibrant_spectra(I, calibrants.loc[i, "mz"], randomlist, dist)

        # compute the metrics for bulk calibrant accuracies
        calibrants = collect_calibrant_stats(cal_spectra, calibrants, i)

        # plot the spectral data of a calibrant
        plot_calibrant_spectra(cal_spectra,
                               calibrants, i,
                               format_flags,
                               dist, pdf_pages)


    # barplot of the accuracies
    plot_accuracy_barplots(calibrants, pdf_pages)

    # calculate per pixel for nearest loc-max the accuracy
    accuracy_images, pixel_order = collect_accuracy_stats(I, calibrants, dist, format_flags)

    # calculate coverage from accuracy images
    calibrants = collect_calibrant_converage(accuracy_images, calibrants, ppm)

    # make accuracy images
    plot_accuracy_images(I, accuracy_images, calibrants, pixel_order, ppm, x_lims, y_lims, pdf_pages)

    # sumamary with coverage and avg. accuracy in non-zero pixels
    write_calibrant_summary_table(calibrants, pdf_pages)

    pdf_pages.close()
    print("QC sussefully generated at: ", outfile_path+"_calibrant_QC.pdf")


def report_regions_qc(I,  # m2.imzMLReader (passing by ref allows faster computation)
                      outfile_path: str,  # path for output file
                      regionfile_path=False,  # path to tsv file for region annotation
                      ):
    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "_region_QC")

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # parse the regionAnnotations:
    if regionfile_path:
        # readable to get dataFrame, col=0 x , col1= y, clo2= name
        region_table, region_image, nr_regions = parse_regionfile(regionfile_path, "annotation", x_lims, y_lims)
    else:
        region_table, region_image, nr_regions = label_connected_region(
            I)  # in nothing provides, make the conCompAnalysis
        write_region_tsv(region_table, outfile_path)

    # from annotation, get unique names in list()

    # additioally, get unique colors to correspond to the regions (neither white nor black)

    # plot the whole binary image
    image_cropped_binary(I.GetMaskArray()[0],
                         pdf_pages, x_lims, y_lims)

    # Plot the regions as colored blobs
    # 0 as non-recorded pixels, 1 as non-annotated pixels, 2-> end for
    # add numbers written on the pixel centra (with black border and their resp. color fill0)
    image_regions(I.GetMaskArray()[0], region_image,
                  pdf_pages, x_lims, y_lims)

    # intensity boxplot analysis
    # collect the metrics
    image_stats = collect_image_stats(I, ['index_nr', 'tic_nr'])

    # get the data grouped by annotation column:
    names_tic_bp, tic_bp = group_region_stat(region_image, I.GetIndexArray()[0], nr_regions, image_stats, "tic_nr")

    # plot the grouped data in a boxplot
    plot_boxplots(names_tic_bp, tic_bp, pdf_pages)

    # plot the averaged spectra of each region
    plot_regions_average(I, format_flags, region_image, nr_regions, pdf_pages)

    pdf_pages.close()
    print("QC sussefully generated at: ",  outfile_path+"_region_QC.pdf")
