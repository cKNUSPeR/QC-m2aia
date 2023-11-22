from dependencies import *
from visualization import *
from utils import *


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




