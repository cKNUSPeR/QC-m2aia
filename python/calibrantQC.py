
from utils import *
from visualization import *
from cal_utils import *


def report_calibrant_QC(I, # m2.imzMLReader (passing by ref allows faster computation)
                        file_path: str,  # path for output file
                        calfile_path: str,  # path to tsv file for calibrants
                        dist: float # allowed distance to check for signals around thoe. masses
                        ):

    #  read in the calibrants
    cal_names, cal_masses = read_calibrants(calfile_path)

    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(file_path, "calibrant_QC")

    #Make a subsample to test accuracies on
    randomlist = make_subsample(I.GetNumberOfSpectra(), 0.01)

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    # Create the data points for calibrant bulk accuracy cals
    cal_spectra = extract_calibrant_spectra(I, cal_masses, randomlist, dist)

    # compute the metrics for bulk calibrant accuracies
    bulkcal_mask, map_list, wavg_list = collect_calibrant_stats(cal_spectra, cal_masses)

    # plot how the calibrants show up in their spectra
        # differentiante the plotting :
        # 1) with profile or centriod  map&wavg
        #2) only data points + map&wavg
        #2.2) zoom of 150% around both metrics with only data points
        #3) zoom on minimal and maximal data points ()
    plot_calibrant_spectra(cal_names, cal_masses, cal_spectra,
                           bulkcal_mask, map_list, wavg_list,
                           dist, pdf_pages)


    # summarize bulk data:
    # bar charts with data (rel.) of bulk measurements in ppm
        # convert dist to ppm

    # calculate per pixel for nearest loc-max the accuracy

    # make accuracy images

    # calc coverage

    # sumamary with coverage and avg. accuracy in non-zero pixels