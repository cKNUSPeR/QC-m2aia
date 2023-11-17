
from utils import *
from visualization import *
from cal_utils import *


def report_calibrant_qc(I, # m2.imzMLReader (passing by ref allows faster computation)
                        outfile_path: str,  # path for output file
                        calfile_path: str,  # path to tsv file for calibrants
                        dist: float # allowed distance to check for signals around thoe. masses
                        ):

    #  read in the calibrants
    calibrants = read_calibrants(calfile_path)

    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "_calibrant_QC")

    #Make a subsample to test accuracies on
    randomlist = make_subsample(I.GetNumberOfSpectra(), 0.01)

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
    accuracy_images = collect_accuracy_stats(I, calibrants, dist, format_flags)

    # make accuracy images

    # calc coverage

    # sumamary with coverage and avg. accuracy in non-zero pixels

    pdf_pages.close()
    print("QC sussefully generated at: ", outfile_path+"_calibrant_QC.pdf")


if __name__ == "__main__":
    file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\conv_output_centroided.imzML"
    I = m2.ImzMLReader(file_name)
    report_calibrant_qc(I,
                        r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\kidney",
                        r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\calibrants_9AA.csv",
                        0.025)
