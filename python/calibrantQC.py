
from utils import *
from visualization import *
from cal_utils import *

def report_agnostic_QC(file_path: str, calfile_path: str):

    # execute pasing of imzML file
    I = m2.ImzMLReader(file_path)

    #  read in the calibrants
    cal_names, cal_masses = read_calibrants(calfile_path)

    # Create a PDF file to save the figures
    pdf_file_path = file_path[:-6] + "_calibrant_QC.pdf"
    pdf_pages = matplotlib.backends.backend_pdf.PdfPages(pdf_file_path)

    #Make a subsample to test accuracies on
    randomlist = make_subsample(I.GetNumberOfSpectra(), 0.01)

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    cal_spectra = extract_calibrant_spectra(I, cal_masses, randomlist, mz_bin=0.03)



