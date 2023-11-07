
from utils import *
from visualization import *

def report_agnostic_QC(file_name: str):

    # execute pasing of imzML file
    I = m2.ImzMLReader(file_name)

    # Create a PDF file to save the figures
    pdf_file_path = file_name[:-6] + "_output_QC(1_4).pdf"
    pdf_pages = matplotlib.backends.backend_pdf.PdfPages(pdf_file_path)

    # create format flag dict to check formatting of imzML file
    format_flags = evaluate_formats(I.GetSpectrumType())

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    # get a randomlist for computationally intense tasks.
    randomlist = make_subsample(I.GetNumberOfSpectra(), 0.01)

    image_full_binary(I.GetMaskArray()[0],
                           pdf_pages)

    image_cropped_binary(I.GetMaskArray()[0],
                              pdf_pages, x_lims, y_lims)

    image_pixel_index(I.GetIndexArray()[0],
                           pdf_pages, x_lims, y_lims)

    image_stats = collect_image_stats(I)

    # visualize the feature numbers
    plot_feature_number(image_stats[0],image_stats[1],
                        pdf_pages)
    image_feature_number(mask_bad_image(image_stats[0], image_stats[1], I.GetIndexArray()[0]),
                         pdf_pages, x_lims, y_lims)

    # vis the tic metrics
    plot_tic_number(image_stats[0], image_stats[2],
                        pdf_pages)
    image_tic_number(mask_bad_image(image_stats[0], image_stats[2], I.GetIndexArray()[0]),
                         pdf_pages, x_lims, y_lims)

    # vis the mab metrics
    plot_max_abun_number(image_stats[0], image_stats[8],
                    pdf_pages)
    image_max_abun_number(mask_bad_image(image_stats[0], image_stats[8], I.GetIndexArray()[0]),
                     pdf_pages, x_lims, y_lims)

    # vis the median metrics
    plot_median_number(image_stats[0], image_stats[3],
                    pdf_pages)
    image_median_number(mask_bad_image(image_stats[0], image_stats[3], I.GetIndexArray()[0]),
                     pdf_pages, x_lims, y_lims)

    # vis the max intensitsy metrics
    plot_max_int_number(image_stats[0], image_stats[4],
                    pdf_pages)
    image_max_int_number(mask_bad_image(image_stats[0], image_stats[4], I.GetIndexArray()[0]),
                     pdf_pages, x_lims, y_lims)

    # vis the  min intensitsy metrics
    plot_min_int_number(image_stats[0], image_stats[5],
                    pdf_pages)
    image_min_int_number(mask_bad_image(image_stats[0], image_stats[5], I.GetIndexArray()[0]),
                     pdf_pages, x_lims, y_lims)

    # vis the max intensitsy metrics
    plot_max_mz_number(image_stats[0], image_stats[6],
                        pdf_pages)
    image_max_mz_number(mask_bad_image(image_stats[0], image_stats[6], I.GetIndexArray()[0]),
                         pdf_pages, x_lims, y_lims)

    # vis the  min intensitsy metrics
    plot_min_mz_number(image_stats[0], image_stats[7],
                        pdf_pages)
    image_min_mz_number(mask_bad_image(image_stats[0], image_stats[7], I.GetIndexArray()[0]),
                         pdf_pages, x_lims, y_lims)



    pdf_pages.close()
    print("QC sussefully generated at: ", pdf_file_path)
