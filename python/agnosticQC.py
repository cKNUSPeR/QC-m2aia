from visualization import *


def report_agnostic_qc(I,  # m2.imzMLReader (passing by ref allows faster computation)
                       outfile_path: str,  # path for output file
                       ):
    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "_agnostic_QC")

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

    write_summary_table(generate_table_data(I, x_lims, y_lims, image_stats),
                        pdf_pages)

    pdf_pages.close()
    print("QC sussefully generated at: ", outfile_path+"_agnostic_QC.pdf")

if __name__ == "__main__":
    file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\conv_output_centroided.imzML"
    I = m2.ImzMLReader(file_name)
    report_agnostic_qc(I,
                      r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\kidney_w_regions")