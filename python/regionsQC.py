
import m2aia as m2
from region_utils import *
from utils import *
from visualization import *



# if no regions provided, make a connected component analysis
# and perform region control on that instead
# boxplot of intensities:
# avg full spectra vs avg spectra of the def. region


def report_regions_QC(I, # m2.imzMLReader (passing by ref allows faster computation)
                      outfile_path: str,  # path for output file
                      regionfile_path = False,  # path to tsv file for region annotation
                      ):
    # Create a PDF file to save the figures
    pdf_pages = make_pdf_backend(outfile_path, "region_QC")

    # get the image limits to crop and display only relevant parts
    x_lims, y_lims = evaluate_image_corners(I.GetMaskArray()[0])

    # parse the regionAnnotations:
    if regionfile_path:
        # readable to get dataFrame, col=0 x , col1= y, clo2= name
        region_table, region_image, nr_regions = parse_regionfile(regionfile_path, "annotation",  x_lims, y_lims)
    else:
        region_table, region_image, nr_regions = label_connected_region(I)# in nothing provides, make the conCompAnalysis
        write_region_tsv(region_table, outfile_path)

    # from annotation, get unique names in list()


    # additioally, get unique colors to correspond to the regions (neither white nor black)


    # plot the whole binary image
    image_cropped_binary(I.GetMaskArray()[0],
                         pdf_pages, x_lims, y_lims)

    # Plot the regions as colored blobs
    # 0 as non- recorded pixels, 1 as non-annotated pixels, 2-> end for
    # add numbers written on the pixel centra (with black border and their resp. color fill0)
    image_regions(I.GetMaskArray()[0], region_image,
                  pdf_pages, x_lims, y_lims)

    # intensity boxplot analysis
    # collect the metrics
    image_stats = collect_image_stats(I, ['index_nr', 'tic_nr'])

    # get the data grouped by annotation column:
    names_tic_bp, tic_bp = group_region_stat(region_image, I.GetIndexArray()[0], nr_regions, image_stats, "tic_nr")

    print(tic_bp)
    # plot the grouped data in a boxplot
    plot_boxplots(names_tic_bp, tic_bp, pdf_pages)

    plot_region_intesities()

    pdf_pages.close()


if True:
    file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\conv_output_centroided.imzML"
    I = m2.ImzMLReader(file_name)
    report_regions_QC(I, r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\kidney_w_regions",
                         r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\kidney_annotated_regions.tsv")

   # report_regions_QC(I, r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\kidney_wo_regions",


