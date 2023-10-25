from QC_v1_4 import generate_TTF_QC

# specify location of input data
sample_list = ["D:\\wittej\\data\\Mannheim centroided data TMA6.2\\freiburg_metabolite measurement_no normalization_cemos.imzML",
               "D:\\wittej\\data\\ICC_TMA62_jannik_versuch2",
               "D:\\wittej\\data\\TMA4_centroided\\FTICR.imzML"
               
               ]
               
cal =  "D:\\programs\\M2aia\\calibrants_9AA.csv"
            
for sample in sample_list:
    generate_TTF_QC(sample,cal)
    # output data is saved next to input data           
          