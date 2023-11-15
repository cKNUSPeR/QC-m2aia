
import pandas as pd
import numpy as np
import skimage.measure as skim

import m2aia as m2

file_name = r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\conv_output_centroided.imzML"
I = m2.ImzMLReader(file_name)

labeled_image = skim.label(I.GetMaskArray()[0], connectivity=1)

# shape of image array:
rows, cols = labeled_image.shape
print(rows, cols)



# get a meshed grid to make x-y accesible
x_coords, y_coords =
print(x_coords, y_coords)
