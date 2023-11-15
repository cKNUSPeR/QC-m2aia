

from agnosticQC import report_agnostic_QC
from calibrantQC import report_calibrant_QC

#report_agnostic_QC(r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\F13_FFPE_Core2.imzML")



report_calibrant_QC(I,
                    r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\data\exmpl_cont\F13_FFPE_Core2",
                    r"C:\Users\Jannik\Documents\Uni\Master_Biochem\4_Semester\M2aia\calibrants_9AA.csv",
                    0.025)