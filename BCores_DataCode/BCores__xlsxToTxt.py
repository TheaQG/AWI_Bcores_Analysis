import numpy as np
import os
import pandas as pd
from pandas import ExcelWriter
import openpyxl
from os import walk
import xlrd
import csv

"""
    Code to generate .txt files from a number of .xlsx files in directory with a number of sheets each.
"""

# Define old path from where to get .xlsx files from
oldPath = '/home/thea/Desktop/MesterTesen/Data/datasets/B_cores_AWI/AWI_Bcores__Cleaned_xlsx/'
# Define new path to save new .csv files to
newPath ='/home/thea/Documents/KUFysik/MesterTesen/Data/datasets/B_cores_AWI/AWI_Bcores__Cleaned_CSV/'


# Find the names of .xlsx files in old directory
f = []
for (dirpath, dirnames, filenames) in walk(oldPath):
    f.extend(filenames)


# Define function to create .txt file from .xlsx file
def xlsx_to_txt(path, output_txt, input_xlsx, sheetNo):

    # open output .txt
    with open(path + output_txt, 'w') as myTxtFile:
        # define a writer wr
        wr = csv.writer(myTxtFile, delimiter = '\t')

        # open input xlsx file
        with xlrd.open_workbook(input_xlsx) as myXlsxFile:
            # get wanted sheet
            mysheet = myXlsxFile.sheet_by_index(sheetNo)
            # write the rows to the output file
            for rownum in range(mysheet.nrows):
                wr.writerow(mysheet.row_values(rownum))

    return



# Go through each file in directory
for file in f:
    print(file)
    # Get the filename w.o. the extension '.xlsx'
    filename = '.'.join(file.split('.')[:-1]) if '.' in file else file
    # Get all sheet names in file
    sheets = pd.ExcelFile(oldPath + file).sheet_names
    # Go through each sheet name and create a '.txt' file
    for s in range(len(sheets)):
        filename_out = filename + '__' + sheets[s] + '.txt'
        print(s)
        xlsx_to_txt(newPath, filename_out, oldPath + file, s)
