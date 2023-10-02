# analysis_utils.py
#
#
#
#  Contains functions for inputs, outputs, and basic plotting
#

import sys
import numpy as np


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // read_file_list //
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def read_file_list():

    with open('files.in', 'r') as file:
        root = file.readline()
        root = root.rstrip("\n")
        num = file.readline()
        num = int(num)
        filelist = np.array([None] * num, dtype=object)
        for i in range(num):
            j=i+1
            fname = file.readline()
            fname = fname.rstrip("\n")            
            filelist[i] = fname
    return root, num, filelist

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# // read_file_list //
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def print_data_to_file(num, filelist_short, datacube, varnames):


    outfile = "output.txt"

    char_count_array = np.array([len(s) for s in filelist_short])                      
    maxchar = np.amax(char_count_array)
    istr = "i  filenames"
    istr = "{:<{}}".format(istr, maxchar+2) 
    format_number = "{:9.4f}"
    format_string = "{:>9}"
    index = np.where(np.array(datacube[:,0]) != 0.0)[0]

    with open(outfile,"w") as f:
        for i in range(num):

            endstr = ' '
            ii=i+1
            if ii==1:
                # This prints the header information
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", file=f)
                print("CESM ExoCAM diagnostic output using analysis.py", file=f)
                print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~", file=f)
                print(istr, end=endstr,flush=True, file=f)
                for j in index:
                    if (j == np.amax(index)): endstr=('\n')
                    print(format_string.format(varnames[j]), flush=True,  end=endstr, file=f)
                endstr =' '
            # This prints file index and name    
            if ii < 10:   istr2 = str(ii) + '   ' + filelist_short[i]
            if ii >= 10:  istr2 = str(ii) + '  '  + filelist_short[i]
            if ii >= 100: istr2 = str(ii) + ' '   + filelist_short[i]
            istr2 = "{:<{}}".format(istr2, maxchar+4)
            print(istr2, end=endstr, flush=True, file=f)            
            for j in index:                
                # This prints the output data
                if (j == np.amax(index)): endstr=('\n')
                print(format_number.format(datacube[j,i]), end=endstr, flush=True,  file=f)
    print("output data written to ", outfile)
