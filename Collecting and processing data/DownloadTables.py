# =====================================================================
# Script that downloads tables from 2MASS and AllWISE catalogs using
# IRSA's Simple Cone Search. The center of the cones are given by the
# user. The program also searchs for equal rows.
# ---------------------------------------------------------------------
# Author: Pablo Gómez Nicolás
# Year: 2018/19
# ---------------------------------------------------------------------
# Documentation: https://irsa.ipac.caltech.edu/docs/vo_scs.html
# ---------------------------------------------------------------------
# Input (when asked by the program):
#   * File, in ASCII format, that contains the position (right
#     ascension and declination) of the search region centers. The
#     right ascension must be in the first column and the declination
#     in the second column.
#
# Output:
#   * The tables from 2MASS and AllWISE catalogs corresponding to the
#     positions in the file, in FITS format. The cone radius is 15'.
# =====================================================================

# Useful packages
import sys
import numpy as np
import subprocess
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab

# START OF THE SCRIPT
print('\nSCRIPT RUNNING...\n')

# Filename with the positions needed
print('Please write the name of the file that contains the positions of the '
      'search region centers (with extension)') #CoordinateFile.txt
filename = input();

# Reading the data
print('Trying to read data file...')
try:
    data = np.loadtxt(filename);
    print('Data file read successfully')
except:
    print('Warning: The input file does not exist. Please check the name')
    exit(1)

# Counting the repeated positions
print('Counting repeated rows...')
n_repeatedrows = 0;
for i in range(len(data)):
    for j in range(len(data)-i-1):
        if data[i,0] == data [i+j+1,0]:
            if data[i,1] == data[i+j+1,1]:
                print('The rows ' + str(i+1) + ' and ' + str(i+j+2) + 
                      ' are equal')
                n_repeatedrows = n_repeatedrows + 1;
print('There were found ' + str(n_repeatedrows) + ' repeated rows')

# Downloading the files
print('Downloading the files...')
for i in range(len(data)):
    bash2MASS = ('wget -O 2MASSTable' + str(i+1) + '.tbl'
                 ' "https://irsa.ipac.caltech.edu/SCS?table=fp_psc&RA=' +
                 str(data[i,0]) + '&DEC=' + str(data[i,1]) + 
                 '&SR=0.25&format=fits"');
    bashAllWISE = ('wget -O AllWISETable' + str(i+1) + '.tbl'
                   ' "https://irsa.ipac.caltech.edu/SCS?table=allwise_p3as_psd'
                   '&RA=' + str(data[i,0]) + '&DEC=' + str(data[i,1]) + 
                   '&SR=0.25&format=fits"');
    subprocess.run(bash2MASS, shell=True)
    subprocess.run(bashAllWISE, shell=True)

print('Files downloaded successfully')

print('\nSCRIPT FINISHED\n')
