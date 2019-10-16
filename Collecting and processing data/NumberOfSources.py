# =====================================================================
# Script that takes the downloaded tables and analyses the number of
# sources found in each of the cones. You can also choose an upper and
# lower threshold for the number of sources in 2MASS and AllWISE
# regions. The program informs about the number of the fields that 
# surpass the thresholds.
# ---------------------------------------------------------------------
# Author: Pablo Gómez Nicolás
# Year: 2018/19
# ---------------------------------------------------------------------
# Documentation: http://www.star.bris.ac.uk/~mbt/stilts/
# ---------------------------------------------------------------------
# Input (when asked by the program):
#   * File, in ASCII format, that contains the position (right
#     ascension and declination) of the search region centers (just for
#     counting the number of sources).
#   * Tables obtained from DownloadTables.py script.
#   * Thresholds for the number of sources in each cone, for 2MASS and
#     AllWISE.
#   * Stilts is needed!
#
# Output:
#   * Files with information with the size of the tables, in ASCII
#     format.
#   * Number (labels) of tables that surpass the thresholds. Shown in
#     the terminal.
#   * Histograms of the number of sources in the fields, for 2MASS and
#     AllWISE, before and after the cut, in PNG and EPS format
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

print('Type 0 if you don\'t have created the files with the number of sources.')
print('If the files are already created, type anything else')
ind = input();
if ind == '0':
    # Filename with the positions needed
    print('Please write the name of the file that contains the positions '
          'of the search region centers (with extension)') #CoordinateFile.txt
    filename = input();

    # Reading the data
    print('Trying to read data file...')
    try:
        data = np.loadtxt(filename);
        print('Data file read succesfuly')
    except:
        print('Warning: The input file does not exist. Please check the name')
        exit(1)
    ndata = len(data); #Number of fields studied
    
    # Counting the size of each table
    print('Creating files with size data...')
    # Neccesary to create files with the size information
    f_2MASS = open('Sizes2MASS.txt', 'w');
    f_AllWISE = open('SizesAllWISE.txt', 'w');
    # Writing the files
    for i in range(ndata):
        bash2MASS = ('echo "$(sh stilts tpipe in=2MASSTable' + str(i+1) +
                     '.tbl omode=count)" >> Sizes2MASS.txt');
        bashAllWISE = ('echo "$(sh stilts tpipe in=AllWISETable' + str(i+1) +
                       '.tbl omode=count)" >> SizesAllWISE.txt');
        subprocess.run(bash2MASS, shell=True)
        subprocess.run(bashAllWISE, shell=True)
    f_2MASS.close()
    f_AllWISE.close()
    print('Files created successfully')

# Plotting histograms
print('Plotting the histograms...')
sizes2MASS = np.loadtxt('Sizes2MASS.txt', usecols=3);
sizesAllWISE = np.loadtxt('SizesAllWISE.txt', usecols=3);
# 2MASS
plt.hist(sizes2MASS, log=True, histtype='step', facecolor=[0.5,0,0], 
         ec=[0.5,0,0], bins=range(0,7250,250)) #Bins limits chosen manually
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plotfile = 'NumberOfSources2MASS_Dirty';
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()
# AllWISE
plt.hist(sizesAllWISE, log=True, histtype='step', facecolor=[0,0,0.5],
         ec=[0,0,0.5], bins=range(600,5600,200)) #Bins limits chosen manually
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plotfile = 'NumberOfSourcesAllWISE_Dirty';
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()

# Selecting the thresholds
print('Histograms created; you can open them in the current folder')
print('Choose a maximum for the number of sources in 2MASS cones')
max_2MASS = input(); max_2MASS = int(max_2MASS); #4000
print('Choose a minimum for the number of sources in 2MASS cones')
min_2MASS = input(); min_2MASS = int(min_2MASS); #0
print('Choose a maximum for the number of sources in AllWISE cones')
max_AllWISE = input(); max_AllWISE = int(max_AllWISE); #5000
print('Choose a minimum for the number of sources in AllWISE cones')
min_AllWISE = input(); min_AllWISE = int(min_AllWISE); #1500

# Plots with vertical lines
print('Plotting the histograms...')
# 2MASS
plt.hist(sizes2MASS, log=True, histtype='step', facecolor=[0.75,0,0],
         ec=[0.75,0,0], bins=range(0,7250,250)) #Bins limits chosen manually
plt.axvline(x=max_2MASS, c='k', linewidth=1, linestyle='--')
plt.axvline(x=min_2MASS, c='k', linewidth=1, linestyle='--')
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plotfile = 'NumberOfSources2MASS_Dirty';
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()
# AllWISE
plt.hist(sizesAllWISE, log=True, histtype='step', facecolor=[0,0,0.5],
         ec=[0,0,0.5], bins=range(600,5600,200)) #Bins limits chosen manually
plt.axvline(x=max_AllWISE, c='k', linewidth=1, linestyle='--')
plt.axvline(x=min_AllWISE, c='k', linewidth=1, linestyle='--')
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plotfile = 'NumberOfSourcesAllWISE_Dirty';
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()

# Updating the 2MASS array
indmax_2MASS = sizes2MASS <= max_2MASS;
indmin_2MASS = sizes2MASS >= min_2MASS;
ind_2MASS = indmax_2MASS * indmin_2MASS;
print('Labels of the tables that are not in the wanted range (2MASS):')
ndata = len(sizes2MASS);
for i in range(ndata):
    if not ind_2MASS[i]:
        print(i+1)
sizes2MASS = sizes2MASS[ind_2MASS];
# Updating the AllWISE array
indmax_AllWISE = sizesAllWISE <= max_AllWISE;
indmin_AllWISE = sizesAllWISE >= min_AllWISE;
ind_AllWISE = indmax_AllWISE * indmin_AllWISE;
print('Labels of the tables that are not in the wanted range (AllWISE):')
for i in range(ndata):
    if not ind_AllWISE[i]:
        print(i+1)
sizesAllWISE = sizesAllWISE[ind_AllWISE];

# Plotting new histograms
print('Plotting the new histograms...')
# 2MASS
plt.hist(sizes2MASS, log=True, histtype='step', facecolor=[0.75,0,0],
         ec=[0.75,0,0], bins=range(0,2875,125)) #Bins limits chosen manually
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plotfile = 'NumberOfSources2MASS_Clean';
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()
# AllWISE
plt.hist(sizesAllWISE, log=True, histtype='step', facecolor=[0, 0, 0.5],
         ec=[0,0,0.5], bins=range(1750,5125,125)) #Bins limits chosen manually
plt.xlabel('Number of sources', fontsize=18)
plt.ylabel('Frequency', fontsize=18)
plt.xticks(fontsize=16)
plt.yticks(fontsize=16)
plotfile = 'NumberOfSourcesAllWISE_Clean';
plt.subplots_adjust(left=0.15, right=0.95, top=0.95, bottom=0.15)
plt.savefig(plotfile+'.png')
plt.savefig(plotfile+'.eps')
plt.close()

print('New histograms plotted successfully')

print('\nSCRIPT FINISHED\n')
