# =====================================================================
# Script that creates the table with the crosscorrelated objects from
# 2MASS and AllWISE catalogs. In order to do that, it takes the three
# tables obtained from xMatch and uses the information contained in
# them.
# ---------------------------------------------------------------------
# Author: Pablo Gómez Nicolás
# Year: 2018/19
# ---------------------------------------------------------------------
# Input (when asked by the program):
#   * Minimum probability of association (number from 0 to 1).
#   * Quality flags chosen for 2MASS and AllWISE catalogs. The possible
#     choices are A, AB and ABC.
#   * Files obtained as output from the xMatch service, in FITS format.
#   * 2MASS and AllWISE files that contains the photometric data, in
#     FITS format.
#   * It is neccesary to create the folder to save the data manually.
#
# Output:
#   * One table with all the data that will be used to apply machine
#     learning techniques, in ASCII and FITS format.
# =====================================================================

# Useful packages
import sys
import numpy as np
from astropy.io import fits
import subprocess

# START OF THE SCRIPT
print('\nSCRIPT RUNNING...\n')

# Data input - probability
print('Please write the minimum probability of association (from 0 y 1)')
probQuality = input(); #0.69133
probQuality = float(probQuality);
# Data input - flags
print('Please write the quality flags admitted for 2MASS and AllWISE '
      'respectively (A, AB, ABC)');
f2MASS = input(); fAllWISE = input(); #ABC,ABC
# We create two lists with all the possible values for the flags
n = len(f2MASS);
listflag2MASS = [];
for i in range(n):
    for j in range(n):
        for k in range(n):
            listflag2MASS.append(f2MASS[i] + f2MASS[j] + f2MASS[k]);

n = len(fAllWISE);
listflagAllWISE = [];
for i in range(n):
    for j in range(n):
        for k in range(n):
            listflagAllWISE.append(fAllWISE[i] + fAllWISE[j] + fAllWISE[k]);

# We create the root for saving files
root = 'Prob' + str(probQuality) + '_Flags' + f2MASS + '_' + fAllWISE;

# Reading the data
print('Trying to read data files...')
try:
    hdulist1 = fits.open('XMatchSet1.tbl');
    xMatch1 = hdulist1[1].data;
    hdulist2 = fits.open('XMatchSet2.tbl');
    xMatch2 = hdulist2[1].data;
    hdulist3 = fits.open('XMatchSet3.tbl');
    xMatch3 = hdulist3[1].data;
    xMatch = [xMatch1, xMatch2, xMatch3];
except:
    print('Warning: Please check that you have your output files from the '
          'xMatch service in the folder with the names "XMatchClean_.fits')
    exit(1)

try:
    hdulist1 = fits.open('2MASSSet1.tbl');
    MASS1 = hdulist1[1].data;
    hdulist2 = fits.open('2MASSSet2.tbl');
    MASS2 = hdulist2[1].data;
    hdulist3 = fits.open('2MASSSet3.tbl');
    MASS3 = hdulist3[1].data;
    MASS = [MASS1, MASS2, MASS3];
    hdulist1 = fits.open('AllWISESet1.tbl');
    AllWISE1 = hdulist1[1].data;
    hdulist2 = fits.open('AllWISESet2.tbl');
    AllWISE2 = hdulist2[1].data;
    hdulist3 = fits.open('AllWISESet3.tbl');
    AllWISE3 = hdulist3[1].data;
    AllWISE = [AllWISE1, AllWISE2, AllWISE3];
    print('Data file read successfuly')
except:
    print('Warning: Please check that you have your files with the photometric '
          'data in the folder with the names "2MASSTable_Set_.fits and '
          'AllWISETable_Set_.fits')
    exit(1)

# We work with the three sets separately and then we join the results
final2MASS_names = [];
finalJ = [];
finalsigJ = [];
finalH = [];
finalsigH = [];
finalK = [];
finalsigK = [];
finalAllWISE_names = [];
finalW1 = [];
finalsigW1 = [];
finalW2 = [];
finalsigW2 = [];
finalW3 = [];
finalsigW3 = [];
for i in range(3):
    print('WORKING WITH SET ' + str(i+1))
    # Columns of interest in the xMatch output
    probabilities = xMatch[i]['proba_AB'];
    MASS_sources = xMatch[i]['2mass_designation'];
    AllWISE_sources = xMatch[i]['allwise_designation'];
    # We take only the rows with identified data
    i_cross = xMatch[i]['nPos'] > 1;
    probabilities = probabilities[i_cross];
    MASS_sources = MASS_sources[i_cross];
    AllWISE_sources = AllWISE_sources[i_cross];
    n = len(MASS_sources); 
    print('Number of matches in the set ' + str(i+1) +': ' + str(n))
    
    # Aplying quality filters
    print('Cleansing data with the quality filters...')
    # We create the list with the input data
    MASS_names = MASS[i]['designation'];
    J = MASS[i]['j_m'];
    sigJ = MASS[i]['j_msigcom'];
    H = MASS[i]['h_m'];
    sigH = MASS[i]['h_msigcom'];
    K = MASS[i]['k_m'];
    sigK = MASS[i]['k_msigcom'];
    MASS_quality = MASS[i]['ph_qual'];
    MASS_use = MASS[i]['use_src'];
    AllWISE_names = AllWISE[i]['designation'];
    W1 = AllWISE[i]['w1mpro'];
    sigW1 = AllWISE[i]['w1sigmpro'];
    W2 = AllWISE[i]['w2mpro'];
    sigW2 = AllWISE[i]['w2sigmpro'];
    W3 = AllWISE[i]['w3mpro'];
    sigW3 = AllWISE[i]['w3sigmpro'];
    AllWISE_quality = AllWISE[i]['ph_qual'];
    AllWISE_quality = AllWISE_quality.astype('<U3') #Only W123 information
    AllWISE_use = AllWISE[i]['use_src'];
    
    # We take only the sources with prob>=probQuality
    i_prob = probabilities>=probQuality;
    MASS_sources = MASS_sources[i_prob];
    AllWISE_sources = AllWISE_sources[i_prob];

    # We take only the sources that are crossmatched
    i_2MASS = np.isin(MASS_names,MASS_sources);
    i_AllWISE = np.isin(AllWISE_names,AllWISE_sources);
    MASS_names = MASS_names[i_2MASS];
    J = J[i_2MASS];
    sigJ = sigJ[i_2MASS];
    H = H[i_2MASS];
    sigH = sigH[i_2MASS];
    K = K[i_2MASS];
    sigK = sigK[i_2MASS];
    MASS_quality = MASS_quality[i_2MASS];
    MASS_use = MASS_use[i_2MASS];
    AllWISE_names = AllWISE_names[i_AllWISE];
    W1 = W1[i_AllWISE];
    sigW1 = sigW1[i_AllWISE];
    W2 = W2[i_AllWISE];
    sigW2 = sigW2[i_AllWISE];
    W3 = W3[i_AllWISE];
    sigW3 = sigW3[i_AllWISE];
    AllWISE_quality = AllWISE_quality[i_AllWISE];
    AllWISE_use = AllWISE_use[i_AllWISE];
    
    n = len(MASS_names);
    print('Number of sources that surpass the minimum probability: ' + str(n))
    
    # We order the data so we can work without problems
    i_2MASS = [];
    i_AllWISE = [];
    for l in range(len(MASS_sources)):
        j = np.where(MASS_names == MASS_sources[l])[0][0];
        k = np.where(AllWISE_names == AllWISE_sources[l])[0][0];
        i_2MASS.append(j);
        i_AllWISE.append(k);
    MASS_names = MASS_names[i_2MASS];
    J = J[i_2MASS];
    sigJ = sigJ[i_2MASS];
    H = H[i_2MASS];
    sigH = sigH[i_2MASS];
    K = K[i_2MASS];
    sigK = sigK[i_2MASS];
    MASS_quality = MASS_quality[i_2MASS];
    MASS_use = MASS_use[i_2MASS];
    AllWISE_names = AllWISE_names[i_AllWISE];
    W1 = W1[i_AllWISE];
    sigW1 = sigW1[i_AllWISE];
    W2 = W2[i_AllWISE];
    sigW2 = sigW2[i_AllWISE];
    W3 = W3[i_AllWISE];
    sigW3 = sigW3[i_AllWISE];
    AllWISE_quality = AllWISE_quality[i_AllWISE];
    AllWISE_use = AllWISE_use[i_AllWISE];
    
    # We clean the data using the quality flags
    i_2MASS = np.isin(MASS_quality, listflag2MASS);
    i_AllWISE = np.isin(AllWISE_quality, listflagAllWISE);
    ind = i_2MASS * i_AllWISE;
    MASS_names = MASS_names[ind];
    J = J[ind];
    sigJ = sigJ[ind];
    H = H[ind];
    sigH = sigH[ind];
    K = K[ind];
    sigK = sigK[ind];
    MASS_use = MASS_use[ind];
    AllWISE_names = AllWISE_names[ind];
    W1 = W1[ind];
    sigW1 = sigW1[ind];
    W2 = W2[ind];
    sigW2 = sigW2[ind];
    W3 = W3[ind];
    sigW3 = sigW3[ind];
    AllWISE_use = AllWISE_use[ind];

    i_2MASS = MASS_use > 0;
    i_AllWISE = AllWISE_use > 0;
    ind = i_2MASS * i_AllWISE;
    MASS_names = MASS_names[ind];
    J = J[ind];
    sigJ = sigJ[ind];
    H = H[ind];
    sigH = sigH[ind];
    K = K[ind];
    sigK = sigK[ind];
    AllWISE_names = AllWISE_names[ind];
    W1 = W1[ind];
    sigW1 = sigW1[ind];
    W2 = W2[ind];
    sigW2 = sigW2[ind];
    W3 = W3[ind];
    sigW3 = sigW3[ind];

    print('Quality filters applied')
    n = len(MASS_names);
    print('Number of sources that satisfy the quality filters: ' + str(n))
    # We put the data into the total data list
    final2MASS_names = np.concatenate((final2MASS_names, MASS_names));
    finalJ = np.concatenate((finalJ, J));
    finalsigJ = np.concatenate((finalsigJ, sigJ));
    finalH = np.concatenate((finalH, H));
    finalsigH = np.concatenate((finalsigH, sigH));
    finalK = np.concatenate((finalK, K));
    finalsigK = np.concatenate((finalsigK, sigK));
    finalAllWISE_names = np.concatenate((finalAllWISE_names, AllWISE_names));
    finalW1 = np.concatenate((finalW1, W1));
    finalsigW1 = np.concatenate((finalsigW1, sigW1));
    finalW2 = np.concatenate((finalW2, W2));
    finalsigW2 = np.concatenate((finalsigW2, sigW2));
    finalW3 = np.concatenate((finalW3, W3));
    finalsigW3 = np.concatenate((finalsigW3, sigW3));

# We create the final dataset and save it in ASCII format
print('Saving dataset in ASCII format...')
JHKW123 = np.c_[final2MASS_names, finalAllWISE_names, finalJ,finalsigJ, finalH, 
                finalsigH, finalK, finalsigK, finalW1, finalsigW1, finalW2, 
                finalsigW2, finalW3, finalsigW3];
rootfile = 'Folder_' + root + '/' + root + '_JHKW123'
np.savetxt(rootfile + '.txt', JHKW123, header=('id_2MASS\tid_AllWISE\tJ\tsigJ\t'
           'H\tsigH\tK\tsigK\tW1\tsigW1\tW2\tsigW2\tW3\tsigW3'), delimiter='\t',
           fmt='%s');
print('File saved successfully, check your folder')

# Saving file in FITS format
print('Saving dataset in FITS format...')
bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + rootfile + '.txt '
             'out=' + rootfile + '.fits');
subprocess.run(bashorder, shell=True)
print('Data file saved successfully, check your folder')

\print('\nSCRIPT FINISHED\n')
