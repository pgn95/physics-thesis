# =====================================================================
# Script that takes a data set with photometric information and applies
# different machine learning techniques. The data set must be in the
# directory of work. These techniques apply dimensionality reduction.
# It also transforms the data with labels.
# ---------------------------------------------------------------------
# Author: Pablo Gómez Nicolás
# Year: 2018/19
# ---------------------------------------------------------------------
# Documentation:
#     http://www.star.bris.ac.uk/~mbt/stilts/sun256/sun256.html
#     https://docs.astropy.org/en/stable/index.html
#     https://scikit-learn.org/stable/index.html
# ---------------------------------------------------------------------
# Input (when asked by the program):
#   * Minimum probability of association (number from 0 to 1) that
#     filtered the data set.
#   * Quality flags chosen for 2MASS and AllWISE catalogs. The possible
#     choices are A, AB and ABC.
#   * File obtained as output from the ExtractTable.py script, in FITS
#     format.
#   * Number that indicates the program to execute.
#   * The sdsstmasswisegood_corralra.fit archive has to be in the same
#     folder as this script.
#
# Output:
# 1.* A triangular figure that shows all the plots, in 2D, for all the
#     filters, in PNG and EPS formats.
#   * When asked by the user, a plot in 2D with one of the plots that
#     are in the triangular figure, in PNG and EPS formats.
#
# 2.* Data that characterises the PCA model created (variance, mean and
#     components). The script shows this data in the shell.
#   * A triangular figure that shows all the plots, in 2D, for all the
#     components, in PNG and EPS formats.
#   * When asked by the user, a plot in 2D with one of the plots that
#     are in the triangular figure, in PNG and EPS formats.
#   * Two tables, in ASCII and FITS formats, with all the data obtained
#     after aplying the PCA technique.
#   * Two plots, in PNG and EPS format, with the eigenvalues for the
#     PCA decomposition and the cumulative sum of the variances,
#     normalised to the unity.
#
# 3.* A triangular figure that shows all the plots, in 2D, for all the
#     components, in PNG and EPS formats.
#   * When asked by the user, a plot in 2D with one of the plots that
#     are in the triangular figure, in PNG and EPS formats.
# =====================================================================

# Useful packages
import sys
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn.decomposition import PCA
from sklearn.manifold import LocallyLinearEmbedding
from astroML.plotting import MultiAxes
from astropy.io import fits

# START OF THE SCRIPT
print('\nSCRIPT RUNNING...\n')

# Data input
print('Warning: it is neccesary to have executed the ExtractTable.py script '
      'before and keep the table obtained as output in FITS format')

print('Please write the minimum probability of association (from 0 y 1)')
probQuality = input(); probQuality = float(probQuality);

print('Please write the quality flags admitted for 2MASS and AllWISE '
      'respectively (A, AB, ABC)')
f2MASS = input(); fAllWISE = input();

# Choosing the task to perform
print('Please write the number corresponding to the task you want to perform')
print('1: Triangular representation with 2D plots of filters (JHKW123)')
print('2: Dimensionality reduction using Principal Component Analysis (PCA)')
print('3: Dimensionality reduction using Locally Linear Embedding (LLE)')
number = input();

# Roots to load/save files
root = 'Folder_Prob' + str(probQuality) + '_Flags' + f2MASS + '_' + fAllWISE;
root_file = 'Prob' + str(probQuality) + '_Flags' + f2MASS + '_' + fAllWISE;

# Reading the data
print('Trying to read data files...')
try:
    hdulist = fits.open(root + '/' + root_file + '_JHKW123.fits');
    data = hdulist[1].data;
    print('Data file read successfuly')
except:
    print('Warning: Please check that you have your files in the correct '
          'folder in FITS format')
    exit(1)

try:
    hdulist = fits.open('sdsstmasswisegood_corralra.fit');
    labeleddata = hdulist[1].data;
    print('Data file with labels read successfuly')
except:
    print('Warning: Please check that you have the "sdsstmasswisegood_'
          'corralra.fit" archive in the correct folder')
    exit(1)

# Creating data array
dataset = np.c_[data['J'], data['H'], data['K'], 
                data['W1'], data['W2'], data['W3']];
labeleddataset = np.c_[labeleddata['j'], labeleddata['h'], labeleddata['k'],
                       labeleddata['w1mpro'], labeleddata['w2mpro'],
                       labeleddata['w3mpro']];

# Solving problems with the subClass data from labeleddata
subclass = labeleddata['subClass'];
labelclass = labeleddata['class'];
for i in range(len(subclass)):
    if subclass[i] == '':
        subclass[i] = 'NULL';
    if subclass[i] == 'STARFORMING BROADLINE':
        subclass[i] = 'SFB';
    if subclass[i] == 'STARBURST BROADLINE':
        subclass[i] = 'SBB';
    if subclass[i] == 'AGN BROADLINE':
        subclass[i] = 'AGNB';
    if labelclass[i] == 'STAR':
        subclass[i] = 'NULL';

# TRIANGULAR REPRESENTATION
def func_triang_mag():
    print('\nTRIANGULAR REPRESENTATION WITH 2D PLOTS OF FILTERS (JHKW123)\n')
    # Full data plot
    fig = plt.figure(figsize=(8,8));
    labels = ['J', 'H', 'K', 'W1', 'W2', 'W3'];
    ax = MultiAxes(6, fig=fig, hspace=0, wspace=0);
    ax.scatter(dataset, s=1, color=[0.75,0,0], marker='o', alpha=0.05)
    ax.set_labels(labels)
    plt.title('Magnitudes\nFull data', fontsize=10)
    plotfile = root + '/JHKW123/' + root_file + '_JHKW123';
    fig.savefig(plotfile + '.png')
    fig.savefig(plotfile + '.eps')
    plt.close(fig)
    print('Triangular representation finished, check your JHKW123 folder')
    
    # Individual plots
    print('If you want to make a close-up plot, write 0. If this is not the '
          'case, write anything else')
    ind = input();
    while ind == '0':
        print('Write the magnitudes you want to plot')
        mag1 = input(); mag2 = input();
        # These plots will be done using STILTS
        bashorder = ('sh stilts plot2plane xpix=600 ypix=450 xlabel=' + mag1 +
                     ' ylabel=' + mag2 + ' texttype=latex fontsize=32 legend='
                     'false layer=mark xcrowd=0.5 in=' + root + '/' +
                     root_file + '_JHKW123.fits x=' + mag1 + ' minor=false y=' +
                     mag2 + ' shading=auto size=0 omode=out out=' + plotfile +
                     '_' + mag1 + mag2)
        subprocess.run(bashorder + '.png ofmt=png', shell=True)
        subprocess.run(bashorder + '.eps ofmt=eps', shell=True)
        print('Close-up finished, check your JHKW123 folder')
        print('If you want to make another close-up plot, write 0. If this is '
              'not the case, write anything else')
        ind = input();

    print('\nFINISHED THE REPRESENTATIONS OF THE DATASET\n')

# PCA
def func_pca():
    print('\nDIMENSIONALITY REDUCTION: PCA\n')
    # Creating the model using the photometric data
    print('Fitting the model...')
    pca = PCA(n_components=6);
    pca.fit(dataset);
    print('PCA model created successfully')
    
    # Model parameters
    mean = pca.mean_;
    comps = pca.components_;
    var = pca.explained_variance_;
    cum_var = pca.explained_variance_ratio_;
    cum_var = np.cumsum(cum_var);
    print('MEAN OF THE ORIGINAL DATA:\n', mean)
    print('PARAMETERS THAT CARACTERISE EACH COMPONENT :\n ', comps)
    print('VARIANCES EXPLAINED WITH EACH COMPONENT: \n', var)
    print('CUMULATIVE VARIANCE EXPLAINED WITH EACH COMPONENT: \n', cum_var)
    
    # Variance plots
    print('Plotting the variance and cumulative variance per component...')
    x = [1, 2, 3, 4, 5, 6];
    plt.plot(x, var,'-o', MarkerSize=2.5, Color=[0.75,0,0])
    plt.xlabel('Component number', fontsize=18)
    plt.ylabel('Variance', fontsize=18)
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.15, right=0.925, top=0.925, bottom=0.125)
    plt.tick_params(top=True, right=True, direction='inout')
    plt.savefig(root + '/PCA/VarianceExplained.eps')
    plt.savefig(root + '/PCA/VarianceExplained.png')
    plt.close()

    plt.plot(x, cum_var,'-o', MarkerSize=2.5, Color=[0.75,0,0])
    plt.xlabel('Component number', fontsize=18)
    plt.ylabel('Cumulative variance', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.axhline(y=0.9973, c='k', linewidth=1, linestyle='--')
    plt.subplots_adjust(left=0.15, right=0.925, top=0.925, bottom=0.125)
    plt.tick_params(top=True, right=True, direction='inout')
    plt.savefig(root + '/PCA/CumulativeVarianceExplained.eps')
    plt.savefig(root + '/PCA/CumulativeVarianceExplained.png')
    plt.close()

    # Adjusting the data to the model
    print('Adjusting the data to the model created...')
    proj0 = pca.transform(dataset);
    proj1 = pca.transform(labeleddataset);
    
    # Full data plot
    fig = plt.figure(figsize=(8,8));
    labels = ['PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6'];
    ax = MultiAxes(6, fig=fig, hspace=0, wspace=0);
    ax.scatter(proj0, s=1, color=[0.75,0,0], marker='o', alpha=0.05)
    ax.set_labels(labels)
    plt.title('PCA\nFull data', fontsize=10)
    plotfile = root + '/PCA/' + root_file + '_PCA';
    fig.savefig(plotfile + '.png')
    fig.savefig(plotfile + '.eps')
    plt.close(fig)
    print('Triangular representation finished, check your PCA folder')
    
    # Saving data in ASCII format
    print('Saving obtained data from PCA in ASCII format...')
    dataheading = 'id_2MASS\tid_AllWISE\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6';
    np.savetxt(plotfile + '.txt', np.c_[data['id_2MASS'], data['id_AllWISE'], 
               proj0], header=dataheading, delimiter='\t', fmt='%s')
    np.savetxt(plotfile + '_labeled.txt', np.c_[proj1, labeleddata['z'],
               labeleddata['class'], subclass],
               header=dataheading[20:]+'\tz\tclass\tsubClass', delimiter='\t',
               fmt='%s')
    print('Data file saved successfully, check your PCA folder')
    
    # Saving data in FITS format
    print('Saving obtained data from PCA in FITS format...')
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile + '.txt '
                 'out=' + plotfile + '.fits');
    subprocess.run(bashorder, shell=True)
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile + 
                 '_labeled.txt out=' + plotfile + '_labeled.fits');
    subprocess.run(bashorder, shell=True)
    print('Data file saved successfully, check your PCA folder')

    # Individual plots
    print('If you want to make a close-up plot, write 0. If this is not the '
          'case, write anything else')
    ind = input();
    while ind == '0':
        print('Write the components you want to plot')
        mag1 = input(); mag2 = input();
        # These plots will be done using STILTS
        bashorder = ('sh stilts plot2plane xpix=600 ypix=450 xlabel=' + mag1 +
                     ' ylabel=' + mag2 + ' texttype=latex fontsize=32 legend='
                     'false layer=mark in=' + plotfile + '.fits x=' + mag1 +
                     ' y=' + mag2 + ' shading=auto size=0 omode=out minor=false'
                     ' out=' + plotfile + '_' + mag1 + mag2);
        subprocess.run(bashorder + '.png ofmt=png', shell=True)
        subprocess.run(bashorder + '.eps ofmt=eps', shell=True)
        print('Close-up finished, check your PCA folder')
        print('If you want to make another close-up plot, write 0. If this is '
              'not the case, write anything else')
        ind = input();
    
    print('\nPCA TECHNIQUE APPLIED\n');

# LLE
def func_lle():
    print('\nDIMENSIONALITY REDUCTION: LLE\n')
    k = 50; #Number of neighbours used to perform LLE, chosen empirically
    # Creating the model using the photometric data
    print('Fitting the model...')
    embedding = LocallyLinearEmbedding(n_components=6, n_neighbors=k,
                                       eigen_solver='arpack');
    embedding.fit(dataset)
    print('LLE model created successfully')

    # Adjusting the data to the model
    print('Adjusting the data to the model created...')
    proj0 = embedding.transform(dataset);
    proj1 = embedding.transform(labeleddataset);
    
    # Full data plot
    fig = plt.figure(figsize=(8,8));
    labels = ['LL1', 'LL2', 'LL3', 'LL4', 'LL5', 'LL6'];
    ax = MultiAxes(6, fig=fig, hspace=0, wspace=0);
    ax.scatter(proj0, s=1, color=[0.75,0,0], marker='o', alpha=0.05)
    ax.set_labels(labels)
    plt.title('LLE\nFull data', fontsize=10)
    plotfile = root + '/LLE/' + root_file + '_LLE';
    fig.savefig(plotfile + '.png')
    fig.savefig(plotfile + '.eps')
    plt.close(fig)
    print('Triangular representation finished, check your LLE folder')
    
    # Saving data in ASCII format
    print('Saving obtained data from LLE in ASCII format...')
    dataheading = 'id_2MASS\tid_AllWISE\tLL1\tLL2\tLL3\tLL4\tLL5\tLL6';
    np.savetxt(plotfile + '.txt', np.c_[data['id_2MASS'], data['id_AllWISE'], 
               proj0], header=dataheading, delimiter='\t', fmt='%s')
    np.savetxt(plotfile + '_labeled.txt', np.c_[proj1, labeleddata['z'],
               labeleddata['class'], subclass],
               header=dataheading[20:]+'\tz\tclass\tsubClass', delimiter='\t',
               fmt='%s')
    print('Data file saved successfully, check your LLE folder')
    
    # Saving data in FITS format
    print('Saving obtained data from LLE in FITS format...')
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile + '.txt '
                 'out=' + plotfile + '.fits');
    subprocess.run(bashorder, shell=True)
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile + 
                 '_labeled.txt out=' + plotfile + '_labeled.fits');
    subprocess.run(bashorder, shell=True)
    print('Data file saved successfully, check your LLE folder')

    # Individual plots
    print('If you want to make a close-up plot, write 0. If this is not the '
          'case, write anything else')
    ind = input();
    while ind == '0':
        print('Write the components you want to plot')
        mag1 = input(); mag2 = input();
        # These plots will be done using STILTS
        bashorder = ('sh stilts plot2plane xpix=600 ypix=450 xlabel=' + mag1 +
                     ' ylabel=' + mag2 + ' texttype=latex fontsize=32 legend='
                     'false layer=mark in=' + plotfile + '.fits x=' + mag1 +
                     ' y=' + mag2 + ' shading=auto size=0 omode=out minor=false'
                     ' out=' + plotfile + '_' + mag1 + mag2);
        subprocess.run(bashorder + '.png ofmt=png', shell=True)
        subprocess.run(bashorder + '.eps ofmt=eps', shell=True)
        print('Close-up finished, check your LLE folder')
        print('If you want to make another close-up plot, write 0. If this is '
              'not the case, write anything else')
        ind = input();
    
    print('\nLLE TECHNIQUE APPLIED\n');

# Different options to execute
options = {1 : func_triang_mag,
           2 : func_pca,
           3 : func_lle,
}

# We execute the desired function
options[int(number)]()

