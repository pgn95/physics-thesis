# =====================================================================
# Script that takes a data set with photometric information and applies
# different machine learning techniques. The data set must be in the
# directory of work. These techniques apply clustering. It also applies
# the techniques to the labeled data set.
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
#   * Files obtained as output from the DimensionalityReduction.py
#     script, in FITS format.
#   * Number that indicates the input set of data (PCA, LLE).
#   * Number that indicates the dimension of the data you want to use.
#   * Number that indicates the program to execute.
#
# Output:
# 1.* For each K from 2 to 10, the positions of the centers and the sum
#     of the square distances to the centroids. The script shows this
#     data in the shell.
#   * For each K from 2 to 10, a table, in ASCII and FITS format, with
#     the data used and the labels obtained with K-Means method.
#   * A plot, in PNG and EPS format, that shows the square distances to
#     the centroids as a function of K.
# 2.* Data that characterises the PCA model created (variance, mean and
#     components). The script shows this data in the shell.
#   * A triangular figure that shows all the plots, in 2D, for all
#     the components, in PNG and EPS formats.
#   * When asked by the user, a plot in 2D with one of the plots that
#     are in the triangular figure, in PNG and EPS formats.
#   * A table, in ASCII and FITS format, with all the data obtained
#     after aplying the PCA technique.
#   * Two plots, in PNG and EPS format, with the eigenvalues for the
#     PCA decomposition and the cumulative sum of the variances,
#     normalised to the unity.
# =====================================================================

# Useful packages
import sys
import subprocess
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from sklearn.cluster import KMeans
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.mixture import GaussianMixture
from astroML.plotting import MultiAxes
from astropy.io import fits

# START OF THE SCRIPT
print('\nSCRIPT RUNNING...\n')

# Data input
print('Warning: it is neccesary to have executed the DimensionalityReduction.py'
      ' script before and keep the tables obtained as output in FITS format')

print('Please write the minimum probability of association (from 0 y 1)')
probQuality = input(); probQuality = float(probQuality);

print('Please write the quality flags admitted for 2MASS and AllWISE '
      'respectively (A, AB, ABC)')
f2MASS = input(); fAllWISE = input();
root = 'Folder_Prob' + str(probQuality) + '_Flags' + f2MASS + '_' + fAllWISE;

# Choosing the data to read
print('Please write the number corresponding to the dataset you want to use')
print('1: PCA dataset')
print('2: LLE dataset')
datainput = input();

# We load the dataset chosen
if datainput == '1':
    root_file = ('Prob' + str(probQuality) + '_Flags' + f2MASS + '_' +
                 fAllWISE + '_PCA');
    hdulist = fits.open(root + '/PCA/' + root_file + '.fits');
    data = hdulist[1].data;
    dataset = np.c_[data['PC1'], data['PC2'], data['PC3'], data['PC4'],
                    data['PC5'], data['PC6']];
    hdulist = fits.open(root + '/PCA/' + root_file + '_labeled.fits');
    labeleddata = hdulist[1].data;
    labeleddataset = np.c_[labeleddata['PC1'], labeleddata['PC2'],
                           labeleddata['PC3'], labeleddata['PC4'],
                           labeleddata['PC5'], labeleddata['PC6']];
    dataheading = 'id_2MASS\tid_AllWISE\tPC1\tPC2\tPC3\tPC4\tPC5\tPC6\tGroup';
    root_folder = '/PCA/';
elif datainput == '2':
    root_file = ('Prob' + str(probQuality) + '_Flags' + f2MASS + '_' +
                 fAllWISE + '_LLE');
    hdulist = fits.open(root + '/LLE/' + root_file + '.fits');
    data = hdulist[1].data;
    dataset = np.c_[data['LL1'], data['LL2'], data['LL3'], data['LL4'],
                    data['LL5'], data['LL6']];
    hdulist = fits.open(root + '/LLE/' + root_file + '_labeled.fits');
    labeleddata = hdulist[1].data;
    labeleddataset = np.c_[labeleddata['LL1'], labeleddata['LL2'],
                           labeleddata['LL3'], labeleddata['LL4'],
                           labeleddata['LL5'], labeleddata['LL6']];
    dataheading = 'id_2MASS\tid_AllWISE\tLL1\tLL2\tLL3\tLL4\tLL5\tLL6\tGroup';
    root_folder = '/LLE/';
else:
    print('Warning: Please indicate which dataset you want to use')
    exit(1)

# Choosing the number of components to use
print('Please write the number of components, from 3 to 6, you want to use')
n = input();

dataused = dataset[:,range(int(n))];
labeleddataused = labeleddataset[:,range(int(n))];

# Choosing the task to perform
print('Please write the number corresponding to the task you want to perform')
print('1: Clustering using K-Means (KMeans)')
print('2: Clustering using Mean Shift (MeanShift)')
print('3: Clustering using Gaussian Mixture Model (GMM)')
number = input();

# K-MEANS
def func_kmeans():
    print('\nCLUSTERING: K-MEANS\n')
    K = [2,3,4,5,6,7,8,9,10]; #Possible choices for the number of clusters
    square_error = []; #List to keep the errors for each value in K
    
    # Creating the model using the photometric data (for each K value)
    for i in range(len(K)):
        print('Fitting the model for K=%i...' % K[i])
        kmeans = KMeans(n_clusters=K[i], random_state=42);
        kmeans.fit(dataused);
        print('K-Means problem solved for K=%i' % K[i])
        
        # Model parameters
        centers = kmeans.cluster_centers_;
        square_error.append(kmeans.inertia_);
        print('POSITIONS OF THE CENTERS:\n', centers)
        print('SUM OF SQUARE DISTANCES TO THE CENTROIDS:\n ', square_error[i])
        
        # Adjusting the data to the model
        print('Obtaining the labels for K=%i...' % K[i])
        labels = kmeans.labels_;
        labeledlabels = kmeans.predict(labeleddataused);
        
        # Saving data in ASCII format
        plotfile = (root + '/KMeans/' + n + root_folder + root_file +
                    '_KMeans_K=' + str(K[i]));
        print('Saving data with labels in ASCII format...');
        np.savetxt(plotfile + '.txt', np.c_[data['id_2MASS'],
                   data['id_AllWISE'], dataset, labels], header=dataheading,
                   delimiter='\t', fmt='%s')
        np.savetxt(plotfile + '_labeled.txt', np.c_[labeleddataset,
                   labeledlabels, labeleddata['z'], labeleddata['class'],
                   labeleddata['subClass']],
                   header=dataheading[20:]+'\tz\tclass\tsubClass',
                   delimiter='\t', fmt='%s')
        print('Data file saved successfully, check your KMeans folder')
        
        # Saving data in FITS format
        print('Saving data with labels in FITS format...');
        bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                     '.txt out=' + plotfile + '.fits');
        subprocess.run(bashorder, shell=True)
        bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                     '_labeled.txt out=' + plotfile + '_labeled.fits');
        subprocess.run(bashorder, shell=True)
        print('Data file saved successfully, check your KMeans folder')

        # Saving centers information
        plotfile = (root + '/KMeans/' + n + root_folder + 'Centers_K=' + 
                    str(K[i]));
        print('Saving centers position...');
        np.savetxt(plotfile + '.txt', centers, header=dataheading[20:-6],
                   delimiter='\t', fmt='%s')
        print('Centers positions saved successfully, check your KMeans folder')
    
    # Residue plots
    print('Plotting the sum of the square distance to the centroids...')
    plt.plot(K, square_error,'-o', MarkerSize=2.5, Color=[0.75,0,0])
    plt.xlabel('K number', fontsize=18)
    plt.ylabel('Square distance', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.15, right=0.925, top=0.925, bottom=0.125)
    plt.tick_params(top=True, right=True, direction='inout')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0),
                         useMathText=True)
    plt.savefig(root + '/KMeans/' + n + root_folder + 'SquareError.eps')
    plt.savefig(root + '/KMeans/' + n + root_folder +'SquareError.png')
    plt.close()
    
    print('\nK-MEANS TECHNIQUE APPLIED\n');

# MEAN SHIFT
def func_meanshift():
    print('\nCLUSTERING: MEAN SHIFT\n')
    #Creating the model using photometric data
    print('Fitting the model...')
    print('If you want to choose the bandwidth, please write it. Else, write 0')
    bandwidth = input();
    if bandwidth == '0':
        bandwidth = estimate_bandwidth(dataused);
    else:
        bandwidth = float(bandwidth);
    print('Using a bandwidth of', bandwidth)
    ms = MeanShift(bandwidth=bandwidth);
    ms.fit(dataused);
    print('Mean Shift problem solved')
    
    # Model parameters
    centers = ms.cluster_centers_;
    print('POSITIONS OF THE CENTERS:\n', centers)
    
    # Adjusting the data to the model
    print('Obtaining the labels...')
    labels = ms.labels_;
    labeledlabels = ms.predict(labeleddataused);
    
    # Saving data in ASCII format
    plotfile = (root + '/MeanShift/' + n + root_folder + root_file +
                '_MeanShift');
    print('Saving data with labels in ASCII format...');
    np.savetxt(plotfile + '.txt', np.c_[data['id_2MASS'],
               data['id_AllWISE'], dataset, labels], header=dataheading,
               delimiter='\t', fmt='%s')
    np.savetxt(plotfile + '_labeled.txt', np.c_[labeleddataset, 
               labeledlabels, labeleddata['z'], labeleddata['class'],
               labeleddata['subClass']],
               header=dataheading[20:]+'\tz\tclass\tsubClass',
               delimiter='\t', fmt='%s')
    print('Data file saved successfully, check your MeanShift folder')
        
    # Saving data in FITS format
    print('Saving data with labels in FITS format...');
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                 '.txt out=' + plotfile + '.fits');
    subprocess.run(bashorder, shell=True)
    bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                 '_labeled.txt out=' + plotfile + '_labeled.fits');
    subprocess.run(bashorder, shell=True)
    print('Data file saved successfully, check your KMeans folder')

    # Saving centers information
    plotfile = root + '/MeanShift/' + n + root_folder + 'Centers';
    print('Saving centers position...');
    np.savetxt(plotfile + '.txt', centers, header=dataheading[20:-6],
               delimiter='\t', fmt='%s')
    print('Centers positions saved successfully, check your MeanShift folder')
    
    print('\nMEAN SHIFT TECHNIQUE APPLIED\n');

# GMM
def func_gmm():
    print('\nCLUSTERING: GMM\n')
    K = [2,3,4,5,6,7,8,9,10]; #Possible choices for the number of clusters
    aic = []; #List to keep the AIC for each value in K
    bic = []; #List to keep the BIC for each value in K
    
    # Creating the model using the photometric data (for each K value)
    for i in range(len(K)):
        print('Fitting the model for K=%i...' % K[i])
        gmm = GaussianMixture(n_components=K[i], random_state=42);
        gmm.fit(dataused)
        print('GMM problem solved for K=%i' % K[i])
        
        # Model parameters
        means = gmm.means_;
        covariances = gmm.covariances_;
        aic.append(gmm.aic(dataused))
        bic.append(gmm.bic(dataused))
        print('POSITIONS OF THE MEANS:\n', means)
        print('COVARIANCES FOR EACH GROUP:\n', covariances)
        print('AKAIKE INFORMATION CRITERION:\n ', aic[i])
        print('BAYESIAN INFORMATION CRITERION:\n ', bic[i])
        
        # Adjusting the data to the model
        print('Obtaining the labels for K=%i...' % K[i])
        labels = gmm.predict(dataused);
        labeledlabels = gmm.predict(labeleddataused);
        
        # Saving data in ASCII format
        plotfile = (root + '/GMM/' + n + root_folder + root_file +
                    '_GMM_K=' + str(K[i]));
        print('Saving data with labels in ASCII format...');
        np.savetxt(plotfile + '.txt', np.c_[data['id_2MASS'],
                   data['id_AllWISE'], dataset, labels], header=dataheading,
                   delimiter='\t', fmt='%s')
        np.savetxt(plotfile + '_labeled.txt', np.c_[labeleddataset,
                   labeledlabels, labeleddata['z'], labeleddata['class'],
                   labeleddata['subClass']],
                   header=dataheading[20:]+'\tz\tclass\tsubClass',
                   delimiter='\t', fmt='%s')
        print('Data file saved successfully, check your GMM folder')
        
        # Saving data in FITS format
        print('Saving data with labels in FITS format...');
        bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                     '.txt out=' + plotfile + '.fits');
        subprocess.run(bashorder, shell=True)
        bashorder = ('sh stilts tcopy ifmt=ascii ofmt=fits in=' + plotfile +
                     '_labeled.txt out=' + plotfile + '_labeled.fits');
        subprocess.run(bashorder, shell=True)
        print('Data file saved successfully, check your GMM folder')
    
    # Residue plots
    print('Plotting the AIC and BIC data for each value of K...')
    plt.plot(K, aic,'-o', MarkerSize=2.5, Color=[0.75,0,0])
    plt.plot(K, bic,'--o', MarkerSize=2.5, Color=[0,0,0.5])
    plt.xlabel('K number', fontsize=18)
    plt.ylabel('Criteria', fontsize=18)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.175, right=0.925, top=0.925, bottom=0.125)
    plt.legend(['AIC','BIC'], fontsize=16, loc=0, handlelength=0.7);
    plt.tick_params(top=True, right=True, direction='inout')
    plt.ticklabel_format(axis='y', style='sci', scilimits=(0,0),
                         useMathText=True)
    plt.savefig(root + '/GMM/' + n + root_folder + 'GMMCriteria.eps')
    plt.savefig(root + '/GMM/' + n + root_folder +'GMMCriteria.png')
    plt.close()
    
    print('\nGMM TECHNIQUE APPLIED\n');


# Different options to execute
options = {1 : func_kmeans,
           2 : func_meanshift,
           3 : func_gmm,
}

# We execute the desired function
options[int(number)]()

