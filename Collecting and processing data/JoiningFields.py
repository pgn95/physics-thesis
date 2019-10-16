# =====================================================================
# Script that takes a table with positions and joins the tables that
# overlap in the sky. The script is complex, and each part perform a
# different task. Each part should be executed following the logical
# order.
# ---------------------------------------------------------------------
# Author: Pablo Gómez Nicolás
# Year: 2018/19
# ---------------------------------------------------------------------
# Documentation: https://matplotlib.org/tutorials/index.html
#                https://docs.astropy.org/en/stable/
#             http://www.star.bris.ac.uk/~mbt/stilts/sun256/sun256.html
# ---------------------------------------------------------------------
# Input (when asked by the program):
#   * File, in FITS format, that contains the position (right ascension
#     and declination) of the search region centers. The header of the
#     right ascension column must be RA, and DEC for declination. It
#     also contains in the first column (FOV) the number of the row.
#   * Stilts is needed!
#   * get_area.pl and MOCorder2pixel.topcat are needed for obtaining
#     the areas!
#
# Output:
#   * Table, in FITS format, with the information of the groups.
#
# 0.* Tables, in FITS format, with the coordinates for each group and
#     individual field.
#
# 1.* MOC archives that contain the area information for each group and
#     individual field.
#
# 2.* Tables, in ASCII format, with the area for each individual field
#     and each group.
#
# 3.* Joined 2MASS and AllWISE tables, in FITS format, for each group.
#
# 4.* Tables, in ASCII format, with the number of sources for each
#     individual field and each group, both in 2MASS and AllWISE.
#
# 5.* Tables, in ASCII format, with the density for each individual
#     field and each group.
#
# 6.* Plots, in PNG and EPS formats, with the relation between the
#     2MASS and AllWISE densities, before and after cleansing the data.
#
# 7.* Joined 2MASS and AllWISE tables, in FITS format, for each big 
#     group.
#
# 8.* The area of each group, written in the terminal.
# =====================================================================

# Useful packages
import sys
import numpy as np
import subprocess
from astropy.io import fits
import matplotlib.pyplot as plt

# START OF THE SCRIPT
print('\nSCRIPT RUNNING...\n')

# Filename with the positions needed
print('Please write the name of the file that contains the positions of the '
      'search region centers (with extension)')
filename = input();

# The program is divided in different parts
print('Type 0 if you have not created the separated coordinate files')
print('Type 1 if you have created the separated coordinate files but not the '
      'MOC files')
print('Type 2 if you have created the MOC files but not the areas files')
print('Type 3 if you have created the area files but not joined the tables')
print('Type 4 if you have created the joined tables but not the files with the '
      'number of sources')
print('Type 5 if you have both tables with areas and sizes and you want to '
      'create the sources density table')
print('Type 6 if you have the density table and want to create the plots')
print('Type 7 if you want to join all the data in different groups')
print('Type 8 if you want to know the area of each group')
ind = input();

# Creating a new table with the group division
print('Trying to create the new table with groups...')
try:
    order = ('sh stilts tmatch1 in=' + filename + ' matcher=sky values=\'RA '
             'DEC\' params=\'1800\' action=identify ocmd=\'sort "GroupID"\' '
             'ocmd=\'addcol Single "NULL_GroupSize"\' ocmd=\'addcol Double '
             '"GroupSize==2"\' ocmd=\'addcol Triple "GroupSize==3"\' ocmd='
             '\'addcol Multiple "GroupSize>3"\' ocmd=\'addcol Group '
             '"GroupSize>1"\' out=CoordinateList_Groups.fits');
    subprocess.run(order, shell=True)
except:
    print('Warning: Check that the file is in the path and that the format '
          'is correct')
    exit(1)

hdulist = fits.open('CoordinateList_Groups.fits');
data = hdulist[1].data;
ind_individual = data['Single']; #True when the row is not in a bigger group
indfov = data['FOV'][ind_individual];
ind_group = data['Group']; #True when the row is in a bigger group
groupname = data['GroupID'][ind_group];
groupname = list(set(list(groupname))); #Removing the repeated numbers

if ind == '0':
    # Getting the tables with the coordinates for each field, both individual
    # and groups
    # INDIVIDUAL FIELDS
    print('Creating tables with position information for each group...')
    for i in range(len(indfov)):
        order = ('sh stilts tpipe in=CoordinateList_Groups.fits omode=out cmd='
                 '\'select "FOV==' + str(indfov[i]) + '"\' out=Coordinates_'
                 'Individual' + str(indfov[i]) + '.fits');
        subprocess.run(order, shell=True)
    # GROUPED FIELDS
    for i in range(len(groupname)):
        order = ('sh stilts tpipe in=CoordinateList_Groups.fits omode=out cmd='
                 '\'select "GroupID==' + str(groupname[i]) + '"\' out='
                 'Coordinates_Group' + str(groupname[i]) + '.fits');
        subprocess.run(order, shell=True)
    print('Separated coordinate files created successfully')
    
elif ind == '1':
    # Getting MOC archives, describing the field area
    print('Getting MOC archives...')
    # INDIVIDUAL FIELDS
    for i in range(len(indfov)):
        order = ('sh stilts pixfoot in=Coordinates_Individual' + str(indfov[i])
                 + '.fits ifmt=fits ra=RA dec=DEC radius=0.25 order=16 '
                 'mocfmt=json out=Individual_' + str(indfov[i]) + '.moc');
        subprocess.run(order, shell=True)
    # GROUPED FIELDS
    for i in range(len(groupname)):
        order = ('sh stilts pixfoot in=Coordinates_Group' + str(groupname[i]) +
                 '.fits ifmt=fits ra=RA dec=DEC radius=0.25 order=16 '
                 'mocfmt=json out=Group_' + str(groupname[i]) + '.moc');
        subprocess.run(order, shell=True)
    print('MOC files created successfully')

elif ind == '2':
    # Getting the files with the area information, in ASCII format
    print('Creating the area archives for the group and individual fields...')
    # INDIVIDUAL FIELDS
    order = 'echo "#No.\tArea" >> IndividualAreas.txt';
    subprocess.run(order, shell=True)
    for i in range(len(indfov)):
        order = ('echo "' + str(indfov[i]) + '\t$(perl get_area.pl Individual_' 
                 + str(indfov[i]) + '.moc)" >> IndividualAreas.txt');
        subprocess.run(order, shell=True)
    # GROUP FIELDS
    order = 'echo "#No.\tArea" >> GroupAreas.txt';
    subprocess.run(order, shell=True)
    for i in range(len(groupname)):
        order = ('echo "' + str(groupname[i]) + '\t$(perl get_area.pl Group_' +
                 str(groupname[i]) + '.moc)" >> GroupAreas.txt');
        subprocess.run(order, shell=True)
    print('Area files created successfully')

elif ind == '3':
    # Joining the 2MASS and AllWISE tables in the groups obtained
    print('Joining the individual tables in groups...')
    for i in range(len(groupname)):
        groupind = data['GroupID'];
        groupfov = data['FOV'][groupind == groupname[i]];
        order = 'sh stilts tcat ifmt=fits ';
        for j in range(len(groupfov)):
            order = order + 'in=2MASSTable' + str(groupfov[j]) + '.tbl ';
        order = (order + 'ofmt=fits out=2MASSTableGroup' + str(groupname[i]) +
                 '.tbl');
        subprocess.run(order, shell=True)
        order = ('sh stilts tmatch1 matcher=exact values=designation action='
                 'keep1 ifmt=fits in=2MASSTableGroup' + str(groupname[i]) +
                 '.tbl ofmt=fits out=2MASSTableGroup' + str(groupname[i]) +
                 '.tbl');
        subprocess.run(order, shell=True)
        order = 'sh stilts tcat ifmt=fits ';
        for j in range(len(groupfov)):
            order = order + 'in=AllWISETable' + str(groupfov[j]) + '.tbl ';
        order = (order + 'ofmt=fits out=AllWISETableGroup' + str(groupname[i]) +
                 '.tbl');
        subprocess.run(order, shell=True)
        order = ('sh stilts tmatch1 matcher=exact values=designation action='
                 'keep1 ifmt=fits in=AllWISETableGroup' + str(groupname[i]) + 
                 '.tbl ofmt=fits out=AllWISETableGroup' + str(groupname[i]) +
                 '.tbl');
        subprocess.run(order, shell=True)
    print('Group files created successfully')

elif ind == '4':
    # Creating the files with the number of sources
    print('Creating the files with the number of sources...')
    # INDIVIDUAL FIELDS
    order = 'echo "#No.\tSize" >> Individual2MASSSizes.txt';
    subprocess.run(order, shell=True)
    order = 'echo "#No.\tSize" >> IndividualAllWISESizes.txt';
    subprocess.run(order, shell=True)
    for i in range(len(indfov)):
        order = ('echo "' + str(indfov[i]) + '\t$(sh stilts tpipe in=2MASSTable'
                 + str(indfov[i]) + '.tbl omode=count)" >> Individual2MASSSizes'
                 '.txt');
        subprocess.run(order, shell=True)
        order = ('echo "' + str(indfov[i]) + '\t$(sh stilts tpipe in='
                 'AllWISETable' + str(indfov[i]) + '.tbl omode=count)" >> '
                 'IndividualAllWISESizes.txt');
        subprocess.run(order, shell=True)
    # GROUP FIELDS
    order = 'echo "#No.\tSize" >> Group2MASSSizes.txt';
    subprocess.run(order, shell=True)
    order = 'echo "#No.\tSize" >> GroupAllWISESizes.txt';
    subprocess.run(order, shell=True)
    for i in range(len(groupname)):
        order = ('echo "' + str(groupname[i]) + '\t$(sh stilts tpipe in='
                 '2MASSTableGroup' + str(groupname[i]) + '.tbl omode=count)" >>'
                 ' Group2MASSSizes.txt');
        subprocess.run(order, shell=True)
        order = ('echo "' + str(groupname[i]) + '\t$(sh stilts tpipe in='
                 'AllWISETableGroup' + str(groupname[i]) + '.tbl omode=count)" '
                 '>> GroupAllWISESizes.txt');
        subprocess.run(order, shell=True)
    print('Sizes files created successfully')

elif ind == '5':
    # Loading the files with the information
    print('Loading all the necessary information...')
    data_indareas = np.loadtxt('IndividualAreas.txt');
    data_groupareas = np.loadtxt('GroupAreas.txt');
    # We order the areas
    fov = data_indareas[:,0];
    indareas = data_indareas[:,1];
    indareas = indareas[np.argsort(fov)];
    fov = data_groupareas[:,0];
    groupareas = data_groupareas[:,1];
    groupareas = groupareas[np.argsort(fov)];
    # Number of sources in each field
    data_2mass = np.loadtxt('Individual2MASSSizes.txt', usecols=(0,4));
    data_allwise = np.loadtxt('IndividualAllWISESizes.txt', usecols=(0,4));
    fov = data_2mass[:,0];
    indsizes_2mass = data_2mass[:,1];
    indsizes_2mass = indsizes_2mass[np.argsort(fov)];
    fov = data_allwise[:,0];
    indsizes_allwise = data_allwise[:,1];
    indsizes_allwise = indsizes_allwise[np.argsort(fov)];
    data_2mass = np.loadtxt('Group2MASSSizes.txt', usecols=(0,4));
    data_allwise = np.loadtxt('GroupAllWISESizes.txt', usecols=(0,4));
    fov = data_2mass[:,0];
    groupsizes_2mass = data_2mass[:,1];
    groupsizes_2mass = groupsizes_2mass[np.argsort(fov)];
    fov = data_allwise[:,0];
    groupsizes_allwise = data_allwise[:,1];
    groupsizes_allwise = groupsizes_allwise[np.argsort(fov)];
    # Saving data in ASCII format
    print('Saving density data table...')
    dataheading = 'No.\t2MASS\tAllWISE';
    np.savetxt('IndividualDensity.txt', np.c_[np.sort(data_indareas[:,0]), 
               indsizes_2mass/indareas, indsizes_allwise/indareas],
               header=dataheading, delimiter='\t', fmt='%s')
    np.savetxt('GroupDensity.txt', np.c_[np.sort(data_groupareas[:,0]), 
               groupsizes_2mass/groupareas, groupsizes_allwise/groupareas],
               header=dataheading, delimiter='\t', fmt='%s')
    print('Data file saved successfully, check your folder')

elif ind == '6':
    # Creating the plots
    print('Plotting the 2MASS and AllWISE density...')
    inddensity = np.loadtxt('IndividualDensity.txt');
    groupdensity = np.loadtxt('GroupDensity.txt');
    density_2mass = np.concatenate((inddensity[:,1], groupdensity[:,1]));
    density_allwise = np.concatenate((inddensity[:,2], groupdensity[:,2]));
    plt.plot(density_2mass, density_allwise, 'ro', MarkerSize=2.5,
             MarkerEdgeColor=[0.75,0,0], MarkerFaceColor=[0.75,0,0])
    plt.xlabel('2MASS Density', fontsize=18)
    plt.ylabel('AllWISE Density', fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.15)
    plt.savefig('Density_Dirty.eps')
    plt.savefig('Density_Dirty.png')
    plt.close()
    print('All plots created successfully')
    print('Choose the AllWISE threshold to eliminate the spureous points')
    threshold = input(); threshold = int(threshold); #11300
    # Plotting the circles in the original figure
    circles_2mass = density_2mass[density_allwise < threshold];
    circles_allwise = density_allwise[density_allwise < threshold];
    plt.plot(density_2mass, density_allwise, 'ro', MarkerSize=2.5,
             MarkerEdgeColor=[0.75,0,0], MarkerFaceColor=[0.75,0,0])
    plt.plot(circles_2mass, circles_allwise, 'ro', MarkerSize=15,
             MarkerEdgeColor=[0,0,0.5], MarkerFaceColor='none')
    plt.xlabel(r'$\rho_{2MASS}/(Sources\cdotº^{-2})$', fontsize=18)
    plt.ylabel(r'$\rho_{AllWISE}/(Sources\cdotº^{-2})$', fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.15)
    plt.savefig('Density_Dirty.eps')
    plt.savefig('Density_Dirty.png')
    plt.close()
    density_2mass = density_2mass[density_allwise >= threshold];
    density_allwise = density_allwise[density_allwise >= threshold];
    plt.plot(density_2mass, density_allwise, 'ro', MarkerSize=2.5,
             MarkerEdgeColor=[0.75,0,0], MarkerFaceColor=[0.75,0,0])
    plt.xlabel(r'$\rho_{2MASS}/(Sources\cdotº^{-2})$', fontsize=18)
    plt.ylabel(r'$\rho_{AllWISE}/(Sources\cdotº^{-2})$', fontsize=18)
    plt.xscale('log')
    plt.yscale('log')
    # Deciding how to divide the points in different groups
    print('If you want to plot a limit to divide the points in groups, please '
          'write 0. If not, write anything else')
    indicator = input();
    while indicator == '0':
        print('Write the 2MASS threshold')
        threshold = input(); threshold = int(threshold); #2200, 4000
        plt.axvline(x=threshold, c='k', linestyle='--')
        print('If you want to plot another limit to divide the points in '
              'groups, please write 0. If not, write anything else')
        indicator = input();
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.subplots_adjust(left=0.175, right=0.95, top=0.95, bottom=0.15)
    plt.savefig('Density_Clean.eps')
    plt.savefig('Density_Clean.png')
    plt.close()
    print('All plots created successfully')

elif ind == '7':
    # Joining the tables in the groups
    threshold = [];
    print('Please choose the 2MASS points of separation. To finish the input, '
          'write 0')
    thresholdnum = input(); #2200, 4000
    while thresholdnum != '0':
        threshold.append(int(thresholdnum));
        thresholdnum = input();
    threshold = np.sort(threshold);
    print('Choose the AllWISE threshold to eliminate the spureous points')
    allwisethreshold = input(); allwisethreshold = int(allwisethreshold); #11300
    indiv = np.loadtxt('IndividualDensity.txt');
    group = np.loadtxt('GroupDensity.txt');
    bashorder2mass = 'sh stilts tcat ifmt=fits ';
    bashorderallwise = bashorder2mass;
    for i in range(len(indiv)):
        if indiv[i,1] <= threshold[0] and indiv[i,2] >= allwisethreshold:
            bashorder2mass = (bashorder2mass + 'in=2MASSTable' +
                              str(int(indiv[i,0])) + '.tbl ');
            bashorderallwise = (bashorderallwise + 'in=AllWISETable' +
                              str(int(indiv[i,0])) + '.tbl ');
    for i in range(len(group)):
        if group[i,1] <= threshold[0] and group[i,2] >= allwisethreshold:
            bashorder2mass = (bashorder2mass + 'in=2MASSTableGroup' +
                              str(int(group[i,0])) + '.tbl ');
            bashorderallwise = (bashorderallwise + 'in=AllWISETableGroup' +
                              str(int(group[i,0])) + '.tbl ');
    bashorder2mass = bashorder2mass + 'ofmt=fits out=2MASSSet1.tbl';
    bashorderallwise = bashorderallwise + 'ofmt=fits out=AllWISESet1.tbl';
    subprocess.run(bashorder2mass, shell=True)
    subprocess.run(bashorderallwise, shell=True)
    for j in range(len(threshold)-1):
        bashorder2mass = 'sh stilts tcat ifmt=fits ';
        bashorderallwise = bashorder2mass;
        for i in range(len(indiv)):
            if (indiv[i,1] > threshold[j] and indiv[i,1] <= threshold[j+1] and
            indiv[i,2] >= allwisethreshold):
                bashorder2mass = (bashorder2mass + 'in=2MASSTable' +
                                  str(int(indiv[i,0])) + '.tbl ');
                bashorderallwise = (bashorderallwise + 'in=AllWISETable' +
                                  str(int(indiv[i,0])) + '.tbl ');
        for i in range(len(group)):
            if (group[i,1] > threshold[j] and group[i,1] <= threshold[j+1] and
            group[i,2] >= allwisethreshold):
                bashorder2mass = (bashorder2mass + 'in=2MASSTableGroup' +
                                  str(int(group[i,0])) + '.tbl ');
                bashorderallwise = (bashorderallwise + 'in=AllWISETableGroup' +
                                  str(int(group[i,0])) + '.tbl ');
        bashorder2mass = (bashorder2mass + 'ofmt=fits out=2MASSSet' + str(j+2) +
                          '.tbl');
        bashorderallwise = (bashorderallwise + 'ofmt=fits out=AllWISESet' +
                            str(j+2) + '.tbl');
        subprocess.run(bashorder2mass, shell=True)
        subprocess.run(bashorderallwise, shell=True)
    bashorder2mass = 'sh stilts tcat ifmt=fits ';
    bashorderallwise = bashorder2mass;
    for i in range(len(indiv)):
        if indiv[i,1] > threshold[-1] and indiv[i,2] >= allwisethreshold:
            bashorder2mass = (bashorder2mass + 'in=2MASSTable' +
                              str(int(indiv[i,0])) + '.tbl ');
            bashorderallwise = (bashorderallwise + 'in=AllWISETable' +
                              str(int(indiv[i,0])) + '.tbl ');
    for i in range(len(group)):
        if group[i,1] > threshold[-1] and group[i,2] >= allwisethreshold:
            bashorder2mass = (bashorder2mass + 'in=2MASSTableGroup' +
                              str(int(group[i,0])) + '.tbl ');
            bashorderallwise = (bashorderallwise + 'in=AllWISETableGroup' +
                              str(int(group[i,0])) + '.tbl ');
    bashorder2mass = (bashorder2mass + 'ofmt=fits out=2MASSSet' +
                     str(len(threshold)+1) + '.tbl');
    bashorderallwise = (bashorderallwise + 'ofmt=fits out=AllWISESet' +
                     str(len(threshold)+1) + '.tbl');
    subprocess.run(bashorder2mass, shell=True)
    subprocess.run(bashorderallwise, shell=True)

elif ind == '8':
    # Calculating the area of each group, in steradians
    threshold = [];
    print('Please choose the 2MASS points of separation. To finish the input, '
          'write 0')
    thresholdnum = input(); #2200, 4000
    while thresholdnum != '0':
        threshold.append(int(thresholdnum));
        thresholdnum = input();
    threshold = np.sort(threshold);
    print('Choose the AllWISE threshold to eliminate the spureous points')
    allwisethreshold = input(); allwisethreshold = int(allwisethreshold); #11300
    indiv = np.loadtxt('IndividualDensity.txt');
    group = np.loadtxt('GroupDensity.txt');
    indiv = indiv[np.argsort(indiv[:,0])];
    group = group[np.argsort(group[:,0])];
    indivarea = np.loadtxt('IndividualAreas.txt');
    grouparea = np.loadtxt('GroupAreas.txt');
    indivarea = indivarea[np.argsort(indivarea[:,0])];
    grouparea = grouparea[np.argsort(grouparea[:,0])];
    area = 0;
    for i in range(len(indiv)):
        if indiv[i,1] <= threshold[0] and indiv[i,2] >= allwisethreshold:
            area = area + indivarea[i,1];
    for i in range(len(group)):
        if group[i,1] <= threshold[0] and group[i,2] >= allwisethreshold:
            area = area + grouparea[i,1];
    print('Set 1: ', area*(np.pi/180.0)**2) #0.0607481326560753
    for j in range(len(threshold)-1):
        area = 0;
        for i in range(len(indiv)):
            if (indiv[i,1] > threshold[j] and indiv[i,1] <= threshold[j+1] and
            indiv[i,2] >= allwisethreshold):
                area = area + indivarea[i,1];
        for i in range(len(group)):
            if (group[i,1] > threshold[j] and group[i,1] <= threshold[j+1] and
            group[i,2] >= allwisethreshold):
                area = area + grouparea[i,1];
        print('Set ' + str(j+2) + ': ', area*(np.pi/180.0)**2) #0.0427892571047376
    area = 0;
    for i in range(len(indiv)):
        if indiv[i,1] > threshold[-1] and indiv[i,2] >= allwisethreshold:
            area = area + indivarea[i,1];
    for i in range(len(group)):
        if group[i,1] > threshold[-1] and group[i,2] >= allwisethreshold:
            area = area + grouparea[i,1];
    print('Set ' + str(len(threshold)+1) + ': ', area*(np.pi/180.0)**2) 
    #0.01523939830397774

print('\nSCRIPT FINISHED\n')
