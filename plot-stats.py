#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
from plotfunctions import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# Select number of species
species = ["HNCO", "HCN", "HNC", "SIO", "SO", "SO2", "NH3", "CH3OH", "H2O"]

# Get list of folders (assuming same folders between c and j)
files = sorted(os.listdir('output/jshock/data/'))

enhance = []
params = []

for specIndx,specName in enumerate(species):

    # Loop through those folders
    for file in files:
        if not file.startswith("v10"):
            continue
        if file == '.DS_Store':
            continue

        file_name = file[:-4]

        if file == 'results':
            break

        print('Plotting: '+file+'.pdf for '+specName)

        # Pulls the velocity and density from the folder name
        if file[2] == 'n':
            v =  float(file[1:2])*100000
            n = float(file[3:-4])
        else:
            v = float(file[1:3])*100000 # Ensures v is in cm/s
            n = float(file[4:-4])

        #call read_uclchem.
        ctime,cdens,ctemp,cabundances=read_uclchem('output/cshock/data/{0}'.format(file),species)
        jtime,jdens,jtemp,jabundances=read_uclchem('output/jshock/data/{0}'.format(file),species)

        v = int(v/100000) # Convert v back to km/s
        n = int(n)

        ratio,time=enhancementRatio(cabundances[specIndx],ctime)
        enhance.append(ratio)
        params.append([v,n])

        # Call function to determine the maximum abundance
        max_abund,abund_time,max_shock=maxAbundDiff(cabundances[specIndx], jabundances[specIndx], ctime, jtime)
        fig1=plt.figure(1)
        ax_diff=fig1.add_subplot(111)
        if max_shock == "C":
            ax_diff.loglog(max_abund,abund_time,"x",label="$n=$"+str(n))
        elif max_shock == "J":
            ax_diff.loglog(max_abund,abund_time,"+",label="$n=$"+str(n))

    fig1.legend(loc='best')
    ax_diff.set_xlabel('t (yrs)')
    ax_diff.set_ylabel(r'$\Delta \chi_{max}$'+ ' (K)')
    fig1.savefig("plots/"+specName+"_v"+str(v)+"_max_abund_diff.pdf",dpi=200)
    plt.close(1)
