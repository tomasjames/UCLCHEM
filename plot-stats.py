#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
from plotfunctions import *
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from matplotlib.colors import LogNorm
import matplotlib.ticker as plticker
import matplotlib.patches as patches

# Select number of species
species = []
with open('Makerates/outputFiles/species.csv', 'r') as speciesData:
    reader = csv.reader(speciesData)
    for row in reader:
        species.append(row[0])

# Select dust
dust = []
with open('Makerates/inputFiles/uclgrainbasic.csv', 'r') as dustData:
    reader = csv.reader(dustData)
    for row in reader:
        dust.append(row[3])

species = ["H2O"]

# Define the velocities and densities used
velocities = np.linspace(5,15,11)
densities = np.logspace(6,3,4) # Requires reversal because of imshow ignoring origin command

# Predefine the arrays to be used later for the tick marks
densities_ticks = ["$10^{"+str(len(str(int(power)))-1)+"}$" for power in np.logspace(3,6,4)]

# Empty lists to house enhancement factors
enhance_j = np.array([[0]*len(velocities)]*len(densities))
enhance_c = np.array([[0]*len(velocities)]*len(densities))

# Loop through species
for specIndx,specName in enumerate(species):

    # Loop through those parameters defined above
    for nIndx,n in enumerate(densities):
        for vIndx,v in enumerate(velocities):

            print('Plotting: v{0}n{1}.pdf for {2}'.format(int(v),int(n),specName))

            patchTest = False

            #call read_uclchem.
            print([specName])
            ctime,cdens,ctemp,cabundances=read_uclchem('output/cshock/data/v{0}n{1}.dat'.format(int(v),int(n)),[specName])
            jtime,jdens,jtemp,jabundances=read_uclchem('output/jshock/data/v{0}n{1}.dat'.format(int(v),int(n)),[specName])

            # Determine time of sputtering
            t_sputter_j = jtime[next(x[0] for x in enumerate(jtemp) if x[1] > 100)]
            if np.max(ctemp) > 100:
                t_sputter_c = ctime[next(x[0] for x in enumerate(ctemp) if x[1] > 100)]
            else:
                # t_sputter_c = ((-15.38729*v*v*v)+(2069.56962*v*v)-(90272.826991*v)+1686858.54278)/(n*1e15)
                t_sputter_c = ctime[next(x[0] for x in enumerate(ctemp) if x[1] == np.max(ctemp))]
                print("By sputtering time, C: ",t_sputter_c)
            #print("t_sputter_j=",t_sputter_j)
            #print("t_sputter_c=",t_sputter_c)

            # Find location in arrays where this occurs
            sputter_j_index = next(x[0] for x in enumerate(jtime) if x[1] > t_sputter_j)
            sputter_c_index = next(x[0] for x in enumerate(ctime) if x[1] > t_sputter_c)

            print("sputter_j_index=",sputter_j_index)
            print("sputter_c_index=",sputter_c_index)

            # Determine whether the species in question is present on dust (and therefore relevant to sputter/desorb)
            if ("#"+specName in dust):
                print(specName+" is sputtered")
                sputter_j_abund = jabundances[0][sputter_j_index]
                sputter_c_abund = cabundances[0][sputter_c_index]
            else:
                print(specName+" is not sputtered")
                sputter_j_abund = jabundances[0][0]
                sputter_c_abund = cabundances[0][0]

            print("sputter_j_abund=",sputter_j_abund)            
            print("sputter_c_abund=",sputter_c_abund)

            # Block out data if below detectable limit
            if sputter_j_abund < 1.0e-12:
                patchTest = True
                print("patchTest = True for J")
                fig1=plt.figure(1)
                #plt.bar(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,5),np.linspace(float(nIndx)-0.5,float(nIndx)+0,5,5),hatch='x',alpha=0.5)
                #plt.gca().add_patch(plt.Rectangle((float(vIndx)-0.5,float(nIndx)-0.5),1,1,linewidth=0.5,hatch='///',edgecolor='k',alpha=1,fill=False))
                plt.fill(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,2),np.linspace(float(nIndx)-0.5,float(nIndx)+0.5,2),color='r',alpha=0.5)
                plt.fill(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,2),np.linspace(float(nIndx)-0.5,float(nIndx)+0.5,2),color='None',edgecolor='r',hatch='///')
            if sputter_c_abund < 1.0e-12:
                patchTest = True
                print("patchTest = True for C")
                fig2=plt.figure(2)
                #plt.bar(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,5),np.linspace(float(nIndx)-0.5,float(nIndx)+0.5,5),hatch='x',alpha=0.5)
                #plt.gca().add_patch(plt.Rectangle((float(vIndx)-0.5,float(nIndx)-0.5),1,1,linewidth=0.5,hatch='///',edgecolor='k',alpha=1,fill=False))
                plt.fill(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,2),np.linspace(float(nIndx)-0.5,float(nIndx)+0.5,2),color='r',alpha=0.5)
                plt.fill(np.linspace(float(vIndx)-0.5,float(vIndx)+0.5,2),np.linspace(float(nIndx)-0.5,float(nIndx)+0.5,2),color='None',edgecolor='r',hatch='///')

            # Compute enhancement ratios and assign to arrays
            enhance_ratio_j,time_j=enhancementRatio(jabundances[0],sputter_j_abund,jtime)
            enhance_ratio_c,time_c=enhancementRatio(cabundances[0],sputter_c_abund,ctime)

            print('nIndx: {0} \nvIndx: {1}'.format(nIndx,vIndx))
            print('v{0}n{1}.dat: '.format(v,n), enhance_ratio_j)
            enhance_j[nIndx][vIndx] = enhance_ratio_j
            enhance_c[nIndx][vIndx] = enhance_ratio_c
    
    # Normalise colour scales
    vmin = np.min([np.min(enhance_j),np.min(enhance_c)])
    vmax = np.max([np.max(enhance_j),np.max(enhance_c)])

    fig1=plt.figure(1)
    plt.imshow(enhance_j,origin='lower',norm=LogNorm(vmin=vmin, vmax=vmax),cmap='Blues')
    # plt.pcolor(velocities, np.log10(densities), enhance_j, norm=LogNorm(vmin=enhance_j.min(), vmax=enhance_j.max()))
    plt.colorbar(label='$f_{enhance}$',orientation='horizontal')

    # Set the xticks and yticks (including location)
    #plt.xticks(np.linspace(0,len(velocities),len(velocities)+1),velocities,fontsize=8)
    #plt.yticks(np.linspace(0,len(densities),len(densities)+1),densities_ticks,fontsize=8)
    plt.xticks(np.linspace(0,len(velocities),len(velocities)+1),velocities,fontsize=8)
    plt.yticks(np.linspace(0,len(densities),len(densities)+1),densities_ticks,fontsize=8)
    plt.minorticks_off() # Turns the minor ticks off
    plt.xlim([-0.5,len(enhance_j[0])-0.5])
    plt.ylim([-0.5,len(enhance_j[:,0])-0.5])

    if patchTest == True:
        #plt.gca().add_patch(plt.Rectangle((-2,-2),0.01,0.01,linewidth=0.5,alpha=0.5,hatch='x',facecolor='None',edgecolor='k',label=r"$\chi_{0} < 10^{-12}$"))
        #plt.gca().add_patch(plt.Rectangle((-1,-1),0.01,0.01,linewidth=0.5,facecolor='none',edgecolor='black',hatch='///',label=r"$\chi_{0} < 10^{-12}$"))
        plt.fill(np.linspace(-0.01,-0.02,2),np.linspace(-0.01,-0.02,2),color='None',edgecolor='r',alpha=0.5)
        plt.fill(np.linspace(-0.01,-0.02,2),np.linspace(-0.01,-0.02,2),color='r',edgecolor='r',hatch='///',alpha=0.5,label=r"$\chi_{0} < 10^{-12}$")

    # Plot grid lines between tick marks (hacky, but works)
    for x in range(0,len(velocities)):
        plt.plot(np.linspace(x-0.5,x-0.5,len(densities)+1),np.linspace(-0.5,len(densities),len(densities)+1),linestyle='-',color='k',linewidth=0.5)
    for y in range(0,len(densities)):
        plt.plot(np.linspace(-0.5,len(velocities),len(velocities)+1),np.linspace(y-0.5,y-0.5,len(velocities)+1),linestyle='-',color='k',linewidth=0.5)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False         # ticks along the top edge are off
    )

    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False         # ticks along the top edge are off
    )

    plt.xlabel('v (km s$^{-1}$)')
    plt.ylabel('n (cm$^{-3}$)')
    plt.title("J-shock: {0}".format(specName))
    plt.legend(loc='best',numpoints=1,bbox_to_anchor=(0.7,1.0))
    fig1.savefig('plots/enhance/{0}_j_enhance.pdf'.format(specName),dpi=200,bbox_inches='tight')
    plt.close(1)


    fig2=plt.figure(2)
    plt.imshow(enhance_c,origin='lower',norm=LogNorm(vmin=vmin, vmax=vmax),cmap='Blues')
    plt.colorbar(label='$f_{enhance}$',orientation='horizontal')

    # Plot grid lines between tick marks (hacky, but works)
    for x in range(0,len(velocities)):
        plt.plot(np.linspace(x-0.5,x-0.5,len(densities)+1),np.linspace(-0.5,len(densities),len(densities)+1),linestyle='-',color='k',linewidth=0.5)
    for y in range(0,len(densities)):
        plt.plot(np.linspace(-0.5,len(velocities),len(velocities)+1),np.linspace(y-0.5,y-0.5,len(velocities)+1),linestyle='-',color='k',linewidth=0.5)

    plt.tick_params(
        axis='x',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False         # ticks along the top edge are off
    )

    plt.tick_params(
        axis='y',          # changes apply to the x-axis
        which='minor',      # both major and minor ticks are affected
        bottom=False,      # ticks along the bottom edge are off
        top=False         # ticks along the top edge are off
    )

    # Set the xticks and yticks (including location)
    #plt.xticks(np.linspace(0.5,len(velocities)+0.5,len(velocities)+1),velocities,fontsize=8)
    #plt.yticks(np.linspace(0.5,len(densities)+0.5,len(densities)+1),densities_ticks,fontsize=8)
    plt.xticks(np.linspace(0,len(velocities),len(velocities)+1),velocities,fontsize=8) # Plots major ticks
    plt.yticks(np.linspace(0,len(densities),len(densities)+1),densities_ticks,fontsize=8) # Plots major ticks

    plt.minorticks_off() # Turns the minor ticks off
    plt.xlim([-0.5,len(enhance_j[0])-0.5])
    plt.ylim([-0.5,len(enhance_j[:,0])-0.5])

    if patchTest == True:
        #plt.gca().add_patch(plt.Rectangle((-2,-2),0.01,0.01,linewidth=0.5,hatch='x',facecolor='None',edgecolor='k',label=r"$\chi_{0} < 10^{-12}$"))
        #plt.gca().add_patch(plt.Rectangle((-1,-1),0.01,0.01,linewidth=0.5,facecolor='none',edgecolor='black',hatch='///',label=r"$\chi_{0} < 10^{-12}$"))
        plt.fill(np.linspace(-0.01,-0.02,2),np.linspace(-0.01,-0.02,2),color='None',edgecolor='r',alpha=0.5)
        plt.fill(np.linspace(-0.01,-0.02,2),np.linspace(-0.01,-0.02,2),color='r',edgecolor='r',hatch='///',alpha=0.5,label=r"$\chi_{0} < 10^{-12}$")

    plt.xlabel('v (km s$^{-1}$)')
    plt.ylabel('n (cm$^{-3}$)')
    plt.title("C-shock: {0}".format(specName))
    plt.legend(loc='best',numpoints=1,bbox_to_anchor=(0.7,1.0))
    fig2.savefig('plots/enhance/{0}_c_enhance.pdf'.format(specName),dpi=200,bbox_inches='tight')
    plt.close(2)

    
