#!/usr/bin/env python
# -*- coding: utf-8 -*-

import glob
import os
from plotfunctions import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# Get list of folders (assuming same folders between c and j)
files = sorted(os.listdir('output/jshock/data/'))
#files = []
#for file in os.listdir('output/cshock/data/'):
#    if file.startswith('v10'):
#        files.append(file)
#print(files)

# pick species, any number is fine
species = ["HNCO", "HCN", "HNC", "SIO", "E-", "SO", "SO2", "NH3", "CH3OH", "H2O", "#CH3OH"]
# species = ["#CH3OH"]
#species = ["CS", "OH"]
#species = ["HNCO", "HCN", "HNC"]
#species = ["SIO"]

# Loop through those folders
for file in files:
    if file == '.DS_Store':
        continue

    file_name = file[:-4]

    if file == 'results':
        break

    print('Plotting: '+file+'.pdf')

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

    # Determine shock distance
    d_cshock = [v*(t*(60*60*24*365)) for t in ctime]
    d_jshock = [v*(t*(60*60*24*365)) for t in jtime]

    v = int(v/100000) # Convert v back to km/s
    n = int(n)

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.loglog(ctime,ctemp)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('T (K)')
    fig.savefig('plots/'+file_name+"ctemp.pdf",dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    # ax.loglog(jtime,jtemp)
    ax.loglog(jtime,jtemp)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('T (K)')
    fig.savefig('plots/'+file_name+"jtemp.pdf",dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.loglog(ctime,cdens)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('n cm$^{-3}$')
    fig.savefig('plots/'+file_name+'cdens.pdf',dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.loglog(jtime,jdens)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('$n$ cm$^{-3}$')
    fig.savefig('plots/'+file_name+'jdens.pdf',dpi=200)
    plt.close()    

    for specIndx,specName in enumerate(species):

        if specName == "H2O":
            specNameLabel = "H_{2}O"
        elif specName == "H2":
            specNameLabel = "H_{2}"
        elif specName == "O2":
            specNameLabel = "O_{2}"
        elif specName == "SIO":
            specNameLabel = "SiO"
        elif specName == "SI":
            specNameLabel = "Si"
        elif specName == "NH3":
            specNameLabel = "NH_{3}"
        elif specName == "CH3OH":
            specNameLabel = "CH_{3}OH"
        elif specName == "#CH3OH":
            specNameLabel = "CH_{3}OH"
        else:
            specNameLabel = specName

        fig=plt.figure(figsize=(7,5), dpi=600)
        ax1=fig.add_subplot(311,xscale='log',yscale='log')

        if specNameLabel[0] == "#":
            ax1.loglog(ctime,cabundances[specIndx],label='C-shock: '+'#'+r'$'+specNameLabel+'$',linewidth='1.0')
            ax1.loglog(jtime,jabundances[specIndx],label='J-shock: '+'#'+r'$'+specNameLabel+'$',linewidth='1.0')
        else:
            ax1.loglog(ctime,cabundances[specIndx],label='C-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')
            ax1.loglog(jtime,jabundances[specIndx],label='J-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')

        #ax.plot(ctime,1e-13*np.ones(len(ctime)),linestyle="--")

        ax1.set_ylabel("X$_{Species}$")
        ax1.axvspan(2000, 4000, color='red', alpha=0.2, label="B2 age range", hatch='//')
        if np.min(jabundances[specIndx]) < 1e-13:        
            ax1.set_ylim([1e-13, np.max(jabundances[specIndx])*10])
        ax1.legend(loc='best')
        plt.setp(ax1.get_xticklabels(), visible=False)

        ax2=plt.subplot(312, sharex=ax1)	
        ax2.loglog(ctime,cdens,label="C-shock: $n$",linewidth='1.0')
        ax2.loglog(jtime,jdens,label="J-shock: $n$",linewidth='1.0')
        ax2.axvspan(2e3, 4e3, color='red', alpha=0.2, hatch='//')
        ax2.set_ylabel('$n$ (cm$^{-3}$)')
        ax2.legend(loc='upper left')
        plt.setp(ax2.get_xticklabels(), visible=False)

        ax3=plt.subplot(313, sharex=ax2)
        ax3.loglog(ctime,ctemp,label="C-shock: $T$",linewidth='1.0')
        ax3.loglog(jtime,jtemp,label="J-shock: $T$",linewidth='1.0')
        ax3.axvspan(2000, 4000, color='red', alpha=0.2, hatch='//')
        ax3.set_xlabel('t (yrs)')
        ax3.set_ylabel('T (K)')
        ax3.legend(loc='upper left')

        fig.suptitle("Shock propogating at v="+str(v)+" km s$^{-1}$ through medium of \n pre-shock density n="+str(n)+" cm$^{-3}$")
        # ax.set_xlim([1e-1,1e5])

        fig.savefig('plots/combined/'+specName+"_"+file_name+"_combined.pdf",bbox_inches='tight')
        plt.close()

        fig=plt.figure(dpi=600)
        ax1=fig.add_subplot(111,xscale='log',yscale='log')

        if specNameLabel[0] == "#":
            ax1.loglog(ctime,cabundances[specIndx],label='C-shock: '+'#'+r'$'+specNameLabel+'$',linewidth='1.0')
            ax1.loglog(jtime,jabundances[specIndx],label='J-shock: '+'#'+r'$'+specNameLabel+'$',linewidth='1.0')
        else:
            ax1.loglog(ctime,cabundances[specIndx],label='C-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')
            ax1.loglog(jtime,jabundances[specIndx],label='J-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')

#        ax1.loglog(ctime,cabundances[specIndx],label='C-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')
#        ax1.loglog(jtime,jabundances[specIndx],label='J-shock: '+r'$'+specNameLabel+'$',linewidth='1.0')

        #ax.plot(ctime,1e-13*np.ones(len(ctime)),linestyle="--")

        ax1.set_ylabel("X$_{Species}$")
        ax1.set_xlabel("t (yrs)")
        ax1.set_ylim([1e-13, np.max(jabundances[specIndx])*10])
        ax1.axvspan(2000, 4000, color='red', alpha=0.2, label="B2 age range", hatch='//')
        ax1.legend(loc='best')
        fig.suptitle("Fractional abundance for a shock propogating at v="+str(v)+" km s$^{-1}$ \n through medium of pre-shock density n="+str(n)+" cm$^{-3}$")

        fig.savefig('plots/'+specName+"_"+file_name+".pdf",bbox_inches="tight")
        plt.close()

    '''
    fig,ax=plt.subplots()

    for specIndx,specName in enumerate(species):
        if specName == "H2O":
            specName = "H_{2}O"
        elif specName == "H2":
            specName = "H_{2}"
        elif specName == "O2":
            specName = "O_{2}"
        elif specName == "SIO":
            specName = "SiO"
        elif specName == "SI":
            specName = "Si"
        elif specName == "NH3":
            specName = "NH_{3}"
        elif specName == "CH3OH":
            specName = "CH_{3}OH"
        # ax.loglog(ctime,cabundances[specIndx],label='C-shock: '+specName,linewidth='1.0')
        ax.loglog(jtime,jabundances[specIndx],label='J-shock: $'+specName+'$',linewidth=1,zorder=2)
        ax.set_xlabel('t (yrs)')
        ax.set_ylabel("X$_{Species}$")

        ax2=ax.twinx()
        ax2.loglog(jtime, jtemp, label=r"$T$", color='r',linestyle='--',linewidth=0.6,zorder=1,alpha=0.5)
        ax2.set_xlabel('t (yrs)')
        ax2.set_ylabel('Temperature (K)')
        # fontP = FontProperties()
        # fontP.set_size('small')
        ax.set_ylim([1e-13,1e-2])

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc='best')
        # legend = ax.legend(h1+h2, l1+l2, loc='best',frameon=False,ncol=len(species)+1,bbox_to_anchor=(1.0, -0.2))
        # legend.get_frame().set_facecolor('#ffffff')
        # legend.set_zorder(20)

    ax.set_title("J-shock propogating at v="+str(v)+" km/s through \n medium of pre-shock density n="+str(n)+" cm$^{-3}$")
    fig.savefig('plots/jcombined'+file_name+".png",dpi=400,bbox_inches='tight')
    plt.close()

    
    fig,ax=plt.subplots()

    for specIndx,specName in enumerate(species):
        # ax.loglog(ctime,cabundances[specIndx],label='C-shock: '+specName,linewidth='1.0')
        if specName == "H2O":
            specName = "H_{2}O"
        elif specName == "H2":
            specName = "H_{2}"
        elif specName == "O2":
            specName = "O_{2}"
        elif specName == "SIO":
            specName = "SiO"
        elif specName == "SI":
            specName = "Si"
        elif specName == "NH3":
            specName = "NH_{3}"
        ax.loglog(ctime,cabundances[specIndx],label='C-shock: $'+specName+'$',linewidth=1,zorder=2)
        ax.set_xlabel('t (yrs)')
        ax.set_ylabel("X$_{Species}$")

        ax2=ax.twinx()
        ax2.loglog(ctime, ctemp, label=r"$T$", color='r',linestyle='--',linewidth=0.6,zorder=1,alpha=0.5)
        ax2.set_xlabel('t (yrs)')
        ax2.set_ylabel('Temperature (K)')

        ax.set_ylim([1e-13,1e-3])

        h1, l1 = ax.get_legend_handles_labels()
        h2, l2 = ax2.get_legend_handles_labels()
        ax.legend(h1+h2, l1+l2, loc='best')
        # legend = ax.legend(h1+h2, l1+l2, loc='best',frameon=False,ncol=len(species)+1,bbox_to_anchor=(1.0, -0.2))
        # legend.get_frame().set_facecolor('#ffffff')
        # legend.set_zorder(20)

    ax.set_title("C-shock propogating at v="+str(v)+" km/s through \n medium of pre-shock density n="+str(n)+" cm$^{-3}$")
    fig.savefig('plots/ccombined'+file_name+".png",dpi=400,bbox_inches='tight')
    plt.close()
    '''

    for specIndx,specName in enumerate(species):

        # Call function to determine the maximum abundance
        max_abund,abund_time,max_shock=maxAbundDiff(cabundances[specIndx], jabundances[specIndx], ctime, jtime)
        fig1=plt.figure(1)
        ax_diff=fig.add_subplot(111)
        if max_shock == "C":
            ax_diff.loglog(max_abund,abund_time,"bx")
        elif max_shock == "J":
            ax_diff.loglog(max_abund,abund_time,"r+")

    ax_diff.set_xlabel('t (yrs)')
    ax_diff.set_ylabel('T (K)')
    ax_diff.axvspan(2000, 4000, color='red', alpha=0.2, label="B2 age range", hatch='//')
    fig1.savefig("plots/"+specName+"_max_abund_diff",dpi=200)
    plt.close(1)
