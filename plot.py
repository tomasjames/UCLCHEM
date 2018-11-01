#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
from plotfunctions import *
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use('TkAgg')

# Get list of folders (assuming same folders between c and j)
files = sorted(os.listdir('output/cshock/data/'))

# pick species, any number is fine
species = ["H2", "HCN", "HNCO"]

# Loop through those folders
for file in files:
    if file == '.DS_Store':
        continue

    file_name = file[:-4]

    if file == 'results':
        break

    print('Plotting: '+file+'.png')

    # Pulls the velocity and density from the folder name
    v = float(file[1:3])*100000 # Ensures v is in cm/s
    n = float(file[4:-4])

    #call read_uclchem.
    ctime,cdens,ctemp,cabundances=read_uclchem('output/cshock/data/{0}'.format(file),species)
    jtime,jdens,jtemp,jabundances=read_uclchem('output/jshock/data/{0}'.format(file),species)

    # if file == "v30n10000":
        # indx_lower = ctime.index(40)
        # indx_upper = ctime.index(50)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(40)
        # indx_upper = jtime.index(50)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(0.9)
        # indx_upper = jtime.index(1.0)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
    # elif file == "v30n100000":
        # indx_lower = ctime.index(3.0)
        # indx_upper = ctime.index(5.0)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(3.0)
        # indx_upper = jtime.index(5.0)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
    # elif file == "v40n10000":
        # indx_lower = ctime.index(35)
        # indx_upper = ctime.index(45)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(35)
        # indx_upper = jtime.index(45)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
    # elif file == "v40n100000":
        # indx_lower = ctime.index(3.5)
        # indx_upper = ctime.index(4.5)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(3.5)
        # indx_upper = jtime.index(4.5)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
    # elif file == "v50n10000":
        # indx_lower = ctime.index(40)
        # indx_upper = ctime.index(45)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(40)
        # indx_upper = jtime.index(45)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
    # elif file == "v50n100000":
        # indx_lower = ctime.index(3.5)
        # indx_upper = ctime.index(4.5)
        # del ctime[indx_lower:indx_upper]
        # del cdens[indx_lower:indx_upper]
        # del ctemp[indx_lower:indx_upper]
        # for abund in cabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(3.5)
        # indx_upper = jtime.index(4.5)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]
#
        # indx_lower = jtime.index(0.65)
        # indx_upper = jtime.index(0.7)
        # del jtime[indx_lower:indx_upper]
        # del jdens[indx_lower:indx_upper]
        # del jtemp[indx_lower:indx_upper]
        # for abund in jabundances:
            # del abund[indx_lower:indx_upper]

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
    fig.savefig('plots/'+file_name+"ctemp.png",dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    # ax.loglog(jtime,jtemp)
    ax.loglog(jtime,jtemp)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('T (K)')
    fig.savefig('plots/'+file_name+"jtemp.png",dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.loglog(ctime,cdens)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('n cm$^{-3}$')
    fig.savefig('plots/'+file_name+'cdens.png',dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)
    ax.loglog(jtime,jdens)
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel('$n$ cm$^{-3}$')
    fig.savefig('plots/'+file_name+'jdens.png',dpi=200)
    plt.close()

    fig=plt.figure()
    ax=fig.add_subplot(111)

    for specIndx,specName in enumerate(species):
        ax.loglog(ctime,cabundances[specIndx],label='C-shock: '+specName,linewidth='1.0')
        ax.loglog(jtime,jabundances[specIndx],label='J-shock: '+specName,linewidth='1.0')

    ax.plot(ctime,1e-13*np.ones(len(ctime)),linestyle="--")

    ax.legend(loc='best')
    ax.set_xlabel('t (yrs)')
    ax.set_ylabel("X$_{Species}$")
    ax.set_title("Shock propogating at v="+str(v)+" km/s through \n medium of pre-shock density n="+str(n)+" cm$^{-3}$")
    # ax.set_ylim([1e-13,1e-1])

    fig.savefig('plots/'+file_name+".png",dpi=200,bbox_inches='tight')
    plt.close()


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
