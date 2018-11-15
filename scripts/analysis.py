from analysisfunctions import analyseChemistry

from __future__ import print_function

import numpy as np

from plotfunctions import *

################################################
#User Inputs Go Here
################################################

speciesNames=["HCN","H2O"]
resultFile="../output/jshock/data/v5n1000.dat"
outFile="v5n1000"
reactionFile="../src/reactions.csv"
speciesFile="../src/species.csv"
parameterFile="../src/defaultparameters.f90"

################################################################################################
#NO CHANGES REQUIRED BELOW THIS LINE
################################################################################################

for speciesName in speciesNames:
	times,specAbundances,destructions,formations=analyseChemistry(resultFile,outFile,reactionFile,speciesFile,speciesName,parameterFile)

    fig,ax=plt.subplots()
    # ax.plot(times,specAbundances[0],color="black")
    ax.plot(times,destructions,color="red")
    ax.plot(times,formations,color="green")
    ax.set_xlabel('t (s)')
    # ax.set_ylabel(' (s)')
    ax.set_title(speciesName)
    for time in plotTimes:
        ax.axvline(time)
    ax.set_yscale('log')
    ax.set_ylim(1e-25,1e-6)
    plt.savefig(outFile)
    plt.close()
'''



	keepReactions=[]
	for reaction in reactions:
		if speciesName in reaction:
			keepReactions.append(reaction)


for i in range(0,len(species)):
	if '#' in species[i]:
		cloud['mantle']=cloud['mantle']+abundances[i]
	if species=='H':
		cloud['h']=abundances[i]


changes,reacIndxs=getChanges("CO",species,masses,abundances,network,cloud)#speciesName,species,abundances,network,temp,dens

A=zip(changes,reacIndxs)
A.sort()
changes,reacIndxs=zip(*A)
changes=np.asarray(changes)
totalDestruct=sum(changes[np.where(changes<0)])
totalProd=sum(changes[np.where(changes>0)])
print totalProd,totalDestruct

for i in range(0,5): 
	print reactions[reacIndxs[i]], "{0:.1f}%".format(100.0*changes[i]/totalDestruct)
for i in range(-1,-6,-1): 
	print reactions[reacIndxs[i]], "{0:.1f}%".format(100.0*changes[i]/totalProd)'''
