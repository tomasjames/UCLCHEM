from analysisfunctions import analyseChemistry

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