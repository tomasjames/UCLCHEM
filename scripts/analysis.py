from analysisfunctions import analyseChemistry
from glob import glob
import numpy as np
from plotfunctions import *

################################################
#User Inputs Go Here
################################################

speciesNames=["HCN","HNCO","SO","SO2","H2O"]
resultFiles=glob("../output/jshock/data/*")
#resultFile="../output/jshock/data/v5n1000.dat"
#outFile="v5n1000"
reactionFile="../src/reactions.csv"
speciesFile="../src/species.csv"
parameterFile="../src/defaultparameters.f90"

################################################################################################
#NO CHANGES REQUIRED BELOW THIS LINE
################################################################################################

for resultFile in resultFiles:
    print("Analysing: "+resultFile)
    for speciesName in speciesNames:
        outFile = resultFile[22:-4]
        print(outFile)
        print("Reactions leading to the formation/destruction of: "+speciesName)
        times,specAbundances,destructions,formations=analyseChemistry(resultFile,outFile,reactionFile,speciesFile,speciesName,parameterFile)

        '''
        fig,ax=plt.subplots()
        #ax.plot(times,specAbundances[0],color="black")
        ax.loglog(times,destructions,label="Destruction",color="red")
        ax.loglog(times,formations,label="Formation",color="green")
        ax.set_xlabel('t (s)')
        ax.set_ylabel('Abundance')
        ax.set_title(speciesName+': '+outFile)
        #for time in times:
            #ax.axvline(time)
        #ax.set_yscale('log')
        ax.set_ylim(1e-25,1e-6)
        plt.savefig(speciesName+'_'+outFile)
        plt.close()
        '''