from plotfunctions import *
import matplotlib.pyplot as plt
import os

def analyseChemistry(resultFile,outFile,reactionFile,speciesFile,speciesName,parameterFile):
    species,masses=np.loadtxt(speciesFile,usecols=[0,1],dtype=str,skiprows=1,unpack=True,delimiter=',',comments="%")
    grains=[i for i,item in enumerate(species) if "#" in item]
    network=getNetwork(reactionFile,speciesName)
    cloud=getParameters(parameterFile)
    cloud['fr']=0.0
    times,dens,temp,specAbundances=read_uclchem(resultFile,[speciesName])

    print(len(times))

    times = times[1001::10]
    dens = dens[1001::10]
    temp = temp[1001::10]
    specAbundances = specAbundances[1001::10]

    oldMostForms=[]
    oldMostDestructs=[]
    plotTimes=[]
    oldTotalChange=0.0
    destructions=[]
    formations=[]
    for time in times:
        #print(time)
        timeStep,cloud,species,abundances=readTimestep(resultFile,time,cloud)
        abundances=np.asarray(abundances)
        cloud['mantle']=sum(abundances[grains])
        changes,reacIndxs=getChanges(speciesName,species,masses,abundances,network,cloud)#speciesName,species,abundances,network,temp,dens
        print(changes)
        print(reacIndxs)

        A=list(zip(changes,reacIndxs))
        A.sort()
        print(A)
        changes,reacIndxs=zip(*A)
        changes=np.asarray(changes)

        totalDestruct=sum(changes[np.where(changes<0)])
        destructions.append(-totalDestruct)
        totalProd=sum(changes[np.where(changes>0)])
        formations.append(totalProd)

        totalChange=sum(changes)
        mostForms=[]
        form=0.0
        i=-1
        while form < 0.99*totalProd:
            mostForms.append(reacIndxs[i])
            form+=changes[i]
            i-=1

        mostDestructs=[]    
        j=0
        destruct=0.0
        while abs(destruct) < 0.99*abs(totalDestruct):
            mostDestructs.append(reacIndxs[j])
            destruct+=changes[j]
            j+=1

        if set(oldMostDestructs)!=set(mostDestructs) or set(oldMostForms) !=set(mostForms):
            oldMostDestructs=mostDestructs[:]
            oldMostForms=mostForms[:]
            with open('analysis/'+speciesName+'_'+outFile+'_analysis.txt', 'a') as output_file:
                output_file.write("\n***************************\nNew Important Reactions At: {0:.2e}\n".format(time))
                output_file.write("\n")
                output_file.write("Formation = {0:.2e} from:".format(totalProd))
                output_file.write("\n")
                for k in range(-1,i,-1):
                    outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=network[reacIndxs[k]])
                    outString+=": {0:.2f}%".format(float(changes[k]/totalProd)*100)
                    output_file.write(outString)
                    output_file.write("\n")

                output_file.write("\nDestruction = {0:.2e} from:".format(totalDestruct))
                output_file.write("\n")
                for k in range(0,j):
                    outString="{x[0]} + {x[1]} -> {x[3]} + {x[4]}".format(x=network[reacIndxs[k]])
                    outString+=": {0:.2f}%".format(float(changes[k]/totalDestruct)*100)
                    output_file.write(outString)
                    output_file.write("\n")
                plotTimes.append(time)
                oldTotalChange=totalChange
                
        output_file.close()
        #resultFile.close()
    return times,abundances,destructions,formations

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
