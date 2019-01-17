#A set of functions for working with UCLCHEM outputs
# adding "from plotfunctions import * to any python script in scripts/ will allow their use"
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import csv
from scipy.interpolate import interp1d

#to read from uclchem's full output send the filename of the output and a list of species names to read_uclchem()
#the list for example would be ["H","C","CO","#CO"]
#the return is a list of times and list of lists of abundances
#call time,abundances=read_uclchem("output-full",["H","C","CO","#CO"])
# abudnance[0] would be a list of "H" abundances
def read_uclchem(filename,species):
    a=open(filename).read()
    a=a.split('\n')
    #there are 68 lines per time step in the UCL_CHEM output file.
    lines=68
    timesteps=len(a)/lines

    abunds=[]
    for spec in species:
        abunds.append([])
    time=[]
    dens=[]
    temp=[]
    #so now do the following until end of file
    with open(filename) as file:
        for line in file:
            bits=line.split()
            #find time line
            if  'age' in bits:
                time.append(float(bits[-2].replace('D','E')))
            #read another line for dens
            if 'density' in bits:
                densi=float(bits[-2].replace('D','E'))
                if densi==0.0:
                    densi=1e-10
                dens.append(densi)
            if "temp" in bits:
                tempi=float(bits[-2].replace('D','E'))
                temp.append(tempi)
            #then read until we hit abundances
            if bits .count('=')>=1:
                for specIndx,specName in enumerate(species):
                    if specName in bits:
                        #print(specName," is in bits")
                        abunds[specIndx].append(float(bits[2+bits.index(specName)].replace('D','E')))

    return time,dens,temp,abunds

def readTimestep(filename,chosenTime,cloud):
    a=open(filename).read()
    a=a.split('\n')
    #there are 68 lines per time step in the UCL_CHEM output file.
    lines=68
    timesteps=len(a)/lines

    abunds=[]
    species=[]
    specPosArray=[0,3,6,9]
    #so now do the following until end of file
    readFlag=False
    with open(filename) as file:
        for line in file:
            bits=line.split()        
            #find time line
            if  'age' in bits:
                if chosenTime==float(bits[-2].replace('D','E')):
                    time=chosenTime
                    readFlag=True
                else:
                    if readFlag:
                        break
                #read another line for dens
            if readFlag:
                if 'density' in bits:
                    densi=float(bits[-2].replace('D','E'))
                    if densi==0.0:
                        densi=1e-10
                    cloud['density']=densi
                elif 'temperature' in bits:
                   cloud['temp']=float(bits[-2].replace('D','E'))
                elif 'extinction' in bits:
                   cloud['av']=float(bits[-2].replace('D','E')) 
                elif 'cosmic' in bits:
                   cloud['zeta']=float(bits[-4].replace('D','E'))
                #then read until we hit abundances
                if bits .count('=')>3:
                    for specPos in specPosArray:
                        species.append(bits[specPos])
                        abunds.append(float(bits[2+specPos].replace('D','E')))
                        if bits[specPos]=="H":
                            cloud['h']=float(bits[2+specPos].replace('D','E'))
    return time,cloud,species,abunds 

def write_cols(filename,times,dens,abundances):
    f=open(filename,"wb")
    for timeIndx,time in enumerate(times):
        outString="{0:.3e} {1:.3e}".format(time,dens[timeIndx])
        for i in range(0,len(abundances)):
            outString+=" {0:.3e}".format(abundances[i][timeIndx])
        outString+="\n"
        f.write(outString)
    f.close()

#send a  list of species names and their abundances in  a list of lists
#with a list of times for each abundance point
#same as the species input and time/abundance output from read_uclchem
#optionally send an output filename to save the plot
#return ax,figure for further manipulation
def plot_species(species,times,abundances,plotFile=None):
    fig=plt.figure()
    ax=fig.add_subplot(111)
    colours=make_colours(len(species))

    for specIndx,specName in enumerate(species):
        ax.plot(times,abundances[specIndx],color=colours.next(),label=specName)

    ax.legend(loc=4,fontsize='small')

    ax.set_xlabel('Time / years')
    ax.set_ylabel("X$_{Species}$")

    ax.set_yscale('log')

    if plotFile is not None:
        fig.savefig(plotFile)
    return ax,fig
    


def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

def make_colours(n):
    return iter(cm.rainbow(np.linspace(0, 1, n)))

def formatSpecies(speciesName):
    speciesName=speciesName.upper()
    speciesName=speciesName.replace("+","$^+$")
    speciesName=speciesName.replace("2","$_2$")
    speciesName=speciesName.replace("3","$_3$")
    speciesName=speciesName.replace("4","$_4$")
    return speciesName


#########################################################################################################
#Analsis Functions
#########################################################################################################

def getRate(reactype,mass,alpha,beta,gamma,cloud):
    if reactype=='CRP':
        rate = alpha*cloud['zeta']

    elif reactype == 'PHOTON':
        #I ignored self-shielding
        #print alpha
        rate = alpha*np.exp(-gamma*cloud['av'])*cloud['radfield']/1.7

    elif reactype=='CRPHOT':
        rate=alpha*gamma*1.0/(1.0-cloud['omega'])*cloud['zeta']*(cloud['temp']/300.0)**beta
   
    elif reactype=='FREEZE':             
        if (cloud['evap'] != 0 or cloud['fr'] == 0.0):
            rate=0.0
        else:    
            if (beta==0.0 ): 
                rate=4.57e4*alpha*np.sqrt(cloud['temp']/mass)*cloud['grainArea']*cloud['fr']
            else:
                cion=1.0+16.71e-4/(cloud['radg']*cloud['temp'])
                rate=4.57e4*alpha*np.sqrt(cloud['temp']/mass)*cloud['grainArea']*cloud['fr']*cion
                    
    elif reactype=='DESOH2':
            if (cloud['desorb'] and cloud['h2desorb']):
                if (gamma < cloud['ebmaxh2'] and cloud['mantle']>1e-30):
                    rate = cloud['epsilon']*1.0e-17*np.sqrt(cloud['temp'])*cloud['h']*1.0/cloud['mantle']
                else:
                    rate=0.0
            else:
                rate = 0.0
            
    elif reactype=='DESCR':
            if (cloud['desorb'] and cloud['crdesorb']):
                if (cloud['mantle'] > 1e-30 and gamma < cloud['ebmaxcr']): 
                    rate = 4.0*np.pi*cloud['zeta']*1.64e-4*(cloud['grainArea'])*(1.0/cloud['mantle'])*cloud['phi']
                else:
                    rate=0.0
            else:
                rate = 0.0
            
    elif reactype=='DEUVCR':
        if (cloud['desorb'] and cloud['uvcr']):
            if( gamma < cloud['ebmaxuvcr'] and cloud['mantle'] > 1.0e-30):
                rate = cloud['grainArea']*cloud['uvy']*4.875e3*cloud['zeta']*(1.0/cloud['mantle'])
                rate = rate *(1.0+(cloud['radfield']/cloud['uvcreff'])*(1.0/cloud['zeta'])*np.exp(-1.8*cloud['av']))
            else:
                rate = 0.0
        else:
            rate=0.0
    else:
        rate = alpha*((cloud['temp']/300.)**beta)*np.exp(-gamma/cloud['temp'])

    return rate

def getChanges(speciesName,species,masses,abundances,network,cloud):
    reacIndxs=[]
    changes=[]
    for i in range(0,len(network)):
        reacIndxs.append(i)
        reacIndx=list(network[i]).index(speciesName)
        rate=getRate(network[i][1],float(masses[species.index(speciesName)]),network[i][-3],network[i][-2],network[i][-1],cloud)
        reactantCount=0
        reactants=network[i][0:3]
        for reactant in reactants:
            if reactant in species:
                indx=species.index(reactant)
                rate*=abundances[indx]
                reactantCount+=1
        if reactantCount>1:
            rate*=cloud['density']
        else:
            if reactants[1]=='FREEZE':
                rate*=cloud['density']
            elif reactants[1]=='DESOH2':
                rate*=cloud['density']
            elif reactants[1]=='electr':
                rate*=cloud['density']
        if reacIndx<3:
            changes.append(-rate)
        else:
            changes.append(rate)
    return changes,reacIndxs

def getNetwork(file,speciesName):
    print(file)
    reactions=np.loadtxt("../src/reactions.csv",dtype=str,skiprows=1,delimiter=',',usecols=
        [0,1,2,3,4,5,6],comments="%")
    alpha,beta,gamma=np.loadtxt(file,usecols=[7,8,9],unpack=True,skiprows=1,delimiter=',',comments="%")
    network=[]
    for i,reaction in enumerate(reactions):
        if speciesName in reaction:
            keep=list(reaction)
            keep.extend([alpha[i],beta[i],gamma[i]])
            network.append(keep)

    return network

def getParameters(file):
    cloud={'temp':10,
    'density':1e2,
    'zeta':1.0,
    'av':1.0,
    'radfield':1.0,
    'fr':1.0,
    'evap':0,
    'mantle':0.0,

    'desorb':True,
    'crdesorb':True,
    'uvcr':True,
    'h2desorb':True,

    'omega':0.5,
    'phi':1.0e5,
    'grainArea':2.4e-22,
    'radg':1.0e-5,
    'ebmaxcr':1.21e3,
    'ebmaxh2':1.21e3,
    'ebmaxcrf':1.21e3,
    'ebmaxuvcr':1.0e4,
    'uvcreff':1.0e-3,
    'epsilon':0.01,
    'uvy':0.1,
    'h':0.0}

    with open(file,"r") as inFile:
        for line in inFile.readlines():
            for bits in line.split(';'):
                #print bits
                vals=bits.split('=')
                if len(vals)>1:
                    #print vals[0].lower()
                    if vals[0] in cloud.keys():
                        if 'desorb' in vals[0] or vals[0]=='uvcr':
                            cloud[vals[0]]=(int(vals[1])==1)
                        else:
                            cloud[vals[0]]=float(vals[1].replace('d','e'))
    return cloud 


def averageAbundance(abundances):
    """
    A function to determine the average abundance 
    Inputs:
        abundances (arr): an array of abundances corresponding to abundance at
        as a function of time
    Outputs:
        average (float): the average abundance 
    """

    average = np.mean(abundances)

    return average


def maxAbundDiff(c_abundances, j_abundances, c_times, j_times):
    """
    A function to determine the maximum abundance and the time at which such
    maximum occurs
    Inputs:
        c_abundances (arr): an array of abundances corresponding to C-shock
        abundance as a function of time
        c_abundances (arr): an array of abundances corresponding to J-shock
        abundance as a function of time
        c_times (arr): the times that correspond to the array c_abundances
        j_times (arr): the times that correspond to the array j_abundances
    Outputs:
        max_abund (float): the maximum abundance difference
        abund_time (float): the time at which max_abund occurs at
        max_shock (str): the shock with greater abundance at abund_time
    """

    # Interpolate between points to ensure same number of elements in each array
    if np.size(j_abundances) > np.size(c_abundances):
        # Changes the shape of the array to be compatible with interp1d
        c_abundances = np.array(c_abundances).reshape(np.size(c_abundances)) 
        c_times = np.array(c_times).reshape(np.size(c_times))

        abundance_interp = interp1d(c_times, c_abundances, fill_value='extrapolate')
        c_abundances = abundance_interp(j_times)
        times = j_times
    elif np.size(c_abundances) > np.size(j_abundances):
        # Changes the shape of the array to be compatible with interp1d
        j_abundances = np.array(j_abundances).reshape(np.size(j_abundances)) 
        j_times = np.array(j_times).reshape(np.size(j_times))

        abundance_interp = interp1d(j_times, j_abundances, fill_value='extrapolate')
        c_abundances = abundance_interp(c_times)
        print(np.size(c_abundances))
        times = c_times

    # Determine the difference between each abundance point 
    abundance_diff = abs(j_abundances - c_abundances)

    # Determine which shock produces greater abundance
    if np.max(j_abundances - c_abundances) < 0:
        max_shock = "C"
    elif np.max(j_abundances - c_abundances) > 0:
        max_shock = "J"
    else:
        max_shock = "None"

    # Find the maximum abundance difference and the time that this occurs
    max_abund = np.max(abundance_diff)
    # [time for times, abund in enumerate(abundance_diff) if abund == max_abund]
    # time = times[np.where(abundance_diff == max_abund)[0][0]]
    time = times[list(abundance_diff).index(max_abund)]

    return max_abund, time, max_shock


def enhancementRatio(abundances,times):

    max_abund = np.max(abundances)
    time = times[list(abundances).index(max_abund)]

    enhance_ratio = max_abund/abundances[0]

    return enhance_ratio, time