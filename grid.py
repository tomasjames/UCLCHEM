import uclchem
import numpy as np
from multiprocessing import Pool

def runuclchem(a):
    modelNo,initial_dens,final_dens,v,phase,outFile,startFile=a
    uclchem.uclchem(initial_dens,final_dens,v,phase,outFile,startFile)

print(__name__)
if __name__ == '__main__':
    outputFolder="output/"
    # f=open(outputFolder+"grid.dat","wb")
    
    velocities = np.linspace(30,35,2)
    densities = np.logspace(3,4,2)
    
    modelNo=1
    modelsOne = []
    for n in densities:
        initial_dens=1e2
        final_dens=n
        startFile=outputFolder+"blank.dat"
        outFile=outputFolder+"start/{0}.dat".format(int(n))
        runuclchem([modelNo,initial_dens,final_dens,10,1,outFile,startFile])
        modelNo+=1

    pool=Pool()
    pool.map(runuclchem,modelsOne)
    pool.close()
    pool.join()

    modelNo=1
    modelsTwo=[]
    for n in densities:
        startFile=outputFolder+"start/{0}.dat".format(int(n))
        for v in velocities:
            outFile=outputFolder+"data/v{0}n{1}.dat".format(int(v),int(n))
            
            initial_dens=n
            final_dens=n**2

            modelsTwo.append([modelNo,initial_dens,final_dens,v,2,outFile,startFile])
            modelNo+=1

    pool=Pool()
    pool.map(runuclchem,modelsTwo)
    pool.close()
    pool.join()    