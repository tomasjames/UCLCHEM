import uclchem
import numpy as np
from multiprocessing import Pool

#uclchem(final_dens,max_temp,shock_vel,phase_flag,outFile,startFile)
def runuclchem(a):
    modelNo,initial_dens,final_dens,max_temp,v,phase,outFile,startFile=a
    print("================")
    print("Model # ", modelNo)
    print("Phase ", phase)
    print("v=",v)
    print("n_init=",initial_dens)
    print("n_final=",final_dens)
    print("startFile=",startFile)
    print("outFile=",outFile)
    uclchem.uclchem(initial_dens,final_dens,max_temp,v,phase,outFile,startFile)

print(__name__)
if __name__ == '__main__':
    outputFolder="output/"
    # f=open(outputFolder+"grid.dat","wb")
    
    velocities = np.linspace(30,60,7)
    densities = np.logspace(3,6,4)

    max_temp=2000
    modelNo=1
    for n in densities:
        startFile=outputFolder+"start/{0}.dat".format(int(n))
        for v in velocities:
            outFile=outputFolder+"data/v{0}n{1}.dat".format(int(v),int(n))
            for phase in [1,2]:
                final_dens=n**phase

                if phase==1:
                    initial_dens=1.0e2
                else:
                    initial_dens=n

                runuclchem([modelNo,initial_dens,final_dens,max_temp,v,phase,outFile,startFile])
                modelNo+=1
    # f.close()
    # pool=Pool()
    # pool.map(runuclchem,models)
    # pool.close()
    # pool.join()