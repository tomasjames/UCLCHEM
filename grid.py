import uclchem
import numpy as np

def runuclchem(a):
    modelNo,initial_dens,final_dens,r_out,v,phase,outFile,startFile=a
    uclchem.uclchem(initial_dens,final_dens,r_out,v,phase,outFile,startFile)

print(__name__)
if __name__ == '__main__':
    outputFolder="output/"
    # f=open(outputFolder+"grid.dat","wb")
    
    velocities = np.linspace(5,5,1)
    densities = np.logspace(3,4,2)

    modelNo=1
    modelsOne = []
    for n in densities:
        initial_dens=1e2
        final_dens=n
        if n <= 1e2:
            r_out = 5.3
        elif 1e2 < n <= 1e3:
            r_out = 2.3
        elif 1e3 < n <= 1e4:
            r_out = 0.5
        elif 1e4 < n <= 1e5:
            r_out = 0.05
        else:
            r_out = 0.005
        print("At n={0}, r_out={1}".format(n,r_out))
        startFile=outputFolder+"start/{0}.dat".format(int(n))
        outFile=outputFolder+"start/{0}-out.dat".format(int(n))
        runuclchem([modelNo,initial_dens,final_dens,r_out,10,1,outFile,startFile])
        modelNo+=1

    modelNo=1
    modelsTwo=[]
    for n in densities:
        startFile=outputFolder+"start/{0}.dat".format(int(n))
        if n <= 1e2:
            r_out = 5.3
        elif 1e2 < n <= 1e3:
            r_out = 2.3
        elif 1e3 < n <= 1e4:
            r_out = 0.5
        elif 1e4 < n <= 1e5:
            r_out = 0.05
        else:
            r_out = 0.005
        for v in velocities:
            outFile=outputFolder+"data/v{0}n{1}.dat".format(int(v),int(n))
            
            initial_dens=n
            final_dens=n**2

            # modelsTwo.append([modelNo,initial_dens,final_dens,r_out,v,2,outFile,startFile])
            runuclchem([modelNo,initial_dens,final_dens,r_out,v,2,outFile,startFile])
            modelNo+=1

        pool=Pool()
        pool.map(runuclchem,modelsTwo)
        pool.close()
        pool.join()    
