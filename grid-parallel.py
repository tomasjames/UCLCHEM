import uclchem
import numpy as np
from multiprocessing import Pool


def runuclchem(a):
    (modelNo, initial_dens, final_dens, v, phase, base_Av, r_out, outFile,
        startFile) = a
    uclchem.uclchem(
        initial_dens, final_dens, v, phase, base_Av, r_out, outFile, startFile
    )


print(__name__)
if __name__ == '__main__':
    outputFolder = "output/"

    velocities = np.linspace(35, 35, 1)
    densities = np.logspace(4, 4, 1)

    modelNo = 1
    modelsOne = []
    for n in densities:
        initial_dens = 1e2
        final_dens = n

        base_Av = 1.0
        # 3.086d18 is number of cm in 1 parsec
        r_out = (base_Av*(1.6e21)/(initial_dens))/(3.086e18)

        startFile = outputFolder+"start/{0}-final.dat".format(int(n))
        outFile = outputFolder+"start/{0}-collapse.dat".format(int(n))

        modelsOne.append(
            [modelNo, initial_dens, final_dens, 10, 1, base_Av,
                r_out, outFile, startFile]
        )
        modelNo += 1

    poolPhaseOne = Pool()
    poolPhaseOne.map(runuclchem, modelsOne)
    poolPhaseOne.close()
    poolPhaseOne.join()

    modelNo = 1
    modelsTwo = []
    for n in densities:
        startFile = outputFolder+"start/{0}-final.dat".format(int(n))
        for v in velocities:
            initial_dens = n
            final_dens = n**2

            if 1e2 < initial_dens <= 1e3:
                base_Av = 5.0
            elif 1e3 < initial_dens <= 1e4:
                base_Av = 10.0
            elif 1e4 < initial_dens <= 1e5:
                base_Av = 100.0
            else:
                base_Av = 1000.0

            print("base_Av=",base_Av)

            # 3.086d18 is number of cm in 1 parsec
            r_out = (base_Av*(1.6e21)/(initial_dens))/(3.086e18)

            outFile = outputFolder+"data/v{0}n{1}.dat".format(int(v), int(n))

            modelsTwo.append(
                [modelNo, initial_dens, final_dens, v, 2, base_Av,
                    r_out, outFile, startFile]
            )
            modelNo += 1

    poolPhaseTwo = Pool()
    poolPhaseTwo.map(runuclchem, modelsTwo)
    poolPhaseTwo.close()
    poolPhaseTwo.join()
