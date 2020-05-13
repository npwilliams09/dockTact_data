from parsers.textParser import textParser
from parsers.targetParser import contactParser, distanceParser, parseStrippedRes
from parsers.inputParsers import coordinateParser
import os
import numpy as np

def main():
    '''
    target = "1axi_BD"
    prefix = "./PPI4DOCK/PPI4DOCK_docking_set/"
    strippedRes = parseStrippedRes(prefix + f"{target}/stripped_res.txt")
    chains = target[-2]

    chainFile = f"{prefix}/{target}/{chains[0]}_model_st.pdb"
    output = coordinateParser(chainFile,strippedRes[0])
    '''
    generateSet("./PPI4DOCK/PPI4DOCK_docking_set/")


#takes in a folder and generates 2 inputs and 1 output target (i.e. one training example)
def generateSet(dir):
    i = 0
    files = next(os.walk(dir))
    with open("./output/protList.txt","w") as outList:
        for target in files[1]:
            targDir = dir + target + "/"
            chains = target[-2:]
            strippedRes = parseStrippedRes(f"{targDir}stripped_res.txt")

            chainA = coordinateParser(f"{targDir}{chains[0]}_model_st.pdb",strippedRes[0])
            chainB = coordinateParser(f"{targDir}{chains[1]}_model_st.pdb",strippedRes[1])
            complex = distanceParser(f"{targDir}{target}_st.pdb",strippedRes)

            makeFolder("./output/trainA")
            np.save(f"./output/trainA/{target}_A.npy",chainA)

            makeFolder("./output/trainB")
            np.save(f"./output/trainB/{target}_B.npy", chainB)

            makeFolder("./output/targets")
            np.save(f"./output/targets/{target}_target.npy", complex)

            outList.write(target + "\n")
            print(target)

            i+=1
            if(i>100):
                break



def makeFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == '__main__':
    main()