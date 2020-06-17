from parsers.textParser import textParser
from parsers.targetParser import distanceParser, dis2contact
from parsers.inputParsers import rawChainParser, msa2pssm
from utilities.dataSplit import trainTestSplit
from utilities.dicNormaliser import dicNormaliser
import os

import numpy as np
import gc

def main():
    makeFolder("./output")
    print("Parsing Text Files...")
    protInfo = textParser()

    print("Splitting Data...")
    train,test = trainTestSplit(protInfo)

    print("Load Train Set...")
    trainDic = loadSet(train)

    normCols = ["hsed","hseu","seqId","aligns"]
    normaliser = dicNormaliser(columns=normCols,coords=True)

    print("Normalise Train & Save...")
    trainDic = normaliser.fit_transform(trainDic)
    saveSet(trainDic)
    del trainDic
    gc.collect()

    print("Load Test Set...")
    testDic = loadSet(test)

    print("Normalise Test & Save...")
    testDic = normaliser.transform(testDic)
    saveSet(testDic)
    del testDic
    gc.collect()

    print("Load and Save Targets...")
    complexes = list(protInfo.keys())
    loadTargets(complexes)

def loadSet(setList,breakPoint=10000):
    master = {}
    entries = min(len(setList),breakPoint)
    i = 0

    for prot in setList:
        if prot[:-3] not in master:
            master[prot[:-3]] = {}
        chains = prot[-2:]
        if chains[0] not in master[prot[:-3]]:
            master[prot[:-3]][chains[0]] = featureExtract(prot, chains[0])
        if chains[1] not in master[prot[:-3]]:
            master[prot[:-3]][chains[1]] = featureExtract(prot, chains[1])

        i+=1
        print(f'\r>>Loading Set: {"{:.2f}".format(i / entries * 100)}%', end='')
        if (i>=breakPoint):
            print()
            return master
    print()
    return master

def loadTargets(complexes):
    end = len(complexes)
    i = 0

    makeFolder("./output/distance")
    makeFolder("./output/contact")

    for prot in complexes:
        filename = f"./PPI4DOCK/PPI4DOCK_docking_set/{prot}/{prot}_st.pdb"
        chains = prot[-2:]
        distance = distanceParser(filename,chains)
        contact = dis2contact(distance)
        with open(f"./output/distance/{prot}_dist.npy","wb") as f:
            np.save(f,distance)
        with open(f"./output/contact/{prot}_cont.npy","wb") as c:
            np.save(c,contact)

        i+=1
        print(f'\r>>Saving Targets: {"{:.2f}".format(i / end * 100)}%',end='')
    print()

def saveSet(dic):
    entries = len(dic)
    i = 0

    makeFolder("./output/chains")

    for prot in dic:
        for chain in dic[prot]:
            filename = f"{prot}_{chain}.npy"
            with open("./output/chains/"+filename,'wb') as f:
                np.save(f,dic[prot][chain].values)
        i += 1
        print(f'\r>>Saving Set: {"{:.2f}".format(i / entries * 100)}%',end='')
    print()

def featureExtract(target, chain):
    prefix = "./PPI4DOCK/PPI4DOCK_docking_set/"
    chainFile = f"{prefix}/{target}/{chain}_model_st.pdb"

    msaFile = f"./PPI4DOCK/PPI4DOCK_MSA/{target}/{target[:-2]}{chain}_coMSA.fasta"
    pssm = msa2pssm(msaFile)

    output = rawChainParser(chainFile, chain, pssm)
    return output

def makeFolder(path):
    if not os.path.exists(path):
        os.makedirs(path)

if __name__ == '__main__':
    main()