from parsers.textParser import textParser
from parsers.targetParser import distanceParser, dis2contact
from parsers.inputParsers import rawChainParser, msa2pssm, adjacencyMat
from utilities.dataSplit import trainTestSplit, surfaceResidues
from utilities.dicNormaliser import dicNormaliser
import os
import numpy as np
import gc
import pandas as pd

def main():
    makeFolder("./output")
    print("Parsing Text Files...")
    protInfo = textParser()

    print("Filtering Proteins Over Size Threshold...")
    protInfo = sizeFilter(protInfo,250000)

    print("Splitting Data...")
    train,test = trainTestSplit(protInfo)

    print("Load Train Set...")
    trainDic = loadSet(train)
    loadTargets(train,trainDic)

    normCols = ["hsed","hseu","seqId","aligns",'res_depth','ca_depth']
    normaliser = dicNormaliser(columns=normCols,coords=False)

    print("Normalise Train & Save...")
    trainDic = normaliser.fit_transform(trainDic)
    saveSet(trainDic)
    del trainDic
    gc.collect()

    print("Load Test Set...")
    testDic = loadSet(test)
    loadTargets(train, trainDic)

    print("Normalise Test & Save...")
    testDic = normaliser.transform(testDic)
    saveSet(testDic)
    del testDic
    gc.collect()

    print("Done, Output saved to folder :)")

def loadSet(setList,breakPoint=100000):
    master = {}
    entries = min(len(setList),breakPoint)
    makeFolder("./output/chains")
    i = 0
    surfaceDict = {}

    for prot in setList:
        if prot[:-3] not in master:
            master[prot[:-3]] = {}
        chains = prot[-2:]
        if chains[0] not in master[prot[:-3]]:
            df = featureExtract(prot, chains[0])
            df, ids = surfaceResidues(df)
            adjMat = adjacencyMat(prot,chains[0],ids)
            makeFolder(f"./output/chains/{prot[:-3]}_{chains[0]}")
            with open(f"./output/chains/{prot[:-3]}_{chains[0]}/adjMat.npy", "wb") as f:
                np.save(f, adjMat)
            master[prot[:-3]][chains[0]] = df
        if chains[1] not in master[prot[:-3]]:
            df = featureExtract(prot, chains[1])
            df, ids = surfaceResidues(df)
            adjMat = adjacencyMat(prot, chains[1], ids)
            makeFolder(f"./output/chains/{prot[:-3]}_{chains[1]}")
            with open(f"./output/chains/{prot[:-3]}_{chains[1]}/adjMat.npy", "wb") as f:
                np.save(f, adjMat)
            master[prot[:-3]][chains[1]] = df
        i+=1
        print(f'\r>>Loading Set: {"{:.2f}".format(i / entries * 100)}%', end='')
        if (i>=breakPoint):
            print()
            return master
    print()
    return master

def sizeFilter(protInfo,size=250000):
    newDic = {}
    for x in protInfo.keys():
        combinedAAs = int(protInfo[x]["nb_AA1"])*int(protInfo[x]["nb_AA2"])
        if (combinedAAs < size):
            newDic[x] = protInfo[x]
    return newDic

def loadTargets(complexes, protDic):
    end = len(complexes)
    i = 0

    makeFolder("./output/distance")
    makeFolder("./output/contact")

    for prot in complexes:
        filename = f"./PPI4DOCK/PPI4DOCK_docking_set/{prot}/{prot}_st.pdb"
        chains = prot[-2:]
        seqIDs = []
        for chain in chains:
            df = protDic[prot[:-3]][chain]
            seqIDs.append([int(x) for x in df["seqId"]])

        distance = distanceParser(filename,chains,seqIDs)
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
            path = f"./output/chains/{prot}_{chain}/input.ft"
            data = dic[prot][chain]
            data.to_feather(path)
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