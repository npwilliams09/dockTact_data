from parsers.textParser import textParser
from parsers.targetParser import distanceParser, dis2contact
from parsers.inputParsers import rawChainParser, msa2pssm, adjacencyMat
from utilities.dataSplit import trainTestSplit, surfaceResidues
from utilities.dicNormaliser import dicNormaliser
import os
import numpy as np
import gc
from threading import Thread
from multiprocessing import Pool
import pandas as pd
import traceback

def main():
    Threads = 8

    makeFolder("./output")
    print("Parsing Text Files...")
    protInfo = textParser()

    print("Filtering Proteins Over Size Threshold...")
    protInfo = sizeFilter(protInfo,250000)

    print("Splitting Data...")
    train,test = trainTestSplit(protInfo)

    print("Load Train Set...")
    trainDic = threaded_process_range(train,Threads)
    loadTargets(train,trainDic)

    normCols = ["hsed","hseu","seqId","aligns",'res_depth','ca_depth']
    normaliser = dicNormaliser(columns=normCols,coords=False)

    print("Normalise Train & Save...")
    trainDic = normaliser.fit_transform(trainDic)
    saveSet(trainDic)
    del trainDic
    gc.collect()

    print("Load Test Set...")
    testDic = threaded_process_range(test,Threads)
    loadTargets(train, trainDic)

    print("Normalise Test & Save...")
    testDic = normaliser.transform(testDic)
    saveSet(testDic)
    del testDic
    gc.collect()

    print("Done, Output saved to folder :)")

def loadProt(prot):
    master = {}
    makeFolder(f"./output/{prot}")

    chains = prot[-2:]

    master[chains[0]] = loadChain(prot,chains[0])
    master[chains[1]] = loadChain(prot, chains[1])

    return master

def loadChain(prot,chain):
    df = featureExtract(prot, chain)
    df, ids = surfaceResidues(df)
    adjMat = adjacencyMat(prot, chain, ids)
    with open(f"./output/{prot}/{chain}_adjMat.npy", "wb") as f:
        np.save(f, adjMat)
    return df

def process_range(id_range):
    """process a number of ids, storing the results in a dict"""
    try:
        store = {}
        for i, id in enumerate(id_range):
            store[id] = loadProt(id)
            print(f'\r>>Loading Set: {"{:.2f}".format(i+1 / len(id_range) * 100)}%', end='')
        return store
    except:
        traceback.print_exc()
        raise

def threaded_process_range(id_range, nthreads=1):
    """process the id range in a specified number of threads"""
    pool = Pool(processes=nthreads)

    ids = []
    # create the threads
    for i in range(nthreads):
        ids.append(id_range[i::nthreads])

    results = pool.map(process_range,ids)

    dic = {}
    for r in results:
        dic.update(r)

    return dic

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

    for prot in complexes:
        filename = f"./PPI4DOCK/PPI4DOCK_docking_set/{prot}/{prot}_st.pdb"
        chains = prot[-2:]
        seqIDs = []
        for chain in chains:
            df = protDic[prot][chain]
            seqIDs.append([int(x) for x in df["seqId"]])

        distance = distanceParser(filename,chains,seqIDs)
        contact = dis2contact(distance)
        with open(f"./output/{prot}/dist.npy","wb") as f:
            np.save(f,distance)
        with open(f"./output/{prot}/cont.npy","wb") as c:
            np.save(c,contact)

        i+=1
        print(f'\r>>Saving Targets: {"{:.2f}".format(i / end * 100)}%',end='')
    print()

def saveSet(dic):
    entries = len(dic)
    i = 0

    for prot in dic:
        for chain in dic[prot]:
            path = f"./output/{prot}/{chain}_input.ft"
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