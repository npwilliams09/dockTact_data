from sklearn.model_selection import train_test_split

def trainTestSplit(protInfo):
    x = []
    categories = []
    for key in protInfo:
        x.append(key)
        categories.append(protInfo[key]["category"])

    train, test, train_cat, test_cat = train_test_split(x,categories,stratify=categories,random_state=2020,train_size=0.8,test_size=0.2)

    with open('./output/trainList.txt', 'w') as f: #write train to file
        for item in train:
            f.write("%s\n" % item)

    with open('./output/testList.txt', 'w') as f: #write test to file
        for item in test:
            f.write("%s\n" % item)

    return train,test

def surfaceResidues(df,threshold=10.0):
    df = df[df["res_depth"] < threshold]
    seqIDs = [int(x) for x in list(df["seqId"])]
    return df,seqIDs