import pandas as pd
class dicNormaliser:
    def __init__(self, columns=None):
        self.columns = columns

    def fit(self, dic):
        mins = []
        maxs = []
        for prot in dic:
            for chain in dic[prot]:
                mins.append(dic[prot][chain].min())
                maxs.append(dic[prot][chain].max())
        max = pd.concat(maxs,axis=1,keys=[s.name for s in maxs).max(axis=1)

    def transform(self,dic):
        #
        return 0
