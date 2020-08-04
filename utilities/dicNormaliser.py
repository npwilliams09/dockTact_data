import pandas as pd
class dicNormaliser:
    def __init__(self, columns, coords = False):
        self.columns = columns
        self.coords = coords

    def fit(self, dic):
        mins = []
        maxs = []
        if (self.coords):
            cols = self.columns + ['x','y','z']
        else:
            cols = self.columns
        for prot in dic:
            for chain in dic[prot]:
                mins.append(dic[prot][chain][cols].min())
                maxs.append(dic[prot][chain][cols].max())
        min = pd.concat(mins,axis=1,keys=[s.name for s in mins]).min(axis=1)
        max = pd.concat(maxs,axis=1,keys=[s.name for s in maxs]).max(axis=1)
        self.range = max - min
        self.min = min

    def transform(self,dic):
        for prot in dic:
            for chain in dic[prot]:
                df = dic[prot][chain]
                df[self.columns] = df[self.columns].sub(self.min[self.columns])
                df[self.columns] = df[self.columns].div(self.range[self.columns])

                if (self.coords):
                    coordRange = self.range[['x','y','z']].max()
                    df[['x','y','z']] = df[['x','y','z']]-(coordRange/-2)
                    df[['x','y','z']] = df[['x','y','z']]/(coordRange)
                    df[['x','y','z']] = (df[['x','y','z']] * 2) - 1
                dic[prot][chain] = df
        return dic

    def fit_transform(self,dic):
        self.fit(dic)
        return self.transform(dic)