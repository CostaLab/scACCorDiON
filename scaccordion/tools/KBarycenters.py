import numpy as np 
import ot 
import pandas as pd

class KBarycenters:
    def __init__(self,k=3, max_iters=100,init='random',random_state=42):
        self.k = k
        self.max_iters = max_iters
        self.flabels = None
        self.bary = []
        self.centroids = None
        self.err = []
        self.seed = random_state
        self.iter = 0
        self.ini = None
    
    def fit(self, data, distr,cost,init='random',reg=1e-2):
        np.random.seed(self.seed)
        self._initialize_centroids(data,init)
        self.data = data
        for _ in range(self.max_iters):
            # Assign points to clusters
            self._assign_clusters()
            self.iter+=1
            # Update centroids
            new_centroids = self._update_centroids(data, distr,cost,reg=reg)
            # Check for convergence
            if np.all(self.flabels==np.argmin(new_centroids, axis=1)):
                break
            self.centroids = new_centroids
     
    def _initialize_centroids(self, data, init):
        indices = []
        if init=='random':
            indices = np.random.choice(len(data.columns), self.k, replace=False)
            self.centroids = data.iloc[:,indices]
        elif init == '++':
            indices = [np.random.choice(len(data.columns),1, replace=False)]
            for idi in range(self.k-1):
                for idj in indices:
                    indices.append(np.argmax(data[:,idj]))
            self.centroids = data.iloc[:,indices]
            self.ini = indices

    def _assign_clusters(self):
        self.flabels = np.argmin(self.centroids, axis=1)
       
    def _update_centroids(self, data, distr,cost,reg=1e-2):
            tmpbary = {}
            tmperr = {}
            tmpdist = {}
            for i in range(self.k):
                currg = self.flabels==i
                bary= ot.barycenter(A=distr.loc[:,currg].apply(lambda x: x/sum(x)).to_numpy(),
                                     M=cost,reg=reg,log=True)
                tmperr[i] = np.median(bary[1]['err'])
                tmpbary[i] = np.round(bary[0],4)/np.sum(np.round(bary[0],4))
                tmpdist[i]={}
                for idp in distr.columns:
                    tmpdist[i][idp] = ot.emd2(a=tmpbary[i],
                                              b=distr[idp].to_numpy()/distr[idp].sum(),
                                              M=cost)
            self.bary=tmpbary
            tmpdist=pd.DataFrame.from_dict(tmpdist)
            ## SSE
            self.err.append(np.sum([np.square(self.centroids[i[1]][self.flabels==i[0]]).sum() for i in enumerate(self.centroids.columns)]))
            return tmpdist
    