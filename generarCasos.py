import numpy as np
import pandas as pd


markov_length_levels = np.array([100000, 500000])
neigh_size_levels = np.array([1500, 2000])
step_reduction_levels = np.array([1.0, 1.05])
nRuns_perCase = 10

data = pd.read_csv('seedExp', header=None)
seeds = data.to_numpy()


nCases = step_reduction_levels.size * neigh_size_levels.size * markov_length_levels.size
nTotal_runs = nRuns_perCase * nCases

f=open('input.txt','a')

matrizDatos = np.zeros((nTotal_runs,8))
randomArray = np.random.rand(nTotal_runs,1)

nCaso = 1
nCorrida = 1
for i in range(markov_length_levels.size):
    for j in range(neigh_size_levels.size):
        for k in range(step_reduction_levels.size):
            for m in range(nRuns_perCase):
                matrizDatos[nCorrida-1,:] = [nCorrida, nCaso, m+1, seeds[m], markov_length_levels[i], neigh_size_levels[j], step_reduction_levels[k], randomArray[nCorrida-1]]#corrida
                nCorrida = nCorrida + 1
            nCaso = nCaso + 1


matrizDatos = matrizDatos[matrizDatos[:, 7].argsort()]
np.savetxt(f, matrizDatos[:,0:7],fmt='%1.2i %1.2i %1.2i %1.2i %1.2i %1.2f %1.2f')
