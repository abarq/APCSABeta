import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def eliminar_outliers(matriz):
    nFilas, nColumnas = matriz.shape
    rowToDelete = np.zeros(nFilas)
    n=0
    for i in range(nFilas):
        if matriz[i,1] >= 1:
            rowToDelete[n] = i
            n = n+1
    rowToDelete = rowToDelete[rowToDelete!=0]
    rowToDelete = rowToDelete.astype(int)
    matriz = np.delete(matriz, rowToDelete,0)
    return matriz

data = pd.read_csv('output.txt',sep= ' ', header=None)
matriz = data.to_numpy()

matriz = matriz[matriz[:, 0].argsort()]
print(matriz)

'''nCasos = max(matriz[:,1])
nCorridas_porCaso = max(matriz[:,2])

datosAgrupados = np.zeros((int(nCasos),int(nCorridas_porCaso)))

for i in range(int(nCasos*nCorridas_porCaso)):
    datosAgrupados[int(matriz[i,1])-1,int(matriz[i,2])-1] = matriz[i,3]'''
