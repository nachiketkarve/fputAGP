# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 15:43:14 2024

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("F:\\FPUT\\python\\fput")
import fputAlpha as fput
import matplotlib.animation as ani

N = 8
E0 = 1.0
k0 = 1
alpha = 0.0

A,V,C,FourierComponents = fput.initialize(N)

FileName = "./../data/alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
df = pd.read_csv(FileName)
period = df["period"].values.tolist()[0]
AGPgrad = np.array(df["AGPgrad"].values.tolist())

AGPgradp = np.array([0])
AGPgradq = np.array([0])

AGPgradp = np.append(AGPgradp,AGPgrad[N:2*N])
AGPgradp = np.append(AGPgradp,0)
AGPgradq = np.append(AGPgradq,AGPgrad[0:N])
AGPgradq = np.append(AGPgradq,0)
AGPgradP = np.matmul(FourierComponents,AGPgradp)
AGPgradQ = np.matmul(FourierComponents,AGPgradq)

print(AGPgradP[2])

w1 = fput.frequency(1,N)
w2 = fput.frequency(2,N)

gradP2 = np.sqrt(2.0/(N+1))*E0*(w2**2 - 2.0*w1**2)/(w2**3 - 4.0*w2*w1**2)
print(gradP2)

print(AGPgradQ)