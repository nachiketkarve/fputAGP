# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:51:57 2024

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("/usr3/graduate/nachiket/FPUT/todaFixedPython")
import toda as br
import matplotlib.animation as ani

def TotalEnergy(q,p,alpha):
    N = len(q)-2
    En = 0
    for i in range(N+1):
        En = En + p[i]**2/2.0 + 1.0/(2.0*alpha)**2*(np.exp(2.0*alpha*(q[i+1]-q[i])) - 2.0*alpha*(q[i+1]-q[i]) - 1.0)

    return En

N = 2
E0 = 1.0
k0 = 1
alphaMax = 4
alpha = 0.01
deltaAlpha = 0.01
frames = 500

figName = ".//..//plot//_agpGrad_fputTodaBr_N" + str(int(N)) + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".pdf"

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

A,V,C,FourierComponents = br.initialize(N)

alphas = np.array([])
gradMag = np.array([])
gradMagX = np.array([])
gradMagP = np.array([])

alphas2 = np.array([])
gradMag2 = np.array([])
gradMagX2 = np.array([])
gradMagP2 = np.array([])

while alpha <= alphaMax:
    FileName = "./../data/todaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        period = df["period"].values.tolist()[0]
        Q = np.array(df["Q"].values.tolist()[0:N+2])
        P = np.zeros(N+2)
        q,p = br.FT(Q,P,FourierComponents)
        En = br.TotalEnergy(q,p,alpha)
        AGPgrad = np.array(df["AGPgrad"].values.tolist())/En
        if period < 99.0:
            alphas = np.append(alphas,np.sqrt(En)*alpha)
            gradMag = np.append(gradMag,np.linalg.norm(AGPgrad))
            gradMagX = np.append(gradMagX,np.linalg.norm(AGPgrad[0:N]))
            gradMagP = np.append(gradMagP,np.linalg.norm(AGPgrad[N+1:2*N]))
       

        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha

alpha = 0.0

while alpha <= alphaMax:
    FileName2 = "./../data/alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df2 = pd.read_csv(FileName2)
        period = df2["period"].values.tolist()[0]
        AGPgrad2 = np.array(df2["AGPgrad"].values.tolist())/E0
        if period < 99.0:
            alphas2 = np.append(alphas2,np.sqrt(E0)*alpha)
            gradMag2 = np.append(gradMag2,np.linalg.norm(AGPgrad2))
            gradMagX2 = np.append(gradMagX2,np.linalg.norm(AGPgrad2[0:N]))
            gradMagP2 = np.append(gradMagP2,np.linalg.norm(AGPgrad2[N+1:2*N]))
        

        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha



fig = plt.figure()
plt.plot(alphas,gradMag)
plt.plot(alphas2,gradMag2)
plt.xlabel(r"$\sqrt{E}\alpha$",fontsize=22)
plt.ylabel(r"$|\nabla A_\alpha|/E$",fontsize=22)
plt.grid()
fig.set_size_inches(12,8)
plt.tight_layout()
plt.savefig(figName,dpi=100)
plt.show()
