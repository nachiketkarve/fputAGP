# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 17:34:19 2024

@author: karve
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("/usr3/graduate/nachiket/FPUT/fputAlphaFixedPython")
import fputAlpha as fput
import matplotlib.animation as ani

N = 2
E0 = 1.0
k0 = 1
alphaMax = 4
alpha = 0.0
deltaAlpha = 0.01
frames = 500

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

#A,V,C,FourierComponents = fput.initialize(N)

figName = ".//..//plot//_agpGrad_alphaBr_N" + str(int(N)) + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".pdf"

alphas = np.array([])
gradMag = np.array([])
gradMagX = np.array([])
gradMagP = np.array([])

while alpha <= alphaMax:
    FileName = "./../data/alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        period = df["period"].values.tolist()[0]
        AGPgrad = np.array(df["AGPgrad"].values.tolist())
        if period < 99.0:
            alphas = np.append(alphas,np.sqrt(E0)*alpha)
            gradMag = np.append(gradMag,np.linalg.norm(AGPgrad))
            gradMagX = np.append(gradMagX,np.linalg.norm(AGPgrad[0:N]))
            gradMagP = np.append(gradMagP,np.linalg.norm(AGPgrad[N+1:2*N]))
        
        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha


fig = plt.figure()
plt.plot(alphas,gradMag)
plt.xlabel(r"$\sqrt{E}\alpha$",fontsize=22)
plt.ylabel(r"$|\nabla A_\alpha|$",fontsize=22)
plt.grid()
fig.set_size_inches(12,8)
plt.tight_layout()
plt.savefig(figName,dpi=100)
plt.show()
