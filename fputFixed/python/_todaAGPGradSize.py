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

isSingleCol = 1
colFrac = 0.47

doubleColpt = 246
singleColpt = 510

fig_width_pt = (isSingleCol*singleColpt + (1-isSingleCol)*doubleColpt)
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size = [fig_width*colFrac,fig_width*colFrac]
fig_font = 10
params = {'backend' : 'ps',
          'axes.labelsize': fig_font,
          'font.size': fig_font,
          'legend.fontsize': fig_font,
          'xtick.labelsize': fig_font,
          'ytick.labelsize': fig_font,
          'text.usetex': True,
          'figure.figsize': fig_size,
          'font.family': 'serif',
          'font.serif': 'STIX',
          'mathtext.fontset': 'stix'}

plt.rcParams.update(params)

def TotalEnergy(q,p,alpha):
    N = len(q)-2
    En = 0
    for i in range(N+1):
        En = En + p[i]**2/2.0 + 1.0/(2.0*alpha)**2*(np.exp(2.0*alpha*(q[i+1]-q[i])) - 2.0*alpha*(q[i+1]-q[i]) - 1.0)

    return En

Ns = np.array([2,3,4])
E0 = 1.0
k0 = 1
alphaMax = 2.7
alpha = 0.01
deltaAlpha = 0.01
frames = 500

figName = ".//..//plot//_agpGrad_todaBr" + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".pdf"

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

fig = plt.figure()

for N in Ns:
    A,V,C,FourierComponents = br.initialize(N)

    alphas = np.array([])
    gradMag = np.array([])
    gradMagX = np.array([])
    gradMagP = np.array([])

    alpha = 0.01

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

    plt.plot(alphas,gradMag,label = r"$N = $ %i"%N)


plt.xlabel(r"$\sqrt{E}\alpha$")
plt.ylabel(r"$|\nabla A_\alpha|/E$")
plt.grid()
plt.tight_layout()
plt.legend()
plt.savefig(figName,dpi=100)
plt.show()
