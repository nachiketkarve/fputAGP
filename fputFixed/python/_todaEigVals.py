# -*- coding: utf-8 -*-
"""
Created on Wed Dec  4 12:54:41 2024

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("/usr3/graduate/nachiket/FPUT/todaFixedPython")
import toda as br
import matplotlib.animation as ani

N = 3
E0 = 1.0
k0 = 1
alphaMax = 4
alpha = 0.01
deltaAlpha = 0.01
frames = 500

A,V,C,FourierComponents = br.initialize(N)

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

#A,V,C,FourierComponents = fput.initialize(N)

alphas = np.array([])

gifName = ".//..//plot//_eigVals_todaBr_N" + str(int(N)) + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".gif"

while alpha <= alphaMax:
    FileName = ".//..//data//todaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        period = df["period"].values.tolist()[0]
        EigValReal = np.array(df["EigValReal"].values.tolist())
        EigValImag = np.array(df["EigValImag"].values.tolist())
        #plt.scatter(EigValReal,EigValImag)
        if period < 99.0:
            alphas = np.append(alphas,alpha)
        
        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha

if frames >= len(alphas):
    frames = len(alphas)-1

interval = int(len(alphas)/frames)
alphasReduced = np.array([])
for i in range(len(alphas)):
    if i%interval == 0:
        alphasReduced = np.append(alphasReduced, alphas[i])
    
frames = len(alphasReduced)

#print(alphasReduced)

def animate(i):
    plt.clf()
    alpha = alphasReduced[i]
    FileName = ".//..//data//todaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    EigValReal = np.array(df["EigValReal"].values.tolist())
    EigValImag = np.array(df["EigValImag"].values.tolist())
    Q = np.array(df["Q"].values.tolist()[0:N+2])
    P = np.zeros(N+2)
    q,p = br.FT(Q,P,FourierComponents)
    En = br.TotalEnergy(q,p,alpha)
    plt.scatter(EigValReal,EigValImag,color = 'b')
    plt.plot(xCircle,yCircle,'r')
    plt.plot(-xCircle,yCircle,'r')
    plt.xlim([-1.1,1.1])
    plt.ylim([-1.1,1.1])
    plt.grid()
    plt.title(r"$N$ = %i, $k0$ = %i, $\sqrt{E}\alpha$ = %.3f"%(N,k0,np.sqrt(En)*alpha),fontsize=20)

fig = plt.figure()
fig.set_size_inches(8, 8)
anim = ani.FuncAnimation(fig, animate, frames = int(frames), interval=1, repeat=False, save_count = int(frames))
writergif = ani.PillowWriter(fps=5) 
anim.save(gifName,writer=writergif)

