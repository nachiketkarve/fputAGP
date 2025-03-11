# -*- coding: utf-8 -*-
"""
Created on Sat Aug 17 17:03:10 2024

@author: karve
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
#sys.path.append("F:\\FPUT\\fput")
#import fputAlpha as fput
import matplotlib.animation as ani

N = 2
E0 = 1.0
k0 = 1
alphaMax = 4
alpha = 0.01
deltaAlpha = 0.01
frames = 500

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

#A,V,C,FourierComponents = fput.initialize(N)

alphas = np.array([])

gifName = ".//..//plot//_eigVals_alphaBr_N" + str(int(N)) + "_K" + str(int(k0)) + "_E" + "{:.6f}".format(E0) + ".gif"

while alpha <= alphaMax:
    FileName = ".//..//data//alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
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
    FileName = ".//..//data//alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    EigValReal = np.array(df["EigValReal"].values.tolist())
    EigValImag = np.array(df["EigValImag"].values.tolist())
    plt.scatter(EigValReal,EigValImag,color = 'b')
    plt.plot(xCircle,yCircle,'r')
    plt.plot(-xCircle,yCircle,'r')
    plt.xlim([-1.1,1.1])
    plt.ylim([-1.1,1.1])
    plt.grid()
    plt.title(r"$N$ = %i, $E$ = %.2e, $k0$ = %i, $\alpha$ = %.3f"%(N,E0,k0,alpha),fontsize=20)

fig = plt.figure()
fig.set_size_inches(8, 8)
anim = ani.FuncAnimation(fig, animate, frames = int(frames), interval=1, repeat=False, save_count = int(frames))
writergif = ani.PillowWriter(fps=5) 
anim.save(gifName,writer=writergif)

