# -*- coding: utf-8 -*-
"""
Created on Sun Dec  1 23:08:44 2024

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
E0 = 1
k0 = 1
alphaMax = 0.5
alpha = 0.0
deltaAlpha = 0.01
frames = 500

yCircle = np.linspace(-1,1,100)
xCircle = np.sqrt(1-yCircle**2)

#A,V,C,FourierComponents = fput.initialize(N)

alphas = np.array([])
Es = np.array([])
periods = np.array([])

while alpha <= alphaMax:
    FileName = "./../data/alphaBr-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    try:
        df = pd.read_csv(FileName)
        period = df["period"].values.tolist()[0]
        E = df["E0"].values[0]
        if period < 99.0:
            alphas = np.append(alphas,alpha)
            Es = np.append(Es,np.linalg.norm(E))
            periods = np.append(periods,period)
        
        alpha = alpha + deltaAlpha
    
    except FileNotFoundError:
        alpha = alpha + deltaAlpha


plt.figure()
plt.plot(alphas,Es)
#plt.xscale("log")
#plt.yscale("log")

plt.figure()
plt.plot(alphas,periods)