# -*- coding: utf-8 -*-
"""
Created on Mon Nov  4 14:49:21 2024

@author: karve
"""

# -*- coding: utf-8 -*-
"""
Created on Sun Nov  3 15:15:21 2024

@author: karve
"""

#import sys
#sys.path.append('F://FPUT//python//fput')

import numpy as np
import matplotlib.pyplot as plt
#import fputAlpha as br
import pandas as pd
from scipy.optimize import curve_fit
import matplotlib as mpl


isSingleCol = 1
colFrac = 0.45

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


N = 64
E0 = 1
k0 = 1

markers = ['^','*']
lines = ['-', '-']
colors = ['r', 'k']

cm = mpl.colormaps['viridis']

alphas = np.arange(0.21,0.6,0.01)
#alphas = np.array([0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])

nonLins = np.array([])
lifeTimes = np.array([])
entLifeTimes = np.array([])
entNonLins = np.array([])
metaLifeTimes = np.array([])
metaNonLins = np.array([])

fig1 = plt.figure()

for i in range(len(alphas)):

    alpha = alphas[i]

    try:
        FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
        df = pd.read_csv(FileName)
        t = df['Time'].values
        AGPnorm = df['AGPnorm'].values
        Ent = df['Entropy'].values
        EntAvg = np.zeros(len(t))
        for i in range(1,len(t)):
            EntAvg[i] = np.mean(Ent[0:i])
        
        FileName2 = "./../todaInitAvg/todaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
        df2 = pd.read_csv(FileName2)
        t2 = df2['Time'].values
        Ent2 = df2['Entropy'].values
        Ent2Avg = np.zeros(len(t2))
        for i in range(1,len(t2)):
            Ent2Avg[i] = np.mean(Ent2[0:i])
    
        for j in range(len(t)):
            if AGPnorm[j] > 100.0:
                lifeTimes = np.append(lifeTimes,t[j])
                nonLins = np.append(nonLins,E0*alpha**2)
                break
            
        for j in range(len(t)):
            if Ent[j] > 0.87:
                entLifeTimes = np.append(entLifeTimes,t[j])
                entNonLins = np.append(entNonLins,E0*alpha**2)
                break
        
        idx = np.where(np.divide(np.abs(EntAvg[2:len(t)]-Ent2Avg[2:len(t)]),Ent2Avg[2:len(t)]) > 0.1)[0]
        if len(idx) > 1:
            metaNonLins = np.append(metaNonLins,E0*alpha**2)
            metaLifeTimes = np.append(metaLifeTimes, np.min(t[idx]))

    except FileNotFoundError:
        continue


coeff = np.polyfit(np.log(nonLins),np.log(lifeTimes),1)
print(coeff)

entCoeff = np.polyfit(np.log(entNonLins),np.log(entLifeTimes),1)
print(entCoeff)

plt.figure(fig1.number)
plt.plot(nonLins,lifeTimes, linestyle = lines[0], marker = markers[0], color = colors[0], label = r"$t_c$")
#plt.plot(entNonLins,entLifeTimes, linestyle = lines[1], marker = markers[1], color = colors[1], label = r"$t_{therm}$")
plt.plot(metaNonLins,metaLifeTimes, linestyle = lines[1], marker = markers[1], color = colors[1], label = r"$t_{therm}$")
#plt.plot(nonLins,np.exp(coeff[0]*np.log(nonLins) + coeff[1]),'k',label=r"$T \sim (E\alpha^2)^{%.2f}$"%coeff[0])
#plt.plot(entNonLins,np.exp(entCoeff[0]*np.log(entNonLins) + entCoeff[1]),'k',label=r"$T \sim (E\alpha^2)^{%.2f}$"%entCoeff[0])
#plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$E\alpha^2$")
plt.ylabel(r"$t$")
plt.legend()
plt.grid()
plt.tight_layout(pad = 0.1)
plt.savefig('_lifeTimeEntMeta.pdf',dpi=100)

plt.show()

