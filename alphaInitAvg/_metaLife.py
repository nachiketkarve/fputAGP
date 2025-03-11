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


isSingleCol = 1
colFrac = 0.9

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

#plt.rcParams.update(params)


N = 64
E0 = 1
k0 = 1

alphas = np.arange(0.2,0.6,0.01)
#alphas = np.array([0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])
#alphas = np.array([0.2,0.3,0.4,0.5])

nonLins = np.array([])
tMeta = np.array([])

fig1 = plt.figure()

for alpha in alphas:

    found = 0
    iStart = 0

    FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    t = df['Time'].values
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

    idx = np.where(np.abs(EntAvg-Ent2Avg) > 0.1)[0]
    print(len(idx))
    if len(idx) > 1:
        nonLins = np.append(nonLins,E0*alpha**2)
        tMeta = np.append(tMeta, np.min(t[idx]))

    #plt.figure(fig1.number)
    #plt.plot(t,EntAvg,label="Alpha")
    #plt.plot(t,Ent2Avg,label="Toda")
    #plt.plot(t,EntAvg-Ent2Avg)

plt.figure(fig1.number)
plt.plot(nonLins,tMeta)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$T$")
plt.ylabel(r"$||A_\alpha(T)||^2$")
#plt.ylim([0.1,100])
#plt.xlim([np.min(t),np.max(t)])
plt.grid()
#plt.legend()
plt.tight_layout()
#fig1.set_size_inches(12,8)
#plt.savefig("_AGPNorm.pdf",dpi=100)

#print(tMeta)
plt.show()

