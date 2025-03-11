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
htFrac = 0.8

doubleColpt = 246
singleColpt = 510

fig_width_pt = (isSingleCol*singleColpt + (1-isSingleCol)*doubleColpt)
inches_per_pt = 1.0/72.27
golden_mean = (np.sqrt(5)-1.0)/2.0
fig_width = fig_width_pt*inches_per_pt
fig_height = fig_width*golden_mean
fig_size = [fig_width*colFrac,fig_width*colFrac*htFrac]
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


N = 32
E0 = 1
k0 = 1

#cm = mpl.colormaps['viridis']

alphas = np.array([0.2,0.3,0.4])
#alphas = np.array([1.2,1.3,2])

clrs = ['b','r','g']
lines = ['-','-','-']
mrkrs = ['o','^','*']

nonLins = np.array([])
lifeTimes = np.array([])
fig1 = plt.figure()

for i in range(len(alphas)):

    alpha = alphas[i]

    try:

        FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
        df = pd.read_csv(FileName)
        t = df['Time'].values
        AGPnorm = df['AGPnorm'].values

        AGPnorm = AGPnorm[np.where(t <= 1e6)]
        t = t[np.where(t <= 1e6)]

        tRed = np.array([])
        AGPnormRed = np.array([])
        for j in range(0,len(t),int(len(t)/20)):
            tRed = np.append(tRed,t[j])
            AGPnormRed = np.append(AGPnormRed,AGPnorm[j])

        plt.figure(fig1.number)
        plt.plot(tRed,AGPnormRed,label=r"$E\alpha^2$ = %.2f"%(E0*alpha**2),marker=mrkrs[i], linestyle=lines[i], color=clrs[i])
        #plt.plot(t,AGPnorm,label=r"$E\beta$ = %.2e"%(E0*beta),color=s_m.to_rgba(i))

    except FileNotFoundError:
        continue

plt.figure(fig1.number)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$T$")
plt.ylabel(r"$\chi_\alpha(T)$")
#plt.ylim([0.1,100])
plt.xticks([0,0.2e6,0.4e6,0.6e6,0.8e6,1e6])
plt.ticklabel_format(style='sci',axis='x',scilimits=(0,0))
plt.ticklabel_format(style='sci',axis='y',scilimits=(0,0))
plt.xlim([np.min(t),1e6])
plt.grid()
plt.legend()
plt.tight_layout(pad=0.1)
#fig1.set_size_inches(12,8)
plt.savefig("_AGPNorm.pdf",dpi=100)
plt.savefig("_AGPNorm.eps",dpi=100)



plt.show()
