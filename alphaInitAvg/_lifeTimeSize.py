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


Ns = np.array([32, 64, 128, 256])
E0 = 1
k0 = 1

markers = ['o','^','*','v']
lines = ['-','-','-','-']
colors = ['b','r','g','purple']
thresholds = [50,50,10,10]

cm = mpl.colormaps['viridis']

alphas = np.arange(0.2,0.6,0.01)
#alphas = np.array([0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])

fig1 = plt.figure()

for iN in range(len(Ns)):

    N = Ns[iN]
    nonLins = np.array([])
    lifeTimes = np.array([])

    for i in range(len(alphas)):

        alpha = alphas[i]

        try:
            FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
            df = pd.read_csv(FileName)
            t = df['Time'].values
            AGPnorm = df['AGPnorm'].values

            for j in range(len(t)):
                if AGPnorm[j]/N > thresholds[iN]:
                    lifeTimes = np.append(lifeTimes,t[j])
                    nonLins = np.append(nonLins,E0*alpha**2)
                    break

        except FileNotFoundError:
            continue

    coeff = np.polyfit(np.log(nonLins),np.log(lifeTimes),1)
    print(coeff)

    plt.figure(fig1.number)
    plt.plot(nonLins,lifeTimes,label=r"$N = $ %i"%N, linestyle = lines[iN], marker = markers[iN], color = colors[iN])
    #plt.plot(nonLins,np.exp(coeff[0]*np.log(nonLins) + coeff[1]),'k',label=r"$T \sim (E\alpha^2)^{%.2f}$"%coeff[0])

plt.figure(fig1.number)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$E\alpha^2$")
plt.ylabel(r"$t_c$")
plt.legend()
plt.grid()
plt.tight_layout(pad=0.1)
plt.savefig('_lifeTimeSize.pdf',dpi=100)
plt.savefig('_lifeTimeSize.eps',dpi=100)


plt.show()

