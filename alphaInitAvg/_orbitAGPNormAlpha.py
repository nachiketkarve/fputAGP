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


N = 32
E0 = 1
k0 = 1

#alphas = np.arange(0.3,0.4,0.01)
#alphas = np.array([0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])
alphas = np.array([0.5])

fig1 = plt.figure()

for alpha in alphas:

    found = 0
    iStart = 0

    FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    t = df['Time'].values
    AGPnorm = df['AGPnorm'].values
    
    AGPnormGradient = np.gradient(AGPnorm,t[1]-t[0])
    
    tFit = t[np.where((t > 10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**5))]
    AGPnormGradientFit = AGPnormGradient[np.where((t > 10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**5))]
    AGPnormFit = AGPnorm[np.where((t > 10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**5))]


    if len(tFit) > 2:
        p = np.polyfit(np.log(tFit), np.log(AGPnormGradientFit), 1)
        n = p[0] + 1
        a = np.exp(p[1])/n
        b = AGPnormFit[0] - a*tFit[0]**n
        plt.figure(fig1.number)
        plt.plot(tFit,a*tFit**n + b,'k--')
        print(n)

    plt.figure(fig1.number)
    plt.plot(t,AGPnorm,label=r"$E\alpha^2$ = %.2e"%(E0*alpha**2))

plt.figure(fig1.number)
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$T$")
plt.ylabel(r"$||A_\alpha(T)||^2$")
#plt.ylim([0.1,100])
plt.xlim([np.min(t),np.max(t)])
plt.grid()
plt.legend()
plt.tight_layout()
#fig1.set_size_inches(12,8)
#plt.savefig("_AGPNorm.pdf",dpi=100)

plt.show()

