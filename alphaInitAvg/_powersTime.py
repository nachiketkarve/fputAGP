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
alpha = 0.15

FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
df = pd.read_csv(FileName)
t = df['Time'].values
AGPnorm = df['AGPnorm'].values
        
AGPnormGradient = np.gradient(AGPnorm,t[1]-t[0])
tFit = t[np.where((t > 1*10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**(-10)))]
AGPnormGradientFit = AGPnormGradient[np.where((t > 1*10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**(-10)))]
AGPnormFit = AGPnorm[np.where((t > 1*10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**(-10)))]

powers = np.array([])

for i in range(len(tFit)):

    if i <= 20:
        powers = np.append(powers,0)
    else:
        p = np.polyfit(np.log(tFit[0:i]), np.log(AGPnormGradientFit[0:i]), 1)
        n = p[0] + 1
        #p = np.polyfit(np.log(t[1:i]), np.log(AGPnorm[1:i]), 1)
        #n = p[0]
        powers = np.append(powers,n)
        
fig1 = plt.figure()
plt.plot(tFit,powers)
#plt.xlabel(r"$E\beta$")
#plt.ylabel(r"$\gamma$")
plt.grid()
plt.tight_layout()
#fig1.set_size_inches(12,8)
plt.savefig("_powers.pdf",dpi=100)

fig2 = plt.figure()
plt.plot(t,AGPnormGradient)
#plt.xlabel(r"$E\beta$")
#plt.ylabel(r"$\gamma$")
plt.grid()
plt.tight_layout()
#fig1.set_size_inches(12,8)
#plt.savefig("_powers.pdf",dpi=100)




plt.show()
