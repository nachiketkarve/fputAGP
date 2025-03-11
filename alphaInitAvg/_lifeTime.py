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

cm = mpl.colormaps['viridis']

alphas = np.arange(0.1,0.2,0.001)
#alphas = np.array([0.3,0.31,0.32,0.33,0.34,0.35,0.36,0.37,0.38,0.39,0.4])

nonLins = np.array([])
lifeTimes = np.array([])
powers = np.array([])
nonLinsPower = np.array([])
entLifeTimes = np.array([])
entNonLins = np.array([])

fig1,ax1 = plt.subplots()
fig2 = plt.figure()
fig3 = plt.figure()
fig4,ax4 = plt.subplots()

colors = cm(np.linspace(0,1,len(alphas)))
norm = mpl.colors.Normalize(vmin = np.min(E0*alphas**2),vmax = np.max(E0*alphas**2))
s_m = mpl.cm.ScalarMappable(cmap=cm,norm=norm)
s_m.set_array([])

for i in range(len(alphas)):

    alpha = alphas[i]

    found = 0
    iStart = 0

    FileName = "alphaInit-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + "-A" + "{:.6f}".format(alpha) + ".csv"
    df = pd.read_csv(FileName)
    t = df['Time'].values
    AGPnorm = df['AGPnorm'].values
    #Ent = df['Entropy'].values
    Ent = np.zeros(len(AGPnorm))

    AGPnorm = AGPnorm/N
    
    AGPnormGradient = np.gradient(AGPnorm,t[1]-t[0])

    tFit = t[np.where((t > 10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**5))]
    AGPnormGradientFit = AGPnormGradient[np.where((t > 10**5) & (AGPnormGradient > 0) & (AGPnorm > 10**5))]

    if len(tFit) > 2:
        p = np.polyfit(np.log(tFit), np.log(AGPnormGradientFit), 1)
        n = p[0] + 1
        powers = np.append(powers,n)
        nonLinsPower = np.append(nonLinsPower,E0*alpha**2)
    
    #for j in range(len(t)):
        #if t[j] > 10**5 and iStart == 0:
            #p = np.polyfit(np.log(t[j+1:len(t)]-t[j]), np.log(np.abs(AGPnorm[j+1:len(t)]-AGPnorm[j])),1)
            #p = np.polyfit(np.log(t[j:len(t)]),np.log(np.abs(AGPnormGradient[j:len(t)])),1)
            #n = p[0] + 1

            #p = np.polyfit(np.log(t[j:len(t)]), np.log(AGPnorm[j:len(t)]), 1)
            #n = p[0]

            #powers = np.append(powers,n)
            #nonLinsPower = np.append(nonLinsPower,E0*alpha**2)
            #break

    for j in range(len(t)):
        if AGPnorm[j] > 100.0:
            lifeTimes = np.append(lifeTimes,t[j])
            nonLins = np.append(nonLins,E0*alpha**2)
            break
            
    for j in range(len(t)):
        if Ent[j] > 0.85:
            entLifeTimes = np.append(entLifeTimes,t[j])
            entNonLins = np.append(entNonLins,E0*alpha**2)
            break

    plt.figure(fig1.number)
    plt.plot(t,AGPnorm,label=r"$E\alpha^2$ = %.2e"%(E0*alpha**2),color=colors[i])

    plt.figure(fig4.number)
    plt.plot(t,Ent,label=r"$E\alpha^2$ = %.2e"%(E0*alpha**2),color=colors[i])


plt.figure(fig1.number)
cb = plt.colorbar(s_m,ax=ax1,label=r"$E\alpha^2$")
#cb.ax.set_label(r"$E\alpha^2$")
plt.xscale("log")
plt.yscale("log")
plt.xlabel(r"$T$")
plt.ylabel(r"$||A_\alpha(T)||^2$")
#plt.ylim([0.1,100])
plt.xlim([np.min(t),np.max(t)])
plt.grid()
#plt.legend()
plt.tight_layout()
#fig1.set_size_inches(12,8)
plt.savefig("_AGPNorm.pdf",dpi=100)

coeff = np.polyfit(np.log(nonLins),np.log(lifeTimes),1)
print(coeff)

#entCoeff = np.polyfit(np.log(entNonLins),np.log(entLifeTimes),1)
#print(entCoeff)

plt.figure(fig2.number)
plt.scatter(nonLins,lifeTimes)
#plt.scatter(entNonLins,entLifeTimes)
plt.plot(nonLins,np.exp(coeff[0]*np.log(nonLins) + coeff[1]),'k',label=r"$T \sim (E\alpha^2)^{%.2f}$"%coeff[0])
#plt.plot(entNonLins,np.exp(entCoeff[0]*np.log(entNonLins) + entCoeff[1]),'k',label=r"$T \sim (E\alpha^2)^{%.2f}$"%entCoeff[0])
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$E\alpha^2$")
plt.ylabel(r"$T$")
plt.legend()
plt.grid()
plt.tight_layout()
plt.savefig('_lifeTime.pdf',dpi=100)

plt.figure(fig3.number)
plt.plot(nonLinsPower,powers)
plt.xlabel(r"$E\alpha^2$")
plt.ylabel(r"$n$")
plt.grid()
plt.tight_layout()
#fig1.set_size_inches(12,8)
plt.savefig("_powers.pdf",dpi=100)

plt.figure(fig4.number)
cb = plt.colorbar(s_m,ax=ax4,label=r"$E\alpha^2$")
#cb.ax.set_label(r"$E\alpha^2$")
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel(r"$T$")
plt.ylabel(r"$S(T)$")
#plt.ylim([0.1,100])
plt.xlim([np.min(t),np.max(t)])
plt.grid()
#plt.legend()
plt.tight_layout()
#fig1.set_size_inches(12,8)
plt.savefig("_entropy.pdf",dpi=100)


plt.show()

