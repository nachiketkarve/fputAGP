# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 16:15:18 2024

@author: karve
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys
sys.path.append("F:\\FPUT\\python\\fput")
import fputAlpha as fput
import matplotlib.animation as ani

N = 2
E0 = 1.0
k0 = 1

tmax = 200.0
deltaT = 0.001
frames = 500

A,V,C,FourierComponents = fput.initialize(N)

FileName = "./../data/alphaAd-N" + str(int(N)) + "-K" + str(int(k0)) + "-E" + "{:.6f}".format(E0) + ".csv"
df = pd.read_csv(FileName)
time = np.array(df["Time"].values.tolist())
nonLin = np.array(df["NonLin"].values.tolist())
Energy = np.array(df["Energy"].values.tolist())

FileName = "./../data/endstate.csv"
df = pd.read_csv(FileName)
q = np.array(df["q"].values.tolist())
p = np.array(df["p"].values.tolist())

t,qArray,pArray,EArray,_,_ = fput.Evolve(q,p,nonLin[len(nonLin)-1],tmax,deltaT,frames)

spectrum = np.linspace(0,N+1,N+2)

fig = plt.figure()
Y,X = np.meshgrid(t,spectrum[1:N+1])
cs1 = plt.pcolormesh(X,Y,pArray[1:N+1,:],cmap='plasma',rasterized=True)
cbar1 = plt.colorbar(cs1)
cbar1.set_ticks([-0.1,0,0.1])
cbar1.ax.tick_params(rotation=90)
#cbar1.set_label(label='$p_n$')
plt.xlabel('$n$')
plt.ylabel('$t$')
plt.yticks(rotation=90)
plt.minorticks_on()
plt.tight_layout()

fig2 = plt.figure()

plt.plot(nonLin,Energy)

pMax = np.amax(pArray)
pMin = np.amin(pArray)

def animate(i):
    plt.clf()
    plt.plot(spectrum[0:N+1],pArray[0:N+1,i])
    plt.xlabel(r"$n$",fontsize=12)
    plt.ylabel(r"$p_n$",fontsize=12)
    plt.title(r'$N =$ %i, $E =$ %.1f, $\alpha =$ %.1f'%(N,E0,nonLin[len(nonLin)-1]) + '\n' + r'$t =$ %.2f'%(tmax*i/frames),fontsize=12)
    plt.ylim([pMin*1.1,pMax*1.1])
    plt.xlim([spectrum[0],spectrum[N]])
    plt.grid()
    plt.tight_layout()
    
fig = plt.figure()
fig.set_size_inches(12,8)
anim = ani.FuncAnimation(fig, animate, frames=frames,  interval=1, repeat=False)
writergif = ani.PillowWriter(fps=30) 
anim.save('fput_alpha.gif',writer=writergif) 