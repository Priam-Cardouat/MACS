# -*- coding: utf-8 -*-
"""
Created on Sun Feb 14 15:39:49 2021

@author: Priam CARDOUAT
"""

import numpy as np
import matplotlib.pyplot as plt

A=np.array([[0,1,0],[0,0,1],[-1,-1,-1]])


def exo2_2(t0,tN,h):
    t=np.arange(t0,tN,h)
    x=np.array([0,-0.25,0.5])
    res=[x[0]]
    for n in range(1,len(t)):
        x+= h*((A @ x) + np.array([0,0,np.cos(t[n])]))
        res.append(x[0])
    return res

def z(t):
    return -0.25*t*np.cos(t)+0.25*t*np.sin(t)

t=np.arange(0,20,0.1)
x=np.array([z(i) for i in t])
y=np.array(exo2_2(0,20,0.1))
plt.plot(t,y,label='résolution numérique par la méthode de Euler')
plt.plot(t,x,c='red', label='solution exacte')
plt.legend()
plt.show()
