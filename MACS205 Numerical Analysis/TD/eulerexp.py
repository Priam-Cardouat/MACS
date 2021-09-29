# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:40:54 2021

@author: Priam CARDOUAT
"""


import numpy as np
import matplotlib.pyplot as plt

def f_exp(Y, t, Iext):
    v,w = Y[0],Y[1]
    A = v - (v**3)/3 - w + Iext(t)
    B = 1/11*(v + 0.95 - 0.4*w)
    return [A,B]

def euler_exp(f, Y0, Iext, debut, fin, pas):
    nb_pas = int((fin-debut)//pas)
    T = [debut]
    Y = [Y0]
    for i in range(nb_pas):
        t = debut + pas*i
        lastY = Y[-1]
        pas_f = f(Y[-1], t, Iext)
        newY = [lastY[0] + pas * pas_f[0], lastY[1] + pas * pas_f[1]]
        T.append(t+pas)
        Y.append(newY)
    return T, Y

T=300

def Iext(temps):
    #Bibi voit la carotte au bout du tiers du temps jusqu'aux deux tiers
    if temps < T/3:
        return 0.59
    if temps < 2*T/3:
        return 0.6
    else:
        return 0.59
    
def affichage_exp(debut, fin, pas):
    T, Y = euler_exp(f_exp, [0,0], Iext, debut, fin, pas)
    Y = np.array(Y)
    v = (Y.T)[0]
    w = (Y.T)[1]
    plt.plot(T,v, 'b', label = "v")
    plt.plot(T,w,'r', label = "w")
    I = [Iext(t) for t in T]
    plt.plot(T,I, '--g', label = "Iext")
    plt.legend()
    plt.show()
    
affichage_exp(0,300,0.1)