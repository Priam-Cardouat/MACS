# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 16:01:01 2021

@author: Priam CARDOUAT
"""

import numpy as np
import matplotlib.pyplot as plt

def Iext(t):
    if t<5:
        return 0.59
    else : 
        return 0.6

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
        newY = Y[-1] +[pas*f(Y[-1], t, Iext)[0],pas*f(Y[-1], t, Iext)[1]]
        T.append(t+pas)
        Y.append(newY)
    return T, Y
print(euler_exp(f_exp, [0,0], Iext,0,300,.1))
    
def affichage_exp(debut, fin, pas):
    T, Y = euler_exp(f_exp, [0,0], Iext, debut, fin, pas)
    Y = np.array(Y)
    print(Y.shape)
    v = (Y.T)[0]
    w = (Y.T)[1]
    plt.plot(T,v, 'b', label = "v")
    plt.plot(T,w,'r', label = "w")
    I = [Iext(t) for t in T]
    plt.plot(T,I, '--g', label = "Iext")
    plt.legend()
    plt.show()
    
#affichage_exp(0, 300, 0.1)