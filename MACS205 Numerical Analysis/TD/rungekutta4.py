# -*- coding: utf-8 -*-
"""
Created on Mon Mar  8 23:33:06 2021

@author: Priam CARDOUAT
"""

import numpy as np
import matplotlib.pyplot as plt

a=0.95
b=0.4
tau=11
T=300
n=1000
h=T/n
t=[0+k*h for k in range(n+1)]

def RungeKutta4(T,n,y0,boolean):
    result = np.zeros((n+1, len(y0)))
    result[0]=y0
    for i in range(len(t)-1):
        pk1 = f(result[i], t[i], boolean, i)
        pk2 = f(result[i] + pk1 * h / 2, t[i] + h / 2, boolean, i)
        pk3 = f(result[i] + pk2 * h / 2, t[i] + h / 2, boolean, i)
        pk4 = f(result[i] + pk3 * h, t[i] + h, boolean, i)
        result[i+1] = result[i] + (h / 6) * (pk1 + 2*pk2 + 2*pk3 + pk4)
    return result

#La fonction f va gérer les 2 neurones selon le cas
def f (y, temps, boolean, i):
    if boolean == True :
        #neurone 1, entrée Iext
        return np.array([y[0]-((y[0]**3)/3)-y[1]+I(temps),(y[0]+a-b*y[1])/tau])
    else:
        #neurone 2, entrée = sortie du 1er neurone
        return np.array([y[0]-((y[0]**3)/3)-y[1]+resultat[i][0],(y[0]+a-b*y[1])/tau])
        

def I(temps):
    #Bibi voit la carotte au bout de la moitié du temps
    if temps < T/3:
        return 0.59
    if temps < 2*T/3:
        return 0.6
    else:
        return 0.59
    


resultat=RungeKutta4(T,n,np.array([0,0]), True)
"""neurone2=RungeKutta4(T,n,np.array([0,0]), False)
plt.plot(t, resultat[:, 0], label = "v neurone 1 rk4")
plt.plot(t, neurone2[:, 0], label = "v neurone 2 rk4")
plt.plot(t, [I(k) for k in t], label = "Tension nerf optique")
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.show()"""