# -*- coding: utf-8 -*-
"""
Created on Wed Mar 17 15:22:35 2021

@author: Priam CARDOUAT
"""

import numpy as np
import matplotlib.pyplot as plt
from rungekutta4 import RungeKutta4

tau = 11
a = 0.95
b = 0.4
T=300
def f (y,t):
  return(np.array([y[0]-((y[0]**3)/3)-y[1]+Iext(t),(y[0]+a-b*y[1])/tau]))

def RungeKutta2(T,n,y0):
  RK=np.zeros((n,2))
  RK[0]=y0
  h=T/n
  t=[0+k*h for k in range (0,n)]
  for i in range (n-1):
    val1=f(RK[i],t[i])
    val2=f(h*val1+RK[i],t[i+1])
    RK[i+1]=RK[i]+h*(0.5*val1+0.5*val2)  
  return t,RK 

def Iext(temps):
    #Bibi voit la carotte au bout du tiers du temps jusqu'aux deux tiers
    if temps < T/3:
        return 0.59
    if temps < 2*T/3:
        return 0.6
    else:
        return 0.59

t,Z=RungeKutta2(300,10000,np.array([0,0]))
n=10000
RK=np.zeros((n,2))
RK[0]=np.array([0,0])
h=T/n
for i in range (n-1):
    val1=(np.array([RK[i][0]-((RK[i][0]**3)/3)-RK[i][1]+Z[i][0],(RK[i][0]+a-b*RK[i][1])/tau]))
    y=h*val1+RK[i]
    val2=(np.array([y[0]-((y[0]**3)/3)-y[1]+Z[i][0],(y[0]+a-b*y[1])/tau]))
    RK[i+1]=RK[i]+h*(0.5*val1+0.5*val2)
plt.plot(t,Z[:,0], label = 'v neurone 1')
plt.plot(t,RK[:,0], label = 'v neurone 2')
plt.plot(t,[Iext(t[i]) for i in range(len(t))], label='tension ext')
n=1000
h=T/n
t=[0+k*h for k in range(n+1)]
resultat=RungeKutta4(T,n,np.array([0,0]), True) 
neurone2=RungeKutta4(T,n,np.array([0,0]), False)
plt.plot(t, resultat[:, 0], label = "v neurone 1 rk4")
plt.plot(t, neurone2[:, 0], label = "v neurone 2 rk4")
plt.legend(loc='best')
plt.xlabel('t')
plt.grid()
plt.legend()
