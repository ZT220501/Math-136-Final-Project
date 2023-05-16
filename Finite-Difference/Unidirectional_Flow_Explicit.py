# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 20:21:44 2023

@author: zheng
"""


'''
Code for running explicit finite difference method for solving the unidirecitonal flow
This code is a modification of the code in C. Pozrikidis's book
'''


'''
TODO: Replace pressure gradient to be non-constant
'''


import numpy as np
import matplotlib.pyplot as plt


'''
Initialze Values
'''
h = 1
mu = 0.6 
rho = 5.5
N = 32
dpdx = -2
gx = 0
V1 = 0
V2 = 0

#The diffusion number
#The method is stable iff alpha<=0.5
al = 0.49                       #Alpha value                  
nstep = 20000                   #Step numbers

#Prepare
nu = mu/rho                     #kinematic viscosity
Dy = h/N
Dt = al*Dy*Dy/nu



'''
Grid and initial condition
'''
y = []
u = []
for i in range(N+1):
    y.append(i*Dy)
    u.append(0)
    
u[0] = V1
u[N] = V2

'''
Run Iterations
'''
t = 0
for step in range(nstep):
    t += Dt
    unew = [V1]
    for i in range(1, N):
        unew.append(al*u[i-1] + (1-2*al)*u[i] + al*u[i+1] + Dt*(-dpdx/rho+gx))
    unew.append(V2)
    u = unew
    plt.plot(u, y, 'o-')
    plt.axis([0, 0.5, 0, h])
    plt.pause(0.01)
        



































