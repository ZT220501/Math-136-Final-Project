# -*- coding: utf-8 -*-
"""
Created on Sat May 13 16:05:08 2023

@author: zheng
"""


'''
Code for running implicit finite difference method for solving the unidirecitonal flow
This code is a modification of the code in C. Pozrikidis's book
'''

'''
TODO: Replace pressure gradient to be non-constant
TODO: Change the tridiagonal solver to others such as the Arnoldi/Lanczos iteration
'''


import numpy as np
import matplotlib.pyplot as plt

'''
Initial Values
'''
h = 1.0
mu = 0.2
rho = 0.5
N = 32
dpdx = -2.0
gx = 0
V1 = 0 
V2 = 0

#The diffusion number
alpha = 0.4;               
nstep = 20000;                          #Number of steps
nu = mu/rho;
Dy = h/N;
Dt = alpha*Dy*Dy/nu;

'''
Grid and initial condition
'''
y = []
u = []
for i in range(N+1):    
    y.append(i*Dy)
    u.append(0)
    
u[0] = V1;
u[N] = V2;
'''
Formulate the tridiagonal projection matrix
atr is the diagonal line of the coefficient matrix
btr is the superdiagonal line of the coefficient matrix
ctr is the subdiagonal line of the coefficient matrix
'''
atr = []
btr = []
ctr = []
for i in range(N-1):
    atr.append(1.0 + 2*alpha)
    btr.append(-alpha)
    ctr.append(-alpha)
    
#Tridiagonal Solver
'''
Thomas algorithm for a tridiagonal system
n: system size
a,b,c: diagonal, superdiagonal, and subdiagonal elements
s: right-hand side
'''
def thomas(n,a,b,c,s):

    
    '''
    % Reduction to upper bidiagonal
    '''
    d = [b[0]/a[0]]
    y = [s[0]/a[0]]

    for i in range(n-2):
        j = i + 1
        den = a[j] - c[j] * d[i]
        d.append(b[j]/den)
        y.append((s[j]-c[j]*y[i])/den)

    den = a[n-1]-c[n-1]*d[n-2]
    y.append((s[n-1]-c[n-1]*y[n-2])/den)
    
    
    '''
    Back substitution
    '''
    x = [y[n-1]]
    for i in range(n-2, -1, -1):
        x.insert(0, y[i] - d[i] * x[0])
        
    return x

'''
Run iterations
'''
t = 0
for step in range(nstep):
    s = []
    for i in range(N-1):
        s.append(u[i+1] + Dt * (-dpdx/rho + gx))
    t = t + Dt
    u[0] = V1
    u[N] = V2
    s[0] = s[0] + alpha * u[0]
    s[N-2] = s[N-2] + alpha * u[N]
    sol = thomas(N-1,atr,btr,ctr,s);
    for i in range(1, N):
        u[i] = sol[i-1]
        
    plt.plot(u, y, 'o-')
    plt.axis([0, 0.5, 0, h])
    plt.pause(0.01)
    






    
    
    
    
    
    
    