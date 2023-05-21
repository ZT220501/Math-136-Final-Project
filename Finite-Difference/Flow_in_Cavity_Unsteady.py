# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:58:17 2023

@author: zheng
"""


'''
Code for running explicit finite difference method for solving the flow in cavity.
This code is a modification of the code in C. Pozrikidis's book
'''

import numpy as np
import matplotlib.pyplot as plt
import copy
from matplotlib.animation import FuncAnimation
from matplotlib import cm

import pandas as pd


'''
Parameters
'''
Vlid = 1.0                                  #lid velocity
Lx = 2                                    #cavity dimensions
Ly = 1                               
T = 10                                      #maximum time

Nx = 32
Ny = 16                                     #grid size
Nt = 500
visc = 0.01                                 #viscosity
rho = 1.0                                   #density
relax = 0.5                                 #relaxation parameter
Niteri = 5                                  #number of inner iterations
Niterg = 100                                #number of global iterations
vort_init = 0.0                             #initial vorticity


'''
Preparation
'''
Dx = Lx/Nx
Dy = Ly/Ny
Dt = T/Nt
Dx2 = 2.0*Dx
Dy2 = 2.0*Dy
Dxs = Dx*Dx
Dys = Dy*Dy
beta = Dxs/Dys 
beta1 = 2.0*(1.0+beta)
nu = visc/rho                               #kinematic viscosity
tol = 1e-8

'''
Generate the grid
Initialize stream function and vorticity
'''
x = np.zeros((Nx+1, Ny+1)).tolist()
y = np.zeros((Nx+1, Ny+1)).tolist()
psi = np.zeros((Nx+1, Ny+1)).tolist()
vort = np.zeros((Nx+1, Ny+1)).tolist()
for i in range(Nx+1):
    for j in range(Ny+1):
        x[i][j] = i*Dx
        y[i][j] = j*Dy
        psi[i][j] = 0.0                      #stream function
        vort[i][j] = -vort_init
        
ux = np.zeros((Nx+1, Ny+1))
uy = np.zeros((Nx+1, Ny+1))
#Set the initial speed for the grid along the lid
for i in range(Nx+1):
    ux[i][Ny] = Vlid

psi_total = [np.array(psi).T]
vort_total = [np.array(vort).T]
ux_total = [np.array(ux).T]
uy_total = [np.array(uy).T]


for time_step in range(Nt):
    psi = copy.deepcopy(psi_total[-1].T)
    vort = copy.deepcopy(vort_total[-1].T)
    ux = copy.deepcopy(ux_total[-1].T)
    uy = copy.deepcopy(uy_total[-1].T)
    '''
    Obtain vorticity using the definition
    '''
    #Top and bottoem wall
    for i in range(Nx+1):
        vort[i][0] = (-4*ux[i][1]+ux[i][2])/(2*Dy)
        vort[i][Ny] = (-3*Vlid+4*ux[i][Ny-1]-ux[i][Ny-2])/(2*Dy)
        
    #Left and right wall
    for i in range(Ny+1):
        vort[0][i] = (4*uy[1][i]-uy[2][i])/(2*Dx)
        vort[Nx][i] = (-4*uy[Nx-1][i]+uy[Nx-2][i])/(2*Dx)
     
    
    '''
    Obtain the interior vorticity after a timestep forward
    '''
    for i in range(1, Nx):
        for j in range(1, Ny):
            G = -ux[i][j]*(vort[i+1][j]-vort[i-1][j])/(2*Dx)-uy[i][j]*(vort[i][j+1]-vort[i][j-1])/(2*Dy)+nu*((vort[i+1][j]-2*vort[i][j]+vort[i-1][j])/Dxs+(vort[i][j+1]-2*vort[i][j]+vort[i][j-1])/Dys)
            vort[i][j] = vort[i][j] + Dt*G
    
    
    '''
    Solve Poisson's equation for the stream function
    '''
    for interior_itr in range(Niteri):
        for j in range(1, Ny):
            for i in range(1, Nx):
                res = (psi[i+1][j]+psi[i-1][j]+ beta*psi[i][j+1]+beta*psi[i][j-1]+Dxs*vort[i][j])/beta1-psi[i][j]
                psi[i][j] = psi[i][j] + relax*res
    
    '''
    Update the velocity components in the interior grids
    '''
    for i in range(1, Nx):
        for j in range(1, Ny):
            ux[i][j] = (psi[i][j+1]-psi[i][j-1])/(2*Dy)
            uy[i][j] = -(psi[i+1][j]-psi[i-1][j])/(2*Dx)
            
    psi_total.append(np.array(psi).T)
    vort_total.append(np.array(vort).T)
    ux_total.append(np.array(ux).T)
    uy_total.append(np.array(uy).T)
    

'''
Do the result animation
'''
xgr = []
ygr = []
for i in range(Nx+1):
    xgr.append(Dx*i)
for j in range(Ny+1):
    ygr.append(Dy*j)
    
X, Y = np.meshgrid(np.array(xgr), np.array(ygr))
    
for t in range(Nt):
    plt.quiver(X, Y, ux_total[t], uy_total[t], color='blue')
    plt.contour(xgr,ygr,np.array(psi_total[t]).tolist(), 8)
    plt.grid()
    plt.show()
    plt.pause(0.01)














