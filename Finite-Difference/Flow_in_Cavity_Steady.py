# -*- coding: utf-8 -*-
"""
Created on Tue May 16 18:38:41 2023

@author: zheng
"""

'''
Code for running explicit finite difference method for solving the flow in cavity.
This code is a modification of the code in C. Pozrikidis's book
'''

import numpy as np
import matplotlib.pyplot as plt
import copy
import pandas as pd

'''
Parameters
'''
Vlid = 1.0                                  #lid velocity
Lx = 2.0                                    #cavity dimensions
Ly = 1.0                               

Nx = 32
Ny = 16                                     #grid size
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
        
    


'''
Global iterations
'''
for iter in range(Niterg):
    #We need to make save to be a copy of the vort, otherwise changing vort 
    #will automatically result in changes of save!!!
    save = copy.deepcopy(vort)
    '''
    Jacobi updating of the stream function at the interior nodes
    '''
    for interior_itr in range(Niteri):
        for j in range(1, Ny):
            for i in range(1, Nx):
                res = (psi[i+1][j]+psi[i-1][j]+ beta*psi[i][j+1]+beta*psi[i][j-1]+Dxs*vort[i][j])/beta1-psi[i][j]
                psi[i][j] = psi[i][j] + relax*res
    
    
    '''
    Compute the vorticity at boundary grid points using the velocity boundary conditions
    (lower-order boundary conditions are commented out)
    '''
    
    '''
    Top and bottom walls
    '''
    for i in range(1, Nx):
        #The following two lines are lower order boundary conditions
        #vort[i][0] = 2.0*(psi[i][0]-psi[i][1])/Dys
        #vort[i][Ny] = 2.0*(psi[i][Ny]-psi[i][Ny-1])/Dys-2.0*Vlid/Dys
        
        vort[i][0] = (7.0*psi[i][0]-8.0*psi[i][1]+psi[i][2])/(2.0*Dys)
        vort[i][Ny] = (7.0*psi[i][Ny]-8.0*psi[i][Ny-1]+psi[i][Ny-2])/(2.0*Dys)-3.0*Vlid/Dy


    '''
    left and right walls
    '''
    for j in range(1, Ny):
        #The following two lines are lower order boundary conditions
        #vort[0][j] = 2.0*(psi[0][j]-psi[1][j])/Dxs
        #vort[Nx][j] = 2.0*(psi[Nx][j]-psi[Nx-1][j])/Dxs
        
        vort[0][j] = (7.0*psi[0][j]-8.0*psi[1][j]+psi[2][j])/(2.0*Dxs)
        vort[Nx][j] = (7.0*psi[Nx][j]-8.0*psi[Nx-1][j]+psi[Nx-2][j])/(2.0*Dxs)


    '''
    Compute the velocity at the interior grid points by central differences
    '''
    ux = np.zeros((Nx+1, Ny+1)).tolist()
    uy = np.zeros((Nx+1, Ny+1)).tolist()
    for j in range(1, Ny):
        for i in range(1, Nx):
            ux[i][j] = (psi[i][j+1]-psi[i][j-1])/Dy2
            uy[i][j] = - (psi[i+1][j]-psi[i-1][j])/Dx2


    '''
    Iterate on Poisson’s equation for the vorticity
    '''
    source = np.zeros((Nx+1, Ny+1)).tolist()
    for iteri in range(Niteri):
        for j in range(1, Ny):
            for i in range(1, Nx):
                source[i][j] = ux[i][j]*(vort[i+1][j]-vort[i-1][j])/Dx2+ uy[i][j]*(vort[i][j+1]-vort[i][j-1])/Dy2
                source[i][j] = -source[i][j]/nu
                res = (vort[i+1][j]+vort[i-1][j] + beta*vort[i][j+1]+beta*vort[i][j-1]+Dxs*source[i][j])/beta1-vort[i][j]
                vort[i][j] = vort[i][j] + relax * res
    
    
    '''
    Check the error tolerance
    '''
    cormax = 0.0
    for i in range(Nx+1):
        for j in range(Ny+1):
            res = abs(vort[i][j]-save[i][j])
            print(res)
            if res > cormax:
                cormax = res
        
    #End the iteration if the relative error between iterations is small enough
    if cormax < tol:
        break



'''
Plot the results
'''
xgr = []
ygr = []
for i in range(Nx+1):
    xgr.append(Dx*i)
for j in range(Ny+1):
    ygr.append(Dy*j)
    


plt.contour(xgr,ygr,np.array(psi).T.tolist(), 32)
# =============================================================================
# plt.xlabel('x','fontsize',15)
# plt.ylabel('y','fontsize',15)
# plt.zlabel('\psi','fontsize',15)
# =============================================================================
plt.axis([0, Lx, 0, Ly])












