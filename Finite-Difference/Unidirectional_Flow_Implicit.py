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
import matplotlib.animation as animation
import copy
import time

'''
Initial Values
'''
h = 1.0
mu = 0.6
rho = 0.5
N = 32
dpdx = -2.0
gx = 0

#The diffusion number
alpha = 0.4               
nstep = 2000                            #Number of steps
nu = mu/rho
Dy = h/N
Dt = alpha*Dy*Dy/nu


'''
Velocity for upper and lower bound
'''
oscillation_const = 10
def V1(t):
    #return 0.5
    return np.sin(oscillation_const*t)

def V2(t):
    #return 1
    return np.cos(oscillation_const*t)

'''
Grid and initial condition
'''
y = []
u = []
for i in range(N+1):    
    y.append(i*Dy)
    u.append(0)
    
u[0] = V1(0)
u[N] = V2(0)
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
u_total = [u]

t = 0
start = time.time()
for step in range(nstep):
    u = u_total[-1]
    s = []
    for i in range(N-1):
        s.append(u[i+1] + Dt * (-dpdx/rho + gx))
    t = t + Dt
    
    u_new = copy.deepcopy(u)
    u_new[0] = V1(t)
    u_new[N] = V2(t)
    s[0] = s[0] + alpha * u_new[0]
    s[N-2] = s[N-2] + alpha * u_new[N]
    sol = thomas(N-1,atr,btr,ctr,s);
    for i in range(1, N):
        u_new[i] = sol[i-1]
        
    u_total.append(copy.deepcopy(u_new))
        
print("Implicit finite difference running time: " + str(time.time()-start) + " seconds.")

'''
Plot the result animation
'''
# creating a blank window
# for the animation
fig = plt.figure()
axis = plt.axes(xlim =(-5, 5),
                ylim =(0, h))
axis.set_xlabel('u')
axis.set_ylabel('y')
 
line, = axis.plot([], [], lw = 2)

 
# animation function
def animate(i):

    line.set_data(u_total[i], y)
    line.set_marker('o')
     
    return line,
 
# calling the animation function    
anim = animation.FuncAnimation(fig, animate,
                            frames = nstep,
                            interval = 20,
                            repeat = False)
 
# saves the animation in our desktop
anim.save('test.gif', writer = 'ffmpeg', fps = 30)




print("I've finished computing.")

'''
Test by comparing to the analytic solution
'''
#Analytic solution to the constant case
# =============================================================================
# def u_analytic(t,y):
#     C = -dpdx/rho + gx
#     r = 0
#     for n in range(1,500):
#       if n % 2 == 1:
#         p1 = 4*C/(n**3*(np.pi)**3*nu) - 4*C/(n**3*(np.pi)**3*nu)*np.exp(-n**2*(np.pi**2)*nu*t)
#         r += p1 * np.sin(n*np.pi*y)
#     return r 
# =============================================================================
# =============================================================================
# def u_analytic(t,y):
#     C = -dpdx/rho + gx
#     k = 10
#     gamma = (1-y) * np.sin(k*t) + y * np.cos(k*t)
#     r = 0
#     for n in range(1,500):
#       den = n*np.pi*(n**4*(np.pi)**4*nu**2+k**2)
#       if n % 2 == 1:
#         s1 = 2*k*np.cos(k*t)*(n**2*(np.pi)**2*nu+k) / den
#         s2 = 2*k*np.sin(k*t)*((n**2)*(np.pi)**2*nu-k) / den
#         s3 = 4*C/(n**3*(np.pi)**3*nu)
#         p1 = -s1 + s2 + s3
#         p2 = -2/(n*np.pi) + 2*k*(n**2*(np.pi)**2*nu+k) / den - s3
#         p3 = p2 * np.exp(-n**2*(np.pi)**2*nu*t)
#         u = (p1 + p3) * np.sin(n*np.pi*y)
#         r += u
#       else:
#         s1 = 2*k*np.cos(k*t)*(n**2*(np.pi)**2*nu-k) / den
#         s2 = 2*k*np.sin(k*t)*(n**2*(np.pi)**2*nu+k) / den
#         p1 = -s1 - s2
#         p2 = 2/(n*np.pi) + 2*k*(n**2*(np.pi)**2*nu-k) / den
#         p3 = p2 * np.exp(-n**2*(np.pi)**2*nu*t)
#         u = (p1 + p3) * np.sin(n*np.pi*y)
#         r += u
#     result = gamma + r
#     return result
# =============================================================================
def u_analytic(t,y):
    C = -dpdx/rho + gx
    r = 0
    gamma = 0.5 + 0.5*y
    for n in range(1,500):
      if n % 2 == 1:
        p = 4*C/(n**3*(np.pi)**3*nu) - (3/(n*np.pi)+4*C/(n**3*(np.pi)**3*nu))*np.exp(-n**2*(np.pi)**2*nu*t)
        r += p * np.sin(n*np.pi*y)
      else:
        p = np.exp(-n**2*(np.pi)**2*nu*t)/(n*np.pi)
        r += p * np.sin(n*np.pi*y)
    return r + gamma

u_anal = []
for i in range(nstep+1):
    u_actual_store = []
    for j in range(N+1):
        u_actual_store.append(u_analytic(i*Dt, j*Dy))
    u_anal.append(u_actual_store)


err = []
for i in range(len(u_total)):
    u_total[i] = np.array(u_total[i])
    u_anal[i] = np.array(u_anal[i])
    diff = u_total[i] - u_anal[i]
    d = np.sqrt(Dy * np.dot(diff, diff))
    err.append(d)
err = np.array(err)
print("Maximum error:", np.amax(err))
print("Minimum error:", np.amin(err))


    






    
    
    
    
    
    
    