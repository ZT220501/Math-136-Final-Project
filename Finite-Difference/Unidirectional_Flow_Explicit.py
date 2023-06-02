# -*- coding: utf-8 -*-
"""
Created on Sun Apr  9 20:21:44 2023

@author: zheng
"""


'''
Code for running explicit finite difference method for solving the unidirecitonal flow
This code is a modification of the code in C. Pozrikidis's book
'''





import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import time
import copy


'''
Initialze Values
'''
h = 1
mu = 0.6 
rho = 0.5
#rho = 20
N = 32
dpdx = -2
#gx=0, represents that the channel is completely horizontal
gx = 0
V1 = 0
V2 = 0

#The diffusion number
#The method is stable iff alpha<=0.5
al = 0.49                           #Alpha value                  
nstep = 2000                        #Step numbers

#Prepare
nu = mu/rho                         #kinematic viscosity
Dy = h/N
Dt = al*Dy*Dy/nu


'''
Velocity for upper and lower bound
'''
oscillation_const = 100
def V1(t):
    return 0
    #return np.sin(oscillation_const*t)

def V2(t):
    return 0
    #return np.cos(oscillation_const*t)



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
Run Iterations
'''
u_total = [u]
t = 0
#Track for running time
start = time.time()
for step in range(nstep):
    t += Dt
    u = u_total[-1]
    unew = [V1(t)]
    for i in range(1, N):
        unew.append(al*u[i-1] + (1-2*al)*u[i] + al*u[i+1] + Dt*(-dpdx/rho+gx))
    unew.append(V2(t))
    u_total.append(copy.deepcopy(unew))

print("Explicit finite difference running time: " + str(time.time()-start) + " seconds.")




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
def u_analytic(t,y):
    C = -dpdx/rho + gx
    r = 0
    for n in range(1,500):
      if n % 2 == 1:
        p1 = 4*C/(n**3*(np.pi)**3*nu) - 4*C/(n**3*(np.pi)**3*nu)*np.exp(-n**2*(np.pi**2)*nu*t)
        r += p1 * np.sin(n*np.pi*y)
    return r 
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
# =============================================================================
# def u_analytic(t,y):
#     C = -dpdx/rho + gx
#     r = 0
#     gamma = 0.5 + 0.5*y
#     for n in range(1,500):
#       if n % 2 == 1:
#         p = 4*C/(n**3*(np.pi)**3*nu) - (3/(n*np.pi)+4*C/(n**3*(np.pi)**3*nu))*np.exp(-n**2*(np.pi)**2*nu*t)
#         r += p * np.sin(n*np.pi*y)
#       else:
#         p = np.exp(-n**2*(np.pi)**2*nu*t)/(n*np.pi)
#         r += p * np.sin(n*np.pi*y)
#     return r + gamma
# =============================================================================

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





































