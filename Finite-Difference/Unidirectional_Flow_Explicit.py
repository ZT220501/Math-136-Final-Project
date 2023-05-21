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


'''
Initialze Values
'''
h = 1
mu = 0.6 
#rho = 5.5
rho = 20
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
oscillation_const = 1
def V1(t):
    return np.sin(oscillation_const*t)

def V2(t):
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
Run Iterations
'''
u_total = [u]
t = 0
for step in range(nstep):
    t += Dt
    u = u_total[-1]
    unew = [V1(t)]
    for i in range(1, N):
        unew.append(al*u[i-1] + (1-2*al)*u[i] + al*u[i+1] + Dt*(-dpdx/rho+gx))
    unew.append(V2(t))
    u_total.append(unew)



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
 
# what will our line dataset
# contain?
def init():
    line.set_data([], [])
    return line,
 
# initializing empty values
# for x and y co-ordinates
xdata, ydata = [], []
 
# animation function
def animate(i):

    line.set_data(u_total[i], y)
    line.set_marker('o')
     
    return line,
 
# calling the animation function    
anim = animation.FuncAnimation(fig, animate,
                            init_func = init,
                            frames = nstep,
                            interval = 20,
                            blit = True)
 
# saves the animation in our desktop
anim.save('uniflow_explicit_oscillation1_low_viscosity.gif', writer = 'ffmpeg', fps = 30)





































