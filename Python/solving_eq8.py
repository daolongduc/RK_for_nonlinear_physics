# -*- coding: utf-8 -*-
"""
This code is for Section 3, paper DOI: ---
The code solves the following differential equation:
    d(dx/dt)/dt + c/m* dx/dt*abs(dx/dt) + k/m*x -f0*sin(w*t) = 0
This differential equation represents the forced vibration of an sdof
system with quadratic damping. Where:
    m = mass
    k = stiffenss
    c = damping coefficient
    x = displacement
    t = time
    f0*sin(w*t) = applied harmonic force with amplitude f0 and angular
        frequency w
The system starts moving at initial displacement x0 and velocity v0
"""
import numpy as np # import nympy library for math functions
from scipy.integrate import solve_ivp # import solve_ipv function for solving an initial value problem for a system of ODEs.
import matplotlib.pyplot as plt # import ploting library
#=============================================
# Inputs:
m = 2.0
k = 50.0
c = 2.0*0.02*np.sqrt(k*m) # the formula is to generate a nice value of c. You can put the value of your problem here.
f0 = 5.0
w = 1.5*np.sqrt(k/m) # the formula is to generate a nice value of w. You can put the value of your problem here.
x0 = 0.0
v0 = 0.0
t0 = 0.0
tf = 20.0*2.0*np.pi*np.sqrt(m/k) #the formula is to generate a nice value of tf. You can put your desired value here.
#=============================================
# Define the system of ODEs:
def myode(t, y, m, k, c, f0, w):
    dxdt = y[1]
    dvdt = f0*np.sin(w*t)-c/m*y[1]*abs(y[1])-k/m*y[0]
    return (dxdt,dvdt)
#==============================================
# Process:
t_span = [t0,tf] # time span to be integrated
y0 = [x0,v0] # initial conditions
dtmax = 2.0*np.pi*np.sqrt(m/k)/100.0 # max step size of t
res = solve_ivp(myode, t_span, y0, method='RK45',
                args=(m,k,c,f0,w),
                rtol=1.0e-6, atol=1.0e-9, max_step = dtmax) # solve the odes
t = res.t # get time vector
x = res.y[0] # get x (the first unknown function of the odes)
v = res.y[1] # get v (the second unknown function of the odes)
plt.plot(t,x) # plot x vs t
plt.xlim([t0, tf]) # set x-limits
plt.grid(visible=True) # add grid
plt.xlabel('t') # add label to the horizontal axis of the plot
plt.ylabel('x') # add label to the vertical axis of the plot
plt.show() # show the plot