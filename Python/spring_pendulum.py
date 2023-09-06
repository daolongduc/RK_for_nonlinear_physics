# -*- coding: utf-8 -*-
"""
This code is for Section 4.3, paper DOI: ---
This code analyzes the movement of a spring pendulum subjected to a pivot movement.
"""
import numpy as np # import nympy library for math functions
from scipy.integrate import solve_ivp # import solve_ipv function for solving an initial value problem for a system of ODEs.
import matplotlib.pyplot as plt # import ploting library
#=============================================
# Inputs:
g = 9.81 # acceleration of gravity
L = 0.5 # length of the spring at rest
k = 50 # stiffness of the spring
m = 0.1 # mass of the bob
c1 = 0.02 # linear damping coefficient
c2 = 0.02 # squared damping coefficient
xs0 = 0 # amplitude of the pivot movement. The pivot movement follows a sin function: xs = xs0*sin(Omega*t)
Omega = 5 # angular frequency of the pivot movement
x0 = 0 # initial elongation of the spring
v0 = 0 # initial velocity along the spring
theta0 = np.pi/2 # initial angle of the pendulum
omega0 = 0 # initial angular velocity of the pendulum
t0 = 0 # start time
tf = 1.7 # end time
#=============================================
# Define the system of ODEs:
def myode(t, y, g, L, k, m, c1, c2, xs0, Omega):
    xsdot = Omega*xs0*np.cos(Omega*t)
    xsddot = -Omega**2*xs0*np.sin(Omega*t);
    x, v, theta, omega = y
    dxdt = v
    dvdt = g*np.cos(theta)-k/m*x-(c1+c2*abs(v+xsdot*np.sin(theta)))/m*(v+xsdot*np.sin(theta))-xsddot*np.sin(theta)
    dthetadt = omega
    R = L+x
    domegadt = -g/R*np.sin(theta)-(c1+c2*abs(R*omega+xsdot*np.cos(theta)))/m*(omega+xsdot*np.cos(theta)/R)-xsddot/R*np.cos(theta)
    return [dxdt,dvdt,dthetadt,domegadt]
#==============================================
# Process:
t_span = [t0, tf] # time span to be integrated
y0 = [x0,v0,theta0, omega0] # initial condition
res = solve_ivp(myode, t_span, y0, method='RK45',
                args=(g, L, k, m, c1, c2, xs0, Omega),
                rtol=1.0e-6, atol=1.0e-9) # solve the odes
t = res.t # get time vector
x = res.y[0] # get x (the first unknown function of the odes)
theta = res.y[2] # get theta (the third unknown function of the odes)
#==============================================
# Plot:
fig = plt.figure()
plt.grid(visible=True) # add grid
plt.plot(t,x,'k-',linewidth=1) # plot x
plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
plt.ylabel(r'$x \; (m)$') # add label to the vertical axis of the plot
plt.xlim([t0,tf])
plt.show() # show the plot
fig = plt.figure()
plt.grid(visible=True) # add grid
plt.plot(t,theta,'k-',linewidth=1) # plot displaced angle
plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
plt.ylabel(r'$\theta \; (rad)$') # add label to the vertical axis of the plot
plt.xlim([t0,tf])
plt.show() # show the plot
