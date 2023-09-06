# -*- coding: utf-8 -*-
"""
This code is for Section 4.2, paper DOI: ---
This code analyzes the movement of a mass sliding on a friction cylinder surface.
"""
import numpy as np # import nympy library for math functions
from scipy.integrate import solve_ivp # import solve_ipv function for solving an initial value problem for a system of ODEs.
import matplotlib.pyplot as plt # import ploting library
#=============================================
# Inputs:
g = 9.81 # gravitational acceleration
m = 0.2 # mass of the sliding object
r = 0.5 # inner radius of the cylinder
mu = 0.5 # friction coefficient of the surface
theta0 = np.pi/2 # initial angle
omega0 = 0 # initial angular velocity
t0 = 0 # start time
tf = 1 # end time
abstol = 1.0e-6 # absolute tolerance
#=============================================
# Define the system of ODEs:
def myode(t, y, g, r, mu, abstol):
    theta = y[0]
    omega = y[1]
    if abs(omega)>=abstol:
        dtheta_dt = omega
        domega_dt = -mu*(omega**2+g/r*np.cos(theta))*np.sign(omega)-g/r*np.sin(theta)
    else:
        dtheta_dt = 0
        if abs(np.tan(theta)) <= mu:
            domega_dt = 0
        else:
            domega_dt = g/r*(mu*np.cos(abs(theta))-np.sin(abs(theta)))*np.sign(theta);
    return [dtheta_dt,domega_dt]
#==============================================
# Process:
t_span = [t0, tf] # time span to be integrated
y0 = [theta0, omega0] # initial condition
res = solve_ivp(myode, t_span, y0, method='RK45',
                args=(g, r, mu, abstol),
                rtol=1.0e-6, atol=abstol) # solve the odes
t = res.t # get time vector
theta = res.y[0] # get theta (the first unknown function of the odes)
omega = res.y[1] # get omega (the second unknown function of the odes)
N = m*np.square(omega)*r+m*g*np.cos(theta); # normal force
#==============================================
# Plot:
fig = plt.figure()
plt.grid(visible=True) # add grid
plt.plot(t,theta,'k-',linewidth=1) # plot angle
plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
plt.ylabel(r'$\theta \; (rad)$') # add label to the vertical axis of the plot
plt.xlim([t0,tf])
plt.show() # show the plot
fig = plt.figure()
plt.grid(visible=True) # add grid
plt.plot(t,omega,'k-',linewidth=1) # plot angular velocity
plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
plt.ylabel(r'$\dot{\theta} \; (rad/s)$') # add label to the vertical axis of the plot
plt.xlim([t0,tf])
plt.show() # show the plot
fig = plt.figure()
plt.grid(visible=True) # add grid
plt.plot(t,N,'k-',linewidth=1);
plt.plot(t,m*g*np.cos(theta),'m--',linewidth=1);
plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
plt.ylabel(r'$N \; (N)$') # add label to the vertical axis of the plot
plt.legend(['dynamic','static']);