# -*- coding: utf-8 -*-
"""
This code is for Section 4.1, paper DOI: ---
The code solves the following differential equation:
d(d(theta)/dt)/dt +
(c1+c2*abs(xsddot*cos(theta)+L*d(theta)/dt))*(xsdot*cos(theta)+L*d(theta)/dt)/(m*L)
+ g/L*sin(theta) + xsddot/L*cos(theta) = 0
The above equation is the governing differential equation of a nonlinear
pendulum subjected to a pivot motion, where:
  theta = displaced angle
  L = length of the pendulum
  m = mass of the bob
  c1, c2 = damping coefficients
  xs = xs0*sin(Omega*t) = pivot motion
  xsdot, xsddot = velocity and acceleration of the pivot motion
  g = gravitational acceleration
The system starts moving at initial angle theta0 and angular velocity
omega0 (= d(theta)/dt at t=0)
"""
import numpy as np # import nympy library for math functions
from scipy.integrate import solve_ivp # import solve_ipv function for solving an initial value problem for a system of ODEs.
import matplotlib.pyplot as plt # import ploting library
#=============================================
# Inputs:
g = 9.81; # gravitational acceleration
m = 0.2; # mass of the bob
L = 0.2; # length of the pendulum
c1 = 0.0; # linear damping coefficient
c2 = 0.0; # square damping coefficient
xs0 = 0.0; # pivot moving amplitude, xs = xs0*sin(Omega*t)
Omega = 5; # angular frequency of the support motion.
theta0 = np.pi/2; # initial angle
omega0 = 0; # initial angular velocity
t0 = 0; # start time
tf = 3; # end time
#=============================================
# Define the system of ODEs:
def myode(t, y, g, m, L, c1, c2, xs0, Omega):
    xsdot = Omega*xs0*np.cos(Omega*t)
    xsddot = -Omega**2*xs0*np.sin(Omega*t)
    theta = y[0]
    omega = y[1]
    dtheta_dt = omega
    domega_dt = -(c1+c2*abs(xsdot*np.cos(theta)+L*omega))*(xsdot*np.cos(theta)+L*omega)/m/L - g/L*np.sin(theta)-xsddot/L*np.cos(theta)
    return [dtheta_dt,domega_dt]
#==============================================
# Analyze for responses
def analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf):
    t_span = [t0, tf] # time span to be integrated
    y0 = [theta0, omega0] # initial condition
    dtmax = 2.0*np.pi*np.sqrt(L/g)/100.0 # max step size of t
    res = solve_ivp(myode, t_span, y0, method='RK45',
                    args=(g, m, L, c1, c2, xs0, Omega),
                    rtol=1.0e-6, atol=1.0e-9, max_step = dtmax) # solve the odes
    t = res.t # get time vector
    theta = res.y[0] # get theta (the first unknown function of the odes)
    omega = res.y[1] # get omega (the second unknown function of the odes)
    return [t,theta,omega]
#==============================================
# Process:
# nonlinear solution:
t, theta_non, omega_non = analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf)
# linear solution:
omega_n = np.sqrt(g/L) # natural angular frequency
A = theta0;
B = omega0/omega_n;
theta_lin = A*np.cos(omega_n*t)+B*np.sin(omega_n*t)

plt.plot(t,theta_non,'m:',linewidth=2) # plot theta_non vs t
plt.plot(t,theta_lin,'k-',linewidth=1) # plot theta_lin vs t
plt.xlim([t0, tf]) # set x-limits
plt.grid(visible=True) # add grid
plt.xlabel('t') # add label to the horizontal axis of the plot
plt.ylabel(r'$\theta$') # add label to the vertical axis of the plot
plt.legend(['nonlinear','linear'],loc='upper right')
plt.show() # show the plot