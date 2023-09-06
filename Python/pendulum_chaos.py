# -*- coding: utf-8 -*-
"""
This code is for Section 4.1, paper DOI: ---
This code demonstrates the chaos of a nonlinear pendulum whose
governing differential equation is:
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
g = 9.81 # gravitational acceleration
m = 0.2 # mass of the bob
L = 0.2 # length of the pendulum
c1 = 0.03 # linear damping coefficient
c2 = 0.03 # square damping coefficient
xs0List = [0.075, 0.085] # pivot moving amplitudes, xs = xs0*sin(Omega*t)
Omega = 14
theta0List = [0.1,0.1001] # initial angles
omega0 = 0 # initial angular velocity
nAnal = 200 # number of cycles to be analyzed
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
T = 2*np.pi/Omega
t0 = 0
tf = nAnal*T
linetype = ['m:','k-']
linewidth = [1.5, 0.5]
for xs0 in xs0List:
    fig = plt.figure()
    plt.grid(visible=True) # add grid
    cnt = 0
    for theta0 in theta0List:
        t, theta, omega = analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf)
        plt.plot(t,theta,linetype[cnt],linewidth=linewidth[cnt]) # plot theta_non vs t
        cnt +=1
    plt.xlabel(r'$t \; (s)$') # add label to the horizontal axis of the plot
    plt.ylabel(r'$\theta \; (rad)$') # add label to the vertical axis of the plot
    plt.xlim([t0,tf])
    plt.show() # show the plot