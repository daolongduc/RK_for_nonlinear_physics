%% This code is for Section 4.1, paper DOI: ---
% The code solves the following differential equation:
% d(d(theta)/dt)/dt +
% (c1+c2*abs(xsddot*cos(theta)+L*d(theta)/dt))*(xsdot*cos(theta)+L*d(theta)/dt)/(m*L)
% + g/L*sin(theta) + xsddot/L*cos(theta) = 0
% The above equation is the governing differential equation of a nonlinear
% pendulum subjected to a pivot motion, where:
%   theta = displaced angle
%   L = length of the pendulum
%   m = mass of the bob
%   c1, c2 = damping coefficients
%   xs = xs0*sin(Omega*t) = pivot motion
%   xsdot, xsddot = velocity and acceleration of the pivot motion
%   g = gravitational acceleration
% The system starts moving at initial angle theta0 and angular velocity
% omega0 (= d(theta)/dt at t=0)
%%
function section_41_nonlinear_pendulum_free_vibration() % you can remove this line and line 48 in version 18 or later.
%% Inputs
g = 9.81; % gravitational acceleration
m = 0.2; % mass of the bob
L = 0.2; % length of the pendulum
c1 = 0.0; % linear damping coefficient
c2 = 0.0; % square damping coefficient
xs0 = 0.0; % pivot moving amplitude, xs = xs0*sin(Omega*t)
Omega = 5;
theta0 = pi/2; % initial angle
omega0 = 0; % initial angular velocity
t0 = 0; % start time
tf = 3; % end time
%% Process
% nonlinear solution:
[t, theta_non, ~] = analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf);
% linear solution:
omega_n = sqrt(g/L); % natural angular frequency
A = theta0;
B = omega0/omega_n;
theta_lin = A*cos(omega_n*t)+B*sin(omega_n*t);
% plot comparision
figure; % create a new figure
hold on; % keep the old objects before plotting a new one
grid on; % add grid to the figure
plot(t,theta_non,'m:','linewidth',1.5); % plot nonlinear solution
plot(t,theta_lin,'k-','linewidth',0.5); % plot linear solution
xlabel('$t (s)$','interpreter','latex'); % add label to the horizontal axis
ylabel('$\theta (rad)$','interpreter','latex'); % add label to the vertical axis
legend('nonlinear','linear'); % add legend
end % you can remove this line and line 18 in version 18 or later.
%% Define the system of odes:
function dydt = myode(t, y, g, m, L, c1, c2, xs0, Omega) % define the system of odes
xsdot = Omega*xs0*cos(Omega*t);
xsddot = -Omega^2*xs0*sin(Omega*t);
theta = y(1);
omega = y(2);
dtheta_dt = omega;
domega_dt = -(c1+c2*abs(xsdot*cos(theta)+L*omega))*(xsdot*cos(theta)+L*omega)/m/L - g/L*sin(theta)-xsddot/L*cos(theta);
dydt = [dtheta_dt;domega_dt];
end
%% Analyze for responses
function [t, theta, omega] = analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf)
ode_fun = @(t, y) myode(t, y, g, m, L, c1, c2, xs0, Omega); % get the function handle of the ode function (defined below)
t_span = [t0, tf]; % time span to be integrated
y0 = [theta0; omega0]; % initial condition
opts = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-9); % define tolerance
[t, y] = ode45(ode_fun, t_span, y0, opts); % solve the ode
theta = y(:,1); % displaced angle
omega = y(:,2); % angular velocity
end
