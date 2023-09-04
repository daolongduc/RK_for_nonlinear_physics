%% This code is for Section 4.1, paper DOI: ---
% This code demonstrates the chaos of a nonlinear pendulum whose
% governing differential equation is:
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
function section_41_nonlinear_pendulum_chaos() % you can remove this line and line 51 in version 18 or later.
%% Inputs
g = 9.81; % gravitational acceleration
m = 0.2; % mass of the bob
L = 0.2; % length of the pendulum
c1 = 0.03; % linear damping coefficient
c2 = 0.03; % square damping coefficient
xs0List = [0.075,0.085]; % pivot moving amplitude, xs = xs0*sin(Omega*t)
Omega = 14;
theta0List = [0.1, 0.1001]; % initial angle
omega0 = 0; % initial angular velocity
nAnal = 200; % number of cycles to be analyzed.
%% Process
T = 2*pi/Omega;
t0 = 0;
tf = nAnal*T;
linetype = {'m:','k-'};
linewidth = [1.5,0.5];
for xs0 = xs0List
    figure;
    hold on;
    grid on;
    cnt = 1;
    for theta0 = theta0List
        [t, theta, ~] = analyze(g, m, L, c1, c2, xs0, Omega, theta0, omega0, t0, tf);
        plot(t,theta,linetype{cnt},'linewidth',linewidth(cnt));
        cnt = cnt + 1;
    end
    xlabel('$t (s)$','interpreter','latex');
    ylabel('$\theta (rad)$','interpreter','latex');
    legend(['\theta(0)=',num2str(theta0List(1)),' rad'],['\theta(0)=',num2str(theta0List(2)),' rad']);
end
end % you can remove this line and line 19 in version 18 or later.
%% Define the system of odes:
function dydt = myode(t, y, g, m, L, c1, c2, xs0, Omega) % define the system of ode
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
