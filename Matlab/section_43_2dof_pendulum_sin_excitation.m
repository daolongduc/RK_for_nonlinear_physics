%% This code is for Section 4.3, paper DOI: ---
% This code analyzes the movement of a spring pendulum subjected to a pivot
% movement.
%%
function section_43_2dof_pendulum_sin_excitation() % you can remove this line and line 40 in version 18 or later.
%% Inputs
g = 9.81; % acceleration of gravity
L = 0.5; % length of the spring at rest
k = 50; % stiffness of the spring
m = 0.1; % mass of the bob
c1 = 0.02; % linear damping coefficient
c2 = 0.02; % squared damping coefficient
xs0 = 0; % amplitude of the pivot movement. The pivot movement follows a sin function: xs = xs0*sin(Omega*t)
Omega = 5; % angular frequency of the pivot movement
x0 = 0; % initial elongation of the spring
v0 = 0; % initial velocity along the spring
theta0 = pi/2; % initial angle of the pendulum
omega0 = 0; % initial angular velocity of the pendulum
t0 = 0; % start time
tf = 1.7; % end time
%% Process
ode_fun = @(t, y) myode(t, y, g, L, k, m, c1, c2, xs0, Omega); % get the function handle of the ode function (defined below)
t_span = [t0, tf]; % time span to be integrated
y0 = [x0,v0,theta0,omega0]; % initial condition
opts = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-9); % define tolerance
[t, y] = ode45(ode_fun, t_span, y0, opts); % solve the ode
%% Plot
figure; % create a new figure
plot(t,y(:,1),'k-','linewidth',1); % plot x
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
xlabel('$t (s)$','interpreter','latex'); % add label to the horizontal axis
ylabel('$x (m)$','interpreter','latex'); % add label to the vertical axis
figure; % create a new figure
plot(t,y(:,3),'k-','linewidth',1); % plot x
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
xlabel('$t (s)$','interpreter','latex'); % add label to the horizontal axis
ylabel('$\theta (rad)$','interpreter','latex'); % add label to the vertical axis
end % you can remove this line and line 5 in version 18 or later.
%%
function dydt = myode(t, y, g, L, k, m, c1, c2, xs0, Omega) % define the system of ode
xsdot = Omega*xs0*cos(Omega*t);
xsddot = -Omega^2*xs0*sin(Omega*t);
x = y(1); v=y(2); theta=y(3); omega=y(4);
dxdt = v;
dvdt = g*cos(theta)-k/m*x-(c1+c2*abs(v+xsdot*sin(theta)))/m*(v+xsdot*sin(theta))-xsddot*sin(theta);
dthetadt = omega;
R = L+x;
domegadt = -g/R*sin(theta)-(c1+c2*abs(R*omega+xsdot*cos(theta)))/m*(omega+xsdot*cos(theta)/R)-xsddot/R*cos(theta);
dydt = [dxdt;dvdt;dthetadt;domegadt];
end