%% This code is for Section 3, paper DOI: ---
% The code solves the following differential equation:
% d(dx/dt)/dt + c/m* dx/dt*abs(dx/dt) + k/m*x -f0*sin(w*t) = 0
% This differential equation represents the forced vibration of an sdof
% system with quadratic damping. Where:
%   m = mass
%   k = stiffenss
%   c = damping coefficient
%   x = displacement
%   t = time
%   f0*sin(w*t) = applied harmonic force with amplitude f0 and angular
%       frequency w
% The system starts moving at initial displacement x0 and velocity v0
%%
function section_3_solving_eq8() % you can remove this line and line 39 in version 18 or later.
%% Inputs
m = 2;
k = 50;
c = 2*0.02*sqrt(k*m); % the formula is to generate a nice value of c. You can put any value here.
f0 = 5;
w = 1.5*sqrt(k/m); % the formula is to generate a nice value of w. You can put any value here.
x0 = 0;
v0 = 0;
t0 = 0;
tf = 20*2*pi*sqrt(m/k); % the formula is to generate a nice value of tf. You can put any value here.
%% Process
ode_fun = @(t, y) myode(t, y, m, k, c, f0, w); % get the function handle of the ode function (defined below)
t_span = [t0, tf]; % time span to be integrated
y0 = [x0; v0]; % initial condition
opts = odeset('RelTol', 1.0e-6, 'AbsTol', 1.0e-9); % define tolerance
[t, y] = ode45(ode_fun, t_span, y0, opts); % solve the ode
x = y(:,1); % get displacement history
%% Plot
figure; % create a new figure
plot(t,x); % plot x
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
xlabel('t'); % add label to the horizontal axis
ylabel('x'); % add label to the vertical axis
end % you can remove this line and line 15 in version 18 or later.
%%
function dydt = myode(t, y, m, k, c, f0, w) % define the system of ode
x = y(1);
v = y(2);
dydt = [v;
        f0*sin(w*t)-c/m*v*abs(v)-k/m*x];
end
