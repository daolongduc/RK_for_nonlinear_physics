%% This code is for Section 4.2, paper DOI: ---
% This code analyzes the movement of a mass sliding on a friction cylinder
% surface.
%%
function section_42_sliding_on_cylinder_surface() % you can remove this line and line 53 in version 18 or later.
%% Inputs
g = 9.81; % gravitational acceleration
m = 0.2; % mass of the sliding object
r = 0.5; % inner radius of the cylinder
mu = 0.5; % friction coefficient of the surface
theta0 = pi/2; % initial angle
omega0 = 0; % initial angular velocity
t0 = 0; % start time
tf = 1; % end time
abstol = 1.0e-6; % absolute tolerance
%% Process
ode_fun = @(t, y) myode(t, y, g, r, mu, abstol); % get the function handle of the ode function (defined below)
t_span = [t0, tf]; % time span to be integrated
y0 = [theta0; omega0]; % initial condition
opts = odeset('RelTol',1.0e-6,'AbsTol', abstol); % define tolerance
[t, y] = ode45(ode_fun, t_span, y0, opts); % solve the ode
theta = y(:,1);
omega = y(:,2);
%% Plot
figure; % create a new figure
plot(t,theta,'k-','linewidth',1); % plot angle
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
xlabel('t (s)','interpreter','latex'); % add label to the horizontal axis
ylabel('$\theta (rad)$','interpreter','latex'); % add label to the vertical axis
figure; % create a new figure
plot(t,omega,'k-','linewidth',1); % plot x
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
xlabel('t (s)','interpreter','latex'); % add label to the horizontal axis
ylabel('$\dot{\theta} (rad/s)$','interpreter','latex'); % add label to the vertical axis
figure;
hold on; % keep what have been plotted before addiing/plotting new objects
grid on; % add grid
N = m*omega.^2*r+m*g*cos(theta);
plot(t,N,'k-','linewidth',1);
plot(t,m*g*cos(theta),'m--','linewidth',1);
xlabel('t (s)','interpreter','latex'); % add label to the horizontal axis
ylabel('$N (N)$','interpreter','latex'); % add label to the vertical axis
legend('dynamic','static');
Wf = 0; % Work of friction force
Ff = mu*N;
for iChay=1:length(Ff)-1
    Wf = Wf + (Ff(iChay)+Ff(iChay+1))/2*r*abs(theta(iChay+1)-theta(iChay));
end
disp(['Work dissipated by friction=', num2str(Wf)]);
disp(['Potential energy difference=',num2str(m*g*r*cos(theta(end)))]);
end % you can remove this line and line 5 in version 18 or later.
%%
function dydt = myode(t, y, g, r, mu, abstol) % define the system of ode
theta = y(1);
omega = y(2);
if abs(omega)>=abstol
    dtheta_dt = omega;
    domega_dt = -mu*(omega^2+g/r*cos(theta))*sign(omega)-g/r*sin(theta);
else
    dtheta_dt = 0;
    if abs(tan(theta)) <= mu
        domega_dt = 0;
    else
        domega_dt = g/r*(mu*cos(abs(theta))-sin(abs(theta)))*sign(theta);
    end
end
dydt = [dtheta_dt;domega_dt];
end
