%% Final Project:
matlabrc; clc; cla
addpath(genpath('src'))


%% Setup:
dt = 1;
duration = 5*93*60; 
tspan = dt:dt:duration;

% Physical Constants:
mu = 3.986004418*1e14; %(m^3/s^2) Earth Standard Gravitational Parameter


%% Simulate Truth Trajectories:
X_iss = zeros(6,length(tspan));
X_dragon = zeros(6,length(tspan));
r_hill = zeros(3,length(tspan));
v_hill = zeros(3,length(tspan));

% Initial angle:
angle = 130; %(deg)

% Define ISS:
a = (6378.137+415)*1000; %(m) Semi-major axis of the ISS
e = 0; %() Eccentricity
i = 52.1; %(deg) Inclination
[r,v] = kep2rv(a,e,i,0,0,angle,mu,true);
X_iss(:,1) = [r;v];


% Define Dragon:
[r,v] = kep2rv(a,0.0001,i+.001,0,0,angle+.01,mu,true);
X_dragon(:,1) = [r;v];


% Initial control:
u = [0;0;0];


% Simulate:
for ii = 1:length(tspan)    
    % Convert to hill frame:
    [r_hill(:,ii),v_hill(:,ii)] = eci2hill(X_iss(1:3,ii), X_iss(4:6,ii),...
                                           X_dragon(1:3,ii), X_dragon(4:6,ii));
   
    % Propagate inertial orbits:
    X_iss(:,ii+1)    = rk4(@orbitalDynamics,dt,X_iss(:,ii),mu,zeros(3,1));
    X_dragon(:,ii+1) = rk4(@orbitalDynamics,dt,X_dragon(:,ii),mu,u);
end


%% Simulate the Linearized Dynamics:

% Initial setup:
X = zeros(6,length(tspan));
X(:,1) = [r_hill(:,1);v_hill(:,1)];

% Simulate trajectory:
for ii = 1:length(tspan)-1
    X(:,ii+1) = rk4(@cweq,dt,X(:,ii),mu,a,u);
end

% Plot the results as a comparison:
% plot3(r_hill(1,:),r_hill(2,:),r_hill(3,:),'b'); hold on
% plot3(X(1,:),X(2,:),X(3,:),'--r'); hold on
% axis equal
% grid on
% rotate3d on
% 
% figure()
% subplot(3,1,1)
%     plot(r_hill(1,:) - X(1,:));
% subplot(3,1,2)
%     plot(r_hill(2,:) - X(2,:));
% subplot(3,1,3)
%     plot(r_hill(3,:) - X(3,:));
    
    
%% Calculate LQR Gain:
% Calculate mean motion:
n = sqrt(mu/a^3);

% Define linear state space model:
A = [  0   0   0    1    0    0;
       0   0   0    0    1    0;
       0   0   0    0    0    1;
     3*n^2 0   0    0   2*n   0;
       0   0   0  -2*n   0    0;
       0   0 -n^2   0    0    0];

B = [0 0 0;
     0 0 0;
     0 0 0;
     1 0 0;
     0 1 0;
     0 0 1];
 
Q1 = diag([1e0*ones(1,3), 1e1*ones(1,3)]);
R1 = 1e8*eye(3);
Q2 = diag([1e3*ones(1,3), 1e3*ones(1,3)]);
R2 = 1e6*eye(3);

[K1,~,~] = lqr(A,B,Q1,R1);
[K2,~,~] = lqr(A,B,Q2,R2);

%% Simulate with LQR Control:
% Initial setup:
X = zeros(6,length(tspan));
X(:,1) = [r_hill(:,1);v_hill(:,1)];
x_targ1 = [500; 0; 0; 
           0;   0; 0];
x_targ2 = [0;0;0;
           0;0;0];
      
sig_r = 1;
sig_v = 0.1;

stable_count = 0;
dock = false;

% Simulate trajectory:
for ii = 1:length(tspan)-1
    % Perform estimation:
    x_hat = X(:,ii) + [sig_r*randn(3,1); sig_v*randn(3,1)];
    
    % Determine if stable around first waypoint:
    if norm(x_hat - x_targ1) < 3*sig_r
        stable_count = stable_count+1;
    end
    
    if stable_count > 120
        dock = true;
    end
    
    if dock
        x_targ = x_targ2;
        K = K2;
    else
        x_targ = x_targ1;
        K = K1;
    end
    
    % Calculate control input:
    u = -K*(x_hat - x_targ);
    X(:,ii+1) = rk4(@cweq,dt,X(:,ii),mu,a,u);
end

plot3(r_hill(1,:),r_hill(2,:),r_hill(3,:),'b'); hold on
plot3(X(1,:),X(2,:),X(3,:),'--r')
% plot3(X(1,1),X(2,1),X(3,1),'.k','MarkerSize',20)
xlabel('Radial')
ylabel('In-track')
zlabel('Cross-track')
axis equal
grid on

%% Animate:
un = plot3(r_hill(1,1),r_hill(2,1),r_hill(3,1),'.b','MarkerSize',20); hold on
tr = plot3(X(1,1),X(2,1),X(3,1),'.r','MarkerSize',20);

for ii = 1:length(tspan)
    set(un,'XData',r_hill(1,ii),'YData',r_hill(2,ii),'ZData',r_hill(3,ii))
    set(tr,'XData',X(1,ii),'YData',X(2,ii),'ZData',X(3,ii))
    drawnow
    pause(1/30)
end