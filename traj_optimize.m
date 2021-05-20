%% Final Project:
matlabrc; clc; cla
addpath(genpath('src'))


%% Setup:
dt = 1;
dt_u1 = 2*60;
dt_u2 = 30;
duration = 93*60; 
tspan = dt:dt:duration;

% Physical Constants:
mu = 3.986004418*1e14; %(m^3/s^2) Earth Standard Gravitational Parameter


%% Simulate Truth Trajectories:
X_iss = zeros(6,length(tspan));
X_dragon = zeros(6,length(tspan));

% Initial angle:
angle = -90; %(deg)

% Define ISS:
a = (6378.137+415)*1000; %(m) Semi-major axis of the ISS
e = 0; %() Eccentricity
i = 52.1; %(deg) Inclination
[r,v] = kep2rv(a,e,i,0,0,angle,mu,true);
X_iss(:,1) = [r;v];

% Define Dragon:
[r,v] = kep2rv(a,0.0001,i+.0001,0,.001,angle-.0115,mu,true);
X_dragon(:,1) = [r;v];

% Convert to hill frame:
[r_hill,v_hill] = eci2hill(X_iss(1:3,1), X_iss(4:6,1),...
                           X_dragon(1:3,1), X_dragon(4:6,1));

% Initial setup:
X = zeros(6,length(tspan));
X(:,1) = [r_hill;v_hill];


%% Generate Optimal Trajectory:
terminal_1 = [600;0;0;
              0;0;0]; % Get into position for docking
terminal_2 = [0;0;0; 
              0;0;0]; % Perform actual docking maneuver

N1 = round(length(tspan)*(3/4));
N2 = round(length(tspan)*(1/4));

u1_0 = zeros(3*round(N1/dt_u1),1);

u2_0 = zeros(3*round(N2/dt_u2),1);

A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = [];
options = optimoptions(@fmincon,'Display','iter-detailed',...
                                'MaxFunctionEvaluations',10000,...
                                'OptimalityTolerance', 1e-4);
                            
% Optimize the first segment:
X0_1 = X(:,1);
u1 = fmincon(@(u) cost1(u,X0_1,mu,dt,dt_u1,N1,a,terminal_1),u1_0,A,b,Aeq,beq,lb,ub,nonlcon,options);
save('data/u1.mat','u1','dt','dt_u1','N1','terminal_1','X0_1')

% Optimize the second segment:
X0_2 = terminal_1;
u2 = fmincon(@(u) cost2(u,X0_2,mu,dt,dt_u2,N2,a,terminal_2),u2_0,A,b,Aeq,beq,lb,ub,nonlcon,options);
save('data/u2.mat','u2','dt','dt_u2','N2','terminal_2','X0_2')

%% Simulate first segment with Control:
load('data/u1.mat')

u_in = reshape(u1,3,[]);

% Number of runs to complete a burn
n = dt_u1/dt;

% Initial setup:
X_u1 = zeros(6,N1);
X_u1(:,1) = X0_1;
kk = 1;
for ii = 1:N1
    if mod(ii,n) == 0
        u_in_vec = u_in(:,kk);
        kk = kk+1;
    else
        u_in_vec = [0;0;0];
    end

    % Simulate dynamics:
    X_u1(:,ii+1) = rk4(@cweq,dt,X_u1(:,ii),mu,a,u_in_vec);
end

%% Simulate second segment with control:
load('data/u2.mat')
u_in = reshape(u2,3,[]);

% Number of runs to complete a burn
n = dt_u2/dt;

% Initial setup:
X_u2 = zeros(6,length(N2));
X_u2(:,1) = X0_2;
kk = 1;
for ii = 1:N2
    if mod(ii,n) == 0
        u_in_vec = u_in(:,kk);
        kk = kk+1;
    else
        u_in_vec = [0;0;0];
    end

    % Simulate dynamics:
    X_u2(:,ii+1) = rk4(@cweq,dt,X_u2(:,ii),mu,a,u_in_vec);
end
    
%% Simulate un-controlled:
X = zeros(6,N1);
X(:,1) = [r_hill;v_hill];
for ii = 1:N1
    X(:,ii+1) = rk4(@cweq,dt,X(:,ii),mu,a,[0;0;0]);
end

%% Plot the Results:
cla;
plot3(X(1,:),X(2,:),X(3,:),'b'); hold on
plot3(X_u1(1,:),X_u1(2,:),X_u1(3,:),'--r')
plot3(X_u2(1,:),X_u2(2,:),X_u2(3,:),'--r')
xlabel('Radial')
ylabel('In-track')
zlabel('Cross-track')
axis equal
grid on

save('data/traj.mat','X_u1','X_u2','dt','a','mu')