%% Final Project:
matlabrc; clc; cla
addpath(genpath('src'))


%% Setup:
dt = 1;
duration = 5*93*60; 
tspan = dt:dt:duration;

% Physical Constants:
mu = 3.986004418*1e14; %(m^3/s^2) Earth Standard Gravitational Parameter

% Load in the desired trajectory:
load('data/u1.mat')
load('data/u2.mat')
load('data/traj.mat')

traj = [X_u1, X_u2];


%% Track via LQR:
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
 
Q = diag([1e3*ones(1,3), 1e-1*ones(1,3)]);
R = 1e9*eye(3);

[K,~,~] = lqr(A,B,Q,R);

% Initial setup:
N = length(X_u1) + length(X_u2);
traj_lqr = zeros(6,N);
traj_lqr(:,1) = X_u1(:,1);
      
sig_r = 0.5;
sig_v = 0.01;

stable_count = 0;
dock = false;
u_lqr = zeros(3,N-1);

% Simulate trajectory:
for ii = 1:N-1
    % Perform estimation:
    x_hat = traj_lqr(:,ii) + [sig_r*randn(3,1); sig_v*randn(3,1)];

    % Calculate control input:
    u_lqr(:,ii) = -K*(x_hat - traj(:,ii+1));
    traj_lqr(:,ii+1) = rk4(@cweq,dt,traj_lqr(:,ii),mu,a,u_lqr(:,ii));
end

% Total fuel consumption:
lqr_sum = sum(sum(abs(u_lqr)));
u1 = reshape(u1,3,[]);
u2 = reshape(u2,3,[]);
u_opt = [u1, u2];
opt_sum = sum(sum(abs(u_opt)));

disp(lqr_sum)
disp(opt_sum)

%% Plot Results:
figure(1)
plot3(traj_lqr(1,:),traj_lqr(2,:),traj_lqr(3,:),'b'); hold on
grid on
axis equal
xlabel('Radial (m)')
ylabel('In-Track (m)')
zlabel('Cross-Track (m)')


% Calculate necessary time vectors:
t_plt = dt:dt:dt*N;
t_opt = [dt_u1:dt_u1:(dt_u1*length(u1)), (dt_u1*length(u1))+(dt_u2:dt_u2:(dt_u2*length(u2)))];
t_lqr = t_plt(1:end-1);

figure(2)
subplot(3,2,1)
    cla
    plot(t_plt,traj(1,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(1,:),'--r');
    xlim([0 inf])
    ylabel('Radial Pos (m)')
    title('Position History')
    legend('Optimal Trajectory','LQR Tracked Trajectory')
subplot(3,2,2)
    cla
    title('Velocity History')
    plot(t_plt,traj(4,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(4,:),'--r');
    xlim([0 inf])
    ylabel('Radial Vel (m/s)')
subplot(3,2,3)
    cla
    plot(t_plt,traj(2,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(2,:),'--r');
    xlim([0 inf])
    ylabel('In-Track Vel (m)')
subplot(3,2,4)
    cla
    plot(t_plt,traj(5,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(5,:),'--r');
    xlim([0 inf])
    ylabel('In-Track Vel (m/s)')
subplot(3,2,5)
    cla
    plot(t_plt,traj(3,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(3,:),'--r');
    xlim([0 inf])
    ylabel('Cross-Track Pos (m)')
    xlabel('Time (sec)')
subplot(3,2,6)
    cla
    plot(t_plt,traj(6,:),'b'); hold on; grid on
    plot(t_plt,traj_lqr(6,:),'--r');
    xlim([0 inf])
    ylabel('Cross-Track Vel (m/s)')
    xlabel('Time (sec)')
    
figure(3)
subplot(3,1,1)
    cla
    plot(t_opt,u_opt(1,:),'b'); hold on; grid on;
    plot(t_lqr,u_lqr(1,:),'--r');
    legend('Forward Optimal','LQR')
    title('Control Input Accelerations')
    ylabel('Accel_x (m/s^2)')
    xlim([0 inf])
subplot(3,1,2)
    cla
    plot(t_opt,u_opt(2,:),'b'); hold on; grid on;
    plot(t_lqr,u_lqr(2,:),'--r');
    ylabel('Accel_y (m/s^2)')
    xlim([0 inf])
subplot(3,1,3)
    cla
    plot(t_opt,u_opt(3,:),'b'); hold on; grid on;
    plot(t_lqr,u_lqr(3,:),'--r');
    ylabel('Accel_z (m/s^2)')
    xlim([0 inf])
    xlabel('Time (sec)')
    
figure(4)
subplot(3,1,1)
    cla
    plot(traj(1,:) - traj_lqr(1,:),'b'); hold on; grid on
    title('Tracking Error')
    ylabel('Error_x (m)')
    xlim([0 inf])
subplot(3,1,2)
    cla
    plot(traj(2,:) - traj_lqr(2,:),'b'); hold on; grid on
    title('Tracking Error')
    ylabel('Error_y (m)')
    xlim([0 inf])
subplot(3,1,3)
    cla
    plot(traj(3,:) - traj_lqr(3,:),'b'); hold on; grid on
    title('Tracking Error')
    ylabel('Error_z (m)')
    xlim([0 inf])
    xlabel('Time (sec)')