% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022
% clc; clear;
% close all


%% import nominal data
% load('novDataFull3.mat')
% load('mayDataFull_00B.mat')
% load('mayDataFull_10B.mat')

lla = ecef2lla([-2700404.467 -4292605.260  3855137.600]); % main antenna middle of durand
rng(6)

%% Extended Kalman Filter
%% simulation parameters
dt = 1;  T = length(dateTvec);  t = 0:dt:T;
nx = 3;   % three state variables
mx = 4;

% November Data Full3
Q = diag([6.6290e-12, 5.0559e-12, 0.3448]);  % process covariance
R = diag([0.5957, 0.1188, 3.3734]);  % measurement covariance
% November Data Short1
Q = 0.1*diag([5.14E-12, 4.18E-12, 0.2979]);  % process covariance
R = 0.1*diag([0.466, 0.1188, 3.6448]);  % measurement covariance

% November Data 3 hours
Q = 0.1*diag([5.14E-12, 4.18E-12, 0.2979]);  % process covariance
R = 0.1*diag([0.466, 1.4184, 0.1188, 4.4435]);  % measurement covariance



%% initialize all the variables
x = zeros(nx,length(t)-1);    % true states
y = zeros(mx,length(t)-1);     % measurement
x_ekf = zeros(nx, length(t)-1);
sig_ekf = zeros(nx,nx,length(t)-1);    % EKF
m = zeros(mx-nx,length(t)-1);

lla = [37.4269, -122.1733, lla(3)];
% initialize belief
x(:,1) = 1000*lla;
x_ekf(:,1) = x(:,1); % known initial state
sig_ekf(:,:,1) = diag([1E-3 1E-3 1]); % 1E-3*eye(3); % strong belief

%% functions
f = @(x)[x];   % system model
g = @(x,m)[x;m];    % measurement model
Ja = @(x)[1 0 0;
    0 1 0;
    0 0 1];               % jacobian for EKF
Jc = @(x)[1 0 0;
    0 1 0;
    0 0 1;
    0 0 0];   % jacobian for EKF
%% simulation actual state
% for i = 2:length(t)-1
%     x(:,i) =  1000*[latdata(i); londata(i); heightData(i)];  % actual state
%     y(:,i) = g(x(:,i),m(:,i)) + mvnrnd(zeros(1,mx),R)';     % actual measurement
% end

y = [M_L1; M_L2; agcCntdata(1,:); agcCntdata(2,:)];

%% EKF
tic;
for i = 2:length(t)-1
    % predict step
    A = Ja( x_ekf(:,i-1) );
    x_ekf(:,i) = f( x_ekf(:,i-1) );
    sig_ekf(:,:,i) = A * sig_ekf(:,:,i - 1) * A' + Q;
    % update step
    mrsInno(:,i) = y(:,i) - g( x_ekf(:,i), m(:,i) );      %measurement innovation
    C = Jc( x_ekf(:,i) );   %measurement matrix
    K = sig_ekf(:,:,i) * C' * inv(C * sig_ekf(:,:,i) * C' + R);   %Kalman Gain
    x_ekf(:,i) = x_ekf(:,i)  + K * mrsInno(:,i);
    sig_ekf(:,:,i) = (eye(3) - K * C) * sig_ekf(:,:,i);
    if mod(i,10000) == 0
        disp(i)
    end
end
t_EKF = toc;
x_ekf = x_ekf/1000;







