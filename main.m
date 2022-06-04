% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022
% main script to run the following scripts:
%   processNominal.m - to find noise values for 
%       nominal conditions (Q, R matrices)
%   runEKF.m - runs EKF on position data using receiver measurements,

clc; clear;
close all

% These import the following .mat file:
%    novDataFull3.mat XXXXX
% 
% load('novData3hours.mat')
load('nomANDjam.mat')

%% run nominal variance analysis
displayStats = false;
displayHists = false;
processNominal
%% run the EKF over the data series
% runEKF % most vanilla version 3x1 measurement vector
% runEKF2 % more complex, 8x1 measurement vector
runEKF3 % same as runEKF2 but without prRes, 7x1 measurement vector
% runEKF4 % like runEKF, 4x1 meas vector, (2 CN0max, 2 AGC)
runEKF5 % same as runEKF3 but with GLRT

% meanlat = mean(latdata(3600:3600*2));
% meanlon = mean(londata(3600:3600*2));

%% plotting details
plot_lonlat = false; % plot all lon/lat data
plot_lonlat_ss = false; % plot lon/lat data after steady state estimates
plot_lonlat_geo = true; % plot lon/lat on google maps
plot_states = false; % plot x, y, height over time
    include_covariance = false; % to include cov on state plots
    trange = 1:3600;
plot_states4 = false; % plot x, y, height over time for meas = 8
    include_covariance = false; % to include cov on state plots
    trange = 1:3600;
plot_innov = true; % plot innovation vector values over time

plot_it_all

%% Analyze the innovation vector
analyzeInnovation






 



% 


disp('Done')