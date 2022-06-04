% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022
% clc; clear;
% close all

%% import data
% load('novDataFull3.mat')
% shortened for nov, 15 hours to 30 hours
trange = 1:3*3600;

% load('mayDataFull_10B.mat')
% % shortened for may to 19 hours
% trange = 1:19*3600;



lla = ecef2lla([-2700404.467 -4292605.260  3855137.600]); % main antenna middle of durand

% [M_L1,~] = max(cn0rawx(:,:,1)); % M is the max value in each column, I is its index in each column
% [M_L2,~] = max(cn0rawx(:,:,2)); % M is the max value in each column, I is its index in each column
% Note: after observation M_L1 is basically the same as M, which is
% expected since the NAV_SAT cn0 data is simply the L1 part

agcCntdata = (100/8191).*agcCntdata;

%{
%% compute R based on measurement noise
% measurements are:
% M - max CN0 over the time series
% AGCcnt - separate to only include L1
% prRes - pseudorange residuals

figure
subplot(1,3,1)
% scatter(dateTvec,M); hold on
plot(dateTvec,M_L1,'-')
% scatter(dateTvec,M_L2,'.')
title('CN_{0_{max}}')
xlabel('Datetime'); ylabel('CN_0[dB Hz]')

subplot(1,3,2)
plot(dateTvec, agcCntdata)
title('AGC (%)')
xlabel('Datetime'); ylabel('AGC [%]')

subplot(1,3,3)
plot(dateTvec, prResdata)
title('Pseudorange Residuals')
xlabel('Datetime'); ylabel('idk')

MaxCN0_L1_mean = mean(M_L1);
MaxCN0_L1_var = var(M_L1);

MaxCN0_L2_mean = mean(M_L2);
MaxCN0_L2_var = var(M_L2);

AGC_L1_mean = mean(agcCntdata(1,:));
AGC_L1_var = var(agcCntdata(1,:));

AGC_L2_mean = mean(agcCntdata(2,:));
AGC_L2_var = var(agcCntdata(2,:));

prResdata_mean = mean(prResdata(~isnan(prResdata)),1);
prResdata_var = var(prResdata(~isnan(prResdata)),0,1);

%%
MaxCN0_L1_mean
MaxCN0_L1_var
MaxCN0_L2_mean
MaxCN0_L2_var
AGC_L1_mean
AGC_L1_var
AGC_L2_mean
AGC_L2_var
prResdata_mean
prResdata_var

%% compute Q based on process noise

figure
scatter(londata,latdata,'.')
xlabel('Longitude'); ylabel('Latitude');

lon_var = var(londata);
lat_var = var(latdata);
height_var = var(heightData);

%%
lon_var
lat_var
height_var

%}



MaxCN0_L1_mean = mean(M_L1(trange));
MaxCN0_L1_var = var(M_L1(trange));

MaxCN0_L2_mean = mean(M_L2(trange));
MaxCN0_L2_var = var(M_L2(trange));

AGC_L1_mean = mean(agcCntdata(1,trange));
AGC_L1_var = var(agcCntdata(1,trange));

AGC_L2_mean = mean(agcCntdata(2,trange));
AGC_L2_var = var(agcCntdata(2,trange));

% prResdatatemp = prResdata(:,trange);
% prResdata_mean = mean(prResdatatemp(~isnan(prResdatatemp)),1);
% prResdata_var = var(prResdatatemp(~isnan(prResdatatemp)),0,1);
% prResdata(isnan(prResdata)) = 0;


%% compute Q based on process noise
%
% figure
% scatter(londata,latdata,'.')
% xlabel('Longitude'); ylabel('Latitude');

lon_var = var(londata(trange));
lat_var = var(latdata(trange));
height_var = var(heightData(trange));

if displayHists
% MaxCN0_L1_mean = mean(M_L1(trange));
% MaxCN0_L1_var = var(M_L1(trange));
figure
histfit(M_L1(trange));
title('MaxCN0_{L1}')

figure
histfit(M_L2(trange));
title('MaxCN0_{L2}')

figure
histfit(agcCntdata(1,trange));
title('AGC_{L1}')

figure
histfit(agcCntdata(2,trange));
title('AGC_{L2}')
end

%%
if displayStats
    MaxCN0_L1_mean
    MaxCN0_L1_var
    MaxCN0_L2_mean
    MaxCN0_L2_var
    AGC_L1_mean
    AGC_L1_var
    AGC_L2_mean
    AGC_L2_var
%     prResdata_mean
%     prResdata_var

    lon_var
    lat_var
    height_var
end









