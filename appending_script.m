% Nicolas San Miguel - May 2022
% AA 273 - Spring 2022

clc; clear;
close all



load('jam2append.mat')
epochs2 = epochs;
dateTvec2 = dateTvec;
M2 = M;
I2 = I;
cn0data2 = cn0data;
svIddata2 = svIddata;
prResdata2 = prResdata;
azimdata2 = azimdata;
elevdata2 = elevdata;
agcCntdata2 = agcCntdata;
ofsIdata2 = ofsIdata;
magIdata2 = magIdata;
ofsQdata2 = ofsQdata;
magQdata2 = magQdata;
cn0rawx2 = cn0rawx;
svIdrawx2 = svIdrawx;
prStdrawx2 = prStdrawx;
londata2 = londata;
latdata2 = latdata;
heightData2 = heightData;
hAccData2 = hAccData;
pDOPdata2 = pDOPdata;
vAccData2 = vAccData;
numSVdata2 = numSVdata;

load('novData3hours.mat')
epochs3 = epochs;
dateTvec3 = dateTvec;
M3 = M;
I3 = I;
cn0data3 = cn0data;
svIddata3 = svIddata;
prResdata3 = prResdata;
azimdata3 = azimdata;
elevdata3 = elevdata;
agcCntdata3 = agcCntdata;
ofsIdata3 = ofsIdata;
magIdata3 = magIdata;
ofsQdata3 = ofsQdata;
magQdata3 = magQdata;
cn0rawx3 = cn0rawx;
svIdrawx3 = svIdrawx;
prStdrawx3 = prStdrawx;
londata3 = londata;
latdata3 = latdata;
heightData3 = heightData;
hAccData3 = hAccData;
pDOPdata3 = pDOPdata;
vAccData3 = vAccData;
numSVdata3 = numSVdata;

% concatenate
epochs = [epochs3 epochs2];
dateTvec = [dateTvec3 dateTvec2];

epp = 1:length(dateTvec2);
% newdatetime = dateTvec3(end)*ones(1,length(dateTvec2)) + epp;

newstart = dateTvec3+1/(24*3600);
newend = dateTvec3+1/(24*3600);
for i = 1:length(dateTvec2)
    dateTvec3 = [dateTvec3 dateTvec3(end)+1/(24*3600)];
end

dateTvec = dateTvec3;
% M = M;
% I = I;
% cn0data = [cn0data3 cn0data2];

% svIddata = svIddata;
% prResdata = prResdata;
% azimdata = azimdata;
% elevdata = elevdata;
agcCntdata = [agcCntdata3 agcCntdata2];
% ofsIdata = ofsIdata;
% magIdata = magIdata;
% ofsQdata = ofsQdata;
% magQdata = magQdata;
% cn0rawx = cn0rawx;
[M_L1_3,~] = max(cn0rawx3(:,:,1)); % M is the max value in each column, I is its index in each column
[M_L2_3,~] = max(cn0rawx3(:,:,2)); % M is the max value in each column, I is its index in each column
[M_L1_2,~] = max(cn0rawx2(:,:,1)); % M is the max value in each column, I is its index in each column
[M_L2_2,~] = max(cn0rawx2(:,:,2)); % M is the max value in each column, I is its index in each column
M_L1 = [M_L1_3 M_L1_2];
M_L2 = [M_L2_3 M_L2_2];



% svIdrawx = svIdrawx;
% prStdrawx = prStdrawx;
londata = [londata3 londata2];
latdata = [latdata3 latdata2];
heightData = [heightData3 heightData2];
% hAccData = hAccData;
% pDOPdata = pDOPdata;
% vAccData = vAccData;
% numSVdata = numSVdata;
%% savemat
save('nomANDjam','epochs','dateTvec',...
    'agcCntdata','M_L1','M_L2',...
    'cn0rawx',...
    'londata','latdata','heightData')




