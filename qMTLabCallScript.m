% Sensitivity Analysis for real data, variable R1obs map
clear; clc;

% From participant folder
Protocol = load('~/Documents/CSProject/qMTLab-master/SPGR/Parameters/Protocol_QMRIMS_Jan4_2018_gbw200_Perscribed_gauss.mat');
FitOpt = load('~/Documents/CSProject/qMTLab-master/SPGR/Parameters/Fit_QMRIMS_KrFree.mat');
load('rMaskNoSkull256.mat');
load('rB0map256.mat');
load('rB1map256.mat');
load('rR1map256.mat');  
data.Mask = double(Mask);
data.B0map = double(B0map);
data.B1map = double(B1map);
data.R1map = double(R1map);

% Test 1
load('MTdata_4_V1.mat');
data.MTdata = double(MTdata);
QMRIMS_011_4_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_4_V1','QMRIMS_011_4_V1');

% Test 2
load('MTdata_4_V2.mat');
data.MTdata = double(MTdata);
QMRIMS_011_4_V2 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_4_V2','QMRIMS_011_4_V2');

% Test 3
load('MTdata_4_V3.mat');
data.MTdata = double(MTdata);
QMRIMS_011_4_V3 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_4_V3','QMRIMS_011_4_V3');

% Test 4
load('MTdata_4_V4.mat');
data.MTdata = double(MTdata);
QMRIMS_011_4_V4 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_4_V4','QMRIMS_011_4_V4');

% Test 1-2
load('MTdata_8_V1.mat');
data.MTdata = double(MTdata);
QMRIMS_011_8_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_8_V1','QMRIMS_011_8_V1');

% Test 2-2
load('MTdata_8_V2.mat');
data.MTdata = double(MTdata);
QMRIMS_011_8_V2 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_8_V2','QMRIMS_011_8_V2');

% Test 3-2
load('MTdata_8_V3.mat');
data.MTdata = double(MTdata);
QMRIMS_011_8_V3 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_8_V3','QMRIMS_011_8_V3');

% Test 4-2
load('MTdata_8_V4.mat');
data.MTdata = double(MTdata);
QMRIMS_011_8_V4 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_8_V4','QMRIMS_011_8_V4');

% Test 1-3
load('MTdata_12_V1.mat');
data.MTdata = double(MTdata);
QMRIMS_011_12_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
save('FitResults/QMRIMS_011_12_V1','QMRIMS_011_12_V1');




% % Next Fit Option
% FitOpt = load('~/Documents/CSProject/qMTLab-master/SPGR/Parameters/Fit_QMRIMS_KrFree.mat');
% 
% % Test 1
% load('MTdata_24_V1.mat');
% data.MTdata = double(MTdata);
% QMRIMS_011_24_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% save('FitResults/QMRIMS_011_24_V1','QMRIMS_011_24_V1');
% 
% % Test 2
% load('MTdata_28_V1.mat');
% data.MTdata = double(MTdata);
% QMRIMS_011_28_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% save('FitResults/QMRIMS_011_28_V1','QMRIMS_011_28_V1');
% 
% % Test 3
% load('MTdata_32_V1.mat');
% data.MTdata = double(MTdata);
% QMRIMS_011_32_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% save('FitResults/QMRIMS_011_32_V1','QMRIMS_011_32_V1');
% 
% % Test 4
% load('MTdata_36_V1.mat');
% data.MTdata = double(MTdata);
% QMRIMS_011_36_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% save('FitResults/QMRIMS_011_36_V1','QMRIMS_011_36_V1');
% 
% % Test 5
% load('MTdata_40_V1.mat');
% data.MTdata = double(MTdata);
% QMRIMS_011_40_V1 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% save('FitResults/QMRIMS_011_40_V1','QMRIMS_011_40_V1');

% % 
% % 
% % data.MTdata = double(MTdata);
% % data.Mask = double(Mask);
% % data.B0map = double(B0map);
% % data.B1map = double(B1map);
% % 
% % % Move to function location
% % cd ~/Documents/CSProject/qMTLab-master/Common
% % 
% % % Call each R1map and fit data
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map1.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map1_krFree = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% % 
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map2.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map2_krFree = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% % 
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map3.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map3_krFree = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% % 
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map4.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map4_krFree = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% % 
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map5.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map5 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
% % 
% % load('~/Documents/CSProject/outputData/2017_10_03_qMRIms_009/Slice30/rR1map6.mat');
% % data.R1map = double(R1map);
% % QMRIMS_009_0_V1_R1map6 = FitData( data, Protocol, FitOpt, 'SPGR', 1);
