% Melany Mclean Sept 21, 2017
% Script to coregister qMT maps and images and to segment GM WM
% in SPM

%% WM GM Segmentation - RUN THIS SECTION ONLY

% This opens a structural DICOM image, segments it and saves as nifti files
% This step FAILS, just continue with script
% REMEMBER to change subject ID
% In the following line enter:
%       path to T1-weighted anatomical DICOM folder, path to output SPM folder

preprocess_structural('/home/melanymclean/Documents/CSProject/sourceData/2018_04_20_qMRIms_013/DICOMS/T1weighted_Anatomical','/home/melanymclean/Documents/CSProject/sourceData/2018_04_20_qMRIms_013/SPMDataTest/')

%% Get a DICOM image and save it as a nifti to get the headers

% IN TERMINAL
% go to the qMT Dicom directory
% Type:
% dcm2nii -o ../../SPMdata IM-0008-0001.dcm - pick any image from qMT folder
% This saves it as a super long name in SPMdata folder - copy the name for next line

%% Save first qMT volume, Mask, R1, B0, and B1 maps as nifti

% Load the qMT volume and maps
% Don't forget to change subject folder and nii filename
qmt_nii = load_nii('~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/20180720_160830PIKELABSMAMqMTs023a1001.nii');
load('~/Documents/CSProject/outputData/2018_07_20_qMRIms_018/B0map.mat', 'B0map');
load('~/Documents/CSProject/outputData/2018_07_20_qMRIms_018/B1map.mat', 'B1map');
load('~/Documents/CSProject/outputData/2018_07_20_qMRIms_018/R1map.mat', 'R1map');
load('~/Documents/CSProject/outputData/2018_07_20_qMRIms_018/mask.mat', 'mask');

% Modify matlab files to fit nifti format
B1map = flipud(rot90(B1map,3));
R1map = flipud(rot90(R1map,3));
B0map = flipud(rot90(B0map,3));
Mask = flipud(rot90(mask,3));

% output filenames
% don't forget to change subject ID
fname_qMT = '~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/MTdata_0_V1_vol1.nii';
fname_mask = '~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/Mask.nii';
fname_B0map = '~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/B0map.nii';
fname_B1map = '~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/B1map.nii';
fname_R1map = '~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata2/R1map.nii';

% Make a nifti file for qMT data
qmt_1 = qmt_nii;
% Save the first volume only
qmt_1.img = qmt_nii.img(:,:,:,1);
qmt_1.hdr.dime.dim = [3,256,256,60,1,1,1,1];
qmt_1.hdr.dime.bitpix = 32;
qmt_1.hdr.dime.datatype = 16;

% Make nifti files for maps
Mask_nii = qmt_1;
B0map_nii = qmt_1;
B1map_nii = qmt_1;
R1map_nii = qmt_1;

% Save matricies in image space
Mask_nii.img = Mask;
R1map_nii.img = R1map; %R1map still fine
B1map_nii.img = B1map;
B0map_nii.img = B0map;

% Save the files
save_nii(qmt_1,fname_qMT);
save_nii(Mask_nii,fname_mask);
save_nii(R1map_nii,fname_R1map);
save_nii(B1map_nii,fname_B1map);
save_nii(B0map_nii,fname_B0map);

%% Coregister all images

% Don't forget to change subject folder
cd ~/Documents/CSProject/sourceData/2018_07_20_qMRIms_018/SPMdata3

% Coregister structural to qMT vol # 1 and bring all other maps
coreg_spm('MTdata_0_V1_vol1','structural.nii',{'c1structural.nii','c2structural.nii','R1map.nii','B1map.nii','B0map.nii','Mask.nii'})
% Registered images saved as 'rFilename.nii'
% This should probably be set up to register each of B0, B1, and R1 (with
% c1 and c2 brought along) separately

%% Convert back to matlab matricies for qMTLab to use

% Get the image data
Mask = load_nii('rrMask.nii');
Mask = Mask.img;
B0map = load_nii('rrB0map.nii');
B0map = B0map.img;
B1map = load_nii('rrB1map.nii');
B1map = B1map.img;
R1map = load_nii('rrR1map.nii');
R1map = R1map.img;
WMmask = load_nii('rc2structural.nii');
WMmask = WMmask.img;
GMmask = load_nii('rc1structural.nii');
GMmask = GMmask.img;

% Put maps in the same ref frame as MTdata
B0map = fliplr(rot90(B0map));
B1map = fliplr(rot90(B1map));
R1map = fliplr(rot90(R1map));
Mask = fliplr(rot90(Mask));
WMmask = fliplr(rot90(WMmask));
GMmask = fliplr(rot90(GMmask));

% Save as correct variables
cd ~/Documents/CSProject/outputData/2018_07_20_qMRIms_018

save('rB0map','B0map')
save('rB1map','B1map')
save('rR1map','R1map')
save('rMask','Mask')
save('rWMmask','WMmask')
save('rGMmask','GMmask')

% Save single slice versions
cd Slice28

Mask = Mask(:,:,28);
B0map = B0map(:,:,28);
B1map = B1map(:,:,28);
R1map = R1map(:,:,28);
WMmask = WMmask(:,:,28);
GMmask = GMmask(:,:,28);

save('rB0map256','B0map')
save('rB1map256','B1map')
save('rR1map256','R1map')
save('rMask256_skull','Mask')
save('rWMmask','WMmask')
save('rGMmask','GMmask')


