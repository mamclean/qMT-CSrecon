%   Create Field maps from qMT protocol DICOMs
%   
%   Ethan Macdonald, 2016
%   Edited by MM
%
%   Make sure ethanCode folder is on path
%
%% Set image display flags

disp( 'Setting Display Flags' ) ;
doDisplaySourceImageData = false ;
doDisplayMasks = true ;
doDisplayMotionRegistrations = true ;
doDisplayIntermediateImages = true ;
doDisplayQMTinputImages = false ;
doDisplayQMTsourceImages = false ;
doDisplayQMTfinalImagesImages = false ;
doDisplaySementations = true ;
doDisplayShowMeasurements = true ;

%% Set up machine and data options

disp( 'Setting Machine and Data Options' ) ;
Machine = 'MELLIN' ;
% Go into path manager amd change file locations
o = PathManager( Machine ) ;
o.subject_directories = { '2018_04_20_qMRIms_013/DICOMS' } ;

o.ME_SPGR_subdir = '/AX_ME_SPGRE_MPRI' ; %B0
o.B1_mapping_subdir = '/MULT_TR_B1' ; %B1
o.DESPOT1_subdir = '/DESPOT1' ; %T1
%o.qMT_subdir = '/qMRIms_010_qMT' ; %qMT

PathToOutput = '~/Documents/CSProject/outputData';
SubjectID = '2018_04_03_qMRIms_012';

% Create Storage Location
cd (PathToOutput);
mkdir(SubjectID);

totalTimeTimer = tic ; 

%% Open Source Data

cellTimer = tic ;

disp( 'Opening source data' ) ; % Dicom data only
objSourceData = SourceData( doDisplaySourceImageData ) ;

objSourceData.zLength = 60 ; % Number of slices
objSourceData.METE = ...
    [ 1.396 3.800 5.664 7.528 9.392 11.256 13.120 14.984 ] ;
objSourceData.qMT_data_table = ...  % Enter MT offsets [TR alpha alphaMT offset(Hz)]
    ...                             % Not used, could remove this
    [ ...
        25 7 142 443
        25 7 142 1088
        25 7 142 2732
        25 7 142 6862
        25 7 142 1735
        25 7 400 443
        25 7 400 1088
        25 7 400 2732
        25 7 400 6862
        25 7 400 1735
    ] ; 
objSourceData.numEchoesSPGR = length( objSourceData.METE ) ;

objSourceData = fOpenMechoSPGR( objSourceData, o, doDisplaySourceImageData ) ; % dicomread
objSourceData = fOpenB1mappingSourceImages( objSourceData, o, ...
    doDisplaySourceImageData ) ;
objSourceData = fOpenDESPOT1images( objSourceData, o, ...
    doDisplaySourceImageData ) ;
objSourceData.zLength = 60 ;
objSourceData.numTimePts = 25 ;
% objSourceData = fOpenQMT( objSourceData, o, doDisplayQMTsourceImages ) ;

totalRunningTime = toc( totalTimeTimer ) ; 
cellRunTime = toc( cellTimer ) ;
disp( [ 'Time to open source data ' num2str( cellRunTime ) ] ) ; 
disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% Masks and Registrations

cellTimer = tic ;

disp( 'Calculating Masks and Registrations' ) ;
PercentThreshold = 0.07 ; % Why is this 7% ?
objMask_MESPGR = Mask( abs( objSourceData.volMESPGR(:,:,:,1) ), ...
    PercentThreshold, doDisplayMasks ) ; % I can't find 'Mask'

PercentThreshold = 0.055 ; % Higher #, less tissue. This is the qMT mask
objMask_DESPOT1 = Mask( objSourceData.volSPGR1, ...
    PercentThreshold, doDisplayMasks ) ;

% objMask.Seed                = [ 128 128 30 ] ;
% objMask.numAngles           = [180*1.2 360*1.2]*3 ;
% objMask.innerSphereSize     = 20 ;
% objMask.RayThreshold        = 0.3 ;
% objMask.maxRayLength        = 256 ;
% objMask.rayResolution       = 0.5 ;
% objMask.erosionLevel        = 6 ;
% objMask.rayFilterKernel     = gausswin(7) * gausswin(7)' ;
% objMask = BrainExtractorMEX( objMask, 1 ) ;
% objMask = FillErrorsInMask( objMask, 1 ) ;
% timeToComputeMasksAndRegistrations = toc ;
% disp( [ 'Time to compute masks and registrations ' ...
%     num2str( timeToComputeMasksAndRegistrations ) ] ) ; 

% Save the mask for use in qMTLab
mask = objMask_DESPOT1.volMask; % Must be named 'Mask' for qmt input
save([PathToOutput, '/', SubjectID, '/mask.mat'], 'mask');


cellRunTime = toc( cellTimer ) ;
totalRunningTime = toc( totalTimeTimer ) ; 
disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% B0 map

cellTimer = tic ;

disp( 'Calculating B0 Map' ) ;
objB0 = B0field( objSourceData ) ; % go into B0field folder/ file
objB0.seedPt = [ 128 128 30 ] ;
objB0 = ApplyMask( objB0, objMask_MESPGR.volMask, 0 ) ; % ApplyMask is function in B0field file - ApplyMask(obj, mask, display)
% objB0 = UnwrapPhase3D( objB0, doDisplayIntermediateImages ) ; % more functions
% objB0 = CalcMagFeildUnipolar( objB0 ) ;

objB0 = ComplexFittingUnipolar( objB0, doDisplayIntermediateImages ) ;
objB0 = displayMagField( objB0 ) ;

% Convert to Hz
B0map = objB0.volB0; % in mT
B0map = B0map.*42.58; % in Hz
% Save the B0map for use in qMTLab
save([PathToOutput, '/', SubjectID, '/B0map.mat'], 'B0map');

time2runCell = toc( cellTimer ) ;
time2CalculateB0maps = toc( totalTimeTimer ) ; 
disp( [ 'Time to calculate B0 maps ' num2str( time2CalculateB0maps ) ] ) ;
disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% B1 map

cellTimer = tic ;

disp( 'Calculating B1 Map' ) ;
objB1map = B1fieldMap( objSourceData );
objB1map = ReconB1Field( objB1map, doDisplayIntermediateImages );

% Filter the map - OPTIONAL
sigma = 2.1233; % Conversion from FWHM = 5 voxels
filter = fspecial('gaussian',ceil(3*sigma),sigma);
for i = 1:60
    objB1map.volB1(:,:,i) = roifilt2(filter,objB1map.volB1(:,:,i),objMask_DESPOT1.volMask(:,:,i));
end

% Save the B1map for use in qMTLab
B1map = objB1map.volB1; % In percent
save([PathToOutput, '/', SubjectID, '/B1map.mat'], 'B1map');

time2runCell = toc( cellTimer ) ;
totalRunningTime = toc( totalTimeTimer ) ; 
disp( [ 'Time to calculate B1 maps ' num2str( time2runCell ) ] ) ;
disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% T1 map

% % Hack - Subject 10 only
% objMask_DESPOT1.volMask = ones(256,256,60);
% % End hack

cellTimer = tic ;

disp( 'Calculating T1 Map' ) ;

% Hack to swap the volumes
vol1 = objSourceData.volSPGR2;
vol2 = objSourceData.volSPGR1;
objSourceData.volSPGR1 = vol1;
objSourceData.volSPGR2 = vol2;
% end of hack

objDespot1 = Despot1( objSourceData, objB1map ) ;
objDespot1 = CalcDespot1( objDespot1, doDisplayIntermediateImages ) ;
objDespot1.volT1 = objDespot1.volT1 .* objMask_DESPOT1.volMask ;
FourDviewer( objDespot1.volT1 ) ;

% Save the R1map for use in qMTLab
T1map = objDespot1.volT1;
T1map = T1map./ 1000; % Convert from ms to s
T1map(T1map ~= T1map) = 0; % Replace NaN with zero
T1map(T1map < 0.001) = 1; % Remove zeros

R1map = 1./T1map;

save([PathToOutput, '/', SubjectID, '/T1map.mat'], 'T1map');
save([PathToOutput, '/', SubjectID, '/R1map.mat'], 'R1map');


time2runCell = toc( cellTimer ) ;
totalRunningTime = toc( totalTimeTimer ) ; 
disp( [ 'Time to calculate T1 maps ' num2str( time2runCell ) ] ) ;
disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% T2* maps
% cellTimer = tic ;
% 
% disp( 'Calculating T2* Map' ) ;
% objT2star = T2star( objSourceData, objMask_MESPGR ) ; 
% objT2star = fitT2star( objT2star, doDisplayIntermediateImages ) ; 
% FourDviewer( objT2star.volT2star ) ; 
% 
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to calculate T2* maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% qMT

% cellTimer = tic ;
% % I'm not sure what this does, but it doesn't calculate qMT maps.
% disp( 'Calculating qMT Maps' ) ;
% objQMT = MagnetizationTransfer( objSourceData ) ; 
% objQMT = MagnetizationTransfer( objQMT, doDisplayQMTfinalImagesImages ) ; 
% 
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to calculate T2* maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% Raw Image Data

% cellTimer = tic ;
% 
% disp( 'Reconstruction of Raw Data' ) ;
% % ADD CODE
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to reconstruct raw maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% Undersampling Raw Data

% cellTimer = tic ;
% 
% disp( 'Undersampling Raw Data' ) ;
% % ADD CODE
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to reconstruct raw maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% Compressed Sensing Reconstructions

% cellTimer = tic ;
% 
% disp( 'Compressed Sensing Reconstructions' ) ;
% % ADD CODE
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to reconstruct raw maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% Segmentations and Measurements

% cellTimer = tic ;
% 
% disp( 'Performing Sementations and Measurements' ) ; 
% 
% time2runCell = toc( cellTimer ) ;

%% Write out Images to Dicom

% cellTimer = tic ;
% 
% disp( 'Writing Dicoms Files for Computed Images from Experiment' ) ; 
% % ADD CODE
% time2runCell = toc( cellTimer ) ;
% totalRunningTime = toc( totalTimeTimer ) ; 
% disp( [ 'Time to reconstruct raw maps ' num2str( time2runCell ) ] ) ;
% disp( [ 'Total Running time ' num2str( totalRunningTime ) ] ) ;

%% total time to perform processing

