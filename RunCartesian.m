% ---------------------------------------------------------------------
%  Melany Mclean
%  July 11, 2016
%
%  Reconstruct p-files to cartesian k-space grid.
%  Make sure orchestra folder is on path.
% ---------------------------------------------------------------------
%
%  HOW TO USE
%  Make sure you are in the directory with all the p-files you want to open.
%
%  Adjust reconstruction settings
% ---------------------------------------------------------------------
clear
numDataChunks = 1;                      % How many matricies would you like to bereak the data into (for memory efficiency)?
settings.SaveDicoms = false;            % Do you want to save DICOM files from reconstructed P-files?
settings.ySqueeze = 0.75;               % What fraction of the x-dimension is covered in y? (eg. 128,96 = 0.75)
numCoils = 32;                          % Enter the number of channels from acquisition coil
virtCoils = 32;                         % Enter the number of 'virtual coils' you want for memory efficiency
numMToffsets = 10;                      % How many MT volumes are thgere?
RefVolLocation = 11;                    % Which volume is the reference volume?
settings.kzTransformedData = true;      % Has the data in the p-file been FFT'd in the kz direction already?
SubSampType = 'STATIC';                 % Choose 'STATIC' or 'DYNAMIC' (only static for now)
PathToOutputDir = '~/Documents/CSProject/CSpreparedData/Data/';
filename = '2018_04_20_qMRIms_013_Kspace_test';

% Display Settings
settings.DisplayFinalImage = true;      % Do you want every full image to pop up after reconstruction?
settings.DisplayFinalSubImage = true;   % Do you want every undersampled image to pop up after reconstruction?
settings.DisplayImageMont = true;       % Do you want a montage of fully sampled z slices for the first P-file?
settings.DisplaySubAxMont = true;       % Do you want a montage of subsampled z slices for the first P-file?
settings.DisplaySubSagMont = true;      % Do you want a montage of subsampled y slices for the first P-file?
settings.DisplaySubKMont = true;        % Do you want a montage of subsampled k-space z slices for the first P-file?
% ---------------------------------------------------------------------

% Create list of P files
% There can be one file for all volumes or each volume separate, both will
% work
files = dir('*.7'); % P files end in .7
pFileNames = {files.name}; % Files in order of time they were made

display('Reading P-files');

% Open files one at a time
for file = 1 : length(pFileNames) % works for volumes saved as multiple p-files
    settings.num = file;
        
    % Call CartesianRecon script with the specified Pfile
    [ fullImage , fullCKspace ] = CartesianRecon(pFileNames(file), settings);
    % fullCKspace organized as ( kx , ky , slice , "phase" , channel , echo )
    
    display('Files Opened: ');
    disp(file/length(pFileNames));
    
    volsPerFile = size(fullCKspace,4);
    
    % Save output variables as 4D and 5D outputs
    output.fullChannelKspace(:,:,:,((file-1)*volsPerFile)+1:file*volsPerFile,:) = fullCKspace; % MAIN DATA ENTRY
    eval(['output.fullImage_' num2str(file) '= fullImage;']);
    
    clear fullImage fullCKspace; % Remove extra variables to avoid confusion
end % file loop

% SAVE DATA

% Specific format for Marc Lebel's code
% [x , y , z , vol , ch]

volsPerChunk = numMToffsets/numDataChunks;

for chunk = 1:numDataChunks;
    
    % Create data matrix
    varname = ['data' num2str(chunk)];
    if RefVolLocation == numMToffsets+1 % Last Volume is reference
        eval([varname '= single(output.fullChannelKspace(:,:,:,((chunk-1)*volsPerChunk)+1:volsPerChunk*chunk,:));']); % volumes 1 to 10 
    else % First volume is reference
        eval([varname '= single(output.fullChannelKspace(:,:,:,((chunk-1)*volsPerChunk)+2:(volsPerChunk*chunk) + 1,:));']);
    end
    
    % Compress number of coils if needed
    if virtCoils ~= numCoils
        eval(['[' varname ',nr] = geometric_coil_compression(' varname ',ceil(virtCoils));'])
    end
    
    % Save data matrix
    if chunk == 1
        save([PathToOutputDir,filename],varname);
    else
        save([PathToOutputDir,filename], varname, '-append');
    end
end

% Save the fully sampled reference Volume
refVol = single(output.fullChannelKspace(:,:,:,RefVolLocation,:));
if virtCoils ~= numCoils
   [refVol,nr] = geometric_coil_compression(refVol,ceil(virtCoils));
end
save([PathToOutputDir,filename], 'refVol', '-append');


