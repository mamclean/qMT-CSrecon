
%   Melany Mclean, 2016
%
%   This script is designed to load k-space data and sampling masks and
%   perform retrospective undersampling and call sparseSENSE code from
%   Marc Lebel. It also stores data in an appropreate format for qMTLab.
%
%   Ensure sparseSENSE code is on path
%   -----------------------------------------------------------------------
clear
SubjectID = '2018_04_20_qMRIms_013';
SampleRate = '4';
VersionNumber = '10';
NumMtOffsets = 10;          % Enter the TOTAL number of mt-offset volumes (not including baseline)
NumProcessingChunks = 1;    % How many processing units do you want the recon broken into?
FillSizeX = 128;            % Enter size of zero filling in kx
FillSizeY = 96;             % Enter size of zero filling in ky
Slices = 64;                % Enter number of slices
ReconWithFullRef = false;   % Do you want the fully-sampled reference volume included in CS recon?
UndersampRef = true;        % Do you want the reference volume to be included and undersampled?
SliceNum = 30;              % Enter the slice you want stored for 2D qMT (after extra slices removed)
                            % NOTE THAT IF THE REFERENCE VOLUME IS NOT INCLUDED IN THE RECON, M NEEDS TO BE SCALED SEPARATELY
%   -----------------------------------------------------------------------

% Name of k-space data file opened with orchestra
datafilename = [SubjectID,'_Kspace.mat'];

%   Set dirictory options
%   Paths to prepared data, masks, and hamming filters
PathToPreparedData = '~/Documents/CSProject/CSpreparedData/Data/';
PathToMask = ['~/Documents/CSProject/CSpreparedData/Masks/',num2str(FillSizeX),'.',num2str(FillSizeY),'.',num2str(Slices),'.',num2str(NumMtOffsets),'/'];
PathToHammingFilter = '~/Documents/CSProject/CSpreparedData/HammingFilters/';
filename = ['CSinputData_',SampleRate,'_V',VersionNumber,'_P',num2str(NumProcessingChunks),'.mat']; % Mask filename
%   Path to reconstructed (output) data
OutputDir = '~/Documents/CSProject/outputData/'; % MT scaled
OutputDirUnscaled = '~/Desktop/Storage/Documents/CSProject/outputData/'; % No MT-scaling. Different location so it doesn't take too much space
%   Set output filename with baseline description
if ~ ReconWithFullRef && UndersampRef
    outputFilename = ['/MTdata_',SampleRate,'_V',VersionNumber];
elseif ReconWithFullRef && ~ UndersampRef
    outputFilename = ['/MTdata_',SampleRate,'_V',VersionNumber,'_FullRef'];
elseif ~ ReconWithFullRef && ~ UndersampRef
    outputFilename = ['/MTdata_',SampleRate,'_V',VersionNumber,'_NonCSref'];
end
%   Set filename of unscaled data
outputFilenameUnscaledData = ['MTunscaled_',SampleRate,'_V',VersionNumber,'.mat'];

% Reconstruct one subset of the data at a time (optional)
for dataChunk = 1:NumProcessingChunks
    
    %   LOAD DATA
  
    % Set data names accordingly
    dataname = ['data',num2str(dataChunk)];
    maskname = ['mask',num2str(dataChunk)];
    numOffsets = NumMtOffsets/NumProcessingChunks;
    % Load files 
    data = load([PathToPreparedData,datafilename],dataname);
    load([PathToPreparedData,datafilename],'refVol'); % refVol is the baseline image
    eval(['data = data.' dataname ';']);
    
    if ReconWithFullRef || UndersampRef
        % Add the ref volume to the CS matrix
        data = cat(4,data,refVol); % Add the ref volume as the last volume in data stack
        clear refVol
    end

    %   Get size
    %   [freq-enc,phase-enc,slice,volume,channel]
    [np,nv,ns,nt,nr] = size(data);

    %   Convert to full 3D k-space
    data = fFastFT(data,3,0);

    if ns ~= 1 % More than 1 slice
        %  Do 3D hamming filter (15 voxels around edges)
        load([PathToHammingFilter,'3DFilter_128*96*64_15.mat'],'ham3');

        for vol = 1:nt
            for ch = 1:nr
                data(:,:,:,vol,ch) = data(:,:,:,vol,ch).*ham3;
            end
        end
        
    else % Only 1 slice
        %  Do 2D hamming filter (15 voxels around edges)
        load([PathToHammingFilter,'2DFilter_128*88_15.mat'],'ham2');

        for vol = 1:nt
            for ch = 1:nr
                data(:,:,1,vol,ch) = data(:,:,1,vol,ch).*ham2;
            end
        end
    end % Done applying filters

    %  Zero fill k-space - THIS MAY CAUSE MEMORY PROBLEMS
    %  Optional - nothing happens if you set fill sizes equal to existing
    %  k-space size
    if FillSizeX ~= np || FillSizeY ~= nv
        zfKspace = single(zeros(FillSizeX,FillSizeY,ns,nt,nr));
        zfKspace((FillSizeX-np)/2+1:np+(FillSizeX-np)/2,(FillSizeY-nv)/2+1:nv+(FillSizeY-nv)/2,:,:,:) = data; %zfKspace(65:192,81:176,:,:,:) = data;
        data = zfKspace; % Gets to 92% mem
        clear zfKspace
    end

    %   MAIN CS RECON
    
    %  Get new size
    [np,nv,ns,nt,nr] = size(data);

    % %   OPTIONAL: perform geometric coil compression to conver large number of real
    % %   coils to fewer virtual coils (12 to 4 virtual coils)
    % [data,nr] = geometric_coil_compression(data,ceil(nr/4));
    % mask = mask(:,:,:,:,1:nr);

    %   Perform image-based fftshift to centre image
    %   (do this in k-space to save forward then inverse transform)
    data = fftshiftF(data,[1 2]); % Throws out of memory error (15 vols, 64 slice)

    %   Estimate coil sensitivity maps from first phase
    S = repmat(FermiFilt([np nv ns],[64 64 64],1/20),[1 1 1 1 nr]) .* data(:,:,:,1,:);
    S = iFastFT(S,[1 2 3],1);
    S = S./repmat(sqrt(sum(abs(S).^2,5)),[1 1 1 1 nr]);
 
    %   Convert k-space to hybrid image (readout) and k-space (phase encodes)
    %   This is a trick to avoid repeat fft's in the fully sampled readout direction 
    data = iFastFT(data,1,1);
    % This step stretches k-space in the y-direction

    %   Load the mask
    if SampleRate == '0'
        mask = single(ones(np,nv,ns,nt,nr));
    else
        load([PathToMask,filename],maskname)
        eval(['mask =' maskname ';']);
        eval(['clear ' maskname]);
        
        if ~ ReconWithFullRef && UndersampRef
            % Repeat the mask  on the first volume unless you made a 31
            % volume mask file
            refmask = mask(:,:,:,1,:);
            mask = cat(4,mask,refmask);
        elseif ReconWithFullRef && ~ UndersampRef
            % Use the FULLY SAMPLED ref vol in the CS recon
            refmask = single(ones(size(mask,1),size(mask,2),size(mask,3),1,nr));
            mask = cat(4,mask,refmask);
        % Otherwise, no refmask needed.
        end 
        
        %  Zero fill the mask
        if FillSizeX ~= np
            % Go from size 128.96 to 256.192
            zfmask = single(zeros(np,nv,ns,nt,nr));
            zfmask(65:192,49:144,:,:,:) = mask;
            mask = zfmask;
        end
        clear zfmask refmask
    end

    %   Pre-apply fftshift to data and mask
    %   This is a trick to avoid needed many fftshifts in the iterative recon
    data = fftshift(fftshift(data,2),3);
    mask = fftshift(fftshift(mask,2),3);

    %   Apply sampling mask and only keep acquired points
    data = data(mask==1); % data size depends on size of mask

    %   Setup recon options
    opt = SPSENSE_optset;
    opt.size = [np nv ns nt nr]; 
    opt.acqSize = [np nv ns nt nr];
    opt.compression_levels = [0 0 0 0 0 0 0];   %   Disable experimental feature. This can be handy for automatically setting the lambda's below
    % SET LAMBDAS
    % Initial Lambda (RML) [0 0.001 0.0001 0 0 0.001 0]
    % Thesis Lambda (MM) [0 0.0005 0.0001 0 0.005 0.005 0]
    if SampleRate == 0
        opt.lambda = [0 0 0 0 0 0 0];           %   Provides Wavelet, total variation, and a PCA constraint
    else
       opt.lambda = [0 0.0005 0.0001 0 0.005 0.005 0];
       % [?, wavelet, TV(1), k-spave filter (TV2), view-sharing, ref image,
       % dynamic gradient]
    end
    opt.Npc_local = 1;                          %   The constraint in the 4th dimention is a global PCA...
    opt.Npc = 2;                                %   ...with 2 components
    opt.FDorder1 = {1,2,3};                     %   Finite differences in all 3 dimensions
    opt.S = S; clear S;                         %   Coil sensitivity operators
    opt.U = find(mask);                         %   Undersampled Fourier operator
    opt.plot = 1;

    clear mask % To save memory

    %   Scale data near unity for consistent relation to regularization factors
    img = opt.iFTU(data,opt); % Fourier Transform?
    img = opt.iPI(img,opt); % Channel combination % This line breaks if sizes are different - MM
    sf = max(abs(img(:)));
    data = data./sf;
    clear img;

    %   Perform SparseSENSE recon
    imgR = SPSENSE_recon(data,opt);

    %   Rescale image
    imgR = imgR * sf; % This is scaled too much - MM
    % imgR is still the acquisition matrix size (128*96*64)

    %  MODIFY DATA FOR QMTLAB
    
    %  Resize Image - from 128.96 rectangle to 128*128 square, also get rid of top/bottom slices
    %  This step is superficial only and could be removed if rectangular images are prefered
    if FillSizeX ~= FillSizeY
        if ns ~= 1 % CHANGE THIS TO ASK FOR # OF EXTRA SLICES - currently assumes 2 extra slices on each side
            sqrImgR = zeros(FillSizeX,FillSizeX,ns-4,nt); 
            sqrImgR(:,(FillSizeX-FillSizeY)/2+1:FillSizeY+(FillSizeX-FillSizeY)/2,:,:) = imgR(:,:,3:ns-2,:); % Fill edges with zeros
        else
            sqrImgR = zeros(FillSizeX,FillSizeX,1,nt); 
            sqrImgR(:,(FillSizeX-FillSizeY)/2+1:FillSizeY+(FillSizeX-FillSizeY)/2,:,:) = imgR(:,:,:,:); % Fill edges with zeros
        end
    else
        sqrImgR = imgR(:,:,3:ns-2,:);
    end
    clear imgR
    
    %  Interpolation (using RESIZEM. Nearest neighbour interpolation. Size change only, no true interpolation)
     if FillSizeX ~= 256
        FinalImage = zeros(256,256,ns-4,nt);  % CHANGE THIS TO ASK FOR # OF EXTRA SLICES
        for vol = 1:nt
            for slice = 1:size(sqrImgR,3)
                FinalImage(:,:,slice,vol) = resizem(sqrImgR(:,:,slice,vol),2);
            end
        end
    else
        FinalImage = sqrImgR;
    end
    clear sqrImgR
    
    %  MT NORMALIZATION
    
    %  Get baseline image
    if ReconWithFullRef || UndersampRef
        M0 = FinalImage(:,:,:,nt); % last volume in stack (always true)
    elseif ~ ReconWithFullRef && ~ UndersampRef % Not tested
       % Add code to recon refVol (with could cobination and resizing) and add it to the image stack
    end
    
    % Perform MT-normalization (ie. divide by baseline image)
    for o = 1:nt-1
        scaledImgR(:,:,:,o) = FinalImage(:,:,:,o)./M0; % complex/complex
    end
        
    %  Store the data
    MTdata(:,:,:,((dataChunk-1)*numOffsets)+1:numOffsets*dataChunk) = abs(scaledImgR);
    
    
end % END OF DATA CHUNK

% Get rid of NaN values
MTdata(MTdata~=MTdata) = 0; % Caused by division of zeros in baseline

%  Save the data
%  3D data
save([OutputDir,SubjectID,outputFilename,'.mat'] , 'MTdata' );
% Save single slice version (2D)
MTdata = MTdata(:,:,SliceNum,:);
save([OutputDir,SubjectID,'/Slice30/',outputFilename,'.mat'] , 'MTdata' ); % Change storage location
% Save 3D unscaled data
save([OutputDirUnscaled,SubjectID,'/MTunscaled/',outputFilenameUnscaledData] , 'FinalImage' );


% Save lower resolution - single slice
% resizem converts this back to exactly what it was before it was scaled to
% 256*256. Could also just remove that section of code
new = zeros(128,128,1,10);
for i = 1:10
    new(:,:,1,i) = resizem(MTdata(:,:,1,i),0.5);
end
MTdata = new;
save([OutputDir,SubjectID,'/Slice30/',outputFilename,'_128.mat'] , 'MTdata' ); % Change storage location



