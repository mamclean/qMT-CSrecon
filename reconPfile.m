% ---------------------------------------------------------------------
% P-FILE READER - qMT Offset Scans
% Melany Mclean
% June 21, 2016
%
% MT-offset scans (saved as p-files) can be opened and reconstructed.
%
% NOT AN ACCURATE RECONSTRUCTION TECHNIQUE. GE uses many other data
% correction functions in its reconstruction. Better to use Orchestra (from
% GE).
% ---------------------------------------------------------------------

% HOW TO USE
% ---------------------------------------------------------------------
% FileNames         is a list of text strings that equals the names of the P-files {'file1','file2','file3'}
% DataDimension     is a vector [ kx , ky , kz , num coils ]
% KImage            is a matrix that is the k-space data
% RealImage         is a matrix that is the reconstructed image


function [ KImage , RealImage ] = reconPfile ( FileNames , DataDimension ) % inputs go here ( input1, input2, ...)

    % Start Timer
    tic;
    disp('Reconstructing P-file');
    
    % Determine the fourth dimension and final image dimensions
    dimension4 = length(FileNames);
    kx = DataDimension(1);
    ky = DataDimension(2);
    kz = DataDimension(3);
    numCoils = DataDimension(4);

    % Allocate empty storage places
    rawDataMatrixCO = zeros( kx , ky , kz , numCoils ,dimension4); % Complex matrix
    rawDataMatrixAB = zeros( kx , ky , kz , numCoils ,dimension4); % Absolute matrix
    imageDataMatrixCO = zeros( max(kx,ky) , max(kx,ky) , kz , numCoils ,dimension4); % Square in-plane %CHANGED
    imageCoilAvgMatrixCO = zeros( max(kx,ky) , max(kx,ky) , kz , dimension4); % CHANGED
    rawCoilAvgMatrixCO = zeros( max(kx,ky) , max(kx,ky) , kz , dimension4); % CHANGED
    
    % Open and store each file
    for imageNum = 1:dimension4

        fileID = fopen ( char(FileNames(imageNum)) );
        disp(fileID); % FOR DEBUGGING
        rawDataVector = fread ( fileID , 'int16' ); % P-files contain int16 values
        fclose ( fileID );   
        
        % Convert the raw data vector into a 5D array (x,y,z,coil,imageNum);
        cursor = 0; % mark place in vector
        for z = 1:kz
            for coil = 1:numCoils
                for y = 1:ky
                    for x = 1:kx
                    % DONT FORGET TO READ BACKWARD EVERY OTHER LINE
                    % Find the complex and absolute components of the voxel
                    if mod(y,2) == 0 % If even, read backwards
%                         realX = rawDataVector( (-2*x) + 0 + cursor + (2*kx) );
%                         imX = rawDataVector( (-2*x) + 1 + cursor + (2*kx) );
%                         radX = sqrt ( realX^2 + imX^2 ); % POSSIBLY REPLACE WITH ABS(...)
                        realX = -1 * rawDataVector((2*x) - 1 + cursor);
                        imX = -1 * rawDataVector(2*x + cursor);
                        radX = sqrt ( realX^2 + imX^2 );
                    else % If odd, read forward
%                         realX = rawDataVector((2*x) - 1 + cursor);
%                         imX = rawDataVector(2*x + cursor);
%                         radX = sqrt ( realX^2 + imX^2 ); % POSSIBLY REPLACE WITH ABS(...)
                        realX = rawDataVector((2*x) - 1 + cursor);
                        imX = rawDataVector(2*x + cursor);
                        radX = sqrt ( realX^2 + imX^2 ); % POSSIBLY REPLACE WITH ABS(...)
                    end
                    
                    % Store the voxel
                    rawDataMatrixAB(x,y,z,coil,imageNum) = radX; % COULD PROBABLY GET RID OF THIS....
                    rawDataMatrixCO(x,y,z,coil,imageNum) = realX + (1i*imX);
                    end
                    cursor = cursor + (2*kx); % We have gone through a block of x values
                end
            end
        end
        disp(cursor); % FOR DEBUGGING
    end
    
    % Perform Zero Filling
    filledDataMatrixCO = zeros(max(kx,ky),max(kx,ky),kz,numCoils,dimension4);
    filledDataMatrixCO(:,17:112,:,:,:) = rawDataMatrixCO;
   
    

    % Calculate Image
    for imageNum = 1:dimension4
        
        % DONT FORGET ABOUT ZERO FILLING!
        % Fourier transform each coil volume
        for coil = 1:numCoils
            imageDataMatrixCO(:,:,:,coil,imageNum) = ifft2(filledDataMatrixCO(:,:,:,coil,imageNum)); % CHANGED TO FILLED MATRIX
            
            % Average the coil signals - step 1, summation
            rawCoilAvgMatrixCO(:,:,:,imageNum) = rawCoilAvgMatrixCO(:,:,:,imageNum) + filledDataMatrixCO(:,:,:,coil,imageNum); % CHANGED TO FILLED
            imageCoilAvgMatrixCO(:,:,:,imageNum) = imageCoilAvgMatrixCO(:,:,:,imageNum) + imageDataMatrixCO(:,:,:,coil,imageNum);
        end
        % Average the coil signals - step 2, division
        rawCoilAvgMatrixCO(:,:,:,imageNum) = rawCoilAvgMatrixCO(:,:,:,imageNum)/numCoils;
        imageCoilAvgMatrixCO(:,:,:,imageNum) = imageCoilAvgMatrixCO(:,:,:,imageNum)/numCoils;
    end
    
    % Set output variables
    RealImage = imageCoilAvgMatrixCO;
    KImage = rawCoilAvgMatrixCO;
    
    % Save output variables to workspace
    assignin('base','realimage',RealImage);
    assignin('base','kspace',KImage);
    
    % Display raw and reconstructed slice; middle slice, middle image.
    figure;
    I = mat2gray(abs(KImage(:,:,round(kz/2),round(dimension4/2))));
    imshow(I)
    title('K-space')
    
    figure;
    J = mat2gray(abs(RealImage(:,:,round(kz/2),round(dimension4/2))));
    imshow(J)
    title('Reconstructed Image')
    
    % DEBUGGING FIGURES
    figure;
    K = mat2gray(rawDataMatrixAB(:,:,round(kz/2),2,1)); % coil 6, image 1
    imshow(K)
    title('Raw Data')
    
    % End Timer
    toc;
   
end