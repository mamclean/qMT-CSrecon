% ---------------------------------------------------------------------
% CREATE UNDERSAMPLING MASKS FOR RETROSPECTIVE CS
%
% Melany Mclean
% April 29, 2016
%
% Four dimensional MT-offset scans can be undersampled in the phase encode
% and slice select directions, as well as the MT-offset direction.
% ---------------------------------------------------------------------

% ENTER IMAGE DATA
clear
% ---------------------------------------------------------------------
UnderSamp = 4;   % Enter RATE of undeersampling
VersionNum = 10;  % Enter the version number for this sampling rate and type
CoreSize = 10;   % Enter PERCENT of samples in the core
FE = 128;        % Enter matrix size in frequency encode direction
PE = 96;         % Enter matrix size in phase encode direction
SS = 64;         % Enter matrix size in slice select direction
MT = 10;         % Enter number of MT-offsets
CO = 32;         % Enter the number of virtual coils % Equal to or a fraction of the number of channels
MaxPD = 3;       % Enter maximum radius of Poisson Disc
Dist = 'POWER';  % Enter distribution type 'POWER', 'POISS', or 'EXP'
Lam = 300;       % Enter lambda of poisson distribution of samples (if using poisson)
%PathToMatrixOutput = '~/Documents/CSProject/Subsample_Patterns/128.96.64.30';  % MAY NOT NEED THIS IN THE FUTURE
PathToMaskOutput = ['~/Documents/CSProject/CSpreparedData/Masks/',num2str(FE),'.',num2str(PE),'.',num2str(SS),'.',num2str(MT),'/']; % Saves the masks used by ReconQMT.m
%PathToDataOutput = '~/Documents/CSProject/UndersampleMasks/StaticPatterns';
%PathToTextOutput = '~/Documents/CSProject/Subsample_Lookup_Tables/128.96.64.10';
% ---------------------------------------------------------------------

% INITIALIZE VARIABLE

% In-Slice values
% Centre of k-space
midX = (SS+(2*MaxPD))/2;
midY = (PE+(2*MaxPD))/2;
slice_size = PE*SS;
% Calculate k-space properties
num_samps = round(0.01 * (100/UnderSamp) * slice_size); % number of samples in the slice
core_r = round(sqrt((num_samps * 0.01 * CoreSize)/pi)); % radius of core
[rr cc] = meshgrid(1:PE+(2*MaxPD),1:SS+(2*MaxPD)); % core template
circle = sqrt((rr-midY).^2+(cc-midX).^2)<=core_r; % core circle
core_samps = size(find(circle==1),1); % number of samples in the core
rand_samps = num_samps - core_samps; % number of samples distributed pseudo randomly

% Volume storage
vol = zeros(SS+(2*MaxPD),PE+(2*MaxPD),MT+(2*MaxPD)); % matrix with extra padding for Poisson Discs
% Prob. Dens. Func.  storage
pdf = zeros(FE,PE,SS,MT);
% Radii storage
lookupRad = zeros(1,rand_samps);
actualRad = zeros(1,rand_samps);
% Poisson Disc Lookup Values
% Pre-defined discs for each "level" of sampling density
disc1 = [1,0;0,-1;-1,0;0,1]; % lookup using disc1(cell,x or y)
disc2 = [1,1;-1,-1;1,-1;-1,1];
disc3 = [0,2;2,0;0,-2;-2,0];
disc4 = [2,1;2,-1;1,-2;-1,-2;-2,-1;-2,1;-1,2;1,2];
disc5 = [2,2;3,0;2,2;0,-3;-2,-2;-3,0;-2,2;0,3];
disc6 = [1,3;3,1;3,-1;1,-3;-1,-3;-3,-1;-3,1;-1,3];
% Can have up to 3 slice radii, only use slice_rad = 2 when disc_rad >= 4.

% Distribution variables
if strcmp(Dist,'POISS')
% Adjust lambda of larger dimension
    if SS == PE
        LamX = Lam;
        LamY = Lam;
        %gamma = 0.25*SS;
    elseif SS > PE
        LamX = (SS/PE)*Lam;
        LamY = Lam;
        %gamma = 0.25*SS;
    else
        LamY = (PE/SS)*Lam;
        LamX = Lam; 
        %gamma = 0.25*PE;
    end
elseif strcmp(Dist,'EXP')
    gammaX = 0.25*PE;
    gammaY = 0.25*SS;
    gamma = 0.5*((PE+SS)/4); % 1/4 the radius length
    
elseif strcmp(Dist,'POWER')
    powerDecay = @(x) (x.^3)/3 - x.^2 + x; % CDF
    pdfFunction = @(x) 1/(x.^2); % PDF
    intList = linspace(2/PE,1,PE/2);
    for i = 0 : (PE/2)-1
        pdf(:,(PE/2)-i,:,:) = pdfFunction(intList(i+1));
        pdf(:,(PE/2)+1+i,:,:) = pdfFunction(intList(i+1));
    end
    pdf(:,(PE/2)-core_r:(PE/2)+1+core_r,:,:) = 1;% Get rid of really high values
    pdf = pdf/max(pdf(:)); % Scale data
    pdf(:,(PE/2)-core_r:(PE/2)+1+core_r,:,:) = 1;% Place the ones back
    
    % Create a list of power distributed radii (radDecay)
    % using power decay (1-r)^d where d = 2 as suggested by Ziljstra
    uniformList = rand(1,10000);
    maxR = sqrt((PE/2)^2 + (SS/2)^2);
    decayList = powerDecay(uniformList);
    phaseDecay = round((PE/2) - (decayList .* 3 .* (PE/2)));
    sliceDecay = round((SS/2) - (decayList .* 3 .* (SS/2)));
    radDecay = maxR - (decayList .* 3 .* maxR);
end


% FOR EACH MT OFFSET "SLICE" 
for mt = 1:MT
    % take the blank slice off the volume "shelf"
    slice = vol(:,:,mt+MaxPD);
   
    % place the fully sampled central core
    slice = slice + circle;
    
    % place the rest of the samples
    s = 1;
    
    % FOR EACH SAMPLE
    while s <= rand_samps
        
        % CHOOSE A LOCATION
        if strcmp(Dist,'POISS')
            % using poisson distribution with lambda distribution parameter
            x = round(poissrnd(LamX)-(LamX-midX)); % shift distribution to middle of "slice"
            y = round(poissrnd(LamY)-(LamY-midY));
            r = sqrt((x-midX)^2+(y-midY)^2);
        elseif strcmp(Dist,'EXP')
            % using exponential distribution with gamma distribution parameter
            x = round(((-1)^randi(2) * exprnd(gammaX)) + midX);
            y = round(((-1)^randi(2) * exprnd(gammaY)) + midY);
            r = sqrt((x-midX)^2+(y-midY)^2);
        elseif strcmp(Dist,'POWER')  
            r = randsample(radDecay,1); % randomly chosen r
            lookupRad(s) = r; % Sample chosen from theoretical distribution

            if round(r) <= 1 
                x = midX;
                y = midY;
            else
                xx = 1.3 * r.*rand(1); % any number between zero and 1.3 r
                yy = 1.3 * r.*rand(1); % note: calculating y from x yields small number rounding errors
                x = (-1)^randi(2) * round(xx) + midX; % positive or negative, shifted to middle
                y = (-1)^randi(2) * round(yy) + midY;
                r = sqrt((x-midX)^2+(y-midY)^2); % actual r
                actualRad(s) = r; % this information is used to ensure sampling distribution not distorted from theoretical
            end 
        end
        
        % CHECK IF IT'S INSIDE PLANE
        if x > MaxPD && x <= (SS+MaxPD) && y > MaxPD && y <= (PE+MaxPD)
           
            % CHECK IF IT'S AVAILABLE
            if slice(x,y) == 0
                % take the spot
                slice(x,y) = 1;
                s = s+1;
                
                % DETERMINE RADIUS OF IT'S DISC
               
                % determine probability of the radius being equal or less
                % than r
                if strcmp(Dist,'POISS')
                    p1 = poisscdf(Lam+r,Lam); % prob less than/equal radius (Area left of centre)
                    p2 = poisscdf(Lam-r,Lam); % prob below radius (Left tail)
                elseif strcmp(Dist,'EXP');
                    p1 = expcdf(r,gamma);
                    p2 = expcdf(-r,gamma);
                elseif strcmp(Dist,'POWER')
                    % Used a normal distribution for this
                    p1 = normcdf(r,0,(maxR/4)); %(x,mu,sig(std dev)) where 3*sig ~ zero prob
                    p2 = normcdf(-r,0,(maxR/4)); % sig = maxR/3 too big!
                end
                prob_r = p1-p2;% The larger the middle is, the wider the radius is
                
                
                % Assign disk size based on probability
                % Values designed for best results with 16x undersample
                %%% Better to adjust the disc size than the prob. value %%%
                if prob_r <= 0.15   % 15th %ile - no disc
                    disc_rad = 0;
                    slice_rad = 0;
                elseif prob_r <= 0.35 % 35th %ile - 1 voxel in plane and 1 betwen planes
                    disc_rad = 1/((100/UnderSamp)/6.25);
                    slice_rad = 1/((100/UnderSamp)/6.25);
                elseif prob_r <= 0.65 % 35th %ile - 2 voxel in plane and 1 betwen planes
                    disc_rad = 2/((100/UnderSamp)/6.25);
                    slice_rad = 1/((100/UnderSamp)/6.25);
                elseif prob_r <= 0.85
                    disc_rad = 3/((100/UnderSamp)/6.25);
                    slice_rad = 2/((100/UnderSamp)/6.25);
                else
                    disc_rad = 4/((100/UnderSamp)/6.25);
                    slice_rad = 3/((100/UnderSamp)/6.25);
                end

                % PLACE A POISSON DISC
                % Flag a "disc" voxel with a value of 2
                % first disc, radius 1
                if disc_rad >= 1
                    for i = 1:length(disc1)
                        slice(x+disc1(i,1),y+disc1(i,2)) = slice(x+disc1(i,1),y+disc1(i,2))+2;
                        % set 
                    end
                end
                % second disc, radius 2
                if disc_rad >= 2
                    for i = 1:length(disc2)
                        slice(x+disc2(i,1),y+disc2(i,2)) = slice(x+disc2(i,1),y+disc2(i,2))+2;
                    end
                end
                % third disc, radius 3
                if disc_rad >= 3
                    for i = 1:length(disc3)
                        slice(x+disc3(i,1),y+disc3(i,2)) = slice(x+disc3(i,1),y+disc3(i,2))+2;
                    end
                end
                
                % PLACE DISK IN OTHER SLICES
                % only stop next slice form being there
                if slice_rad >= 1
                    vol(x,y,mt+MaxPD+1) =  vol(x,y,mt+MaxPD+1)+2;
                end
                % stop next slice from being near
                if slice_rad >= 2
                    for i = 1:length(disc1)
                        vol(x+disc1(i,1),y+disc1(i,2),mt+MaxPD+1) = vol(x+disc1(i,1),y+disc1(i,2),mt+MaxPD+1)+2;
                    end
%                     for i = 1:length(disc2)
%                         vol(x+disc2(i,1),y+disc2(i,2),mt+MaxPD+1) = vol(x+disc2(i,1),y+disc2(i,2),mt+MaxPD+1)+2;
%                     end
                end
                % stop next 2 slices from being there
                if slice_rad >= 3
                    vol(x,y,mt+MaxPD+2) =  vol(x,y,mt+MaxPD+2)+2;
                end
            
            else
                % if point not available, restart while loop
                continue
            end
        else 
            % if out of plane, restart while loop
            continue
        end % end placing sample 
    end % end picking sample 
    
    % STORE THE SLICE
    final_slice = rem(slice,2); % remove poisson disks (even numbers, flagged by "2")
    vol(:,:,mt+MaxPD) = final_slice; 
    disp('Progress');
    disp(mt/MT);
    
end % end 'slice' loop

% SAVE THE FINAL RESULT

% Trim off edges
finalSamplingPattern = vol(MaxPD+1:SS+MaxPD,MaxPD+1:PE+MaxPD,MaxPD+1:MaxPD+MT);

% Save the matrix
% cd (PathToMatrixOutput);
% filename = ['StaticPattern_',num2str(UnderSamp),'_',Dist,'_V', num2str(VersionNum), '.mat'];
%save(filename , 'finalSamplingPattern')

% SAVE THE MATRIX AS A 4D MASK - I don't think this is working ...
% Convert matrix into 4D sampling pattern
    fullSampleMatrix = zeros(FE,PE,SS,MT); % CHANGED from 30
    for v = 1:size(finalSamplingPattern,3) 
        for z = 1:size(finalSamplingPattern,1) 
            for y = 1:size(finalSamplingPattern,2) 
                fullSampleMatrix(:,y,z,v) = finalSamplingPattern(z,y,v); 
                % all the x's are zero, or all the x's are one
            end
        end
    end
    
 mask = single(zeros([size(fullSampleMatrix),CO]));
 for c = 1:CO
     mask(:,:,:,:,c) = logical(fullSampleMatrix); % Same mask on every coil
 end
%  filename = ['Mask_',num2str(UnderSamp),'_V',num2str(VersionNum),'.mat'];
%save([PathToDataOutput,filename], 'mask');

% SAVE THE MATRIX AS 5D MASK IN CSpreparedData/SubjectID

% NOTE: This only works for 1 processing chunk right now!
mask1 = ones(FE,PE,SS,MT,CO);
for c = 1:CO
    mask1(:,:,:,:,c) = fullSampleMatrix;
end
filename = ['CSinputData_',num2str(UnderSamp),'_V',num2str(VersionNum),'_P1.mat'];
save([PathToMaskOutput,filename], 'mask1');

% Save ASCII coordinates
% coord = zeros(num_samps,2,size(finalSamplingPattern,3));
% cd (PathToTextOutput);
% for mtOffst = 1:size(finalSamplingPattern,3)
%     slice = finalSamplingPattern(:,:,mtOffst);
%     [kz,ky] = find(slice == 1);
%     coord(:,1,mtOffst) = ky;
%     coord(:,2,mtOffst) = kz;
%     filedata = coord(:,:,mtOffst); % Save 2 col matrix
%     filename = ['EncodeLookup_',num2str(UnderSamp),'_',num2str(mtOffst),'.txt'];
% %     save(filename,'filedata','-ascii','-tabs'); % Tab Separated
% %     csvwrite(filename,filedata); % Comma Separated
% 
% %   Ethan's double comma format
%     filename = ['EncodeLookup_',num2str(UnderSamp),'.txt'];
%     fid = fopen(filename,'a');
%     for i = 1:num_samps
%         fprintf( fid, '%i,%i,\n', ky(i), kz(i) ); 
%     end
%     fclose(fid);
% end
cd ..;

% % Save the pdf % NEW
% filename = ['PDFStaticPattern_',num2str(UnderSamp),'_',Dist, num2str(Num), '.mat'];
% save(filename , 'pdf')


% DISPLAY THE UNDERSAMPLE PATTERN

% % Figure with dots
% figure;
% % only the first MT "slice"
% spy(vol(:,:,3+MaxPD),'.', 6) %size changed from 5
% grid on

% % Figure with grid lines
% figure;
% hold on
% bar = lines;
% bar(1,:) = [ 1 1 1 ];
% bar(end,:) = [ 0.5 0.5 0.5 ];
% colormap(bar);
% h = pcolor(vol(:,:,3+MaxPD));
% axis equal tight

% Figure with squares
figure;
bar = lines;
bar(1,:) = [ 1 1 1 ];
colormap(bar);
imagesc(finalSamplingPattern(:,:,3)); % third slice % CHANGED
xlabel('Ky') % x-axis label
ylabel('Kz') % y-axis label
axis equal tight

% Figure of all MT slices
figure;
[rowInd,colInd,zInd]=ind2sub(size(finalSamplingPattern),find(finalSamplingPattern)); % CHANGED
zList = zInd(find(zInd<(34))); % get rid of end slices with only poisson disks
scatter3(colInd(1:length(zList)), rowInd(1:length(zList)), zList,'.');
xlabel('Ky') % x-axis label
ylabel('Kz') % y-axis label
zlabel('MT-offset') % z-axis label

% Figure of histograms
if strcmp(Dist,'POWER')
    % remeber not to use poisson discs
    figure;
    h1 = histogram(lookupRad,20);
    hold on
    h2 = histogram(actualRad,20);
    legend('lookup','actual')
end

