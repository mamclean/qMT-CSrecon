% Create Hamming Filters of size X * Y *Z
%
% Melany Mclean
% Oct 21, 2016

%--------------------------------------------------------------------------
X = 128;        % Input x-res
Y = 128;         % Input y-res
Z = 48;         % Input z-res
EdgeSize = 15;  % Input border of k-space to filter
PathToStorageFolder = '~/Documents/CSProject/CSpreparedData/HammingFilters/';
Filename = ['3DFilter_128*128*48_' num2str(EdgeSize) '.mat'];
%--------------------------------------------------------------------------

% Initialize things
sz = (2*EdgeSize)+1;
ham = hamming(sz);
hamx = ones(X,1);
hamy = ones(1,Y);
hamz = ones(Z,1);
ham2 = ones(X,Y);
ham3 = ones(X,Y,Z);

% Add Hamming Values
hamx(1:EdgeSize) = ham(1:EdgeSize);
hamx((X-EdgeSize)+1:X) = ham(EdgeSize+2:sz);

hamy(1,1:EdgeSize) = ham(1:EdgeSize);
hamy(1,(Y-EdgeSize)+1:Y) = ham(EdgeSize+2:sz);

hamz(1:EdgeSize) = ham(1:EdgeSize);
hamz((Z-EdgeSize)+1:Z) = ham(EdgeSize+2:sz);
hamz = permute(hamz, [3,2,1]);

% Build matrix
for i = 1:X
    ham2(i,:) = hamx(i) .* hamy;
    for j = 1:Y
        ham3(i,j,:)=ham2(i,j) .* hamz;
    end
end

save([PathToStorageFolder,Filename],'ham3');
    
    
    
    
    