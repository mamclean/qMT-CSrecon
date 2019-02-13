% ratioVolume=zeros(size(MTdata)-[0 0 0 1]);
% for count = 1 : size( ratioVolume, 4 ) 
%  
%     ratioVolume(:,:,:,count) = ...
%         100*( MTdata(:,:,:,end) - MTdata(:,:,:,count) ) ...
%         ./ MTdata(:,:,:,end) ;
%     
% end
ratioVolume = MTdata;
% ratioVolume( ratioVolume < 0 ) = 0 ; ratioVolume( ratioVolume > 100 ) = 100 ;
%FourDviewer( ratioVolume ) ;
 
sliceNum = 1 ; 

% Comment or uncomment this for the first / other volumes
% clim = [0,1];
% figure;imagesc(ratioVolume(:,:,sliceNum,8),clim);
% mask = roipoly ;
% maskcon = mask ;
% mask = roipoly ;
% maskwat = mask ;
 
%
for count = 1 : size( ratioVolume, 4 )
     
    conValues = squeeze(MTdata(:,:,sliceNum,count)).*maskcon; %Make 2D
    conValues = conValues(:) ; %Make 1D
    conValues= conValues( not( conValues == 0 ) ) ; %Save only non-zero numbers
    maxCon = prctile(conValues,75)+(1.5*iqr(conValues)); %Determine Tukeys fences
    minCon = prctile(conValues,25)-(1.5*iqr(conValues));
    conValues(find(conValues>maxCon))=[]; %Remove Outliers
    conValues(find(conValues<minCon))=[]; 
    
    watValues = squeeze(MTdata(:,:,sliceNum,count)).*maskwat;
    watValues = watValues(:) ;    
    watValues= watValues( not( watValues == 0 ) ) ; 
    maxWat = prctile(watValues,75)+(1.5*iqr(watValues)); %Determine Tukeys fences
    minWat = prctile(watValues,25)-(1.5*iqr(watValues));
    watValues(find(watValues>maxWat))=[]; %Remove Outliers
    watValues(find(watValues<minWat))=[];
 
    %conditionerResponse(count) = sum( conValues ) ./sum(maskcon(:));
    %waterResponse(count) = sum( watValues ) ./sum(maskwat(:));
    conditionerResponse(count) = mean(conValues); % This offset image averaged over ROI
    waterResponse(count) = mean(watValues);
    stdConditionerResponse(count) = std( conValues ) ; % Variance in this ROI
     stdWaterResponse(count) = std( watValues ) ;
    
end
%  %%
% figure;plot(conditionerResponce/conditionerResponce(end));hold all ;plot(waterResponce/waterResponce(end));
% figure;errorbar(1:size( ratioVolume, 4 ),conditionerResponce/conditionerResponce(end),stdConditionerResponce/conditionerResponce(end)) ; hold all ;
% errorbar(1:size( ratioVolume, 4 ),waterResponce/waterResponce(end),stdWaterResponce/waterResponce(end)) ; 
%  
% asymeticZ = (conditionerResponce(16:30)-rot90(conditionerResponce(1:15),2))./conditionerResponce(end) ; 
% figure;plot(asymeticZ);
% figure;plot(waterResponce(16:30));hold all;plot(rot90(waterResponce(1:15)));
% figure;plot(conditionerResponce(16:30));hold all;plot(rot90(conditionerResponce(1:15)));