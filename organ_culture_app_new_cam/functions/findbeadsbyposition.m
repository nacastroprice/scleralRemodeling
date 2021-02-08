
function [centers,BW,radii] = findbeadsbyposition(img,mask_function,varargin)
% Here goes the processing
% img = medfilt2(img);
% img = medfilt2(img);
% img = medfilt2(img);
% img = medfilt2(img);
% Create image mask
[BW,maskedImage] = mask_function(img);

% BW = img;
% T = img < 150;
% BW(T) = 0;
% % cluster what youre interested in
% num_clust = 5;
% [L,Centers] = imsegkmeans(BW,num_clust);
% B = labeloverlay(img,L);
% imshow(B)
% title('Labeled Image')
% BW = L == 2;
% figure;imshow(c)

if nargin == 4
    circleImage = varargin{nargin-3};
    % mask the image
    maskedImage = BW; % Initialize with the entire image.
    maskedImage(~circleImage) = 0; % Zero image outside the circle mask.
    % imshow(maskedImage)
    BW = maskedImage;
end


% BW = imopen(BW,se);

% This is to find all objects and filter them 
obj_foudn = regionprops(BW,'Centroid','Area','Circularity','Solidity','EquivDiameter');
table_obj = struct2table(obj_foudn);


% The break point should go here for debug: (1) uncomment the code at the bottom to plot the detected centers. 
% (2) adjust the values to address the object characteristics found in table_obj; most likely it will require a small circularity or
% solidity adjustment

indexes = table_obj.Area < 2300; % min and max area
table_obj(indexes,:) = [];
indexes = table_obj.Area > 8000;
table_obj(indexes,:) = [];
indexes = table_obj.Circularity < 0.4; % circularity
table_obj(indexes,:) = [];
indexes = table_obj.Solidity < 0.5; % solidity
table_obj(indexes,:) = [];

% bw_only_beads_and_nerve = bwareafilt(BW,nbeads + 1);
% bw_only_beads_and_nerve = bwareafilt(bw_only_beads_and_nerve,nbeads, 'smallest');
% imshow(bw_only_beads_and_nerve)
%[centers,radii,metric] = imfindcircles(BW,[8 20],'ObjectPolarity','bright');

% max_diam = max(table_obj.EquivDiameter);
% max_radii = ceil(max_diam/2);
% se = strel('disk',max_radii);
% dial = imdilate(BW,se);
% 
% % This is to find all objects and filter them 
% obj_foudn = regionprops(dial,'Centroid','Area','Circularity','Solidity','EquivDiameter');
% table_obj = struct2table(obj_foudn);
% 
% 
% % The break point should go here for debug: (1) uncomment the code at the bottom to plot the detected centers. 
% % (2) adjust the values to address the object characteristics found in table_obj; most likely it will require a small circularity or
% % solidity adjustment
% indexes = table_obj.Area < 1000; % min and max area
% table_obj(indexes,:) = [];
% indexes = table_obj.Area > 2000;
% table_obj(indexes,:) = [];
% indexes = table_obj.Circularity < 0.95; % circularity
% table_obj(indexes,:) = [];
% indexes = table_obj.Solidity < 0.75; % solidity
% table_obj(indexes,:) = [];

centers(:,1) = table_obj.Centroid(:,1);
centers(:,2) = table_obj.Centroid(:,2);
radii(:,1) = ceil(table_obj.EquivDiameter(:,1)/2);
% Use this figure to set the parameters above
%                 figure;
%                 imshow(img)
%                 hold on
%                 plot(centers(:,1),centers(:,2),'b*')
%                 hold off

end