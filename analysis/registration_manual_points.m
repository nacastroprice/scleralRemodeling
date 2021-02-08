cpselect(distorted,original,movingPoints,fixedPoints);

% transformation
tform = fitgeotrans(movingPoints,fixedPoints,'nonreflectivesimilarity');
% tform = fitgeotrans(movingPoints,fixedPoints,'affine');

% recover scale and angle
tformInv = invert(tform);
Tinv = tformInv.T;
ss = Tinv(2,1);
sc = Tinv(1,1);
scale_recovered = sqrt(ss*ss + sc*sc);
theta_recovered = atan2(ss,sc)*180/pi;

Roriginal = imref2d(size(original));
recovered = imwarp(distorted,tform,'OutputView',Roriginal);
figure;montage({original,distorted})
figure;montage({original,recovered})
figure;imshowpair(original,distorted);title("original and distorted");
figure;imshowpair(original,recovered);title("original and recovered");
%%%%%%%%%%%%%%%%%%%

[optimizer,metric] = imregconfig('monomodal');
movingRegistered = imregister(distorted, original, 'affine', optimizer, metric);
figure;imshowpair(original,movingRegistered);title("original and movingRegistered");


original_masked = original;
o_mask = roipoly(original_masked);
original_masked(o_mask) = 0;
figure;imshow(original_masked);
distorted_masked = distorted;
d_mask = roipoly(distorted_masked);
distorted_masked(d_mask) = 0;
figure;imshow(distorted_masked);

[movingRegistered_masked,Rreg,tform,Rmoving,Rfixed] = imregister2(distorted_masked, original_masked, 'affine', optimizer, metric);
figure;imshowpair(original_masked,movingRegistered_masked);title("original and movingRegistered_masked");

%% 
[movingRegistered_unmasked,Rreg] = imwarp(distorted,Rmoving,tform,'OutputView',Rfixed, 'SmoothEdges', true);
figure;imshowpair(original,movingRegistered_unmasked);title("original and movingRegistered_unmasked");
%% Feature based registration
%//////////////////////////////////////////////////////////////////////
directory = 'R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\';
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);
exp_folder = strcat(directory,num2str(exp_number));
brx_folder = num2str(brx_number);
exp_dir = exp_folder + "\" + brx_folder + "\";
img_folder = exp_dir + "images";

addpath(exp_dir);
% addpath("R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\organ_culture_app_new_cam\functions");

mask_function = @bead_mask;
% load(strcat(exp_dir, "fixed_roi.mat")); % load roi for masking


% edgeThreshold = 0.8;
% amount = 0.4;
% 
% 
% originalc = localcontrast(original, edgeThreshold, amount);
% distortedc = localcontrast(distorted, edgeThreshold, amount);
originalc = original;
distortedc = distorted;

endPoints = 20;

% set your rois 
figure;imshowpair(originalc,distortedc);
fixed_roi = drawpolygon();

bw_mask = createMask(fixed_roi,original);    
o_roi = fixed_roi;
d_roi = fixed_roi;
% find features
% ptsOriginal  = detectSURFFeatures(originalc,'NumScaleLevels',5);
% ptsDistorted = detectSURFFeatures(distortedc,'NumScaleLevels',5);
ptsOriginal  = detectSURFFeatures(originalc);
ptsDistorted = detectSURFFeatures(distortedc);
% exclude features if not in ROI
ftog = inROI(o_roi,double(ptsOriginal.Location(:,1)),double(ptsOriginal.Location(:,2)));
ftdis = inROI(d_roi,double(ptsDistorted.Location(:,1)),double(ptsDistorted.Location(:,2)));

ptsOriginal(ftog,:) = [];
ptsDistorted(ftdis,:) = [];

ptsOriginal = selectStrongest(ptsOriginal,20);
ptsDistorted = selectStrongest(ptsDistorted,20);

% ptsDistorted.Location = cpcorr(ptsDistorted.Location,ptsOriginal.Location,...
%                               distorted,original);

[featuresOriginal,  validPtsOriginal]  = extractFeatures(originalc,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distortedc, ptsDistorted);





% match features by using their descriptors
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
% retrieve locations in each image
matchedOriginal  = ptsOriginal(indexPairs(:,1));
matchedDistorted = ptsDistorted(indexPairs(:,2));
matchedOriginal = selectStrongest(matchedOriginal,endPoints);
matchedDistorted = selectStrongest(matchedDistorted ,endPoints);

figure;
showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');

% estimate transform
[tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(...
    matchedDistorted, matchedOriginal, 'similarity','MaxNumTrials',5000,'Confidence',99.99);

figure;
showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
title('Matching points (inliers only)');
legend('ptsOriginal','ptsDistorted');

% compute inverse transformation matrix
Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc);
thetaRecovered = atan2(ss,sc)*180/pi;
% tranform the distorted
outputView = imref2d(size(originalc));
recovered  = imwarp(distorted,tform,'OutputView',outputView);

[centerso] = findbeadsbyposition(originalc,mask_function);
[centersd] = findbeadsbyposition(distortedc,mask_function);
[centersr] = findbeadsbyposition(recovered,mask_function);

areao = calcArea(centerso);
aread = calcArea(centersd);
arear = calcArea(centersr);

figure;
title('Area of original, distorted and recovered')
x = categorical({'original','distorted','recovered'});
x = reordercats(x,{'original','distorted','recovered'});
y = (([areao,aread,arear]/areao)-1)*100;
bar(x,y);



figure, imshowpair(original,recovered);
hold on
plot(centerso(:,1),centerso(:,2),'b*')
plot(centersd(:,1),centersd(:,2),'r*')
plot(centersr(:,1),centersr(:,2),'g*')

hold off
%% Feature based registration (projective)
%//////////////////////////////////////////////////////////////////////
directory = 'C:\Users\nacastro\Documents\brx_exp\';
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);
exp_folder = strcat(directory,num2str(exp_number));
brx_folder = num2str(brx_number);
exp_dir = exp_folder + "\" + brx_folder + "\";
img_folder = exp_dir + "images";

addpath(exp_dir);
% addpath("R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\organ_culture_app_new_cam\functions");

mask_function = @bead_mask;
% load(strcat(exp_dir, "fixed_roi.mat")); % load roi for masking
% 
% bw_mask = createMask(fixed_roi,original);         
endPoints = 20;

edgeThreshold = 0.8;
amount = 0.4;


originalc = localcontrast(original, edgeThreshold, amount);
distortedc = localcontrast(distorted, edgeThreshold, amount);
% originalc = original;
% distortedc = distorted;

figure;imshowpair(originalc,distortedc);

% set your rois 
% figure;imshowpair(originalc,imsharpen(original),"montage");
fixed_roi = drawpolygon();

o_roi = fixed_roi;
d_roi = fixed_roi;
% find features
% ptsOriginal  = detectSURFFeatures(originalc,'NumScaleLevels',5);
% ptsDistorted = detectSURFFeatures(distortedc,'NumScaleLevels',5);
ptsOriginal  = detectMSERFeatures(originalc);
ptsDistorted = detectMSERFeatures(distortedc);

% exclude features if not in ROI
ftog = inROI(o_roi,double(ptsOriginal.Location(:,1)),double(ptsOriginal.Location(:,2)));
ftdis = inROI(d_roi,double(ptsDistorted.Location(:,1)),double(ptsDistorted.Location(:,2)));

ptsOriginal(ftog,:) = [];
ptsDistorted(ftdis,:) = [];

% ptsOriginal = selectStrongest(ptsOriginal,40);
% ptsDistorted = selectStrongest(ptsDistorted,40);

% ptsDistorted.Location = cpcorr(ptsDistorted.Location,ptsOriginal.Location,...
%                               distorted,original);

[featuresOriginal,  validPtsOriginal]  = extractFeatures(originalc,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distortedc, ptsDistorted);





% match features by using their descriptors
indexPairs = matchFeatures(featuresOriginal, featuresDistorted);
% retrieve locations in each image
matchedOriginal  = ptsOriginal(indexPairs(:,1));
matchedDistorted = ptsDistorted(indexPairs(:,2));
% matchedOriginal = selectStrongest(matchedOriginal,endPoints);
% matchedDistorted = selectStrongest(matchedDistorted ,endPoints);

figure;
showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
title('Putatively matched points (including outliers)');

% estimate transform
[tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(...
    matchedDistorted, matchedOriginal, 'similarity','MaxNumTrials',5000,'Confidence',99.999);

figure;
showMatchedFeatures(original,distorted,inlierOriginal,inlierDistorted);
title('Matching points (inliers only)');
legend('ptsOriginal','ptsDistorted');

% compute inverse transformation matrix
Tinv  = tform.invert.T;

ss = Tinv(2,1);
sc = Tinv(1,1);
scaleRecovered = sqrt(ss*ss + sc*sc);
thetaRecovered = atan2(ss,sc)*180/pi;
% tranform the distorted
outputView = imref2d(size(original));
recovered  = imwarp(distorted,tform,'OutputView',outputView);

[centerso] = findbeadsbyposition(original,mask_function, bw_mask);
[centersd] = findbeadsbyposition(distorted,mask_function, bw_mask);
[centersr] = findbeadsbyposition(recovered,mask_function, bw_mask);

areao = calcArea(centerso);
aread = calcArea(centersd);
arear = calcArea(centersr);

figure;
title('Area of original, distorted and recovered')
x = categorical({'original','distorted','recovered'});
x = reordercats(x,{'original','distorted','recovered'});
y = (([areao,aread,arear]/areao)-1)*100;
bar(x,y);



figure, imshowpair(original,recovered);
hold on
plot(centerso(:,1),centerso(:,2),'b*')
plot(centersd(:,1),centersd(:,2),'r*')
plot(centersr(:,1),centersr(:,2),'g*')

hold off