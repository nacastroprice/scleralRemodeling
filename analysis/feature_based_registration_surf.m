function [recovered] = feature_based_registration(original,distorted,fixed_roi)
% Registers and two images using SURF Features and similarity transform.
% Automatically detects features and computes transform based on two
% strongest points.
%   Inputs:
%  original = fixed image
%  distorted = moving image
%  roi = roi (will not be used to transform computation)
edgeThreshold = 0.6;
amount = 0.3;
endPoints = 2;

originalc = localcontrast(original, edgeThreshold, amount);
distortedc = localcontrast(distorted, edgeThreshold, amount);
% originalc = original;
% distortedc = distorted;

% set your rois 
% figure;imshow(original);
% o_roi = drawpolygon();
% figure;imshow(distorted);
% d_roi = drawpolygon();
o_roi = fixed_roi;
d_roi = fixed_roi;
% find features
ptsOriginal  = detectSURFFeatures(originalc,'NumScaleLevels',5);
ptsDistorted = detectSURFFeatures(distortedc,'NumScaleLevels',5);
% ptsOriginal  = detectSURFFeatures(originalc);
% ptsDistorted = detectSURFFeatures(distortedc);
% exclude features if not in ROI
ftog = inROI(o_roi,double(ptsOriginal.Location(:,1)),double(ptsOriginal.Location(:,2)));
ftdis = inROI(d_roi,double(ptsDistorted.Location(:,1)),double(ptsDistorted.Location(:,2)));

ptsOriginal(ftog,:) = [];
ptsDistorted(ftdis,:) = [];

ptsOriginal = selectStrongest(ptsOriginal,20);
ptsDistorted = selectStrongest(ptsDistorted,20);

% ptsDistorted.Location = cpcorr(ptsDistorted.Location,ptsOriginal.Location,...
%                               distorted,original);

[featuresOriginal,  validPtsOriginal]  = extractFeatures(original,  ptsOriginal);
[featuresDistorted, validPtsDistorted] = extractFeatures(distorted, ptsDistorted);





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

end

