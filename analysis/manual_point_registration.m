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
load(strcat(exp_dir, "fixed_roi.mat")); % load roi for masking

bw_mask = createMask(fixed_roi,original);





[movingPoints,fixedPoints] = cpselect(distorted,original,'Wait',true);

% use correlation to improve moving point location
movingPointsAdjusted = cpcorr(movingPoints,fixedPoints,distorted,original)

figure; imshow(distorted)
hold on
plot(fixedPoints(:,1),fixedPoints(:,2),'xr') 
plot(movingPoints(:,1),movingPoints(:,2),'xg') 
plot(movingPointsAdjusted(:,1),movingPointsAdjusted(:,2),'xy') 


% transformation
tform = fitgeotrans(movingPointsAdjusted,fixedPoints,'similarity');
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

% get centers and area
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