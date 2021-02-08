img1 = imread('R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\1119\1\images\11191_20112020 170923.tif');
load('R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\1119\1\1119_1_retrieved_from_pics_calib_imhisteq.mat');

img1 = undistortImage(img1,cameraParams);
figure;imshow(img1)
img1_ref = imref2d(size(img1));
sub_img1 = img1;

sub_img1_ref = imref2d(size(sub_img1));
figure;imshow(sub_img1)
% Do translation for middle section
xdisp = 0.05;
y = 0.06;
T = [1 0 0;0 1 0;xdisp y 1]; % in x direction
%T = [x 0 0;0 y 0;0 0 1]; % in x direction
tform = affine2d(T);
    
%[img_trans,img_trans_ref] = imwarp(img1,tform);
sameAsInput = affineOutputView(size(sub_img1),tform,'BoundsStyle','SameAsInput');
sub_img1 = im2double(sub_img1);
[img_trans,img_trans_ref] = imwarp(sub_img1,tform,'OutputView',sameAsInput);
img_trans = im2double(img_trans);

figure;
subplot(1,2,1);
imshow(img1,img1_ref);
subplot(1,2,2);
imshow(img_trans,img_trans_ref)



real_points = 10:17;
ref_points = 1:27;
ref_points(real_points) = [];

 
% get the centers from the last entry and the current entry
%      for t = 1:28
%         centersUINT(t,:) = dtable{1,5+t};
%      end
%          % Get reference points
  [centers1,bw1,radii1] = findbeadsbyposition(sub_img1,mask_function);%,globe_mask
   a = [1:4]'; b = num2str(a); labels = cellstr(b);
    figure;
    imshow(sub_img1)
    hold on
    h = labelpoints(centers1(:,1), centers1(:,2), b, 'N', 0.15);
    hold off
%   radiusr = ceil(mean(radii1));
%   method = 'PhaseCode';
%   sens = 0.9999;
%   edge_thresh = 0.3;
%   [centers1r, radii1r, metric1r] = imfindcircles(sub_img1,radiusr,'Method',method,'Sensitivity',sens,'Edgethreshold',edge_thresh);
 
%   imshow(sub_img1)
%   viscircles(centers1r, radii1r,'EdgeColor','b');
% 
%   centers1p = centers1;
%   crop_factor = 0.53;
%   hshift = 0;
%   vshift = 0;
%   mx = "max";
% % % % %     % Implementation of intensity fitting algorith for subpixel accuracy
%     for k =1:length(centers1)
%         referenced = find_bead_center_intensity_fit(sub_img1,centers1(k,:),radii1(k),"poly22",crop_factor,hshift,vshift,mx);
%         centers1(k,:) = referenced;
%     end

        centers_x_first_ref = centers1(ref_points',1);
        centers_y_first_ref = centers1(ref_points',2);
        centers1_ref = centers1(ref_points',:);
         % Get real trackable points
        centers_x_first_real = centers1(real_points',1);
        centers_y_first_real = centers1(real_points',2);
        centers1_real = centers1(real_points',:);
        area1 = calculate_area_from_coords(centers1_real);
        
  % No get coords for transformed pic
  [centers2,bw2,radii2] = findbeadsbyposition(img_trans,mask_function);%,globe_mask

%   [centers2r, radii2r, metric2r] = imfindcircles(img_trans,radiusr,'Method',method,'Sensitivity',sens,'Edgethreshold',edge_thresh);
  %,'EdgeThreshold',edge_thresh
%   imshow(img_trans)
%   viscircles(centers2r, radii2r,'EdgeColor','b');
  
%   centers2p = centers2;
%       for k =1:length(centers1)
%         referenced = find_bead_center_intensity_fit(img_trans,centers2(k,:),radii2(k),"poly22",crop_factor,hshift,vshift,mx);
%         centers2(k,:) = referenced;
%     end
  
nbeads = 4;
   centersx = centers2(:,1);
   centersy = centers2(:,2);
   sortedcenters = zeros(nbeads,2);
%         Point are named after first pic: the closest points
%         to the last iteration get the name
    for j = 1:nbeads
        [dif, I] = min(sqrt((centers1(j,1)-centersx).^2 + (centers1(j,2)-centersy).^2));
        sortedcenters(j,:) = centers2(I,:);
        centersx(I,:) = [];
        centersy(I,:) = [];
        centers2(I,:)  = [];            
    end

%     
%         centers1(27,:) = [];
%         centers2(27,:) = [];
centers2 = sortedcenters;
    diff_vec = [ones(length(centers2),1)*xdisp,ones(length(centers2),1)*y];
    diff_cent = (centers1-centers2) + diff_vec;
    diff_mean = sqrt(sum(diff_cent.^2)/length(diff_cent))

    
  error = [0.0864,0.0857;0.0515,0.0592;0.0660,0.0921;0.0434,0.0625;0.0682,0.0669];
  error_old_cam = [0.0461,0.0350;0.0252,0.0173;0.0734,0.0356;0.0084,0.0354];  
    
    %%
    
         % Get reference points
        centers_x_second_ref = sortedcenters(ref_points',1);
        centers_y_second_ref = sortedcenters(ref_points',2);
        sortedcenters_ref = sortedcenters(ref_points',:);

         % Get real trackable points
        centers_x_second_real = sortedcenters(real_points',1);
        centers_y_second_real = sortedcenters(real_points',2);
        sortedcenters_real = sortedcenters(real_points',:);

    
        % uci is the displacement in reference x = vector 1xpoints of x'-x
    uci  = centers_x_second_ref - centers_x_first_ref;
    % uvi is the displacement in reference y = vector 1xpoints of y'-y
    vci  = centers_y_second_ref - centers_y_first_ref;
    % ui is the displacement in real x = vector 1xpoints of x'-x
    ui  = centers_x_second_real - centers_x_first_real;
    % vi is the displacement in real y = vector 1xpoints of y'-y
    vi  = centers_y_second_real - centers_y_first_real;
    
    coeff_num = 8;
    num_points_included = length(ref_points);
    % create the x matrix and the d vector 
     d = zeros(num_points_included*2,1);
     x = zeros(num_points_included*2,coeff_num);
     d = [uci; vci];
%      x = [ones(num_points_included,1), centers_x_first_ref, centers_y_first_ref, centers_x_first_ref.^2, centers_x_first_ref.*centers_y_first_ref, zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), (centers_x_second_ref.*(centers_x_second_ref.^2+centers_y_second_ref.^2)-centers_x_first_ref.*(centers_x_first_ref.^2+centers_y_first_ref.^2)); ...
%      zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), centers_x_first_ref, centers_y_first_ref, centers_x_first_ref.*centers_y_first_ref, centers_y_first_ref.^2, (centers_y_second_ref.*(centers_x_second_ref.^2+centers_y_second_ref.^2)-centers_y_first_ref.*(centers_x_first_ref.^2+centers_y_first_ref.^2))];
     x = [ones(num_points_included,1), centers_x_first_ref, centers_y_first_ref, centers_x_first_ref.*centers_y_first_ref, zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1); ...
      zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), centers_x_first_ref, centers_y_first_ref, centers_x_first_ref.*centers_y_first_ref];
    % Solve for the parameters
     psol = (transpose(x) * x)\transpose(x)*d;
     % Assign the parameters
%      a0 = psol(1);
%      a1 = psol(2);
%      a2 = psol(3);
%      a3 = psol(4);
%      a4 = psol(5);
%      b0 = psol(6);
%      b1 = psol(7);
%      b2 = psol(8);
%      b3 = psol(9);
%      b4 = psol(10);
%      k1 = psol(11);
     a0 = psol(1);
     a1 = psol(2);
     a2 = psol(3);
     a3 = 0;
     a4 = psol(4);
     b0 = psol(5);
     b1 = psol(6);
     b2 = psol(7);
     b3 = psol(8); 
     b4 = 0;
     k1 = 0;
     % Correct the initial points with the new displacements
        for i = 1:length(real_points)
            uic(i) = ui(i) - a0 - a1*centers_x_first_real(i) - a2*centers_y_first_real(i)- a3*centers_x_first_real(i)^2 - a4*centers_x_first_real(i)*centers_y_first_real(i) + k1*(centers_x_second_real(i)*(centers_x_second_real(i)^2+centers_y_second_real(i)^2)-centers_x_first_real(i)*(centers_x_first_real(i)^2+centers_y_first_real(i)^2));
            vic(i) = vi(i) - b0 - b1*centers_x_first_real(i) - b2*centers_y_first_real(i)- b3*centers_x_first_real(i)*centers_y_first_real(i) - b4*centers_y_first_real(i)^2 + k1*(centers_y_second_real(i)*(centers_x_second_real(i)^2+centers_y_second_real(i)^2)-centers_y_first_real(i)*(centers_x_first_real(i)^2+centers_y_first_real(i)^2));
            centers_corrected(i,1) = centers_x_first_real(i) + uic(i);
            centers_corrected(i,2) = centers_y_first_real(i) + vic(i);
        end
%        for i = 1:length(ref_points)
%             uicr(i) = uci(i) - a0 - a1*centers_x_first_ref(i) - a2*centers_y_first_ref(i)- a3*centers_x_first_ref(i)^2 - a4*centers_x_first_ref(i)*centers_y_first_ref(i) + k1*(centers_x_second_ref(i)*(centers_x_second_ref(i)^2+centers_y_second_ref(i)^2)-centers_x_first_ref(i)*(centers_x_first_ref(i)^2+centers_y_first_ref(i)^2));
%             vicr(i) = vci(i) - b0 - b1*centers_x_first_ref(i) - b2*centers_y_first_ref(i)- b3*centers_x_first_ref(i)*centers_y_first_ref(i) - b4*centers_y_first_ref(i)^2 + k1*(centers_y_second_ref(i)*(centers_x_second_ref(i)^2+centers_y_second_ref(i)^2)-centers_y_first_ref(i)*(centers_x_first_ref(i)^2+centers_y_first_ref(i)^2));
%             centers_corrected_ref(i,1) = centers_x_first_ref(i) + uicr(i);
%             centers_corrected_ref(i,2) = centers_y_first_ref(i) + vicr(i);
%         end
      area2 = calculate_area_from_coords(centers_corrected);

      centers_theoretical(:,1) = centers1_real(:,1) + xdisp .* ones(length(centers1_real),1);
      centers_theoretical(:,2) = centers1_real(:,2) + ones(length(centers1_real),1).*y;

      figure;
        imshow(img_trans)
        hold on
        plot(centers_corrected(:,1),centers_corrected(:,2),'b*')
        hold on
        plot(sortedcenters_real(:,1),sortedcenters_real(:,2),'r*')
        plot(centers_theoretical(:,1),centers_theoretical(:,2),'g*')
        
        plot(centers1_real(:,1),centers1_real(:,2),'y*')
        diff_cent_trans = centers_theoretical-sortedcenters_real;
        diff_mean_trans = sqrt(sum(diff_cent_trans.^2)/length(diff_cent_trans))
        
        diff_cent_init = centers1_real-centers_corrected;
        diff_corrected_init = sqrt(sum(diff_cent_init.^2)/length(diff_cent_init))
        
areadiff = ((area2-area1)/area1)*100;

%% Plot pixel change in pixels over time 
[r c ] = size(dtable);
figure;title('Change in pixels')
count= 0;
for i = 10:15
    count = count+1;
    diff = pdist2(dtable{:,i},[ones(r,1)* dtable{1,i}(1),ones(r,1)* dtable{1,i}(2)]);
    diffv = diag(diff);
    plot(dtable.dateTaken(1:end),diffv(1:end),'color', rand(1,3))
    hold on
    grid on
    leg{count} = num2str(i);
end
legend(leg);











































