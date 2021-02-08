clear;clc;
%//////////////////////////////////////////////////////////////////////
organ_cult_dir = "R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis";
% directory and file naming schema set up
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);


exp_folder = "\" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = organ_cult_dir + exp_folder + "\" + brx_folder + "\";

tablename_load = exp_dir + num2str(exp_number) + "_" + num2str(brx_number) + "_from_images_nocalib_imhistmatch.mat"; % this is the name of the data table
tablename_new =  exp_dir + num2str(exp_number) + "_" + num2str(brx_number) + "_from_images_nocalib_imhistmatch_table_corrected_dist_model_11_both_real_and_ref_compare_to_previous_correcting_normal_centerorigin_ogcode.mat";
load(tablename_load);
% preprocess the table
    % to delete rows with zero area
    toDelete = dtable.markerArea == 0;
    dtable(toDelete,:) = [];
% Image folder and Image file name set up
img_folder = exp_dir + "images";
fileSchema = strcat(num2str(exp_number),"_", num2str(brx_number)) ;

addpath(exp_dir) % add the exp directory to path to access the masks
mask_function = @bead_mask; % this mask is the one that was saved in the last step of the set_up script
globe_mask = imread(exp_dir + "circle.jpeg"); % this is the mask created in the second step of the set_up.m

% Retrieve all the exp parameters set during setup: exp_info = [exp_number brx_number nimgtavg save_pic_timing timeInterval startDate endDate sensorcap nbeads];
load(exp_dir + "exp_info.mat");
if (exp_number == exp_info{1}) && (brx_number == exp_info{2})
    nimgtavg = exp_info{3}; % how many images do you want to average
    save_pic_timing = exp_info{4};% how often you want images to be saved
    timeInterval = exp_info{5};% how often does it take the images
    startDate = exp_info{6};
    endDate = exp_info{7};
    sensorcap = exp_info{8};% 1 if sensors are connected, 0 if not 
    nbeads = exp_info{9};% how many beads were positioned on the eye 
else
    error('Incorrect experimental parameter') 
end


% Do first image to establish the coordinate order
index_good_data = 1; % when does the good experiment data start
img = imread(strcat(exp_dir,'\','images','\',dtable.imageName{index_good_data}));
[r c] = size(dtable);
column_num_before_centers = 7;

index_of_second_pic = index_good_data + 1;
c_index = c - column_num_before_centers; % this is the column where the center coordinates start

% This is to show the labeled beads in order to input real vs ref
        for t = 1:c_index
            centers12(t,:) = dtable{index_good_data,column_num_before_centers+t};
        end
        a = [1:c_index]'; b = num2str(a); labels = cellstr(b);
        figure;
        imshow(img)
        hold on
        h = labelpoints(centers12(:,1), centers12(:,2), b, 'N', 0.15);
        hold off
        savefig(strcat(exp_dir,'labeled_points.fig'));
    
% This will have to be made more robust in case the bead placements are different in the future    
prompt_real_point_first = 'Enter the bead numbers in the ROI e.j [1,4,5] ';
real_points_one = input(prompt_real_point_first);
% Data point interval for which to fit the distortion model
prompt_interval = 'Enter the data interval to use ';
interval= input(prompt_interval);
% Vectors indicating which points are the real vs ref
points = 1:nbeads;

real_points = real_points_one;

ref_points = 1:nbeads;
ref_points(real_points) = [];

% add a coefficient column
coeff_num = 11; % number of coefficients in the model
coefficients = cell(r,1);
coefficients{1,1} = zeros(coeff_num,1);
coefficients = cell2table(coefficients);
dtable = horzcat(dtable,coefficients);


% Shift origin to center
[img_rows,img_col] = size(img);

for i = 1:r
   for ii = 1:c_index 
      coordinates = dtable{i,column_num_before_centers+ii};
      coordinates(1) = coordinates(1) - img_col/2;
      coordinates(2) = coordinates(2) - img_rows/2;
      dtable{i,column_num_before_centers+ii} = coordinates;
   end
end

dtable_new = dtable;
save(tablename_new,'dtable_new');


%
    % Gets the first valid picture
     for t = 1:c_index
        centers_corrected(t,:) = dtable{index_good_data,column_num_before_centers+t};
     end
     % Calculate area based on those only and save to markerArea
     markerArea = calculate_area_from_coords(centers_corrected(points,:));
     dtable_new.markerArea(index_good_data) = markerArea;
     %Calculate rreal beads area
     markerAreaReal = calculate_area_from_coords(centers_corrected(real_points',:));
     dtable_new.markerAreaReal(index_good_data) = markerAreaReal;
     %Calculate ref beads area
     markerAreaRef = calculate_area_from_coords(centers_corrected(ref_points',:));
     dtable_new.markerAreaRef(index_good_data) = markerAreaRef;

 
  
% Main loop; Gets the coordinates from the relevant rows; Creates ref_point
% first and second vector which it then uses to fit the distortion model
% with; then corrects the real points and ref points with the fitted terms
% and updates the area columns 

for j = index_good_data+interval:interval:r
    
        % get the centers from the last entry and the current entry
         for t = 1:c_index
             % centers of the og table (not corrected)
            centers1(t,:) = dtable{j-interval,column_num_before_centers+t};
%             centers1(t,:) = dtable{index_good_data,column_num_before_centers+t};
            centers2(t,:) = dtable{j,column_num_before_centers+t};
             % centers of the new table (corrected)
%             centers1_c(t,:) = dtable_new{j-interval,column_num_before_centers+t};
% %             centers1_c(t,:) = dtable_new{index_good_data,column_num_before_centers+t};
%             centers2_c(t,:) = dtable_new{j,column_num_before_centers+t};
         end
         
         % Centers expanded of og table
             % Get reference points
            centers_x_first_ref = centers1(ref_points',1);
            centers_y_first_ref = centers1(ref_points',2);
            centers_x_second_ref = centers2(ref_points',1);
            centers_y_second_ref = centers2(ref_points',2);
             % Get real trackable points
            centers_x_first_real = centers1(real_points',1);
            centers_y_first_real = centers1(real_points',2);
            centers_x_second_real = centers2(real_points',1);
            centers_y_second_real = centers2(real_points',2);
         % Centers expanded of corrected table
%              % Get reference points
%             centers_x_first_ref_c = centers1_c(ref_points',1);
%             centers_y_first_ref_c = centers1_c(ref_points',2);
%             centers_x_second_ref_c = centers2_c(ref_points',1);
%             centers_y_second_ref_c = centers2_c(ref_points',2);
%              % Get real trackable points
%             centers_x_first_real_c = centers1_c(real_points',1);
%             centers_y_first_real_c = centers1_c(real_points',2);
%             centers_x_second_real_c = centers2_c(real_points',1);
%             centers_y_second_real_c = centers2_c(real_points',2);

        
            % fit with original points
            coeff_vec = fit_distortion_model([centers_x_first_ref,centers_y_first_ref],[centers_x_second_ref,centers_y_second_ref],coeff_num);
%           coeff_vec = fit_distortion_model([centers_x_first_real,centers_y_first_real],[centers_x_second_real,centers_y_second_real],coeff_num);

         % coefficients calculated from original points --> saved in corrected table
         dtable.coefficients{j} = coeff_vec;
         
         % correct real points from new table
             [centers_real_corrected,disp_vec] = correct_displacement_with_distortion_model_coeffs([centers_x_first_real,centers_y_first_real],[centers_x_second_real,centers_y_second_real],coeff_vec);  
             % Replace the old trackable points with the new ones
                for h = 1:length(centers_real_corrected)
                    dtable{j,column_num_before_centers+real_points(h)} = centers_real_corrected(h,:);
                end
             % Calculate area based on those only and save to markerArea
                markerArea = calculate_area_from_coords(centers_real_corrected);
                dtable.markerAreaReal(j) = markerArea;
                
          % correct ref points from new table
             [centers_ref_corrected,disp_vec_ref] = correct_displacement_with_distortion_model_coeffs([centers_x_first_ref,centers_y_first_ref],[centers_x_second_ref,centers_y_second_ref],coeff_vec);  
             % Replace the old trackable points with the new ones
                for h = 1:length(centers_ref_corrected)
                    dtable{j,column_num_before_centers+ref_points(h)} = centers_ref_corrected(h,:);
                end
             % Calculate area based on those only and save to markerArea
                markerArea = calculate_area_from_coords(centers_ref_corrected);
                dtable.markerAreaRef(j) = markerArea;

    % if interval >1 substitute NaN for all the points not corrected
    if(j~=index_good_data+1)&&(interval>1)
        dtable.markerArea(j-1:-1:j-interval+1) = NaN;
        dtable.markerAreaReal(j-1:-1:j-interval+1) = NaN;
        dtable.markerAreaRef(j-1:-1:j-interval+1) = NaN;
        psol = zeros(11,1);
        for u = j-1:-1:j-interval+1
         dtable.coefficients{u} = psol';
        end
    end
    disp(j)
end
% if interval >1 substitute NaN for all the points not corrected
if(interval>1)
  dtable.markerArea(r-interval-1:r) = NaN;
  dtable.markerAreaReal(j-1:-1:j-interval+1) = NaN;
  dtable.markerAreaRef(j-1:-1:j-interval+1) = NaN;
  for u = r-interval-1:r
     dtable_new.coefficients{u} = psol';
  end
end
save(tablename_new,'dtable');
