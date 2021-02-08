clear;clc;
%//////////////////////////////////////////////////////////////////////
organ_cult_dir = "R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture";
% directory and file naming schema set up
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);
prompt_calib = 'Calib?: ';
calib = input(prompt_calib);


exp_folder = "\" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = organ_cult_dir + exp_folder + "\" + brx_folder + "\";

tablename = exp_dir + num2str(exp_number) + "_" + num2str(brx_number) + "_from_pics_calib_imhistmatch_halfsampledtable.mat"; % this is the name of the data table (can be changed and will create a new table) 
img_folder = exp_dir + "images";
fileSchema = strcat(num2str(exp_number),"_", num2str(brx_number)) ;

addpath(exp_dir) % add the exp directory to path to access the masks
mask_function = @bead_mask; % this mask is the one that was saved in the last step of the set_up script
globe_mask = imread(exp_dir + "circle.jpeg"); % this is the mask created in the second step of the set_up.m
if calib
    load(exp_dir + "cameraParams.mat"); % cameraParameters for calibration
end
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
%
files = dir(img_folder + '\' + '*.tif');
% loop through these
[rows, columns] = size(files);
disp(rows)
% convert to a table for easier handling
dtable = struct2table(files);
folder = dtable.folder{1}; % might have to change this to a string of the directory of images
dtable= removevars(dtable,{'bytes','isdir','folder'});  % remove columns we don't want
dtable.Properties.VariableNames{'name'} = 'imageName';

% here we keep only half o fthe table
dtable = dtable(1:2:end,:);
[rows, columns] = size(dtable);

% Do first image to establish the coordinate order
index_good_data = 28; % when does the good experiment data start
img1 = imread(strcat(folder,'\',dtable.imageName{index_good_data}));
img1 = im2double(img1);
imgref = im2double(imread(strcat(folder,'\',dtable.imageName{index_good_data}))); 

if calib
    img1 = undistortImage(img1,cameraParams);
    imgref = undistortImage(imgref,cameraParams); 
end

% extract centers vector
[centers,bw,radii] = findbeadsbyposition(img1,mask_function,globe_mask); 
% Implementation of intensity fitting algorith for subpixel accuracy
% 11041
% crop_factor = 0.6;
% hshift = 0;
% vshift = 0;
% mx = "min";
% 
% for i =1:length(centers)
%     referenced = find_bead_center_intensity_fit(img1,centers(i,:),radii(i),"poly22",crop_factor,hshift,vshift,mx);
%     centers(i,:) = referenced;
% end
          
a = [1:length(centers)]'; b = num2str(a); labels = cellstr(b);
figure;
imshow(img1)
hold on
h = labelpoints(centers(:,1), centers(:,2), b, 'N', 0.15);
hold off
prompt = 'Enter 1 if the bead location are INCORRECT:';
x = input(prompt);
if (x==1)
    disp('There is a problematic bead location; please check the mask before continuing')
end

% Get which ones are real and ref
% This will have to be made more robust in case the bead placements are different in the future    
prompt_real_point_first = 'Enter the bead numbers in the ROI e.j [1,4,5] ';
real_points_one = input(prompt_real_point_first);

% Vectors indicating which points are the real vs ref
points = 1:nbeads;
real_points = real_points_one;
ref_points = 1:nbeads;
ref_points(real_points) = [];


% calculate area
markerArea = calculate_area_from_coords(centers);
dtable.markerArea(index_good_data) = markerArea;
markerAreaReal = calculate_area_from_coords(centers(real_points',:));
dtable.markerAreaReal(index_good_data) = markerAreaReal;
markerAreaRef = calculate_area_from_coords(centers(ref_points',:));
dtable.markerAreaRef(index_good_data) = markerAreaRef;
initwidth = width(dtable);

% save the data in a table;
centerpad = zeros(height(dtable),2);
for j = 1:nbeads
   namecells(j,1) = strcat('center',string(j)); % name the centers in the order which they came
   dtable = addvars(dtable,centerpad,'NewVariableNames',namecells(j,1)); 
end
for j = 1:nbeads
   dtable{index_good_data,initwidth+j} = centers(j,:); 
end
save(tablename,'dtable');
stable = dtable(index_good_data,:);
ctable = cell(rows,1);
ctable{1,1} = centers;
atable = zeros(rows,1);
etable = zeros(rows,1);
ftable = zeros(rows,1);

%% parfor that deals with all the rest of the analysis
centers = zeros(nbeads,2);
centersx = zeros(nbeads,1);
centersy = zeros(nbeads,1);
sorted_centers = zeros(nbeads,2);
area = zeros(nbeads,1);
% 
% [BW,maskedImage] = mask_function(img1);
% mask = BW;
   
    
parfor i = index_good_data+1:rows
    
    % Do first image to establish the coordinate order
    img = imread(strcat(folder,'\',dtable{i,1}{1}));
    img = im2double(img);

    if calib
        img = undistortImage(img,cameraParams);
    end
%       
%     img = localhistmatch(imgref,img);
%       img = adapthisteq(img,'NumTiles', [10,10],'clipLimit',0.005,'Distribution','rayleigh');
    img = imhistmatch(img,imgref);
    % extract centers vector
    [centers,bw,radii] = findbeadsbyposition(img,mask_function,globe_mask);
    disp(i)
    
% % %     % Implementation of intensity fitting algorith for subpixel accuracy
%     for k =1:length(centers)
%         referenced = find_bead_center_intensity_fit(img,centers(k,:),radii(k),"poly22",crop_factor,hshift,vshift,mx);
%         centers(k,:) = referenced;
%     end
 
%     % order the centers with the first and place them in table
%     % create vector of x coord and y coord
%         centersx = centers(:,1);
%         centersy = centers(:,2);
%         sortedcenters = zeros(nbeads,2);
% %         % Point are named after first pic: the closest points
% %         % to the last iteration get the name
%         for j = 1:nbeads
%             [dif, I] = min(sqrt((stable{1,initwidth+j}(1)-centersx).^2 + (stable{1,initwidth+j}(2)-centersy).^2));
%             sortedcenters(j,:) = centers(I,:);
%             centersx(I,:) = [];
%             centersy(I,:) = [];
%             centers(I,:)  = [];            
%         end

        sorted_centers = sortCoordinateVector(table2cell(stable(1,initwidth+1:initwidth+nbeads)),centers);
%         sorted_centers = sortCoordinateVector(ctable{i-1,1},centers);

% 
%         % calculate area
        markerArea = calculate_area_from_coords(sorted_centers);
        atable(i,1) = markerArea;
        markerAreaReal = calculate_area_from_coords(sorted_centers(real_points',:));
        etable(i,1) = markerAreaReal;
        markerAreaRef = calculate_area_from_coords(sorted_centers(ref_points',:));
        ftable(i,1) = markerAreaRef;


                
        ctable{i,1} = sorted_centers; 
   

end
atable(1:index_good_data) = [];
etable(1:index_good_data) = [];
ftable(1:index_good_data) = [];
ctable(1:index_good_data) = [];
dtable(1:index_good_data-1,:) = [];



dtable{2:end,4} = atable(1:end,1);
dtable{2:end,5} = etable(1:end,1);
dtable{2:end,6} = ftable(1:end,1);
[rows, columns] = size(dtable);


for i = 1:rows-1
    for j = 1:nbeads
       if ~(length(ctable{i,:}) == 0)
        dtable{i+1,initwidth+j} = ctable{i,:}(j,:); 
       end
    end
end

date_string = strcat(num2str(exp_number),num2str(brx_number));
for i = 1:1:rows
        % grab date
        
        gits = regexp(dtable.imageName{i},' ','split');
        str1 = regexp(gits{1},strcat(date_string,'_(\d*)'),'tokens'); % extract the date
        str2 = regexp(gits{2},'(\d*).tif','tokens'); % extract the date

%         strDigits = regexp(dtable.imageName{i},imageSchema,'tokens'); % extract the date
%         d = string(strDigits{1, 1}{1, 1}{1, 1});
%         f = string(strDigits{1, 1}{1, 1}{1, 2});
        timep = [str1{1}{1} ' ' str2{1}{1}];
        timep = datetime(timep,'InputFormat','ddMMyyyy HHmmss');
        dtable.dateTaken(i) = timep;
end

dtable= removevars(dtable,{'date','datenum'}); 
dtable = movevars(dtable,'dateTaken','Before','imageName');
pad = NaN(rows,1);
dtable.pressure(:,1) = pad;
dtable.flow(:,1) = pad;
dtable = movevars(dtable,'pressure','After','markerArea');
dtable = movevars(dtable,'flow','After','pressure');

save(tablename,'dtable');

