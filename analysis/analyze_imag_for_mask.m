clear;clc;
%//////////////////////////////////////////////////////////////////////
organ_cult_dir = "R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture";
% directory and file naming schema set up
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);


exp_folder = "\" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = organ_cult_dir + exp_folder + "\" + brx_folder + "\";

tablename = exp_dir + num2str(exp_number) + "_" + num2str(brx_number) + "retrieved_from_pics.mat"; % this is the name of the data table (can be changed and will create a new table)
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

% Do first image to establish the coordinate order
index_good_data = 1; % when does the good experiment data start
img = imread(strcat(folder,'\',dtable.imageName{index_good_data})); 
% extract centers vector
[centers,bw] = findbeadsbyposition(img,mask_function,globe_mask); 
           
figure;
imshow(img)
hold on
plot(centers(:,1),centers(:,2),'b*')
hold off
prompt = 'Enter 1 if the bead location are INCORRECT:';
x = input(prompt);
if (x==1)
    disp('There is a problematic bead location; please check the mask before continuing')
end


% calculate area
tri = delaunayTriangulation(centers);             
area = zeros(length(tri.ConnectivityList),1);

for j = 1:size(tri.ConnectivityList,1)
    % point coordinate (find out how the points are orderd in csv)
    point1 = centers(tri.ConnectivityList(j,1),:);
    point2 = centers(tri.ConnectivityList(j,2),:);
    point3 = centers(tri.ConnectivityList(j,3),:);
    area(j) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
end

markerArea = sum(area);
dtable.markerArea(1) = markerArea;
initwidth = width(dtable);
% save the data in a table; either existing or
centerpad = zeros(height(dtable),2);
for j = 1:nbeads
   namecells(j,1) = strcat('center',string(j)); % name the centers in the order which they came
   dtable = addvars(dtable,centerpad,'NewVariableNames',namecells(j,1)); 
end
for j = 1:nbeads
   dtable{1,initwidth+j} = centers(j,:); 
end
save(tablename,'dtable');
stable = dtable(1,:);
ctable = cell(rows,1);
atable = zeros(rows,1);
%% parfor that deals with all the rest of the analysis
centers = zeros(nbeads,2);
centersx = zeros(nbeads,1);
centersy = zeros(nbeads,1);
sortedcenters = zeros(nbeads,2);
area = zeros(nbeads,1);

[BW,maskedImage] = mask_function(img);
mask = BW;
count = 0;
   
    
for i = index_good_data:125
    
    % Do first image to establish the coordinate order
    img = imread(strcat(folder,'\',dtable{i,1}{1})); 
    % extract centers vector
    [centers,bw] = findbeadsbyposition(img,mask_function,globe_mask);
    count = count+1;
   i_masks{count} = bw;
    disp(i)
 
    if length(centers) == nbeads
%     % order the centers with the first and place them in table
%     % create vector of x coord and y coord
        centersx = centers(:,1);
        centersy = centers(:,2);
        sortedcenters = zeros(nbeads,2);
%         % Point are named after first pic: the closest points
%         % to the last iteration get the name
        for j = 1:nbeads
            [dif, I] = min(sqrt((stable{1,initwidth+j}(1)-centersx).^2 + (stable{1,initwidth+j}(2)-centersy).^2));
            sortedcenters(j,:) = centers(I,:);
            centersx(I,:) = [];
            centersy(I,:) = [];
            centers(I,:)  = [];
            
        end
% 
%         % calculate area
        tri = delaunayTriangulation(sortedcenters);             
        area = zeros(length(tri.ConnectivityList),1);
% 
        for j = 1:size(tri.ConnectivityList,1)
            % point coordinate (find out how the points are orderd in csv)
            point1 = sortedcenters(tri.ConnectivityList(j,1),:);
            point2 = sortedcenters(tri.ConnectivityList(j,2),:);
            point3 = sortedcenters(tri.ConnectivityList(j,3),:);
            area(j) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
        end
        markerArea = sum(area);
        atable(i,1) = markerArea;
        
        ctable{i,1} = sortedcenters; 
    end

end

dtable{2:end,4} = atable(2:end,1);

for i = 2:rows
    for j = 1:nbeads
       if ~(length(ctable{i,:}) == 0)
        dtable{i,initwidth+j} = ctable{i,:}(j,:); 
       end
    end
end

for i = 1:1:rows
        % grab date
        
        gits = regexp(dtable.imageName{i},' ','split');
        str1 = regexp(gits{1},'10221_(\d*)','tokens'); % extract the date
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

