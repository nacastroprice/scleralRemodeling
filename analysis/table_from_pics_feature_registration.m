clear;clc;
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
tablename = strcat(exp_dir,num2str(exp_number),"_",num2str(brx_number),"_registered_sim.mat"); % the name of the function to create the mask for the specific image set
% picture directory
files = dir(strcat(img_folder,'\*.tif'));
imageSchema = strcat(num2str(exp_number), num2str(brx_number),'_(\d* )(\d*).tif');
mkdir(strcat(exp_dir,"registered_images_sim"));
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
% loop through these
[rows, columns] = size(files);
disp(rows)
% convert to a table for easier handling
dtable = struct2table(files);
dtable= removevars(dtable,{'date','datenum','bytes','isdir'});  % remove columns we don't want
dtable.Properties.VariableNames{'name'} = 'imageName';

for i = 1:1:rows
            % grab date
        strDigits = regexp(dtable.imageName(i),imageSchema,'tokens'); % extract the date
        d = string(strDigits{1, 1}{1, 1}{1, 1});
        f = string(strDigits{1, 1}{1, 1}{1, 2});
        timep = d + f;
        timep = datetime(timep,'InputFormat','ddMMyyyy HHmmss');
        dtable.dateTaken(i) = timep;
        dtable.index(i) = i;
end

dtable = movevars(dtable,'dateTaken','Before','imageName');
% sort on dates, keep the sorting indices 
dtable = sortrows(dtable,1);

if isfile(strcat(exp_dir, "fixed_roi.mat"))
    load(strcat(exp_dir, "fixed_roi.mat")); % load roi for masking
else
    figure;imshowpair(imread(strcat(dtable.folder{i},'\',dtable.imageName{1})),imread(strcat(dtable.folder{i},'\',dtable.imageName{end})));
    fixed_roi = drawpolygon();
    save(strcat(exp_dir, "fixed_roi.mat"),'fixed_roi');

end
% vec1 = 1:156;
% vec2 = 157:30:8283;
% vec3 = 8283:rows;
% 
% vec = [vec1,vec2,vec3]';
% toDelete = ~ismember(dtable.index,vec);
% dtable(toDelete,:) = [];
% [rows, columns] = size(dtable);

for i = 1:50:rows
    disp(i)
     % get data
        % import image
        img = imread(strcat(dtable.folder{i},'\',dtable.imageName{i}));         
        if  (i == 1)  
                og_img = img;
                bw_mask = createMask(fixed_roi,img);         
                % extract centers vector
                [centers] = findbeadsbyposition(img,mask_function, bw_mask);
                % check that the centers calculated are the actual beads
                figure;imshow(img)
                hold on
                plot(centers(:,1),centers(:,2),'b*')
                hold off
                prompt = 'Enter 1 if the bead location are INCORRECT:';
                x = input(prompt);
                if (x==1)
                    disp('There is a problematic bead location; please check the mask before continuing')
                    break
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
                dtable.markerArea(i) = markerArea;
                initwidth = width(dtable);
                % save the data in a table; either existing or
                        centerpad = zeros(height(dtable),2);
                        for j = 1:nbeads
                           namecells(j,1) = strcat('center',string(j)); % name the centers in the order which they came
                           dtable = addvars(dtable,centerpad,'NewVariableNames',namecells(j,1)); 
                        end
                        for j = 1:nbeads
                           dtable{i,initwidth+j} = centers(j,:); 
                        end
                        save(tablename,'dtable');
        else
             
%            original = imread(strcat(dtable.folder{i},'\',dtable.imageName{i-1})); 
           originalc = og_img;
           distortedc = img;
           % registration
%            edgeThreshold = 0.6;
%             amount = 0.3;
% 
%             originalc = localcontrast(original, edgeThreshold, amount);
%             distortedc = localcontrast(distorted, edgeThreshold, amount);

            % set your rois 
            % figure;imshow(original);
            % o_roi = drawpolygon();
            % figure;imshow(distorted);
            % d_roi = drawpolygon();
%             o_roi = fixed_roi;
%             d_roi = fixed_roi;
            % find features
            ptsOriginal  = detectSURFFeatures(originalc);
            ptsDistorted = detectSURFFeatures(distortedc);
            % exclude features if not in ROI
            ftog = inROI(fixed_roi,double(ptsOriginal.Location(:,1)),double(ptsOriginal.Location(:,2)));
            ftdis = inROI(fixed_roi,double(ptsDistorted.Location(:,1)),double(ptsDistorted.Location(:,2)));

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

%             figure;
%             showMatchedFeatures(original,distorted,matchedOriginal,matchedDistorted);
%             title('Putatively matched points (including outliers)');

            % estimate transform
            [tform, inlierDistorted, inlierOriginal] = estimateGeometricTransform(...
                matchedDistorted, matchedOriginal, 'similarity');

%             figure;
%             showMatchedFeatures(originalc,distortedc,inlierOriginal,inlierDistorted);
%             title('Matching points (inliers only)');
%             legend('ptsOriginal','ptsDistorted');

            % compute inverse transformation matrix
            Tinv  = tform.invert.T;

            ss = Tinv(2,1);
            sc = Tinv(1,1);
            scaleRecovered = sqrt(ss*ss + sc*sc);
            thetaRecovered = atan2(ss,sc)*180/pi;
            % tranform the distorted
            outputView = imref2d(size(originalc));
            recovered  = imwarp(distortedc,tform,'OutputView',outputView);
%             figure, imshowpair(original,recovered);
           % save registered image
           imwrite(recovered,strcat(exp_dir,"registered_images_sim\r",dtable.imageName{i}));
           % find the beads
           [centers] = findbeadsbyposition(recovered,mask_function, bw_mask);
           % postprocess center coords found        
                if (length(centers)< nbeads)
                    centers = zeros(5,2);
                end
                
               if (dtable.markerArea(i-1)== 0)
                   h = 1;
                   while dtable.markerArea(i-h)== 0
                       h = h+1;
                   end
               else 
                   h = 1;
               end
           % create vector of x coord and y coord
                        centersx = centers(:,1);
                        centersy = centers(:,2);
                        sortedcenters = zeros(nbeads,2);
                        % Point are named after first pic: the closest points
                        % to the last iteration get the name
                        for j = 1:nbeads
                            [dif, I] = min(sqrt((dtable{i-h,initwidth+j}(1)-centersx).^2 + (dtable{i-h,initwidth+j}(2)-centersy).^2));
                            difm(j,1) = dif;
                            difm(j,2) = I;
                            sortedcenters(j,:) = centers(I,:);
                            centersx(I,:) = [];
                            centersy(I,:) = [];
                            centers(I,:)  = [];
                        end
        
                        % calculate area
                        tri = delaunayTriangulation(sortedcenters);             
                        area = zeros(length(tri.ConnectivityList),1);

                        for j = 1:size(tri.ConnectivityList,1)
                            % point coordinate (find out how the points are orderd in csv)
                            point1 = sortedcenters(tri.ConnectivityList(j,1),:);
                            point2 = sortedcenters(tri.ConnectivityList(j,2),:);
                            point3 = sortedcenters(tri.ConnectivityList(j,3),:);
                            area(j) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
                        end
                        markerArea = sum(area);
                        dtable.markerArea(i) = markerArea;
                        % save the data to the table
                        for j = 1:nbeads
                           dtable{i,initwidth+j} = sortedcenters(j,:); 
                        end
                        save(tablename,'dtable');
       end
end
pressurepad = zeros(height(dtable),1);
dtable = addvars(dtable,pressurepad,pressurepad,'NewVariableNames',{'pressure','flow'}); 
dtable= removevars(dtable,{'folder'});  % remove columns we don't want
dtable = movevars(dtable,'markerArea','After',width(dtable)-2);

save(tablename,'dtable');
%% convert to table and timetable
dtable = rmmissing(dtable);
save(tablename,'dtable');
