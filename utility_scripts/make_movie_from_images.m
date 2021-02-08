%% Make videos from images in directory
clear;clc;
% Change projectdir to directory where desired images live
projectdir = "Z:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\20171\1\images\";
videodir = "Z:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\videos\";
dinfo = dir( fullfile(projectdir, '*.tif') );
dinfo_t = struct2table(dinfo);
clear dinfo
[rows, columns] = size(dinfo_t);
dinfo_t= removevars(dinfo_t,{'folder','date','bytes','isdir','datenum'});

for i = 1:1:rows
    % grab date
    one = split(dinfo_t.name(i),{'_','.'});
    timep = datetime(one{2},'InputFormat','ddMMyyyy HHmmss');
    date_o(i) = [timep];
end

dinfo_t.date_o = date_o';

dinfo_t = sortrows(dinfo_t, 'date_o'); % sort the table by 'DOB'


name = dinfo_t.name(1);
exp_info = strtok(name, '_');
video_name = strcat(videodir,exp_info,'.avi');
video = VideoWriter(video_name); %create the video object
open(video); %open the file for writing
disp(rows)
for ii=1:10:rows
  cla reset % This resets the figure; important so memory doesn;t run out
  clf('reset')
  img = imread(strcat(projectdir,'\',dinfo_t.name{ii})); %read the next image
%   img = med_filt_stack(img,3); % Median filters image n times
  title_string = split(dinfo_t.name{ii},{'_','.'});


    f = gcf;
    imshow(img);
    title(title_string{2});
    set(f,'Position',[50 50 744 700])
    drawnow;
    set(f,"Visible",'off')
   

    F = getframe(f);
%     [X, Map] = frame2im(F);
    
  writeVideo(video,F); %write the image to file
  disp(ii)
end
% If picture frequency is exagerated skip multiple
% for ii=1:15:15000
%   cla reset
%   img = imread(strcat(projectdir,'\',dinfo_t.name{ii})); %read the next image
%   img = med_filt_stack(img,3);
% %   % plot the bead coordinates on pic
% %   centers(1,:) = dtablep.center1(ii,:);
% %   centers(2,:) = dtablep.center2(ii,:);
% %   centers(3,:) = dtablep.center3(ii,:);
% %   centers(4,:) = dtablep.center4(ii,:);
% %   
% title_string = split(dinfo_t.name{ii},{'_','.'});
% 
%   figure('visible','off')
%     imshow(img);
%     title(title_string{2});
% 
%     F = getframe(gcf);
% %     [X, Map] = frame2im(F);
%     
%   writeVideo(video,F); %write the image to file
%   disp(ii)
% end
close(video); %close the file




