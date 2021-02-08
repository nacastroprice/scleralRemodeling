% This script assumes that data tables and directories are in the standar
% way
% 2003 - 
%     1
%     start1 = 28;
%     end1 = end;
%     2
%     start2 = 21
%     end2 = end;
% 2004 - 
%     1
%     start1 = 51;
%     end1 = end;
%     2
%     start2 = 5
%     end2 = end;
% 2025 - 
%     1
%     start1 = 4517;for rate of change = 34685
%     end1 = 'end'; for rate of change = 93882
%     2
%     start2 = 4041;for rate of change = 34225
%     end2 = 'end';for rate of change = 93860
% 2037 - aggrecanase
%     1
%     start1 = 35;
%     end1 = 'end';
%     2
%     start2 = 1;
%     end2 = 'end';
% 19145 - 
%     1
%     start1 = 1;
%     end1 = 'end';
%     2
%     start2 = 24;
%     end2 = 'end';
%  19139 - 
%     1
%     start1 = 20;
%     end1 = 'end';
%     2
%     start2 = 8;
%     end2 = 'end';
%  19118 - 
%     1
%     start1 = 6;
%     end1 = 'end';
%     2
%     start2 = 24;
%     end2 = 'end';
% img_date(1) = "13122019 214540";
% img_date(2) = "13122019 213504";
%  2043 - 
%     1
%     start1 = 3821;for rate of change =29568
%     end1 = 'end';for rate of change =
%     2
%     start2 = 4583;for rate of change =36972
%     end2 = 'end';for rate of change =
% img_date(1) = "02062020 232815";
% img_date(2) = "02062020 133001";
    
    
clear,clc;
experimentnumber = input("What is the number for this experiment:\n",'s');
dir = 'R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\';
projectdir = strcat(dir,experimentnumber,'\');
% savedirect = "Z:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\analysis\demo\";
savedirect = projectdir;
imgdir = strcat(projectdir,'images\');


% Set up the data tables
tables1 = load(strcat(projectdir,'1\',experimentnumber,'_1_from_images_nocalib_imhistmatch.mat'));
tables1 = tables1.dtable;
tables2 = load(strcat(projectdir,'1\',experimentnumber,'_1_from_images_nocalib_imhistmatch_table_corrected_dist_model_11_both_real_and_ref_compare_to_previous_correcting_normal_centerorigin_ogcode.mat'));
tables2 = tables2.dtable;

% to delete rows with zero area
toDelete = tables1.markerArea == 0;
tables1(toDelete,:) = [];
toDelete = tables2.markerArea == 0;
tables2(toDelete,:) = [];
% eliminate duplicates
[C,ia,ic] = unique(tables1(:,1),'rows');
tables1 = tables1(ia,:);
[C,ia,ic] = unique(tables2(:,1),'rows');
tables2 = tables2(ia,:);

variables1 = tables1.Properties.VariableNames;
variables2 = tables2.Properties.VariableNames;

Index1 = find(contains(variables1,'pressure'));
Index2 = find(contains(variables2,'pressure'));

if isempty(Index1)
   pressure =  NaN(height(tables1),1);
   flow =  NaN(height(tables1),1);
   tables1 = addvars(tables1,pressure,'After','markerArea');
   tables1 = addvars(tables1,flow,'After','pressure');

end
if isempty(Index2)
   pressure =  NaN(height(tables2),1);
   flow =  NaN(height(tables2),1);
   tables2 = addvars(tables2,pressure,'After','markerArea');
   tables2 = addvars(tables2,flow,'After','pressure');
end




sensor_cap = 0;
% set up the table array
tables = {tables1;tables2};
clear C tables1 tables2 toDelete ia ic
%%
% This script will plot the area and change in area over time (what time stamp?) (%) of a given dataset
% from the brx and with linear regression extract the slope and RMS
%         img_date(1) = "22102020 171026";% change this to a date that appears in the wanted pictures name (pictures are saved less than data is collected)
%         img_date(2) = "22102020 171026";
% You will have to adjust the indices below to set the time range you want
% to analize
    % set up for dateTakens that you want to look at (default is the entire
    % range)-1
    startdd(1) = tables{1}.dateTaken(1); % first value
    endd(1) = tables{1}.dateTaken(end); % end dateTaken value 
    startdindex(1) = find(tables{1}{:,1}==startdd(1));
    enddindex(1) = find(tables{1}{:,1}==endd(1));
    % set up for dateTakens that you want to look at (default is the entire
    % range)- 2
    startdd(2) = tables{2}.dateTaken(1); % first value
    endd(2) = tables{2}.dateTaken(end); % end dateTaken value 
    startdindex(2) = find(tables{2}{:,1}==startdd(2));
    enddindex(2) = find(tables{2}{:,1}==endd(2));
    variableNames = cell(2);
    
    exp_start = startdd(1); % first value
    day_mark  = exp_start + hours(24);
    [~,ind(1)] = min(abs(datenum(tables{1}{:,1})-datenum(day_mark)));
    [~,ind(2)] = min(abs(datenum(tables{2}{:,1})-datenum(day_mark))); 
    

    
% %%
% for j=1:length(tables)
%    n_col = 7;
%    refbead = 17;
%    xtra_cols = 1;
%    c_index = width(tables{j})- n_col; % this is the number of beads
%    tables{j}.dateTaken.Format = 'ddMMyyyy HHmmss';
%    table = tables{j};
%    for i = 1:c_index - xtra_cols
%         dif = pdist2(table{startdindex(j):enddindex(j),refbead+n_col},table{startdindex(j):enddindex(j),n_col+i},'euclidean');
%         dist_vec(:,i) = diag(dif);              
%    end   
%     % eliminate zero column
%       dist_vec(:,refbead) = [];
%       % normalize to initial distance
%       figure;
%       for b = 1:c_index - 1 - xtra_cols
%         
% 
% %         dist_vec(:,b) = ((dist_vec(:,b)/dist_vec(1,b))-1)* 100;
%         dist_vec(:,b) = dist_vec(:,b) - dist_vec(1,b);
% 
%         plot((tables{j}.dateTaken(startdindex(j):enddindex(j)))', dist_vec(startdindex(j):enddindex(j),b),'color', rand(1,3),'LineWidth',2)
%         % Other aspects of graph
%         if b >= refbead
%             leg{b} = num2str(b+1);
%         else
%             leg{b} = num2str(b);
%         end
%         grid on
%         hold on
%         legend(leg)
%       end
%         title(strcat(experimentnumber,num2str(j),': bead distance relative to bead ',num2str(refbead)))
%         xlabel('Time')
%         ylabel('% change in pixel distance')
%         clear leg
% 
% end

%%
   
    
    figure;
for j = 1: length(tables)
    
    variableNames{j,1} = string(tables{j}.Properties.VariableNames);
    % This figure shows the area over the time of the experiment (then make a
    % function with input dateTakeninit,dateTakenfinal)
    plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.markerArea(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
    % Other aspects of graph
    grid on
    hold on
end
    title(strcat(experimentnumber,'affine: area over time'))
    xlabel('Time')
    ylabel('Area p^2')
    legend('1','2')
%     savefig(strcat(savedirect,experimentnumber,'-','area-over-time-from-pics.fig'));
    
% This figure is the percent change in area over the first viable day
 figure;
for j = 1: length(tables)
   
    tables{j}.areanorm = (tables{j}.markerArea/tables{j}.markerArea(ind(j))-1) * 100;
    % This figure shows the change in area over the time of the experiment
    dur_vector = tables{j}.dateTaken(startdindex(j):ind(j))' - tables{j}.dateTaken(startdindex(j));
    plot(tables{j}.dateTaken(startdindex(j):ind(j))',tables{j}.areanorm(startdindex(j):ind(j)),'color', rand(1,3),'LineWidth',2)
    % Other aspects of graph
    grid on
    hold on

end

    title(strcat(experimentnumber,'affine: Relative area change over the first viable 24h'))
    xlabel('Time')
    ylabel('Area change (% of p^2)')
    legend('1','2')
%     savefig(strcat(savedirect,experimentnumber,'-','area-change-over-first-day-from-pics.fig'));
    
    % This figure is the percent change in area over time
 figure;
for j = 1: length(tables)
    tables{j}.areanorm = (tables{j}.markerArea/tables{j}.markerArea(startdindex(j))-1) * 100;
    tables{j}.areanormreal = (tables{j}.markerAreaReal/tables{j}.markerAreaReal(startdindex(j))-1) * 100;
    tables{j}.areanormref = (tables{j}.markerAreaRef/tables{j}.markerAreaRef(startdindex(j))-1) * 100;

    % This figure shows the change in area over the time of the experiment
    %plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.areanorm(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
%     plot(tables{j}.dateTaken(1:end)',tables{j}.markerAreaReal(1:end),'color', rand(1,3),'LineWidth',2)
%     plot(tables{j}.dateTaken(ind(j):enddindex(j))',tables{j}.areanorm(ind(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)

    % Other aspects of graph
    grid on
    hold on

end
legend('3,5,8','3,5,8,2','3,5,8,2,1','3,5,8,2,1,10','3,5,8,2,1,10,17','3,5,8,2,1,10,17,25,26');
title(strcat(experimentnumber,' Reference Experiment in incubator: Relative area change over time (different bead areas)'))

% figure;
%     plot(tables{3}.dateTaken(startdindex(j):end)',tables{3}.areanorm(startdindex(j):end),'color', rand(1,3),'LineWidth',2)
%     hold on
%     plot(tables{2}.dateTaken(startdindex(j):enddindex(j))',tables{2}.areanorm(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
%     hold on
%     plot(tables{1}.dateTaken(startdindex(j):enddindex(j))',tables{1}.areanorm(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)


    title(strcat(experimentnumber,' Reference Experiment outside incubator: Relative area change over time'))
    xlabel('Time')
    ylabel('Area change (% of p^2)')
    legend('a0a1a2a3a4b0b1b2b3b4k1', 'a0a1a2a4b0b1b2b3k1','a0a1a2a3a4b0b1b2b3b4','a0a1a2a4b0b1b2b3','a0a1a2b0b1b2','beads of interest','reference beads')
    
    figure;
    title(strcat(experimentnumber,' Reference Experiment inside incubator: Relative area change over time'))

    %plot(tables{1}.dateTaken(isfinite(tables{1}.markerArea(1:end))),tables{1}.areanorm(isfinite(tables{1}.markerArea(1:end))))
    grid on
    hold on
    
    sindex = 1;
     plot(tables{1}.dateTaken(sindex:end)',tables{1}.areanormreal(sindex:end),'color', rand(1,3))
     plot(tables{1}.dateTaken(sindex:end)',tables{1}.areanormref(sindex:end),'color', rand(1,3))
     plot(tables{2}.dateTaken(sindex:end)',tables{2}.areanormreal(sindex:end),'color', rand(1,3))
     plot(tables{2}.dateTaken(sindex:end)',tables{2}.areanormref(sindex:end),'color', rand(1,3))
%      plot(tables{2}.dateTaken(1:end)',tables{2}.areanorm(1:end),'color', rand(1,3))
%      plot(tables{1}.dateTaken(1:end)',tables{1}.areanorm(1:end),'color', rand(1,3))
    legend('imhistmatch real','imhistmatch ref');%,'corrected real','corrected ref


     %plot(tables{1}.dateTaken(sindex:end)',tables{1}.markerAreaReal(sindex:end),'color', rand(1,3))

% 
%     plot(tables{3}.dateTaken(1:end),tables{3}.areanormreal(1:end))
%     plot(tables{3}.dateTaken(1:end),tables{3}.areanormref(1:end))
% % 
%     plot(tables{4}.dateTaken(1:end),tables{4}.areanormreal(1:end))
%     plot(tables{4}.dateTaken(1:end),tables{4}.areanormref(1:end))
%     plot(tables{5}.dateTaken(isfinite(tables{5}.markerArea(1:end))),tables{5}.areanormreal(isfinite(tables{5}.markerArea(1:end))))
%     plot(tables{5}.dateTaken(isfinite(tables{5}.markerArea(1:end))),tables{5}.areanormref(isfinite(tables{5}.markerArea(1:end))))
% % 
%     plot(tables{6}.dateTaken(isfinite(tables{6}.markerArea(1:end))),tables{6}.areanormreal(isfinite(tables{6}.markerArea(1:end))))
%     plot(tables{6}.dateTaken(isfinite(tables{6}.markerArea(1:end))),tables{6}.areanormref(isfinite(tables{6}.markerArea(1:end))))


    
%     dtable = tables{j};
%     for i = 1:height(dtable)
%         for j = 1:11
%             coeff(i,j) = dtable.coefficients{i,1}(j);
%             disp(i)
%         end
%     end
%     
%     toDel = coeff > 0.5;
%     coeff(toDel) = NaN;
%     toDel = coeff < -0.5;
%     coeff(toDel) = NaN;
%     
%     for i = 1:11
%         figure;
%         plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),i),'color', rand(1,3),'LineWidth',2)
%         ylim([-0.00005,0.00005])
%         title(num2str(i));
%     end
%     
%      figure;
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),1),'color', rand(1,3),'LineWidth',2)
%      figure;
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),2),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),3),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),4),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),5),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),6),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),7),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),8),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),9),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),10),'color', rand(1,3),'LineWidth',2)
%      plot(dtable.dateTaken(startdindex(j):enddindex(j))',coeff(startdindex(j):enddindex(j),11),'color', rand(1,3),'LineWidth',2)
%      grid on
% 
%     title(strcat(experimentnumber,'Model coefficients over time'))
%     xlabel('Time')
%     ylabel('coefficient value')
%     legend('a0','a1','a2','a3','a4','b0','b1','b2','b3','b4','k1')
% %     
% tables1.areanorm = (tables1.markerArea/tables1.markerArea(startdindex(j))-1) * 100;
% plot(tables1.dateTaken(startdindex(j):enddindex(j))',tables1.areanorm(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)


%legend('original','a-b*((x-x0).^2 + (y-y0).^2)','poly22','poly44','poly22 bisquare circle data','a-b*((x-x0).^2 + (y-y0).^2) bisquare circle data')
% 
tables1 = dtable;
toDelete = tables1.markerArea == 0;
tables1(toDelete,:) = [];
[C,ia,ic] = unique(tables1(:,1),'rows');
tables1 = tables1(ia,:);
n = 4;
tables{n} = tables1;
startdd(n) = tables{n}.dateTaken(1); % first value
endd(n) = tables{n}.dateTaken(end); % end dateTaken value 
startdindex(n) = find(tables{n}{:,1}==startdd(n));
enddindex(n) = find(tables{n}{:,1}==endd(n));

    savefig(strcat(savedirect,experimentnumber,'-','relative-area-change-over-time-comparison-of-intensity-fit-functions-11141.fig'))
%%% specific graph from edited data vectors %%%
% dataforanalysis_2025(:,2) = tables{2}.markerArea(ind(2):enddindex(2));
% dataanalysisnorm_2025(:,2) = ((dataforanalysis_2025(:,2)/dataforanalysis_2025(1,2))-1) * 100;
% figure;
% plot(dataanalysisnorm_2025(:,1),'color', rand(1,3),'LineWidth',2)
% hold on
% plot(dataanalysisnorm_2025(:,2),'color', rand(1,3),'LineWidth',2)
% save('dataanalysisnorm_2025.mat','dataanalysisnorm_2025')
%%% ------------------------------------------------------- %%%

% Fit linear regression and add results to graph
 j = 2; % toggle between 2 and 1
    diffdate = abs(startdindex(j) - enddindex(j))+1;
    dateRange = 1:diffdate;
    X = dateRange;
    y = tables{j}.areanorm(startdindex(j):enddindex(j))';
    ftav = fitlm(X,y);figure;plotAdded(ftav);title(strcat(experimentnumber,num2str(j),': Linear regression'))
    rms = ftav.RMSE;
    r_squared = ftav.Rsquared.Adjusted ;
    slope = ftav.Coefficients{2,1};
    slope_se = ftav.Coefficients{2,2};
    title('original')
%     to clear the outliers
idx = tables{1,1}.areanorm > 0.3;
tables{1,1}(idx,:) = [];

dtable = tables{1};
imgindex = 644;
img = imread(strcat(projectdir,'1\','images','\',dtable.imageName{imgindex}));
[r c] = size(dtable);
c_index = c-3-7;
    for t = 1:c_index
        centers(t,:) = dtable{imgindex,7+t};
    end
    a = [1:c_index]'; b = num2str(a); labels = cellstr(b);
    figure;
    imshow(img)
    hold on
    h = labelpoints(centers(:,1), centers(:,2), b, 'N', 0.15);
    hold off
    savefig(strcat(exp_dir,'labeled_points.fig'));

   %% This portion will extract information about each bead
   % Plot distance relative to other beads
   for j = 1:length(tables)
       n_col = 7;
       c_index = width(tables{j})- n_col; % this is the number of beads
       tables{j}.dateTaken.Format = 'ddMMyyyy HHmmss';
       table = tables{j};
           % Here we will plot an image with labelled points
% %                 img = imread(strcat(projectdir,num2str(j),'\','images','\',experimentnumber,num2str(j),'_',string(tables{j}.dateTaken(startdindex(j))),'.tif'));
%                 img = imread(strcat(projectdir,num2str(j),'\','images','\',experimentnumber,num2str(j),'_',img_date(j),'.tif'));
% 
%                 for t = 1:c_index
%                 centers(t,:) = table{startdindex(j),5+t};
%                 end
%                 a = [1:c_index]'; b = num2str(a); labels = cellstr(b);
%                 figure;
%                 imshow(img)
%                 hold on
%                 h = labelpoints(centers(:,1), centers(:,2), b, 'N', 0.15);
%                 hold off
%                 savefig(strcat(projectdir,experimentnumber,num2str(j),'_labeled_points.fig'));
       for i = 1:c_index
          for h = 1:c_index
%               dif = pdist2(table{startdindex(j):enddindex(j),5+h},table{startdindex(j):enddindex(j),5+i},'euclidean');
              dif = pdist2(table{startdindex(j):10:40000,5+h},table{startdindex(j):10:40000,5+i},'euclidean');

              dist_vec(:,h) = diag(dif);              
          end
              % eliminate zero column
              dist_vec(:,i) = [];
              % normalize to initial distance
              cla reset
              figure;
              for b = 1:c_index - 1
                dist_vec(:,b) = ((dist_vec(:,b)/dist_vec(1,b)) -1)* 100;
                plot((tables{j}.dateTaken(startdindex(j):enddindex(j)))', dist_vec(startdindex(j):enddindex(j),b),'color', rand(1,3),'LineWidth',2)
%                 plot((tables{j}.dateTaken(startdindex(j):10:40000))', dist_vec(:,b),'color', rand(1,3),'LineWidth',2)

                % Other aspects of graph
                grid on
                hold on
              end
                title(strcat(experimentnumber,num2str(j),': bead distance relative to bead ',num2str(i)))
                xlabel('Time')
                ylabel('% change in pixel distance')
                % legend creation
                c_index_vec = 1:c_index;
                idxc = c_index_vec == i;
                c_index_vec(idxc) = [];
                for g = 1:length(c_index_vec)
                    leg{g} = num2str(c_index_vec(g));
                end
                legend(leg);
                savefig(strcat(savedirect,experimentnumber,num2str(j),'-','rel_mvmt_beads_rel',num2str(i),'.fig'));
                clear dist_vec 
       end
   end
  
%% This section has not been tested; it is meant to graph area vs pressure
% Plot pressure over time
% get rid of outliers
% press_limit = [6 6];
% for j = 1: length(tables)
%     for i = 1:height(tables{j})
%         if tables{j}.pressure(i) < press_limit(j)
%            tables{j}.pressurecleaned(i) = NaN; 
%         else
%             tables{j}.pressurecleaned(i) = tables{j}.pressure(i);
%         end
%     end
% end
% 
%    figure;
% for j = 1: length(tables)
%     
%     tables{j}.pressurechange = (tables{j}.pressurecleaned/tables{j}.pressurecleaned(startdindex(j))-1) * 100;
% %     plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.pressurechange(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
%     plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.pressurecleaned(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
% 
%     % Other aspects of graph
%     grid on
%     hold on
% 
% end
% 
% %     title(strcat(imageName,': IOP change over time'))
% %     ylabel('IOP (% mmHg)')
%     title(strcat(imageName,': IOP over time'))
%     ylabel('IOP (mmHg)')
%     xlabel('Time')
%     legend('1','2')
%     savefig(strcat(savedirect,'-','pressure-over-time.fig'));
%     
% % plot pressure and area for each eye    
% for j = 1: length(tables)
%     if j == 1
%        eye = '1';
%     else
%        eye = '2';
%     end
%     figure;
%     plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.areanorm(startdindex(j):enddindex(j)).* 20,'color', rand(1,3),'LineWidth',2)
%     hold on
%     plot(tables{j}.dateTaken(startdindex(j):enddindex(j))',tables{j}.pressurechange(startdindex(j):enddindex(j)),'color', rand(1,3),'LineWidth',2)
%     hold off
%     title(strcat(imageName,eye,': IOP and markerArea change over time'))
%     ylabel('% change')
%     xlabel('Time')
%     legend('area','pressure')
% %     savefig(strcat(savedirect,eye,'-','pressure-area-change-over-time.fig'));
% end
   





