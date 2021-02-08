%% This section will run the main loop
clear,clc,imaqreset;
addpath('./functions')

% directory and file naming schema set up
prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);


exp_folder = "../" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = exp_folder + "/" + brx_folder + "/";

tablename = exp_dir + num2str(exp_number) + "_" + num2str(brx_number) + "bosetest.mat"; % this is the name of the data table (can be changed and will create a new table)
img_folder = exp_dir + "images";
fileSchema = strcat(num2str(exp_number), num2str(brx_number),"_") ;

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

% establishes a connection with the correct raspi and takes a picture
% [img,rpi,cam] = initial_brx_connection(brx_number);
if (brx_number == 1)
    [img,vid] = wired_cam_connection(brx_number);
else
    [img,vid] = wired_cam_connection2(brx_number);
end
% If the sensors are connected open serial connection with arduino
% % % % if (rpi)
% % % %     if (sensorcap)
% % % %         myserialdevice = serialdev(rpi,'/dev/ttyACM0',9600);
% % % %     end
% % % %     maxDist = 1.0e5;
% % % % end% max distance between newly found points and the ones from the previous iteration (not in use)

% Switch statement set up
[state,currentDate] = check_state(startDate,endDate); % This function checks the current date and time and sets the system state
% state = 0 ; current datetime > end date 
% state = 1 ; current datetime < startdate date
% state = 2 ; startdate date < current datetime < end date 

% variable set up
misseddata = 0;
picCount = 0;
extrac_pic_count = 0;
dateTaken = datetime("now");


% timing loop %
while state
    switch state
        case 1
            [state,currentDate] = check_state(startDate,endDate);
        case 2
            if ((datetime('now') - dateTaken) >= timeInterval)
                if ((datetime('now') - dateTaken) > timeInterval)
                    warning('The time in between captures is bigger than you configured')
                end
                tic
                % This takes/averages desired pictures and converts them to
                % grayscale
                extrac_pic_count = extrac_pic_count + 1;
                imgavg = 0;
                icount = 0;
%                     for jj = 1:1:nimgtavg
%                         img = double(snapshot(cam));
%                         start(vid);
%                         while (vid.FramesAvailable < 1)
%                             pause(0.001);
%                         end
                            start(vid);
                            trigger(vid);
                            disp(vid.FramesAcquired)
                            disp(vid.FramesAvailable)
                            memory
                            img = getdata(vid);
                            disp(vid.FramesAvailable)
                            stop(vid);
                            img = im2double(img);


%                         if jj == 1
%                             imgave = img;
%                         else
%                             imgave = imgave + img;
%                         end
%                     end
               %%%  img = rgb2gray(uint8(imgave/nimgtavg));
%                  img = rgb2gray(imgave/nimgtavg); %Grytz mod
%                     img = rgb2gray(img);
                 % bead finding function
                 [centers,bw,radii] = findbeadsbyposition(img,mask_function, globe_mask); 
                 
                 %     % Implementation of intensity fitting algorith for subpixel accuracy
%                 for k =1:length(centers)
%                     referenced = find_bead_center_intensity_fit(img,centers(k,:),radii(k),0.52);
%                     centers(k,:) = referenced;
%                 end

                % update last picture
                dateTaken = datetime('now');
                dateTaken.Format = 'ddMMyyyy HHmmss';
                % Retrieve sensor data
                if (sensorcap)
                    serialData = char(read(myserialdevice,36)); % second output is bit number
                    dataLines=splitlines(serialData);
                    values = split(dataLines{2},',');

                    if (length(values) > 1)        
                        pressure = str2double(values{1});
                        flow = str2double(values{2});
                    else
                       pressure = NaN;
                       flow = NaN;
                    end
                 else
                       pressure = NaN;
                       flow = NaN;  
                 end

                % Naming and saving the data %
                if (length(centers) >= nbeads)
                     disp('All beads detected')
                     picCount = picCount + 1;
                    if  (picCount == 1) 
                        % Naming and saving the image % 
                        imageName = strcat(fileSchema, string(dateTaken), '.tif');
                        Fullfilename = fullfile(img_folder, imageName);
                        imwrite(img, Fullfilename);%Grytz
                        last_pic_saved = datetime('now');
                        last_pic_saved.Format = 'ddMMyyyy HHmmss';
                               
                        % calculate area
                        markerArea = calculate_area_from_coords(centers);
                        
                        % save the data in a table; either existing or
                        % create
                        if exist(tablename)~= 0
                            % load table, append, save
                             load(tablename);

                                 if(sensorcap==1)
                                    Ntablerow = {dateTaken,imageName,markerArea,pressure,flow};
                                    tableconstants = 5;
                                 else
                                    Ntablerow = {dateTaken,imageName,markerArea};
                                    tableconstants = 3;
                                 end
                                 % append
                                for i = 1:nbeads
                                   Ntablerow = [Ntablerow centers(i,:)]; 
                                end
                             dtable = [dtable;Ntablerow];
                             save(tablename,'dtable');
                         else
                            % create table
                                 if(sensorcap==1)
                                    dtable = table(dateTaken,imageName,markerArea,pressure,flow);
                                 else
                                    dtable = table(dateTaken,imageName,markerArea);
                                 end
                             tableconstants = width(dtable);
                                for i = 1:nbeads
                                   namecells(i,1) = strcat('center',string(i)); % name the centers in the order which they came
                                   dtable = addvars(dtable,centers(i,:),'NewVariableNames',namecells(i,1)); 
                                end
                             % append row of data
                                 if(sensorcap==1)
                                    Ntablerow = {dateTaken,imageName,markerArea,pressure,flow};
                                 else
                                    Ntablerow = {dateTaken,imageName,markerArea};
                                 end
                                 % append                                 
                                for i = 1:nbeads
                                   Ntablerow = [Ntablerow centers(i,:)]; 
                                end
                             dtable = [dtable;Ntablerow];
                             % save table
                             save(tablename,'dtable');
                           end
                            % Plot live %
                            if (sensorcap)
                                first_data_point = dateTaken;

                                d = figure('Name', 'Live Data Feed', 'Position', [50 50 1800 600]);

                                A = subplot(1,3,1); % Area Subplot
                                axis
                                xlabel('time');
                                ylabel('Area in pixels squared');
                                A_Line = animatedline('Color',[1.000 0.000 0.000],'LineWidth',2); 

                                P = subplot(1,3,2); % Pressure Subplot
                                xlabel('time');
                                ylabel('Pressure (mmHg)');
                                P_Line = animatedline('Color',[0.812 0.573 0.000],'LineWidth',2); 

                                F = subplot(1,3,3); % Flow Subplot
                                xlabel('time');
                                ylabel('Flow (microl/min)');
                                F_Line = animatedline('Color',[0.812 0.4 0.000],'LineWidth',2);  

                                % An axis array for all data subplot's X axis to be updated and visualised as date time format 
                                All_Axes = [A P F]; % Add all subplot's X axes to this array
                            else
                                first_data_point = dateTaken;
                                d = figure('Name', 'Live Data Feed', 'Position', [50 50 1800 600]);
                                A = subplot(1,1,1); % Area Subplot
                                title(strcat("BRX",num2str(brx_number), "Area"))
                                axis
                                xlabel('time');
                                ylabel('Area in pixels squared');
                                A_Line = animatedline('Color',[1.000 0.000 0.000],'LineWidth',2);

                                % An axis array for all data subplot's X axis to be updated and visualised as date time format 
                                All_Axes = [A]; % Add all subplot's X axes to this array
                            end
                            
                                
                    else
                         picCount = height(dtable)+1;
                          if ((datetime('now') - last_pic_saved) >= save_pic_timing)
                                % Naming and saving the image % 
                                imageName = strcat(fileSchema, string(dateTaken), '.tif');
                                Fullfilename = fullfile(img_folder, imageName);
                                imwrite(img, Fullfilename); %Grytz
                                last_pic_saved = datetime('now');
                                format long
                                time_point = datenum(last_pic_saved);
                                last_pic_saved.Format = 'ddMMyyyy HHmmss';
                          end
                         
                        load(tablename);
                        % create vector of x coord and y coord
                        sortedcenters = zeros(nbeads,2);
                        
                         % handles order if previous entry is zero
                           if (dtable.markerArea(picCount-1) ==0)
                               h = 1;
                               while dtable.markerArea (picCount-h) ==0
                                   h = h+1;
                               end
                           else
                               h = 1;
                           end
                           
                        % Point are named after first pic: the closest points
                        % to the last iteration get the name
                        for x = 1:nbeads
                            centers_previous(x,1) = dtable{picCount-h,tableconstants+x}(1);
                            centers_previous(x,2) = dtable{picCount-h,tableconstants+x}(2);
                        end

                        idx = knnsearch(centers,centers_previous);
                        sortedcenters = centers(idx,:);

                        % calculate area
                        markerArea = calculate_area_from_coords(sortedcenters);

                        % load table, append, save
                         if(sensorcap==1)
                            Ntablerow = {dateTaken,imageName,markerArea,pressure,flow};
                         else
                            Ntablerow = {dateTaken,imageName,markerArea};
                         end
                        % append
                        for i = 1:nbeads
                           Ntablerow = [Ntablerow sortedcenters(i,:)]; 
                        end
                        dtable = [dtable;Ntablerow];
                        save(tablename,'dtable');
                    end
                    
                    if (sensorcap)
                        addpoints(A_Line,datenum(dateTaken),markerArea);
                        addpoints(P_Line,datenum(dateTaken),pressure);
                        addpoints(F_Line,datenum(dateTaken),flow);
                        drawnow
                        
                        N = 5;
                        A.XLim = datenum([first_data_point dateTaken+seconds(N)]);
                        P.XLim = datenum([first_data_point dateTaken+seconds(N)]);
                        F.XLim = datenum([first_data_point dateTaken+seconds(N)]);
                        
                        % ddMMyyyy HHmmss
                        dynamicDateTicks(All_Axes, 'HH:MM:SS') % Uses the array of all X axis
                        drawnow
                    else
                        addpoints(A_Line,datenum(dateTaken),markerArea);
                        drawnow
                        
                        N = 5;
                        A.XLim = datenum([first_data_point dateTaken+seconds(N)]);
                        
                        % ddMMyyyy HHmmss
                        dynamicDateTicks(All_Axes, 'HH:MM:SS') % Uses the array of all X axis
                        drawnow
                    end
                        
                        
                        
                else
                     disp('The beads are not being recognized; check the lighting; check for bubbles or anomalies in the chamber')
                    % if image processing fails save every picture to have more info on what
                    % went wrong
                    imageName = strcat(fileSchema, string(dateTaken), '.tif');
                    Fullfilename = fullfile(img_folder, imageName);
                    imwrite(img, Fullfilename); %Grytz
                    last_pic_saved = datetime('now');
                    format long
                    time_point = datenum(last_pic_saved);
                    last_pic_saved.Format = 'ddMMyyyy HHmmss';

            end
            end
            
            [state,currentDate] = check_state(startDate,endDate);
            
        case 3
        % create timetable from table
        ttable = table2timetable(dtable,'Row','dateTaken');
        state = 0;
    end
end





