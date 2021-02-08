function exp_info = parameter_setup()

prompt_exp_number = 'Enter the exp number for this experiment: ';
exp_number = input(prompt_exp_number);
prompt_brx_number = 'Enter the brx number for this experiment: ';
brx_number = input(prompt_brx_number);

% Experiment parameter set up
prompt_nimgtavg = 'Enter the number of images you want to average (press enter for default(5)):\n';
nimgtavg = input(prompt_nimgtavg);  % how many images do you want to average
    if isempty(nimgtavg)
     nimgtavg = 5; % default
    end
prompt_save_pic_timing = 'Enter how often you want images to be saved in minutes (press enter for default(3):\n';
save_pic_timing = minutes(input(prompt_save_pic_timing));  % how often you want images to be saved
    if isempty(save_pic_timing)
     save_pic_timing = minutes(3); % default
    end
prompt_timeInterval = 'Enter how often you want images to take pictures in seconds (press enter for default(0=fastest):\n';
timeInterval = seconds(input(prompt_timeInterval));  % how often does it take the images 
    if isempty(timeInterval)
     timeInterval = seconds(0); % default
    end    

% Timing Set Up %
% start date
values_start = 1;
while ~(length(values_start) == 6)
    prompt_startDate = 'Enter the start date for the experiment (yyyy,mm,dd,hh,mm,ss):\n';
    start_date = input(prompt_startDate,'s');
    values_start = split(start_date,',');
end
startDate = datetime(str2double(values_start{1}),str2double(values_start{2}),str2double(values_start{3}),str2double(values_start{4}),str2double(values_start{5}),str2double(values_start{6}));
% end date
values_end = 1;
while ~(length(values_end) == 6)
    prompt_endDate = 'Enter the end date for the experiment (yyyy,mm,dd,hh,mm,ss):\n';
    end_date = input(prompt_endDate,'s');
    values_end = split(end_date,',');
end
endDate = datetime(str2double(values_end{1}),str2double(values_end{2}),str2double(values_end{3}),str2double(values_end{4}),str2double(values_end{5}),str2double(values_end{6}));
 
prompt_sensorcap = 'Enter (1) if sensors are connected to the arduino; press Enter if they are not:\n';
sensorcap = input(prompt_sensorcap);  % 1 if sensors are connected, 0 if not 
    if isempty(sensorcap)
     sensorcap = 0; % default
    end 
    
prompt_nbeads = 'Enter the number of beads placed on the sclera (default is 5):\n';
nbeads = input(prompt_nbeads);  % how many beads were positioned on the eye 
    if isempty(nbeads)
     nbeads = 5; % default
    end 




% save exp info
exp_info = {exp_number brx_number nimgtavg save_pic_timing timeInterval startDate endDate sensorcap nbeads};

end

