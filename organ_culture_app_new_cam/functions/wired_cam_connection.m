function [img,vid] = wired_cam_connection(brx_number)
%this function returns the camera object and an initial picture

%   It is crucial for BRX1 to be plugged in to the computer first so that
%   the numbers correspond with the ports
% imaqfind % look at available devices
% info = imaqhwinfo('winvideo',1);
% brx_num
% camera object
if(brx_number==1)
    brx_number = 2;
else
    brx_number = 1;
end
vid = videoinput('tisimaq_r2013_64', brx_number, 'Y800 (2448x2048)');
triggerconfig(vid, 'manual');
src = getselectedsource(vid);

vid.FramesPerTrigger = 1;
vid.TimeOut = 100;
% 
% src.GainAuto = 'Off';
% src.ToneMappingAuto = 'Off';
% src.ExposureAuto = 'Off';
% src.Exposure = 0.02;
% src.Gain = 20;
% src.Gamma = 102;
% src.Sharpness = 5;
% src.Contrast = 0;
% src.Brightness = 0;
% src.Denoise = 0;

src.AutoFunctionsROI = 'Disable';
src.ExposureAuto = 'Off';
src.ToneMappingAuto = 'Off';
src.Brightness = 181;
% src.Exposure = 0.25002;
src.Exposure = 0.05;
src.Denoise = 6;
src.Sharpness = 0;






start(vid);
trigger(vid);
img = getdata(vid);
stop(vid);


end

