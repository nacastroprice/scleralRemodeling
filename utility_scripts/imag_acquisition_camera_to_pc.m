


imaqfind % look at available devices
info = imaqhwinfo('winvideo',1);
% brx_num
% camera object
cam = imaq.VideoDevice('winvideo',1);% Device ID has to do with the order they were plugged in!!!!
% Formats and resolutions
cam.VideoFormat = 'YUY2_2592x1944';
cam.ReturnedColorSpace= 'RGB';
cam.ReturnedDataType = 'double';
cam.ReadAllFrames = 'off';

% set camera properties
cam.DeviceProperties.BacklightCompensation = 'off';
cam.DeviceProperties.Brightness = 20;
cam.DeviceProperties.Contrast = 44;
cam.DeviceProperties.ExposureMode = 'manual';
cam.DeviceProperties.Exposure = -1;
cam.DeviceProperties.FrameRate = '2.0000';
cam.DeviceProperties.Gain = 0;
cam.DeviceProperties.Gamma = 100;
cam.DeviceProperties.Hue = 0;
cam.DeviceProperties.Saturation = 64;
cam.DeviceProperties.Sharpness = 3;
cam.DeviceProperties.WhiteBalanceMode = 'manual';
cam.DeviceProperties.WhiteBalance = 4501;


img = step(cam);
figure;
imshow(img)
