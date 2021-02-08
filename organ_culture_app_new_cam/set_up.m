%% Takes a picture that is then used to create:
%  (1) A mask of the eye globe;
%  (2) An intensity mask for the beads using the 

clear,clc;
addpath('./functions')

% Experiment parameter set up
exp_info = parameter_setup();


exp_number = exp_info{1};
brx_number = exp_info{2};

% create exp folder
exp_folder = "../" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = exp_folder + "/" + brx_folder + "/";
img_folder = exp_dir + "images";
status = mkdir(img_folder); % write error handling code for this

save(exp_dir + "exp_info.mat","exp_info");

% Creates the connection to the raspi; sets desired parameters in the
% camera; outputs a picture
if (brx_number == 1)
    [img,vid] = wired_cam_connection(brx_number);
else
    [img,vid] = wired_cam_connection2(brx_number);
end

% img = rgb2gray(img);
img = im2double(img);
disp('click along the edge of the eye to create the ROI')
globe_mask1 = roipoly(img);
globe_mask = ~globe_mask;
% globe_mask = ~globe_mask;
img(~globe_mask) = 0;
imwrite(globe_mask, exp_dir + "circle.jpeg");

% This section prompts the user to create the bead mask

    disp("(1) Open the image segmenter app;choose 'img' from workspace") 
    disp("(2) Adjust the intensity threshold so that the beads are clearly visible")
    disp("(3) Export as function and save with this name in the experiment directory as 'bead_mask'")
    disp("After the mask is saved copy and paste this command into the command window: centers = verify_bead_location(img,exp_number,brx_number,globe_mask);")
