function cameraParams = cam_params_from_calib(fc,cc,kc,nx,ny)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

intrinsic_matrix = [fc(1) 0 0; 0 fc(2) 0; cc(1) cc(2) 1];
radial = [kc(1) kc(2)];
tangential = [kc(3) kc(4)];
% error = mean(err_std);
imgsize = [ny nx];

cameraParams = cameraParameters('IntrinsicMatrix',intrinsic_matrix,'ImageSize',imgsize,'RadialDistortion',radial,'TangentialDistortion',tangential); 

end

