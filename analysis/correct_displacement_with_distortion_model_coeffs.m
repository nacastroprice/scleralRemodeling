function [coord_vec_corrected,disp_vec] = correct_displacement_with_distortion_model_coeffs(coord_vec_1,coord_vec_2,coeff_vec)
% CORRECT_DISPLACEMENT_WITH_DISTORTION_MODEL - Corrects the displacement between two coordinate vectors with distortion model proposed by (Bing Pan  et al 2014 Meas. Sci. Technol. 25 025001) 
%   INPUT - 
%    coord_vec_1 = nx2 vector of point coordinates in the first position
%    coord_vec_2 = nx2 vector of point coordinates in the second position
%    coeff_vec = nx1 vector of model coefficients
%   OUTPUTS 
%    coord_vec_corrected = coord_vec_2 corrected by the model
%    disp_vec = corrected displacement vector (coord_vec_corrected-coord_vec_1)

if length(coeff_vec)~=11
   error('Error: coefficient vector should be of length 11')
elseif ~(length(coord_vec_1)==length(coord_vec_2))
   error('Error: Coordinate vectors must be of the same length')
end 

% uci is the displacement in x
ui  = coord_vec_2(:,1) - coord_vec_1(:,1);
% uvi is the displacement in y
vi  = coord_vec_2(:,2) - coord_vec_1(:,2);
% unpack coeff_vec
a0 = coeff_vec(1);a1 = coeff_vec(2);a2 = coeff_vec(3);a3 = coeff_vec(4);a4 = coeff_vec(5);b0 = coeff_vec(6);b1 = coeff_vec(7);b2 = coeff_vec(8);b3 = coeff_vec(9);b4 = coeff_vec(10);k1 = coeff_vec(11);
% Correct the initial points with the new displacements
for i = 1:length(coord_vec_1)
    uic(i) = ui(i) - a0 - a1*coord_vec_1(i,1) - a2*coord_vec_1(i,2)- a3*coord_vec_1(i,1)^2 - a4*coord_vec_1(i,1)*coord_vec_1(i,2) + k1*(coord_vec_2(i,1)*(coord_vec_2(i,1)^2+coord_vec_2(i,2)^2)-coord_vec_1(i,1)*(coord_vec_1(i,1)^2+coord_vec_1(i,2)^2));
    vic(i) = vi(i) - b0 - b1*coord_vec_1(i,1) - b2*coord_vec_1(i,2)- b3*coord_vec_1(i,1)*coord_vec_1(i,2) - b4*coord_vec_1(i,2)^2 + k1*(coord_vec_2(i,2)*(coord_vec_2(i,1)^2+coord_vec_2(i,2)^2)-coord_vec_1(i,2)*(coord_vec_1(i,1)^2+coord_vec_1(i,2)^2));
    coord_vec_corrected(i,1) = coord_vec_1(i,1) + uic(i);
    coord_vec_corrected(i,2) = coord_vec_1(i,2) + vic(i);
end
disp_vec = [uic vic];
end

