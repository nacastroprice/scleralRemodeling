function coeff_vec = fit_distortion_model(coord_vec_1,coord_vec_2,coeff_num)
% FIT_DISTORTION_MODEL - Fits camera distortion model proposed by (Bing Pan  et al 2014 Meas. Sci. Technol. 25 025001)
%   INPUTS
%    coord_vec_1 = nx2 vector of point coordinates in the first position
%    coord_vec_2 = nx2 vector of point coordinates in the second position
%    coeff_num = number of coefficients to include in the model
%       6 = a0,a1,a2,b0,b1,b2
%       8 = a0,a1,a2,a4,b0,b1,b2,b3
%       10 = a0,a1,a2,a3,a4,b0,b1,b2,b3,b4
%       odd = coeff_num - 1 + k1
%   OUTPUTS 
%    coeff_vec = coeff_numx1 vector of coefficient values

% check if k1 coeff will be included in model
if mod(coeff_num,2)~=0
    act_coeff_num = coeff_num -1;
    lens = 1;
else
    act_coeff_num = coeff_num;
    lens = 0;
end

% check if coeff_num input is valid
cases = [2,6,8,10];
if ~ismember(act_coeff_num,cases)
    error('Error: Number of coefficient is not valid, the value must be one of these numbers [6,7,8,9,10,11]')
end
% check that both coordinate vector are of equal length
if ~(length(coord_vec_1)==length(coord_vec_2))
    error('Error: Coordinate vectors must be of the same length')
end
num_points_included = (length(coord_vec_1));
% uci is the displacement in x
uci  = coord_vec_2(:,1) - coord_vec_1(:,1);
% uvi is the displacement in y
vci  = coord_vec_2(:,2) - coord_vec_1(:,2);

if (sum(uci) == 0 && sum(vci) == 0)
    warning('Warning: Displacement between points is zero')
    coeff_vec = zeros(coeff_num,1);
else
    coeff_vec = zeros(coeff_num,1);% vector of coefficients
    x = zeros(num_points_included*2,coeff_num);% matrix
    d = [uci; vci];% vector of displacements

    switch act_coeff_num
         case 2
            if lens
                % Matrix setup
                x = [ones(num_points_included,1), zeros(num_points_included,1), (coord_vec_2(:,1).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,1).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2)); ... 
                zeros(num_points_included,1), ones(num_points_included,1), (coord_vec_2(:,2).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,2).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2))];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = 0;coeff_vec(3,1) = 0;coeff_vec(4,1) = 0;coeff_vec(5,1) = 0;coeff_vec(6,1) = psol(2);coeff_vec(7,1) = 0;coeff_vec(8,1) = 0;coeff_vec(9,1) = 0;coeff_vec(10,1) = 0;coeff_vec(11,1) = psol(3);
            else
                % Matrix setup
                x = [ones(num_points_included,1), zeros(num_points_included,1); ... 
                zeros(num_points_included,1), ones(num_points_included,1)];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = 0;coeff_vec(3,1) = 0;coeff_vec(4,1) = 0;coeff_vec(5,1) = 0;coeff_vec(6,1) = psol(2);coeff_vec(7,1) = 0;coeff_vec(8,1) = 0;coeff_vec(9,1) = 0;coeff_vec(10,1) = 0;coeff_vec(11,1) = 0;
            end
        case 6
            if lens
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), (coord_vec_2(:,1).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,1).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2)); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), (coord_vec_2(:,2).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,2).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2))];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = 0;coeff_vec(5,1) = 0;coeff_vec(6,1) = psol(4);coeff_vec(7,1) = psol(5);coeff_vec(8,1) = psol(6);coeff_vec(9,1) = 0;coeff_vec(10,1) = 0;coeff_vec(11,1) = psol(7);
            else
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2)];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = 0;coeff_vec(5,1) = 0;coeff_vec(6,1) = psol(4);coeff_vec(7,1) = psol(5);coeff_vec(8,1) = psol(6);coeff_vec(9,1) = 0;coeff_vec(10,1) = 0;coeff_vec(11,1) = 0;
            end
        case 8
            if lens
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), (coord_vec_2(:,1).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,1).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2)); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2), (coord_vec_2(:,2).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,2).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2))];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = 0;coeff_vec(5,1) = psol(4);coeff_vec(6,1) = psol(5);coeff_vec(7,1) = psol(6);coeff_vec(8,1) = psol(7);coeff_vec(9,1) = psol(8);coeff_vec(10,1) = 0;coeff_vec(11,1) = psol(9);
            else
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2)];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = 0;coeff_vec(5,1) = psol(4);coeff_vec(6,1) = psol(5);coeff_vec(7,1) = psol(6);coeff_vec(8,1) = psol(7);coeff_vec(9,1) = psol(8);coeff_vec(10,1) = 0;coeff_vec(11,1) = 0;

            end
            
        case 10
            if lens
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).^2, coord_vec_1(:,1).*coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), (coord_vec_2(:,1).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,1).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2)); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2), coord_vec_1(:,2).^2, (coord_vec_2(:,2).*(coord_vec_2(:,1).^2+coord_vec_2(:,2).^2)-coord_vec_1(:,2).*(coord_vec_1(:,1).^2+coord_vec_1(:,2).^2))];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = psol(4);coeff_vec(5,1) = psol(5);coeff_vec(6,1) = psol(6);coeff_vec(7,1) = psol(7);coeff_vec(8,1) = psol(8);coeff_vec(9,1) = psol(9);coeff_vec(10,1) = psol(10);coeff_vec(11,1) = psol(11);
            else
                % Matrix setup
                x = [ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).^2, coord_vec_1(:,1).*coord_vec_1(:,2), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1); ... 
                zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), zeros(num_points_included,1), ones(num_points_included,1), coord_vec_1(:,1), coord_vec_1(:,2), coord_vec_1(:,1).*coord_vec_1(:,2), coord_vec_1(:,2).^2];
                % Solve for the coefficients
                psol = (transpose(x) * x)\transpose(x)*d;
                coeff_vec(1,1) = psol(1);coeff_vec(2,1) = psol(2);coeff_vec(3,1) = psol(3);coeff_vec(4,1) = psol(4);coeff_vec(5,1) = psol(5);coeff_vec(6,1) = psol(6);coeff_vec(7,1) = psol(7);coeff_vec(8,1) = psol(8);coeff_vec(9,1) = psol(9);coeff_vec(10,1) = psol(10);coeff_vec(11,1) =0;
                

            end
    end
end
coeff_vec = coeff_vec';
end

