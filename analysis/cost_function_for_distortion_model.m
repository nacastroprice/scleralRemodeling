% show me a picture with the points labelled so I know which points are
% reference and which are not
imgindex = 2;
img = imread(strcat(exp_dir,'\','images','\',dtable.imageName{imgindex}));
[r c] = size(dtable);
column_num_before_centers = 5;
index_of_second_pic = imgindex + 1;
c_index = c - column_num_before_centers;
    for t = 1:c_index
        centers1(t,:) = dtable{imgindex,column_num_before_centers+t};
        centers2(t,:) = dtable{index_of_second_pic,column_num_before_centers+t};
    end
    a = [1:c_index]'; b = num2str(a); labels = cellstr(b);
    figure;
    imshow(img)
    hold on
    h = labelpoints(centers1(:,1), centers1(:,2), b, 'N', 0.15);
    hold off
    savefig(strcat(exp_dir,'labeled_points.fig'));
% pic reference points from pictures
real_points = 10:17;
ref_points = 1:exp_info{end};
ref_points(real_points) = [];

% reference points
centers_x_first_ref = centers1(ref_points',1);
centers_y_first_ref = centers1(ref_points',2);
centers_x_second_ref = centers2(ref_points',1);
centers_y_second_ref = centers2(ref_points',2);
% real trackable points
centers_x_first_real = centers1(real_points',1);
centers_y_first_real = centers1(real_points',2);
centers_x_second_real = centers2(real_points',1);
centers_y_second_real = centers2(real_points',2);
    % do just three points
    num_points_included = length(ref_points);
    centers_x_first_ref(num_points_included+1:end) =  [];
    centers_y_first_ref(num_points_included+1:end) =  [];
    centers_x_second_ref(num_points_included+1:end) =  [];
    centers_y_second_ref(num_points_included+1:end) =  [];
% For 3 reference points
    % uci is the displacement in x = vector 1xpoints of x'-x
    uci  = centers_x_second_ref - centers_x_first_ref;
    % uvi is the displacement in y = vector 1xpoints of y'-y
    vci  = centers_y_second_ref - centers_y_first_ref;
        % uci is the displacement in x = vector 1xpoints of x'-x
    ui  = centers_x_second_real - centers_x_first_real;
    % uvi is the displacement in y = vector 1xpoints of y'-y
    vi  = centers_y_second_real - centers_y_first_real;
% create the x matrix and the d vector 
 d = zeros(num_points_included*2,1);
 x = zeros(num_points_included*2,11);
 for i = 1:num_points_included
     d(i) = uci(i);
     d(num_points_included+i) = vci(i);
     x(i,:) = [1 centers_x_first_ref(i) centers_y_first_ref(i) centers_x_first_ref(i)^2 centers_x_first_ref(i)*centers_y_first_ref(i) 0 0 0 0 0 (centers_x_second_ref(i)*(centers_x_second_ref(i)^2+centers_y_second_ref(i)^2)-centers_x_first_ref(i)*(centers_x_first_ref(i)^2+centers_y_first_ref(i)^2)) ];
     x(num_points_included+i,:) = [0 0 0 0 0 1 centers_x_first_ref(i) centers_y_first_ref(i) centers_x_first_ref(i)*centers_y_first_ref(i) centers_y_first_ref(i)^2 (centers_y_second_ref(i)*(centers_x_second_ref(i)^2+centers_y_second_ref(i)^2)-centers_y_first_ref(i)*(centers_x_first_ref(i)^2+centers_y_first_ref(i)^2))];
 end
 % (centers_x_second_ref(i)*(centers_x_second_ref(i)^2+centers_y_second_ref(i)^2)-centers_x_first_ref(i)*(centers_x_first_ref(i)^2+centers_y_first_ref(i)^2))
 psol = (transpose(x) * x)\transpose(x)*d;
 
 a0 = psol(1);
 a1 = psol(2);
 a2 = psol(3);
 a3 = psol(4);
 a4 = psol(5);
 b0 = psol(6);
 b1 = psol(7);
 b2 = psol(8);
 b3 = psol(9);
 b4 = psol(10);
 k1 = psol(11);

    % fu = a0 + a1*xci + a2yci;
    % fv = b0 + b1*xci +b2yci;
for i = 1:length(real_points)
    uic(i) = ui(i) - a1*centers_x_first_real(i) - a2*centers_y_first_real(i)- a3*centers_x_first_real(i)^2 - a4*centers_x_first_real(i)*centers_y_first_real(i) + k1*(centers_x_second_real(i)*(centers_x_second_real(i)^2+centers_y_second_real(i)^2)-centers_x_first_real(i)*(centers_x_first_real(i)^2+centers_y_first_real(i)^2));
    vic(i) = vi(i) - b1*centers_x_first_real(i) - b2*centers_y_first_real(i)- b3*centers_x_first_real(i)*centers_y_first_real(i) - b4*centers_y_first_real(i)^2 + k1*(centers_y_second_real(i)*(centers_x_second_real(i)^2+centers_y_second_real(i)^2)-centers_y_first_real(i)*(centers_x_first_real(i)^2+centers_y_first_real(i)^2));
    centers_corrected(i,1) = centers_x_first_real(i) + uic(i);
    centers_corrected(i,2) = centers_y_first_real(i) + vic(i);
end
centers_second = [centers_x_second_real centers_y_second_real];
centers_first = [centers_x_first_real centers_y_first_real];
centers_compare = [centers_second centers_first centers_corrected];
disp_compare = [uic;ui';vic;vi'];
% 
%                 figure;
%                 imshow(img)
%                 hold on
%                 plot(centers_x_second_real(:),centers_y_second_real(:),'r*')
%                 hold on
%                 plot(centers_corrected(:,1),centers_corrected(:,2),'b*')
%                 hold off

    % where p = (a0 a1 a2 b0 b1 b2)T; terms in the distortion model equation
    % where d = (uci uvi)T; = (u1 u2 u3 u4 u5 u6 v1 v2 v3 v4 v5 v6)T;
    % displacements of reference points
    
    
% create cost functio with equations above in the form sum((uci -fu)^2 +
% (vci - fv)^2) this is sum of least squares, will solve for the equation
% terms

