%% Script that substitutes markerArea column with new area calculated from desired points
% Requires that you run the first section of area_analysis_compare to load
% the desired tables into the workspace
% USE
% (1)- Run 1st section of post_main_analysis_from_tables with desired exp_number
% (2)- Change j to the desired brx_number 
% (3)- Change the first section of for loop depending on what points you
% want to include: example below

j = 2; % Toggle brx table
clear centers tri point1 point2 point3 area


for i = 1:1: height(tables{j})
    %%%%% EXAMPLE %%%%%
    % nbeads=5; Area wanted from centers 1,3,4
    % Only change the X on right side: .centerX
    % Add more lines if needed for more points
%     centers(1,:) = tables{j}.center1(i,:);
%     centers(2,:) = tables{j}.center3(i,:);
%     centers(3,:) = tables{j}.center4(i,:);
    %%%%% --------------- %%%%%
    
    %%%%% CHANGE THIS %%%%%
%     centers(1,:) = tables{j}.center17(i,:);
%     centers(2,:) = tables{j}.center18(i,:);
%     centers(3,:) = tables{j}.center19(i,:);
%     centers(4,:) = tables{j}.center20(i,:);
%     centers(5,:) = tables{j}.center5(i,:);
%     centers(6,:) = tables{j}.center6(i,:);
%     centers(7,:) = tables{j}.center7(i,:);
%     centers(8,:) = tables{j}.center8(i,:);
%     centers(9,:) = tables{j}.center9(i,:);
%     centers(10,:) = tables{j}.center10(i,:);
%     centers(11,:) = tables{j}.center11(i,:);
%     centers(12,:) = tables{j}.center12(i,:);



%     centers(1,:) = tables{j}.center21(i,:);
%     centers(2,:) = tables{j}.center22(i,:);
%     centers(3,:) = tables{j}.center23(i,:);
%     centers(4,:) = tables{j}.center24(i,:);
%     centers(5,:) = tables{j}.center25(i,:);
%     centers(6,:) = tables{j}.center26(i,:);
%     centers(7,:) = tables{j}.center27(i,:);
%     centers(8,:) = tables{j}.center28(i,:);
%     centers(9,:) = tables{j}.center29(i,:);
%     centers(10,:) = tables{j}.center30(i,:);

%     centers(1,:) = tables{j}.center10(i,:);
%     centers(2,:) = tables{j}.center11(i,:);
%     centers(3,:) = tables{j}.center12(i,:);

    centers(1,:) = tables{j}.center36(i,:);
    centers(2,:) = tables{j}.center37(i,:);
    centers(3,:) = tables{j}.center39(i,:);
    centers(4,:) = tables{j}.center40(i,:);
    centers(5,:) = tables{j}.center42(i,:);
    centers(6,:) = tables{j}.center43(i,:);
%     centers(7,:) = tables{j}.center19(i,:);
%     centers(8,:) = tables{j}.center20(i,:);


    % calculate area
    tri = delaunayTriangulation(centers);             
    area = zeros(length(tri.ConnectivityList),1);

    for h = 1:size(tri.ConnectivityList,1)
        % point coordinate (find out how the points are orderd in csv)
        point1 = centers(tri.ConnectivityList(h,1),:);
        point2 = centers(tri.ConnectivityList(h,2),:);
        point3 = centers(tri.ConnectivityList(h,3),:);
        area(h) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
    end

    markerArea = sum(area);
    tables{j}.markerAreaReal(i) = markerArea;              
end