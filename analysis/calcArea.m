function [area] = calcArea(centers)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
             % calculate area
                tri = delaunayTriangulation(centers);             
                area = zeros(length(tri.ConnectivityList),1);

                for j = 1:size(tri.ConnectivityList,1)
                    % point coordinate (find out how the points are orderd in csv)
                    point1 = centers(tri.ConnectivityList(j,1),:);
                    point2 = centers(tri.ConnectivityList(j,2),:);
                    point3 = centers(tri.ConnectivityList(j,3),:);
                    area(j) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
                end

                markerArea = sum(area);
                area = markerArea;

end

