function total_area = calculate_area_from_coords(coords)
% Calculate_area_from_coords using delaunay triangulation

tri = delaunayTriangulation(coords);             
area = zeros(length(tri.ConnectivityList),1);

for j = 1:size(tri.ConnectivityList,1)
    % Iterate over the triangles found and calculate the area within each
    % one
    point1 = coords(tri.ConnectivityList(j,1),:);
    point2 = coords(tri.ConnectivityList(j,2),:);
    point3 = coords(tri.ConnectivityList(j,3),:);
    area(j) = 1/2*abs(det([point1(1),point1(2),1;point2(1),point2(2),1;point3(1),point3(2),1])); %column vector 1x ntriangles
end

total_area = sum(area);
end

