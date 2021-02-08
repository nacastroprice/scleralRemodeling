function output_cell = sortCoordinateCellArray(ref,test,varargin)
% sortCoodinateCellArray = function takes in two cell arrays of coordinates
% and ouputs a cell array of the size of ref containing coordinates of test
% ordered by minimum distance
%   INPUTS
%       ref = cell array of point coordinates (reference)
%       test = cell array of point coordinates (to sort according to distance to ref)
%       thresh = threshold distance that points must be under (optional)
%           DEFAULT only handles if test is smaller than ref
%   OUTPUT
%       output_cell = cell array of the size of ref containing coordinates of test
%           ordered by minimum distance

% Convert cell arrays to matrices
% ref
    ref_vec = cell2mat(ref)';
    f = ref_vec(2:2:end);
    ref_vec(2:2:end) = [];
    ref_vec(:,2) = f;
    clear f;
% test
    test_vec = cell2mat(test)';
    f = test_vec(2:2:end);
    test_vec(2:2:end) = [];
    test_vec(:,2) = f;
    clear f;
% Compute distance matrix
D = pdist2(ref_vec,test_vec);
[r,c] = size(D);

output_vec = cell(1,length(ref));

if(nargin>2)
    % With a threshold distance we use the following logic
    thresh = varargin{1};
    underThresh = D <= thresh;
    for i=1:r
        for ii=1:c
            if (underThresh(i,ii)==1)
                 % Populate new cell array with the closest points under the 
                 output_vec{1,i}(1) = test_vec(ii,1);
                 output_vec{1,i}(2) = test_vec(ii,2);
                 underThresh(i,ii) = 0;
            end
        end
    end
    
else
    
    % Create vector from D matrix
    nlength = int64(r*c);
    n_dist_vector = zeros(1,nlength);
    const = 0;
    i1 = 1;
    i2 = c;
    for i=1:r
        i1 = i1+const;
        i2 = i2+const;
        n_dist_vector(1,i1:i2) = D(i,1:c);
        const = c;
    end

    % if the new reference vector is bigger than the test vector sort the
    % points in columns based on the ref point that it is closest to
    for j=1:length(ref_vec)
        % grab the index of the minimum element in the vector 
       [diff,I] = min(n_dist_vector);
       % transform the index into a matrix indec
       if (I<=c)
           %output the a vector the size of reference with zeros if diff
           row = 1;
           column = I;
       else
           row = ceil(I/c);
           column = I-((row-1)*c);
       end
       % Populate new cell array with the closest points under the 
         output_vec{1,row}(1) = test_vec(column,1);
         output_vec{1,row}(2) = test_vec(column,2);
         n_dist_vector(I) = NaN;
    end
end
    output_cell = output_vec;
end

