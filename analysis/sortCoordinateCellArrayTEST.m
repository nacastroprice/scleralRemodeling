% sortCoordinateCellArray TEST
% output_cell = sortCoordinateCellArray(ref,test,varargin)

% Test 1: ref.length==test.length no thresh
ref = {[1,2],[34,64],[5,9],[21,9]};
test = {[34,60],[0,2],[19,9],[5,8]};

output_t = {[0,2],[34,60],[5,8],[19,9]};
output_cell = sortCoordinateCellArray(ref,test);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 1 succesfully");
else
    disp("You FAILED Test 1");
end
% Test 2: ref.length==test.length, thresh = 6
ref = {[1,2],[34,64],[5,9],[21,9]};
test = {[34,60],[0,2],[19,9],[5,8]};
thresh = 6;

output_t = {[0,2],[34,60],[5,8],[19,9]};
output_cell = sortCoordinateCellArray(ref,test,thresh);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 2 succesfully");
else
    disp("You FAILED Test 2");
end

% Test 3: ref.length>test.length
ref = {[1,2],[34,64],[5,9],[21,9]};
test = {[34,60],[0,2],[19,9]};

output_t = {[0,2],[34,60],[],[19,9]};
output_cell = sortCoordinateCellArray(ref,test);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 3 succesfully");
else
    disp("You FAILED Test 3");
end
% Test 4: ref.length>test.length, thresh = 6
ref = {[1,2],[34,64],[5,9],[21,9]};
test = {[34,60],[0,2],[19,9]};
thresh = 6;

output_t = {[0,2],[34,60],[],[19,9]};
output_cell = sortCoordinateCellArray(ref,test,thresh);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 4 succesfully");
else
    disp("You FAILED Test 4");
end
% Test 5: ref.length<test.length, thresh=6
ref = {[1,2],[34,64],[21,9]};
test = {[34,60],[0,2],[19,9],[5,8]};
thresh = 6;

output_t = {[0,2],[34,60],[19,9]};
output_cell = sortCoordinateCellArray(ref,test,thresh);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 5 succesfully");
else
    disp("You FAILED Test 5");
end
% Test 6: ref.length<test.length
ref = {[1,2],[34,64],[21,9]};
test = {[34,60],[0,2],[19,9],[5,8]};


output_t = {[0,2],[34,60],[19,9]};
output_cell = sortCoordinateCellArray(ref,test);

err = 0;
for i = 1:length(ref)
    if~(output_cell{i}==output_t{i}) 
        err = err+1;
    end
end
if (err==0)
    disp("You PASSED Test 6 succesfully");
else
    disp("You FAILED Test 6");
end


