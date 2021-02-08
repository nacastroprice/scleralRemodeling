function img = med_filt_stack(imgin,it)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
for i = 1:it
   imgin = medfilt2(imgin);  
end
img = imgin;
end

