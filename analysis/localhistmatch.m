function imgout = localhistmatch(img1,img2)
%LCOALHISTMATCH Performs histmatch on subimages and stitches the image back
%together
%   INPUT-
%       img1 - reference image
%       img2 - image requiring matching
%       cornermat - Matrix of indices for creating the subimages

cornermat = [ 1,1,213,235;
              472,1,744,267;
              533,275,744,480;
              1,307,254,480;
              219,100,510,343
            ];
stichedimage = img2;

for i = 1:size(cornermat,1)
    
    % cropimage
    subimg1 = img1(cornermat(i,2): cornermat(i,4),cornermat(i,1): cornermat(i,3));
    subimg2 = img2(cornermat(i,2): cornermat(i,4),cornermat(i,1): cornermat(i,3));
    
    % match
    matchedsubimage2 = imhistmatch(subimg2,subimg1);
    
    % stitch
    stichedimage(cornermat(i,2): cornermat(i,4),cornermat(i,1): cornermat(i,3)) = matchedsubimage2;
    
    figure; imshowpair(subimg1,subimg2)
end

imgout = stichedimage;
end

