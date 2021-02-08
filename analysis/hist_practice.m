% histogram fun/practice

img = imread('R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\121120\2\images\1211202_11122020 172322.tif');

% show pic
figure;
imshow(img);
title('original image')
% compute histogram
    % assert that it is grayscale
    [rows,columns,channels] = size(img);
    if (channels>1)
        warning('this is not a grayscale image')
    end
    
    hist_vec = zeros(1,256);

    for i=1:rows
       for ii = 1:columns
          int_value = img(i,ii);
          hist_vec(int_value +1) = hist_vec(int_value+1) + 1;
       end
    end

    % plot histogram
    figure;
    bar(0:255,hist_vec)
    title('histogram')
    
 % normalize histogram
 norm_hist = hist_vec/(rows*columns);
 % plot norm histogram
 figure;
 bar(0:255,norm_hist)
 title('normalized histogram')
 
 % CDF
  cdf = zeros(1,256);
  for i = length(norm_hist):-1:1
     cdf(i) = sum(norm_hist(1:i)); 
  end
  % plot norm histogram
     figure;
     bar(0:255,cdf)
     title('cdf')
  % equalize
  nobins= 255;
  eq_norm = round(cdf .* nobins);
    % plot norm histogram
     figure;
     bar(0:255,eq_norm)
     title('eq_norm')
  
     nimg = zeros(rows,columns);
      for i=1:rows
       for ii = 1:columns
          int_value = img(i,ii);
          nimg(i,ii) = eq_norm(int_value+1);
       end
      end
      
      nimg = uint8(nimg);
      imshow(nimg)
      
      % histogram matching
      
      simg = imread('R:\funded_projects\Grytz-R01EY026588-Scleral_remodeling_in_myopia\organ_culture\121120\2\images\1211202_14122020 065743.tif');
      
      cdf2 = (cumsum(imhist(simg)) / numel(simg))';
      
      snimg = zeros(rows,columns);
      cdfmap = zeros(1,256);
      for i=1:length(cdf)
            [~,ind] = min(abs(cdf(i)-cdf2));
            cdfmap(i) = ind-1;
      end
      
      out = uint8(cdfmap(double(img)+1));
      
      
      
      
      %% Break image into chunks, equalize and stitch back together
      
      % (c,r)
      % top-left = (1,1) - (213,235)
      % top-right = (472,start) - (end,267)
      % bottom-right = (533,275) - (end,end)
      % bottom-left = (start,307) -(254,end)
      % middle = (219,113) - (510,343)
      
      
      % cut
        sub_img1 = img(1:235,1:215);
        sub_img2 = simg(1:235,1:215);

        imshowpair(sub_img1, sub_img2,'montage')
      % match
      sub_img1_c = imhistmatch(sub_img2,sub_img1);
    
        imshowpair(sub_img2,sub_img1_c,'montage')
      % stitch
      simgmatched = simg;
      simgmatched(1:235,1:215) = sub_img1_c;
        
      imshow(simgmatched)
  