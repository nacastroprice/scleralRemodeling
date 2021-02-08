function referenced = find_bead_center_intensity_fit(img,center,radi,crop_fact)
% find_bead_center_intensity_fit- Find the center of bead by fitting the
% intensity with a polynomial and solving for the minimum/max
%   The purpose of the function is to achieve subpixel resolution when
%   calculating the center coordinate of the bead fiducial marker

% The function will take a bead center coordinate and estimated radius(from
% findbeadsbyposition) and use that to crop a specific area around the
% previously estimated center. That subimage is fitted with a poynomial
% function and then the absolute minimum/maximum is solved for and the new
% coordinate is returned
% INPUT - 
% center = vector(1,2) of the estimated coordinate of the marker center
% radi = estimated radius of the marker center (regionprops(equivdiam/2))
% crop_factor = fraction of the radius wanted in the cropped image (because of ilumination it could be better to change what area of the marker is used for the fitting)
% img = current image
% OUTPUT -
% referenced = vector(1,2) containing x and y coordinates of the newly
% estimated center

crop_factor = crop_fact; % only half of the radius will be used (so intensity gradient is linear)
cent = ceil(center); % estimated coordinate of the center of marker
radi = ceil(radi* crop_factor); % radius used to crop the marker into the subimage

        % Crop the picture with a given radius:
        % Keep Top/left hand corner pixel for reference
         x0 = cent(1)-radi;
          y0 = cent(2)-radi;
          rect = [x0,y0,radi*2,radi*2];
          cropped = imcrop(img,rect);  
%           figure;imshow(cropped);
          
          xdata = 1:size(cropped,2);
          ydata = 1:size(cropped,1);
          zdata = cropped;
          
          [xgrid,ygrid] = meshgrid(1:size(cropped,2),1:size(cropped,1));
          xcenter = size(cropped,2)/2;
          ycenter = size(cropped,1)/2;
          radius = sqrt((xgrid-xcenter).^2+(ygrid-ycenter).^2);
          idx = radius>radi;
          cropped_mod = double(cropped);
          cropped_mod(idx) = NaN;

        % fit the surface
        % get minimum point
        % https://www.mathworks.com/matlabcentral/answers/310115-how-to-find-local-minima-of-a-function
%         toolboxFile = 'PolyfitnTools.mltbx';
%         installedToolbox = matlab.addons.toolbox.installToolbox(toolboxFile

        [xData, yData, zData] = prepareSurfaceData( xgrid, ygrid, cropped_mod ); % converts grid and image into 3 vectors needed for fitting
        
                          % Set up fittype and options.
                ft = fittype( 'a-b*((x-x0).^2 + (y-y0).^2)', 'independent', {'x', 'y'}, 'dependent', 'z' );
                opts = fitoptions( 'Method', 'NonlinearLeastSquares','Robust','Bisquare');
                opts.Display = 'Off';
                opts.StartPoint = [0.964888535199277 0.157613081677548 0.970592781760616 0.957166948242946];
                [sfit_out, gof] = fit( [xData, yData], zData, ft, opts );
                values = coeffvalues(sfit_out);

                xo = [xData(1),yData(1)];
                f = @(x,y)(values(1)-values(2)*((x-values(3)).^2 + (y-values(4)).^2));
                fun = @(x)f(x(1),x(2));
                [xmin,fval] = fminsearch(fun,xo);
                new_x = xmin(1);
                new_y = xmin(2);
                
%                 ft = fittype( 'poly22');
%                 opts = fitoptions( 'Method', 'LinearLeastSquares','Robust','Bisquare');
% %                 opts.Display = 'Off';
% %                 opts.StartPoint = [0.964888535199277 0.157613081677548 0.970592781760616 0.957166948242946];
%                 [sfit_out, gof] = fit( [xData, yData], zData, ft, opts );
%                 values = coeffvalues(sfit_out);
% 
%                 xo = [xData(1),yData(1)];
%                 f = @(x,y)(values(1) + values(2)*x + values(3)*y + values(4)*x^2 + values(5)*x*y + values(6)*y^2);
%                 fun = @(x)f(x(1),x(2));
%                 [xmin,fval] = fminsearch(fun,xo);
%                 new_x = xmin(1);
%                 new_y = xmin(2);
% 
%         P = polyfitn([xData, yData],zData,2);% Fits a polynomial to the data
%         Ps = polyn2sym(P); % converts to symbolic to be more easily used to solve for the minima/maxima
%         [x1,x2] = solve(gradient(Ps) == 0); % solve for the 0 derivative; will find many solutions
%         
%         cent_new = double([x1,x2]); % convert form sym to double
%         [zdatamin,loc] = min(zData); % find the location of the minimum intensity of the cropped image (where the center will be)
%        
%         % Iterate over the solutions found and get rid of the ones that are
%         % imaginary or have a negative component.
%         % Also, calculate the distance of the viable solutions to the
%         % initial estimate of the marker coordinate (the one that is closes will be the solution we are looking for)
%         for(i=1:size(cent_new,1))
%             if (~isreal(cent_new(i,1)) || ~isreal(cent_new(i,2))|| cent_new(i,1)<0 || cent_new(i,2)<0)
%                 log(i) = 1;
%                 cent_new(i,3) = 2000;
%             else
%                 cent_new(i,3) = sqrt((xData(loc)-cent_new(i,1)).^2 + (yData(loc)-cent_new(i,2)).^2);
% %               cent(i,3) = polyvaln(P,double([x1(i),x2(i)])); % get values at those locations
%                 log(i) = 0;
%             end
%         end
% %         cent(log',:) = []; % these are the candidates for minimum
%         % Locate the solution that is closest to the original estimation of
%         % the marker coordinate and keep that one
%         [diff,I] = min(cent_new(:,3));
%         new_x = cent_new(I,1);
%         new_y = cent_new(I,2);
        % tranform to original reference system using the left hand corner
        % of the cropped rectangle
        referenced_x = new_x - 1 + x0;
        referenced_y = new_y - 1 +y0;
        referenced(1,1) = referenced_x;
        referenced(1,2) = referenced_y;

        % replace previous coordinate
        
        % Plot the initial estimate and the new center coordinate on the
        % original image
%         figure;imshow(img)
%         hold on
%         plot(referenced_x,referenced_y,'r*');
%         plot(center(1),center(2),'go');
end

