function centers = verify_bead_location(img,exp_number,brx_number,globe_mask)
% This displays the detected bead centers in the current experiment to make sure they are correct; 
% If they are not the parameters in the function findbeadsbyposition might have to be adjusted:
% To adjust the parameters open the function; set a breakpoint where indicated  

% globe_mask = num2str(exp_number)+ "_circle_image_"+ num2str(brx_number) +".jpeg";
% mask_function = "mask_for_" + num2str(exp_number) + "_brx" + num2str(brx_number);

exp_folder = "../" + num2str(exp_number);
brx_folder = num2str(brx_number);
exp_dir = exp_folder + "/" + brx_folder + "/";

addpath(exp_dir) % add the exp directory to path to access the masks

[centers] = findbeadsbyposition(img,@bead_mask,globe_mask);

figure;
    imshow(img)
    hold on
    plot(centers(:,1),centers(:,2),'b*')
    hold off
    prompt = 'Enter 1 if the bead location are INCORRECT:';
    x = input(prompt);
    if (x==1)
        disp('There is a problematic bead location; please check the mask before continuing')
    end
    
end

