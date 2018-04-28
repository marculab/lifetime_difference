function [] = handle_ngregistration(im_mov, im_fixed, snr_1, snr_2)

% variables
snr_limit = 15;
close_radius = 3;
mean_w = 3;
% snr mask
noise_mask_1 = snr_1 < snr_limit;
noise_mask_2 = snr_2 < snr_limit;

se = strel('disk', close_radius);
SE = strel('rectangle', [close_radius,close_radius]);
closemask = imclose(noise_mask_1, se);
dilatemask = imdilate(closemask, SE);

im_mov(dilatemask) = 0;

closemask = imclose(noise_mask_2, se);
dilatemask = imdilate(closemask, SE);

im_fixed(dilatemask) = 0;

h = ones(mean_w, mean_w) / (mean_w * mean_w);

im_mov = imfilter(im_mov, h);
im_fixed = imfilter(im_fixed, h);

% save image 1 for cselect
%limit = [0, max(im_mov(:))];
%scale your image to the range [0,1] before saving.
im_mov = im_mov - min(im_mov(:));
im_mov = im_mov ./ max(im_mov(:));
savepath =  'im1_reg.jpg';
[J,~] = gray2ind(im_mov);
imwrite(J, jet, savepath);

% save image 2for cselect
im_fixed = im_fixed - min(im_fixed(:));
im_fixed = im_fixed ./ max(im_fixed(:));
savepath =  'im2_reg.jpg';
[J,~] = gray2ind(im_fixed);
imwrite(J, jet, savepath);


cpselect('im1_reg.jpg', 'im2_reg.jpg');

end