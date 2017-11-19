function [ tform ] = reg_FLIm( im1, im2 )
%REG_FLIM Summary of this function goes here
%   Detailed explanation goes here

im1 = ( im1/max(im1(:))*255);
im2 = ( im2/max(im2(:))*255);
 
[optimizer, metric] = imregconfig('multimodal');
movingRegistered = imregister(uint8(im1/max(im1(:))*255), uint8(im2/max(im2(:))*255), 'affine', optimizer, metric);

for i=0:30
    if isnan(sum(im1(:,end)))
        im1 = im1(:,1:end-1);
    end
    
    if isnan(sum(im2(:,end)))
        im2 = im2(:,1:end-1);
    end
end
im1 = uint8(im1); im2 = uint8(im2);
tform = imregtform([im1,im1], im2, 'translation', optimizer, metric);
end

