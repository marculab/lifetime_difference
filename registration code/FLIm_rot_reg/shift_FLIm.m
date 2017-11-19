function [ im_res ] = shift_FLIm( im1, tform )
%SHIFT_FLIM Summary of this function goes here
%   Detailed explanation goes here

for i=0:30
    if isnan(sum(im1(:,end)))
        im1 = im1(:,1:end-1);
    end
end

shift_ang = mod(tform.T(3,1),size(im1,2));
im_res = [im1(1:end-1,end-round(shift_ang)+1:end),im1(2:end,:)];
im_res = im_res(:,1:size(im1,2));
im_res = [im_res;zeros(1,size(im1,2))];

len = size(im_res, 1);

shift_pb = round(tform.T(3,2));
if shift_pb>0
    pre = zeros( shift_pb, size(im_res,2) );
    im_res = [pre;im_res];
    im_res = im_res(1:len,:);
elseif shift_pb<0
    im_res = im_res(-shift_pb + 1:end,:);
    post = zeros( -shift_pb, size(im_res,2) );
    im_res = [im_res;post];
end

end

