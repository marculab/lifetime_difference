function [ im_res ] = shift_FLIm( im_mov, delta_X, delta_Y )
%SHIFT_FLIM Summary of this function goes here
%   Detailed explanation goes here

for i=0:30
    if isnan(sum(im_mov(:,end)))
        im_mov = im_mov(:,1:end-1);
    end
end

im_res = im_mov;

start_X = 0;
delta_X = round(delta_X);
if delta_X ~= 0
    if delta_X > 0
        start_X = size(im_mov, 1) - delta_X;
    elseif delta_X < 0
        start_X = -delta_X;
    end
    im_res = [im_mov(:, start_X: end) im_mov(:, 1: start_X - 1)];
end

start_Y = 0;
delta_Y = round(delta_Y);
if delta_Y ~= 0
    if delta_Y < 0
        im_res = im_res(- delta_Y: end, :);
        post = zeros(- delta_Y, size(im_res, 2));
        im_res = [im_res; post];
    elseif delta_Y > 0
        pre = zeros(delta_Y, size(im_res, 2));
        im_res = im_res(1: end - delta_Y, :);
        im_res = [pre; im_res];
    end
end

end

