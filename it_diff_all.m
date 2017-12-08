% measure difference of image pairs and compute stats

% add path of registration code
addpath('registration code/FLIm_rot_reg');

% variables
diff_matrices = {};
diff_list = {};
diff_mean_pixel = [];
case_id = {};

show_figs = 1;
mean_w = 3; % smoothing window size
base_channel = 11; % use channel 2's image for registration
erode_radius = 4;
        
% get list of input files (suppose each pair of .mat file is in a seperate folder)
current_folders = dir();
current_folders = current_folders([current_folders.isdir]);
current_folders = current_folders(3:end);

% go through each folder
for i = 1 : length(current_folders) % for each folder
    matfiles = dir(current_folders(i).name);
    matfiles = matfiles(~[matfiles.isdir]);
    img_pair = {};
    for j = 1 : length(matfiles)
        [filepath,name,ext] = fileparts(matfiles(j).name);
        if strcmp(ext, '.mat')
            fullpath = fullfile(matfiles(j).folder, matfiles(j).name);
            load(fullpath); 
            % suppose variable name is AllData
            img_pair = [img_pair; AllData];
        end
    end
    if size(img_pair) ~= 2
        continue;
    end
    
    % variables for this image pair
    diffs = {}; % store difference value of each channel, as 1D array
    outfile = fopen(strcat('stats-',current_folders(i).name, '.csv'), 'wt');
    fprintf(outfile, [',', 'mean,', 'median,', 'std,', 'min,', 'max,\n']);
    rawrun_1 = cell2mat(img_pair(1));
    rawrun_2 = cell2mat(img_pair(2));
    case_id = [case_id, current_folders(i).name];
    
    % use channel 2 to do registration
    im1 = cell2mat(img_pair(1));
    im1 = im1(:,:,base_channel);
    im2 = cell2mat(img_pair(2));
    im2 = im2(:,:,base_channel);
    im1(im1 < 1) = 1;
    im1(im1 > 12) = 12;
    im2(im2 < 1) = 1;
    im2(im2 > 12) = 12;
    tform = reg_FLIm(im2, im1);
    
    % choose the channel (3rd dimension of AllData)
    for channel = 9 : 12
        % plot 4 channels in one figure
        % process image pair
        run_1_orig = rawrun_1(:,:,channel);
        run_2_orig = rawrun_2(:,:,channel);
        run_1 = run_1_orig;
        run_2 = run_2_orig;
        % mask, filter values < 15
        snr_1 = rawrun_1(:,:,channel-4);
        snr_2 = rawrun_2(:,:,channel-4);
        
        % dilate / erode
        noise_mask_1 = (run_1_orig < 1) | (run_1_orig > 12) | (snr_1 < 15);
        noise_mask_2 = (run_2_orig < 1) | (run_2_orig > 12) | (snr_2 < 15);
        se = strel('disk',erode_radius);
        noise_mask_1 = imdilate(noise_mask_1, se);
        noise_mask_2 = imdilate(noise_mask_2, se);
        %noise_mask_1 = imerode(noise_mask_1, se);
        %noise_mask_2 = imerode(noise_mask_2, se);
        
        % displaye noise mask
        %figure, imagesc(noise_mask_1);
        %figure, imagesc(noise_mask_2);
        
        % apply mask to image
        run_1(noise_mask_1) = -1;
        run_2(noise_mask_2) = -1;
        
        % register run2 to run_1
        try
            %tform = reg_FLIm(run_2, run_1);
            run_2_res = shift_FLIm(run_2, tform);
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'\n%s',e.message);
            continue
        end

        % exclude lines where run_2_resare all 0s
        while sum(run_2_res(1, :)) == 0
            run_2_res = run_2_res(2 : end, :);
        end
        run_1 = run_1(size(run_1,1)-size(run_2_res,1)+1:end, :);

        while sum(run_2_res(end, :)) == 0
            run_2_res = run_2_res( 1: end-1, :);
        end
        run_1 = run_1(1:size(run_2_res,1), :);

        % filter nan
        while isnan(sum(run_1(:,end)))
            run_1 = run_1(:,1:end-1);
        end
        if size(run_1, 2) > size(run_2_res, 2)
            run_1 = run_1(:,1:size(run_2_res, 2));
        else
            run_2_res = run_2_res(:,1:size(run_1, 2));
        end
        run_1(run_1 == -1) = NaN;
        run_2_res(run_2_res == -1) = NaN;
        noise_filter = isnan(run_2_res) | isnan(run_1);
        run_1(isnan(run_1)) = NaN;
        run_2_res(isnan(run_2_res)) = NaN;
        
        % smooth using mean
        h = ones(mean_w, mean_w) / (mean_w * mean_w);
        run_1 = imfilter(run_1, h);
        run_2_res = imfilter(run_2_res, h);
        
        % filter noise
        run_1(noise_filter) = NaN;
        run_2_res(noise_filter) = NaN;

        % subtract
        if (size(run_1,2) > size(run_2_res, 2))
            run_1 = run_1(:, 1:size(run_2_res, 2));
        else
            run_2_res = run_2_res(:, 1:size(run_1, 2));
        end
        diff_img = imabsdiff(run_1, run_2_res);
        diff_img(isnan(diff_img)) = 0;

        % sum of difference
        %diff_matrices = [diff_matrices, diff_img];
        %diff_mean_pixel = [diff_mean_pixel, sum(diff_img(:)) / (size(diff_img, 1) * size(diff_img, 2))];
           
        % stats for this pair & channel
        diffs = [diffs, diff_img(:)];
        
        cur_min = min(diff_img(:));
        cur_max = max(diff_img(:));
        cur_mean = mean(diff_img(:));
        cur_med = median(diff_img(:));
        cur_std = std(diff_img(:));
        
        ch = channel-8;
        fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [ch cur_mean cur_med cur_std cur_min cur_max]');
        
        % show images
        if show_figs
            limit = [1, 12];
            %{
            if max(run_1(:)) > max(run_2(:))
                limit = [min(run_1(:)) max(run_1(:))];
            else
                limit = [min(run_2(:)) max(run_2(:))];
            end
            %}
            f = figure('visible', 'off');
            imagesc(run_1);
            colormap(jet), caxis(limit), colorbar;
            set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 7.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 7.5]);
            saveas(f, fullfile('./', current_folders(i).name, strcat('1-',int2str(channel)) ), 'jpeg');
            saveas(f, fullfile('./', current_folders(i).name, strcat('1-',int2str(channel), '.fig') ));
            f = figure('visible', 'off');
            imagesc(run_2_res);
            colormap(jet), caxis(limit), colorbar;
            set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 7.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 7.5]);
            saveas(f, fullfile('./', current_folders(i).name, strcat('2-',int2str(channel)) ), 'jpeg');
            saveas(f, fullfile('./', current_folders(i).name, strcat('2-',int2str(channel), '.fig') ));

            f = figure('visible', 'off');
            imagesc(diff_img);
            %colormap(flipud(gray)), caxis(limit), colorbar;
            c = ceil(max(diff_img(:)));
            if c < 2
                c = 2;
            end
            colormap(jet), caxis([0, c]), colorbar;
            set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 7.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 7.5]);
            saveas(f, fullfile('./', current_folders(i).name, strcat('diff-',int2str(channel)) ), 'jpeg');
            saveas(f, fullfile('./', current_folders(i).name, strcat('diff-',int2str(channel), '.fig') ));
        end
    end
   % stats of this pair in all channel
    diffs_1d = [];
    for i = 1 : size(diffs, 2)
        d = cell2mat(diffs(i));
        diffs_1d = cat(1, diffs_1d, d(:));
    end
    pair_mean = mean(diffs_1d(:));
    pair_med = median(diffs_1d(:));
    pair_std = std(diffs_1d(:));
    pair_min = min(diffs_1d(:));
    pair_max = max(diffs_1d(:));
    fprintf(outfile, 'AllChannels, %f, %f, %f, %f, %f,\n', [pair_mean pair_med pair_std pair_min pair_max]');
    % exclude CH4
    count_1_3 = size(cell2mat(diffs(1)), 1) + size(cell2mat(diffs(2)), 1) + size(cell2mat(diffs(3)), 1);
    diffs_1_3 = diffs_1d(1:count_1_3, :);
    pair_1_3_mean = mean(diffs_1_3);
    pair_1_3_med = median(diffs_1_3);
    pair_1_3_std = std(diffs_1_3);
    pair_1_3_min = min(diffs_1_3);
    pair_1_3_max = max(diffs_1_3);
    fprintf(outfile, 'AllExcludingCH4, %f, %f, %f, %f, %f,\n', [pair_1_3_mean pair_1_3_med pair_1_3_std pair_1_3_min pair_1_3_max]');
    fclose(outfile);
    % add to difflist
    diff_list = [diff_list, {diffs}];

end

% overall stats
diffs_all = [];
diffs_all_1_3 = [];
for i = 1 : size(diff_list, 2)
    diffs = diff_list{1, i};
    diffs_1d = [];
    for j = 1 : size(diffs, 2)
        d = cell2mat(diffs(j));
        diffs_1d = cat(1, diffs_1d, d(:));
    end
    count_1_3 = size(cell2mat(diffs(1)), 1) + size(cell2mat(diffs(2)), 1) + size(cell2mat(diffs(3)), 1);
    diffs_1_3 = diffs_1d(1:count_1_3, :);
    diffs_all = cat(1, diffs_all, diffs_1d(:));
    diffs_all_1_3 = cat(1, diffs_all_1_3, diffs_1_3(:));
end
outfile = fopen(strcat('stats-allcases', '.csv'), 'wt');
fprintf(outfile, [',', 'mean,', 'median,', 'std,', 'min,', 'max,\n']);
fprintf(outfile, 'AllCases, %f, %f, %f, %f, %f,\n', [mean(diffs_all) median(diffs_all) std(diffs_all) min(diffs_all) max(diffs_all)]');
fprintf(outfile, 'AllExcludingCH4, %f, %f, %f, %f, %f,\n', [mean(diffs_all_1_3) median(diffs_all_1_3) std(diffs_all_1_3) min(diffs_all_1_3) max(diffs_all_1_3)]');
fclose(outfile);
save('diff_list.mat', 'diff_list');
    save('case_id.mat', 'case_id');