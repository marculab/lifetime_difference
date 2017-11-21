% measure difference of image pairs and compute stats

% add path of registration code
addpath('registration code/FLIm_rot_reg');

% variables
diff_matrices = {};
diff_list = {};
diff_mean_pixel = [];
case_id = {};

show_figs = 1;
mean_w = 1; % smoothing window size
base_channel = 11; % use channel 2's image for registration

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
    diffs = []; % store difference value of each channel, as 1D array
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
    
        % process image pair
        run_1 = rawrun_1(:,:,channel);
        run_2 = rawrun_2(:,:,channel);
        run_1(run_1 < 1) = 1;
        run_1(run_1 > 12) = 12;
        run_2(run_2 < 1) = 1;
        run_2(run_2 > 12) = 12;
        
        % mask, filter values < 15, replace with 1
        snr_1 = rawrun_1(:,:,channel-4);
        snr_2 = rawrun_2(:,:,channel-4);
        run_1(snr_1<15) = 1;
        run_2(snr_2<15) = 1;
        
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
        
        % smooth using mean
        h = ones(mean_w, mean_w) / (mean_w * mean_w);
        run_1 = imfilter(run_1, h);
        run_2_res = imfilter(run_2_res, h);

        % subtract
        if (size(run_1,2) > size(run_2_res, 2))
            run_1 = run_1(:, 1:size(run_2_res, 2));
        else
            run_2_res = run_2_res(:, 1:size(run_1, 2));
        end
        diff_img = imabsdiff(run_1, run_2_res);

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
            %limit = [1, 12];
            if max(run_1(:)) > max(run_2(:))
                limit = [min(run_1(:)) max(run_1(:))];
            else
                limit = [min(run_2(:)) max(run_2(:))];
            end

            f = figure('visible', 'off');
            imagesc(run_1);
            colormap(jet), caxis(limit), colorbar;
            saveas(f, fullfile('./', current_folders(i).name, strcat('1-',int2str(channel)) ), 'jpeg');

            f = figure('visible', 'off');
            imagesc(run_2_res);
            colormap(jet), caxis(limit), colorbar;
            saveas(f, fullfile('./', current_folders(i).name, strcat('2-',int2str(channel)) ), 'jpeg');

            f = figure('visible', 'off');
            imagesc(diff_img);
            colormap(flipud(gray)), caxis(limit), colorbar;
            saveas(f, fullfile('./', current_folders(i).name, strcat('diff-',int2str(channel)) ), 'jpeg');
        end
    end
   % stats of this pair in all channel
    pair_mean = mean(diffs(:));
    pair_med = median(diffs(:));
    pair_std = std(diffs(:));
    pair_min = min(diffs(:));
    pair_max = max(diffs(:));
    fprintf(outfile, 'AllChannels, %f, %f, %f, %f, %f,\n', [pair_mean pair_med pair_std pair_min pair_max]');
    % exclude CH4
    diffs_1_3 = diffs(:, 1:3);
    pair_1_3_mean = mean(diffs_1_3(:));
    pair_1_3_med = median(diffs_1_3(:));
    pair_1_3_std = std(diffs_1_3(:));
    pair_1_3_min = min(diffs_1_3(:));
    pair_1_3_max = max(diffs_1_3(:));
    fprintf(outfile, 'AllExcludingCH4, %f, %f, %f, %f, %f,\n', [pair_1_3_mean pair_1_3_med pair_1_3_std pair_1_3_min pair_1_3_max]');
    fclose(outfile);
    % add to difflist
    diff_list = [diff_list, diffs];

end

% overall stats
diffs_all = [];
diffs_all_1_3 = [];
for i = 1 : size(diff_list, 2)
    t = cell2mat(diff_list(i));
    t_1_3 = t(:, 1:3);
    diffs_all = cat(1, diffs_all, t(:));
    diffs_all_1_3 = cat(1, diffs_all_1_3, t_1_3(:));
end
outfile = fopen(strcat('stats-allcases', '.csv'), 'wt');
fprintf(outfile, [',', 'mean,', 'median,', 'std,', 'min,', 'max,\n']);
fprintf(outfile, 'AllCases, %f, %f, %f, %f, %f,\n', [mean(diffs_all) median(diffs_all) std(diffs_all) min(diffs_all) max(diffs_all)]');
fprintf(outfile, 'AllExcludingCH4, %f, %f, %f, %f, %f,\n', [mean(diffs_all_1_3) median(diffs_all_1_3) std(diffs_all_1_3) min(diffs_all_1_3) max(diffs_all_1_3)]');
fclose(outfile);
save('diff_list.mat', 'diff_list');
save('case_id.mat', 'case_id');