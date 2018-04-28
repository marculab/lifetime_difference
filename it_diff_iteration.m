% register the lifetime data of 2 runs
% measure pixel wise difference of lifetime pairs and compute stats
% save the difference figures and mat to file system
% saved mat files will be used by stats_figures_all.m

% add path of registration code
addpath('registration code/FLIm_rot_reg');

% get list of input files (suppose each pair of .mat file is in a seperate folder)
data_path = '../All Artery Data';
disp(['root folder: ', data_path]);
all_folder_flag = input(['Process all the cases in ', data_path, ' ? (''y''/''n'')']);
if ~strcmp(all_folder_flag, 'y')
    % folder_to_process = '1698_2';
    folder_to_process = input('Input the folder to process: ');
    disp(['The folder to process: ', folder_to_process]);
    disp(['Will check if ',  folder_to_process, ' is in ', data_path]);
end
current_folders = dir(data_path);
current_folders = current_folders([current_folders.isdir]);
current_folders = current_folders(3:end);

% variables
diff_matrices = {};
diff_list = {};
diff_mean_pixel = [];
case_id = {};

show_figs = 1;
mean_w = 3; % smoothing window size
base_channel = 9; % use channel 2's image for registration
erode_radius = 0;
close_radius = 3;
channel_range = [13 14 15 16];
snr_range = [5 6 7 8];
disp('Register using channel 2''data.');
disp(['Morphology  transformation radius: ', num2str(close_radius)]);
disp(['Mean smoothing window size: ', num2str(mean_w)])

cp_flag = 0; % control point selection, set 1 to be used by cases when the registration code works bad
cp_select = 1; % control point selection

%noise_limit_min = 1;
%noise_limit_max = 12;
noise_limit_min = 0;
noise_limit_max = 12;
disp(['SNR range: [', num2str(noise_limit_min), ', ',num2str(noise_limit_max), ']'])

% go through each folder
for folder_ind = 1 : length(current_folders) % for each folder
    cur_folder = current_folders(folder_ind);
    matfiles = dir(fullfile(data_path, cur_folder.name));
    if ~strcmp(cur_folder.name, folder_to_process)
        continue
    end
    matfiles = matfiles(~[matfiles.isdir]);
    img_pair = {};
    run_code = {};
    for j = 1 : length(matfiles)
        [filepath,name,ext] = fileparts(matfiles(j).name);
        w = strsplit(matfiles(j).name,'_');
        if strcmp(ext, '.mat') && size(w, 2) >= 5
            fullpath = fullfile(matfiles(j).folder, matfiles(j).name);
            load(fullpath); 
            % suppose variable name is AllData
            img_pair = [img_pair; AllData];
            % get file name for storing
            run_code = [run_code; strcat(w(1,2), '_', w(1,3))];
        end
    end
    if size(img_pair) ~= 2
        continue;
    end
    fprintf('%s\n', cur_folder.name);
    % variables for this image pair
    diffs = {}; % store difference value of each channel, as 1D array
    diff_imgs = {}; % store difference value of each channel, as image matrix
    outfile = fopen(fullfile(data_path, cur_folder.name, strcat('stats', '.csv')), 'wt');
    fprintf(outfile, [',', 'mean,', 'median,', 'std,', 'min,', 'max,\n']);
    rawrun_1 = cell2mat(img_pair(1));
    rawrun_2 = cell2mat(img_pair(2));
    case_id = [case_id, cur_folder.name];
    
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
    snr_1 = rawrun_1(:,:,base_channel - 4);
    snr_2 = rawrun_2(:,:,base_channel - 4);
    
    if cp_flag
        if cp_select
            handle_ngregistration(im2, im1, snr_2, snr_1);
        end
        dist = fixedPoints - movingPoints;
        translate_dist = [mean(dist(:, 1)) mean(dist(:, 2))];
    end
    
    % cell array to save figures
    run_1_origs = {};
    run_2_origs = {};
    run_1_figs = {};
    run_2_figs = {};
    run_1_dilates = {};
    run_2_dilates = {};
    diff_figs = {};
    % choose the channel (3rd dimension of AllData)
    for i = 1 : 4
        % plot 4 channels in one figure
        % process image pair
        channel_ind = i;
        run_1_orig = rawrun_1(:,:,channel_range(i));
        run_2_orig = rawrun_2(:,:,channel_range(i));
        run_1 = run_1_orig;
        run_2 = run_2_orig;
        % mask, filter values < 15
        snr_1 = rawrun_1(:,:,snr_range(i));
        snr_2 = rawrun_2(:,:,snr_range(i));
        
        % dilate / erode
        noise_mask_1 = (run_1_orig < noise_limit_min) | (run_1_orig > noise_limit_max) | (snr_1 < 15);
        noise_mask_2 = (run_2_orig < noise_limit_min) | (run_2_orig > noise_limit_max) | (snr_2 < 15);
        %se = strel('disk',erode_radius);
        %noise_mask_1 = imdilate(noise_mask_1, se);
        %noise_mask_2 = imdilate(noise_mask_2, se);
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
            if cp_flag
                run_2_res = shift_ngregistration(run_2, translate_dist(1), translate_dist(2));
            end
        catch e
            fprintf(1,'The identifier was:\n%s',e.identifier);
            fprintf(1,'\n%s',e.message);
            continue
        end

        % exclude lines where run_2_resare all 0s
        while sum(run_2_res(1, :)) == 0
            run_2_res = run_2_res(2 : end, :);
        end
        if size(run_1,1) > size(run_2_res,1)
            run_1 = run_1(size(run_1,1)-size(run_2_res,1)+1:end, :);
        elseif size(run_1,1) < size(run_2_res,1)
            run_2_res = run_2_res(size(run_2_res,1)-size(run_1,1)+1:end, :);
        end

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
        %run_1 = imfilter(run_1, h);
        %run_2_res = imfilter(run_2_res, h);
        %run_1(noise_filter) = 1;
        %run_2_res(noise_filter) = 1;
        I1 = imfilter(run_1, h);
        I2 = imfilter(run_2_res, h);
        falsemask = (isnan(I1) - isnan(run_1));
        I1(falsemask==1) = run_1(falsemask==1);
        falsemask = (isnan(I2) - isnan(run_2_res));
        I2(falsemask==1) = run_2_res(falsemask==1);
        run_1 = I1;
        run_2_res = I2;

        % filter noise
        run_1(noise_filter) = NaN;
        run_2_res(noise_filter) = NaN;

        % dilate / erode / close
        se = strel('disk', close_radius);
        closemask = imclose(noise_filter, se);
        SE = strel('rectangle', [close_radius,close_radius]);
        dilatemask = imdilate(closemask, SE);
        run_1_dilate = run_1;
        run_1_dilate(dilatemask) = NaN;
        run_2_res_dilate = run_2_res;
        run_2_res_dilate(dilatemask) = NaN;
        
        %imshow(dilatemask),figure, imshow(closemask); figure, imshow(noise_filter);
        %imagesc(run_1), figure, imagesc(run_2_res), figure, imagesc(I1), figure, imagesc(I2);
        % subtract
        if (size(run_1,2) > size(run_2_res, 2))
            run_1 = run_1(:, 1:size(run_2_res, 2));
        else
            run_2_res = run_2_res(:, 1:size(run_1, 2));
        end
        
        if (size(run_1_dilate,2) > size(run_2_res_dilate, 2))
            run_1_dilate = run_1_dilate(:, 1:size(run_2_res_dilate, 2));
        else
            run_2_res_dilate = run_2_res_dilate(:, 1:size(run_1_dilate, 2));
        end
        
        %diff_img = imabsdiff(run_1, run_2_res);
        %diff_img(isnan(diff_img)) = 0;
        diff_img = imabsdiff(run_1_dilate, run_2_res_dilate);
        diff_img(isnan(diff_img)) = 0;

        % sum of difference
        %diff_matrices = [diff_matrices, diff_img];
        %diff_mean_pixel = [diff_mean_pixel, sum(diff_img(:)) / (size(diff_img, 1) * size(diff_img, 2))];
           
        % stats for this pair & channel
        diffs = [diffs, diff_img(:)];
        
        run_1_origs = [run_1_origs; run_1_orig];
        run_2_origs = [run_2_origs; run_2_orig];
        run_1_figs = [run_1_figs; run_1];
        run_2_figs = [run_2_figs; run_2_res];
        diff_figs = [diff_figs; diff_img];
        run_1_dilates = [run_1_dilates; run_1_dilate];
        run_2_dilates = [run_2_dilates; run_2_res_dilate];
        
        cur_min = min(diff_img(:));
        cur_max = max(diff_img(:));
        cur_mean = mean(diff_img(:));
        cur_med = median(diff_img(:));
        cur_std = std(diff_img(:));
        
        ch = channel_ind;
        fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [ch cur_mean cur_med cur_std cur_min cur_max]');
        
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
    
    % save diff to mat file
    % diff_list = [diff_list, {diffs}];
    for i = 1 : size(diff_figs, 1)
        imdiff = cell2mat(diff_figs(i));
        savepath = fullfile(data_path, cur_folder.name, char(strcat('diffs_ch', string(i),'.mat')));
        save(savepath, 'imdiff');
    end
    
    % save figs to mat
    savepath = char( fullfile(data_path, cur_folder.name, strcat(run_code(1) ,'_figs.mat')) );
    save(savepath, 'run_1_figs');
    savepath = char( fullfile(data_path, cur_folder.name, strcat(run_code(2) ,'_figs.mat')) );
    save(savepath, 'run_2_figs');
    savepath = char( fullfile(data_path, cur_folder.name, strcat(run_code(1) ,'_dilate_figs.mat')) );
    save(savepath, 'run_1_dilates');
    savepath = char( fullfile(data_path, cur_folder.name, strcat(run_code(2) ,'_dilate_figs.mat')) );
    save(savepath, 'run_2_dilates');
    
    % save subplot figures
    limit = [0, noise_limit_max];
    % run_1_orig
    f = figure('visible', 'off');
    for i = 1 : size(run_1_origs, 1)
        subplot(1,4,i);
        im = cell2mat(run_1_origs(i));
        limit = [0, max(im(:))];
        imagesc(im);
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(1), '-orig' )) );
    saveas(f, savepath, 'jpeg');
    
    
    % run_2_orig
    f = figure('visible', 'off');
    for i = 1 : size(run_2_origs, 1)
        subplot(1,4,i);
        im = cell2mat(run_2_origs(i));
        limit = [0, max(im(:))];
        imagesc(im);
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(2), '-orig' )) );
    saveas(f, savepath, 'jpeg');
    
    
    % run_1
    f = figure('visible', 'off');
    for i = 1 : size(run_1_figs, 1)
        subplot(1,4,i);        
        im = cell2mat(run_1_figs(i));
        limit = [0, max(im(:))];
        imagesc(im);
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(1), '-figs' )) );
    saveas(f, savepath, 'jpeg');
    
    % run_2
    f = figure('visible', 'off');
    for i = 1 : size(run_2_figs, 1)
        subplot(1,4,i);
        im = cell2mat(run_2_figs(i));
        limit = [0, max(im(:))];
        imagesc(im); 
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(2), '-figs' )) );
    saveas(f, savepath, 'jpeg');
    
    
    % run_1_dilate
    f = figure('visible', 'off');
    for i = 1 : size(run_1_dilates, 1)
        subplot(1,4,i);        
        im = cell2mat(run_1_dilates(i));
        lim_m = max(im(:));
        if isnan(lim_m)
            lim_m = 0.1;
        end
        limit = [0, lim_m];
        imagesc(im); 
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(1), '-dilate-figs' )) );
    saveas(f, savepath, 'jpeg');
    
    % run_2_dilate
    f = figure('visible', 'off');
    for i = 1 : size(run_2_dilates, 1)
        subplot(1,4,i);
        im = cell2mat(run_2_dilates(i));
        lim_m = max(im(:));
        if isnan(lim_m)
            lim_m = 0.1;
        end
        limit = [0, lim_m];
        imagesc(im);
        colormap(jet), caxis(limit), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(2), '-dilate-figs' )) );
    saveas(f, savepath, 'jpeg');
    
    % diff
    f = figure('visible', 'off');
    for i = 1 : size(diff_figs, 1)
        subplot(1,4,i);
        im = cell2mat(diff_figs(i));
        imagesc(im);
        c = max(im(:));
        if c == 0
            c = 0.1;
        end
        colormap(jet), caxis([0, c]), colorbar;
    end
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    savepath = cell2mat( fullfile(data_path, cur_folder.name, strcat(run_code(1), '-', run_code(2), '-diffs' )) );
    saveas(f, savepath, 'jpeg');
    
end
