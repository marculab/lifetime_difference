% measure difference of image pairs and compute stats

% add path of registration code
addpath('registration code/FLIm_rot_reg');

% variables
diff_matrices = {};
diff_mean_pixel = [];

show_figs = 1;

% get list of input files (suppose each pair of .mat file is in a seperate folder)
current_folders = dir();
current_folders = current_folders([current_folders.isdir]);
current_folders = current_folders(3:end);

% choose the channel (3rd dimension of AllData)
for channel = 9 : 12
% go through each folder
for i = 1 : length(current_folders)
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
    if size(img_pair) < 2
        continue;
    end
    % process image pair
    run_1 = cell2mat(img_pair(1));
    run_1 = run_1(:,:,channel);
    run_2 = cell2mat(img_pair(2));
    run_2 = run_2(:,:,channel);
    run_1(run_1 < 1) = 1;
    run_1(run_1 > 12) = 12;
    run_2(run_2 < 1) = 1;
    run_2(run_2 > 12) = 12;
    
    % register run2 to run_1
    try
        tform = reg_FLIm(run_2, run_1);
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
    
  
    % smooth using 4*4 mean
    h = ones(4, 4) / 16;
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
    diff_matrices = [diff_matrices, diff_img];
    diff_mean_pixel = [diff_mean_pixel, sum(diff_img(:)) / (size(diff_img, 1) * size(diff_img, 2))];
    
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
    
    % compute stats
    diff_mean = mean(diff_mean_pixel);
    diff_max = max(diff_mean_pixel);
    diff_min = min(diff_mean_pixel);
    diff_std = std(diff_mean_pixel);
    % write to text file
    outfile = fopen(strcat('stats-',int2str(channel), '.txt'), 'wt');
    fprintf(outfile, ['mean\t', 'max\t', 'min\t', 'std', 'r\n']);
    fprintf(outfile, '%f %f %f %f \n', [diff_mean diff_max diff_min diff_std]');
end

end
