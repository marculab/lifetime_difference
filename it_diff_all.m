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
erode_radius = 0;
        
% get list of input files (suppose each pair of .mat file is in a seperate folder)
data_path = '../All Artery Data';
current_folders = dir(data_path);
current_folders = current_folders([current_folders.isdir]);
current_folders = current_folders(3:end);

% go through each folder to load diffimg.mat
for i = 1 : length(current_folders) % for each folder
    matfiles = dir(fullfile(data_path, current_folders(i).name));
    matfiles = matfiles(~[matfiles.isdir]);
    img_pair = {};
    run_code = {};
    for j = 1 : length(matfiles)
        [filepath,name,ext] = fileparts(matfiles(j).name);
        if strcmp(ext, '.mat')
            if matfiles(j).name == 'diffimg.mat'
                fullpath = fullfile(matfiles(j).folder, matfiles(j).name);
                load(fullpath);
                % add to difflist
                continue;
            end
        end
    end
    if size(img_pair) ~= 2
        continue;
    end
    fprintf('%s\n', current_folders(i).name);
    % variables for this image pair
    diffs = {}; % store difference value of each channel, as 1D array
    
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