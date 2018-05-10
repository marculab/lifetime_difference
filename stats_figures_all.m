% After the difference of all the cases are generated, run this script to
% count the statistics of the difference for all the included cases
    % Will load 'diffs_ch1.mat', 'diffs_ch2.mat', 'diffs_ch3.mat',
    % 'diffs_ch4.mat'
% The program will output 2 image files and 1 mat file
% histogramdata.mat
    % channel 4 also included
% Diff Distribution All Cases.jpg
% Diff Distribution All cases All channels (Without CH4).jpg

% get list of input files
data_path = '../All Artery Data';
disp(['root folder: ', data_path]);
current_folders = dir(data_path);
current_folders = current_folders([current_folders.isdir]);
current_folders = current_folders(3:end);

all_case_ch1 = [];
all_case_ch2 = [];
all_case_ch3 = [];
all_case_ch4 = [];
allcase_alldiff_noch4 = [];

for i = 1 : length(current_folders) % for each folder
    diff_ch1 = [];
    diff_ch2 = [];
    diff_ch3 = [];
    diff_ch4 = [];
    max_x = 0;
    cur_folder = current_folders(i);
    matfiles = dir(fullfile(data_path, cur_folder.name));
    matfiles = matfiles(~[matfiles.isdir]);
    matcount = 0;
    for j = 1 : length(matfiles)
        [filepath,name,ext] = fileparts(matfiles(j).name);
        fullpath = fullfile(matfiles(j).folder, matfiles(j).name);
        if strcmp(ext, '.mat') && strcmp(name, 'diffs_ch1')
            disp(['Loading ', fullpath]);
            load(fullpath);
            diff_ch1 = imdiff(:);
            all_case_ch1 = cat(1, all_case_ch1, diff_ch1);
            matcount = matcount + 1;
            if max_x < max(diff_ch1)
                max_x = max(diff_ch1);
            end
        elseif strcmp(ext, '.mat') && strcmp(name, 'diffs_ch2')
            disp(['Loading ', fullpath]);
            load(fullpath);
            diff_ch2 = imdiff(:);
            all_case_ch2 = cat(1, all_case_ch2, diff_ch2);
            matcount = matcount + 1;
            if max_x < max(diff_ch1)
                max_x = max(diff_ch1);
            end
        elseif strcmp(ext, '.mat') && strcmp(name, 'diffs_ch3')
            disp(['Loading ', fullpath]);
            load(fullpath);
            diff_ch3 = imdiff(:);
            all_case_ch3 = cat(1, all_case_ch3, diff_ch3);
            matcount = matcount + 1;
            if max_x < max(diff_ch1)
                max_x = max(diff_ch1);
            end
        elseif strcmp(ext, '.mat') && strcmp(name, 'diffs_ch4')
            disp(['Loading ', fullpath]);
            load(fullpath);
            diff_ch4 = imdiff(:);
            all_case_ch4 = cat(1, all_case_ch4, diff_ch4);
            matcount = matcount + 1;
            if max_x < max(diff_ch1)
                max_x = max(diff_ch1);
            end
        end  
    end
    if matcount ~= 4
        continue;
    end
    fprintf('%s\n', cur_folder.name);
    % histogram
    title_str = strcat('Difference distribution');
    f = figure('visible', 'off');
    title(title_str);
    subplot(1,4,1);
    h = histogram(diff_ch1);
    title('CH1');
    xlim([0 max_x])
    subplot(1,4,2);
    h = histogram(diff_ch2);
    title('CH2');
    xlim([0 max_x])
    subplot(1,4,3);
    h = histogram(diff_ch3);
    title('CH3');
    xlim([0 max_x])
    subplot(1,4,4);
    h = histogram(diff_ch4);
    title('CH4');
    xlim([0 max_x])
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);
    
    savepath = fullfile(data_path, cur_folder.name, 'Diff Distribution');
    saveas(f, savepath, 'jpeg');
    
    % histogram all but ch4
    diff_all_noch4 = cat(1, diff_ch1, diff_ch2, diff_ch3);
    title_str = 'Difference distribution All channels (Without CH4)';
    f = figure('visible', 'off');
    h = histogram(diff_all_noch4);
    xlim([0 max_x])
    savepath = fullfile(data_path, cur_folder.name, 'Diff Distribution All channels (Without CH4)');
    saveas(f, savepath, 'jpeg');
    
    % add to all
    allcase_alldiff_noch4 = cat(1, allcase_alldiff_noch4, diff_all_noch4);
end

% diff for all cases
outfile = fopen(fullfile(data_path, strcat('stats_allcases', '.csv')), 'wt');
fprintf(outfile, [',', 'mean,', 'median,', 'std,', 'min,', 'max,\n']);

allch1_min = min(all_case_ch1);
allch1_max = max(all_case_ch1);
allch1_mean = mean(all_case_ch1);
allch1_med = median(all_case_ch1);
allch1_std = std(all_case_ch1);
fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [1 allch1_mean allch1_med allch1_std allch1_min allch1_max]');

allch2_min = min(all_case_ch2);
allch2_max = max(all_case_ch2);
allch2_mean = mean(all_case_ch2);
allch2_med = median(all_case_ch2);
allch2_std = std(all_case_ch2);
fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [2 allch2_mean allch2_med allch2_std allch2_min allch2_max]');

allch3_min = min(all_case_ch3);
allch3_max = max(all_case_ch3);
allch3_mean = mean(all_case_ch3);
allch3_med = median(all_case_ch3);
allch3_std = std(all_case_ch3);
fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [3 allch3_mean allch3_med allch3_std allch3_min allch3_max]');

allch4_min = min(all_case_ch4);
allch4_max = max(all_case_ch4);
allch4_mean = mean(all_case_ch4);
allch4_med = median(all_case_ch4);
allch4_std = std(all_case_ch4);
fprintf(outfile, 'CH%d, %f, %f, %f, %f, %f,\n', [4 allch4_mean allch4_med allch4_std allch4_min allch4_max]');

allcase_min = min(allcase_alldiff_noch4);
allcase_max = max(allcase_alldiff_noch4);
allcase_mean = mean(allcase_alldiff_noch4);
allcase_med = median(allcase_alldiff_noch4);
allcase_std = std(allcase_alldiff_noch4);
fprintf(outfile, 'CH1~3, %f, %f, %f, %f, %f,\n', [allcase_mean allcase_med allcase_std allcase_min allcase_max]');

% histogram
max_x = max(allcase_alldiff_noch4);
title_str = strcat('Difference distribution');
f = figure('visible', 'off');
title(title_str);
subplot(1,4,1);
h = histogram(all_case_ch1);
title('All cases CH1');
xlim([0 max_x])
subplot(1,4,2);
h = histogram(all_case_ch2);
title('All cases CH2');
xlim([0 max_x])
subplot(1,4,3);
h = histogram(all_case_ch3);
title('All cases CH3');
xlim([0 max_x])
subplot(1,4,4);
h = histogram(all_case_ch4);
title('All cases CH4');
xlim([0 max_x])
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0.04, 0.04, 0.9, 0.9]);

% save histogram data as mat
Disp('Save histogram_data to mat file');
savepath = char(fullfile(data_path, 'histogramdata.mat'));
histogram_data = {all_case_ch1, all_case_ch2, all_case_ch3, all_case_ch4};
save(savepath, 'histogram_data');

savepath = fullfile(data_path, 'Diff Distribution All Cases');
saveas(f, savepath, 'jpeg');

% histogram all but ch4
title_str = 'Difference distribution All channels (Without CH4)';
f = figure('visible', 'off');
h = histogram(allcase_alldiff_noch4);
xlim([0 max_x])
savepath = fullfile(data_path, 'Diff Distribution All cases All channels (Without CH4)');
saveas(f, savepath, 'jpeg');

