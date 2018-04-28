load('diff_list.mat');
load('case_id.mat');

% for each case
for i = 1 : size(diff_list, 2)
diff_case = diff_list{1, i};
diff_case_all = [];
for j = 1 : 4 % for each channel
    diff_case_channel = cell2mat(diff_case(j));
    diff_case_all = cat(1, diff_case_all, diff_case_channel);
    % histogram
    title_str = strcat('Difference distribution: Case ', case_id(i), ' CH', int2str(j));
    f = figure('visible', 'off');
    h = histogram(diff_case_channel);
    title(title_str);
    xlim([0 6.5]) % choose 6.5 because the biggest difference in all cases/channels < 6.5
    % save image
    fname = char(strcat('Case', case_id(i),'-distribution', '-CH', int2str(j)));
    set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 7.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 7.5]);
    saveas(f, fname, 'jpeg');
end
% box plot
% outliers: boxplot draws points as outliers if they are 
% greater than q3 + w × (q3 – q1) or less than q1 – w × (q3 – q1), 
% where w is the maximum whisker length, 
% and q1 and q3 are the 25th and 75th percentiles of the sample data, respectively.
% here we choose w as 3 (extreme Outliers)
f = figure('visible', 'off');
x_labels = {'CH1', 'CH2', 'CH3', 'CH4'};
x = cat(1, diff_case_all, diff_case_all);
g = [ones(1, size(cell2mat(diff_case(:,1)) ,1) ), 2*ones(1, size(cell2mat(diff_case(:,2)) ,1) ), 3*ones(1, size(cell2mat(diff_case(:,3)) ,1) ), 4*ones(1, size(cell2mat(diff_case(:,4)) ,1) ), 5*ones(1, size(diff_case_all,1))];
gcell = num2cell(g);
gcell(g==1) = {'CH1'};
gcell(g==2) = {'CH2'};
gcell(g==3) = {'CH3'};
gcell(g==4) = {'CH4'};
gcell(g==5) = {'Overall'};
boxplot(x, gcell, 'Whisker', 3, 'OutlierSize',5);
title_str = strcat('Box plot: Case ', case_id(i));
title(title_str);
fname = char(strcat('Case', case_id(i),'-box'));
set(gcf, 'Units', 'Inches', 'Position', [0, 0, 4.5, 7.5], 'PaperUnits', 'Inches', 'PaperSize', [4.5, 7.5]);
saveas(f, fname, 'jpeg');

end
% overall