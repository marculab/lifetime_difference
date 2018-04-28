# lifetime_difference

Please make the directory structure as below:

    -Root
        -1697_2
            -20170601_run_003_export_NoRS.mat
            -20170601_run_005_export_NoRS.mat

In each pair's folder only place one pair of mat files (i.e. `20170601_run_003_export_NoRS.mat` and `20170601_run_005_export_NoRS.mat`)

Run `it_diff_iteration.m` to do registration and generate the difference mat files and figures, for each pair (folder). 

* The output files will be stored the same directory as the file loaded

With the mat files ('diffs_ch1.mat', 'diffs_ch2.mat', 'diffs_ch3.mat', 'diffs_ch4.mat') saved, run `stats_figure_all.m` to calculate the distribution of the difference of all the cases included. Output the histogram as figures and mat files.