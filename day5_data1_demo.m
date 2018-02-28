set(0,'DefaultTextInterpreter','none');
load('Data/Subject_1/fMRI_BOLD_FC.mat');
ROI_names=csvimport('Data/Subject_1/freesurfer_regions_68_sort_full.txt');
figure;
plot(ROIts_68(1:10:end,22)-mean(ROIts_68(1:10:end,22))); %Demeaned timeseries of Left Posterior Cingulate Cortex
hold on; 
plot(ROIts_68(1:10:end,56)-mean(ROIts_68(1:10:end,56))); %Demeaned timeseries of Right Posterior Cingulate Cortex
disp(FC_cc_68(22,56));
hold off;
figure;
plot(ROIts_68(1:10:end,39)-mean(ROIts_68(1:10:end,39))); %Demeaned timeseries of Right Entorhinal Cortex
hold on; 
plot(ROIts_68(1:10:end,51)-mean(ROIts_68(1:10:end,51))); %Demeaned timeseries of Right Pars Opercularis 
disp(FC_cc_68(39,51));