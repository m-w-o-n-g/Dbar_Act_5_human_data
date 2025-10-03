%===================================================================================================
%======================================== Specify External Data ====================================
%===================================================================================================
clear all

% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name for .mat file containing EIT data 
datafname = 'Sbj02_2D_16e_24_10_16_12_39_39_93750';  

load([datadir, datafname]);
cur_pattern = cur_pattern(1:16,1:16);
frame_voltage = frame_voltage(1:16,1:16,:);

save([datadir,'modified_16x16_',datafname])




