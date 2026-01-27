%===================================================================================================
% This script trims the number of frames in a multiframe dataset .mat file to a user-specified length.
%
% Authors:            Lydia Lonzarich
% Date Modified:      October, 2025
%
%===================================================================================================
clear all

% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name for .mat file containing EIT data 
datafname = 'Sbj001_93kHz_perf_24_10_15_11_05_09_1';  

load([datadir, datafname]);

% start_frame and end_frame given by Jennifer. (127 frames)
start_frame = 300;
end_frame = 500; 

% update variables that rely on frame number.
EOF_data = EOF_data(:,start_frame:end_frame);
cali_temperature = cali_temperature(:,start_frame:end_frame);
corrected_frame_current = corrected_frame_current(:,start_frame:end_frame);
frame_ECG = frame_ECG(:,:,start_frame:end_frame);
frame_current = frame_current(:,start_frame:end_frame);
frame_voltage = frame_voltage(:,:,start_frame:end_frame);
meas_imp_C = meas_imp_C(:,start_frame:end_frame);
meas_imp_R = meas_imp_R(:,start_frame:end_frame);


save([datadir,'chunk_',datafname])