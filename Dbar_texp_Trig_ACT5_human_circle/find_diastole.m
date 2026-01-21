function diastoles = find_diastole(Vmulti_perf)
%==========================================================================
% This function finds the ventricular mid-systole frames in a multiframe 
% dataset.
%
% Author:              Lydia Lonzarich 
% Date Modified:        December, 2025
%
% Inputs:   Vmulti_perf = the dataset frames that we want to find mid-systole
%                         locations in
%                   
% Outputs:  mid_systoles = the mid_systole indices (aka, the frames where
%                          ventricular mid-systoles occur. 
%
% External files:       None                         
%==========================================================================

[~, ~, num_frames] = size(Vmulti_perf);

cardiac_signal = zeros(num_frames, 1);

for i = 1:num_frames
    V = Vmulti_perf(:, :, i);
    cardiac_signal(i) = norm(V(:), 2);
end

cardiac_signal = detrend(cardiac_signal);

% plot
figure;
plot(cardiac_signal);
title('Extracted Cardiac Signal from EIT Voltages')

% - peaks ~= systole, troughs = diastole
% - more specifically, peaks correspond to full contraction which is exactly what mid-systole represents (b/n start and end of ventricular contraction).
% - find peaks and the indices at which they occur
[~, pks_indices] = findpeaks(-cardiac_signal, 'MinPeakDistance', 10);
% disp(pks_indices)

% return the frames at which start-systoles occur.
diastoles = pks_indices;

