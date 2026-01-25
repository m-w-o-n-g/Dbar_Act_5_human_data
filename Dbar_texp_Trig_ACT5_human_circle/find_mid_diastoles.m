function mid_diastoles = find_mid_diastoles(Vmulti_perf)
%==========================================================================
% This function finds the ventricular mid-diastole frames in a multiframe 
% dataset.
%
% Author:               Lydia Lonzarich 
% Date Modified:        January 2026
%
% Inputs:   Vmulti_perf = the dataset frames that we want to find mid-diasystole
%                         locations in
%                   
% Outputs: mid-diastoles = the mid-diastole indices (aka, the frames where
%                          ventricular mid-diasystoles occur. 
%
% External files:       None                         
%==========================================================================

[~, ~, num_frames] = size(Vmulti_perf);

cardiac_signal = zeros(num_frames, 1);

for i = 1:num_frames
    V = Vmulti_perf(:, :, i);
    cardiac_signal(i) = norm(V(:), 2);
end

% extract cardiac signal.
cardiac_signal = detrend(cardiac_signal);

% OPTIONAL: plot cardiac signal (frame # on x-axis, magnitude of voltage on y-axis). 
figure;
plot(cardiac_signal);
title('Extracted Cardiac Signal from EIT Voltages')

% find peaks (end-systole) and troughs (end-diastole)
prominence = 0.5 * std(cardiac_signal);
[~, pks_indices] = findpeaks(cardiac_signal, 'MinPeakProminence', prominence, 'MinPeakDistance', 15);
[~, trough_indices] = findpeaks(-cardiac_signal, 'MinPeakProminence', prominence, 'MinPeakDistance', 15);

% calculate the mid-diastoles.
mid_diastoles = [];
for j = 1:length(pks_indices)
    curr_peak = pks_indices(j);

    future_troughs = trough_indices(trough_indices > curr_peak);

    if ~isempty(future_troughs)
        next_trough = future_troughs(1);
        
        p_val = cardiac_signal(curr_peak);
        t_val = cardiac_signal(next_trough);
        half_val = (p_val + t_val)/2;

        search_range = curr_peak:next_trough;
        [~, min_idx] = min(abs(cardiac_signal(search_range) - half_val));
        % min_idx = search_range(min_idx);
        mid_diastole = min_idx + curr_peak;
        
        mid_diastoles = [mid_diastoles, mid_diastole];
    end
end


