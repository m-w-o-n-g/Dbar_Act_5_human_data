% a tool to guide frame selection.
% converts raw EIT voltages into GE-expected format
% visualize voltage meshes
% computes a power waveform to help you find a stable section of breathing/heartbeat.
% - power waveforms show you how voltages and currents fluctuate over time,visually 

%% load data (first dataset)
% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name for .mat file containing EIT data 
datafname = 'modified_16x16_Sbj02_2D_16e_24_10_16_12_38_03_93750';  

load([datadir, datafname]);

Vmulti_frames = real(frame_voltage(1:15,1:16,:));  % grab the (real parts of the) voltages from rows 1-15, columns 1-16 ==> size of vmulti is 15x16.

output_directory = './';


%% now, we'll convert each 15x16 frame into a 240x1 vector 
% Need to reshape this to be the same as the GE datasets (num_frames by 992)
[m,n,num_frames]= size(Vmulti_frames); % get the number of frames in the dataset.

% preallocate a 240 x num_frames matrix, where each col = flattened frame
Data_perm = zeros(240,num_frames);

% iterate through each 15x16 frame.
for ii = 1:num_frames
    temp = squeeze(Vmulti_frames(:,:,ii));  % get the frame's 15x16 voltage matrix 
    %return
    temp = temp.';                          % transpose the voltage matrix ==> 16x15
%     size(temp)
%     return
    temp=temp(:);                           % creates a 240x1 vector.
    Data_perm(:,ii) = temp;                 % store this vector as the ii-th col.
end

DataVol0 = Data_perm.';    % transpose data_perm to be num_frames x 240.

[m, n] = size(DataVol0);      % m is the number of frames, n the number of measurements
Slides = m;                  % Slides = total frames of the dataset

Vmulti_t = DataVol0.';                   % Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t, 16, 15, Slides); % convert flattened voltage vectors back into 16x15 slides.


%% Define applied currents - they will be the same every time
current_pat = cur_pattern(1:16,1:15);  % load the current pattern matrix (16 electrodes by 15 patterns). These are in amps;

% [Slides,num_mesh_elts]=size(DataVol0);                  % Slides = total frames of the dataset
% [first_comp,ref_frame] = findRefFrameTest(real(DataVol0),Slides,32);
Vmulti_CP1 = squeeze(real(Vmulti(:, 1, :))); % take all rows (aka, electrodes) in col 1 (aka, current pattern 1) for all frames ==> 16 electrodes x slides frames.


%% plot this as a mesh.
figure
mesh(Vmulti_CP1) 
title('Voltages for CP 1')


%% Compute power waveform
% line 62 gives a single number per frame/slide which represents the overall "strength" ofthe measured voltage pattern for that frame.
% for ventilation: - breathing = large,slow variations in conductivity aslungs go from full of air (low conducitivty) to empty (high conductivity)
%                  - power waveform should look like a slow sinusoidal wave w/ peaks=end of expiration, troughs=end of inspiration.
% 
% for perfusion: - we only want to see cardiac signal = small, fast oscillation and little wiggles on top of breathing curve. 
%                - power waveform should have small peaks=increased and decreased conductivity.


% iterate over each frame
for ii = 1:Slides
  Vframe = squeeze(real(Vmulti(:,:,ii)));     % extract the ii-th frame's voltage matrix (16x15)
  power_wavef_mx(:,:,ii) = current_pat'*Vframe; % should be size L-1 by L-1 by num frames
end

figure
power_wavef = squeeze(power_wavef_mx(1,1,:)); % extract the (1,1) element in each frame ==> a scalar waveform.
plot(power_wavef)                             % plot this
title('Power Waveform First data set')




%%  Second data set
load('Sbj02_2D_16e_24_10_16_12_39_39_93750.mat')
Vmulti_frames = real(frame_voltage(1:15,1:16,:));
output_directory = './';

% Need to reshape this to be the same as the GE datasets (num_frames by 992)
[m,n,num_frames]= size(Vmulti_frames);
Data_perm=zeros(240,num_frames);
for ii=1:num_frames
    temp = squeeze(Vmulti_frames(:,:,ii));  
    %return
    temp = temp.';  % 32 by 31
%     size(temp)
%     return
    temp=temp(:);
    Data_perm(:,ii) = temp;
end

DataVol0 = Data_perm.';  % Num_frames by 992

[m,n]=size(DataVol0);  % m is the number of frames, n the number of measurements
Slides=m;                  % Slides = total frames of the dataset
Vmulti_t = DataVol0.';  %Change to num_mesh_elts by number of frames
Vmulti = reshape(Vmulti_t,16,15,Slides);

%% Define applied currents - they will be the same every time
current_pat = cur_pattern(1:16,1:15);  % These are in amps;

% [Slides,num_mesh_elts]=size(DataVol0);                  % Slides = total frames of the dataset
Vmulti_CP1=squeeze(real(Vmulti(:,1,:)));
figure
mesh(Vmulti_CP1)
title('Voltages for CP 1')

% Compute power waveform
for ii=1:Slides
  Vframe= squeeze(real(Vmulti(:,:,ii)));
  power_wavef_mx(:,:,ii)=current_pat'*Vframe; % should be size L-1 by L-1 by num frames
end
figure
power_wavef = squeeze(power_wavef_mx(1,1,:));
plot(power_wavef)
title('Power Waveform Second data set')
