% ==================================================================================================
% This script generates a D-bar reconstruction "movie" by...
% - loads in gamma_all (a matrix filled with each frame's conductivity reconstruction (aka, gamma)
% - plots each gamma using imagesc

% This script can be used as a general-purpose movie-generating script, as it just requires a filled 
% matrix with gamma reconstruction(s).
%
% Note: this code is a little messy and could be cleaned up some, but it does seem to work fine
%
% Author:            Lydia Lonzarich
% Date Modified:     December 10, 2025
% 
% ==================================================================================================
clear 

% ==================================================================================================
% ================================= Choose What to Plot and Save ===================================
% ==================================================================================================
display_movie_to_screen = 0;
save_movie_as_mat = 0;

% ==================================================================================================
% =================================== load external data ===========================================
% ==================================================================================================
% directory where data is stored.
data_dir = 'gammas/';

% .mat file containing each frame's gamma_real conductivity reconstruction measurement data.
data_fname = 'feb4_2340_gamma_cond_distributions_Sbj002_vent_set';

load([data_dir, data_fname]); % load gamma_all

movie_outdir = 'Dbar_human_recons_movies';                 % directory for the .mat measurements file and .avi recon movie file.
movie_mat_fname = 'feb4_2359_Dbar_movie_sbj002_vent_set';                % filename for output .mat measurements file.
movie_avi_fname = [movie_mat_fname, '.avi'];               % filename for output .avi reconstruction movie.

% create the movie outdir if it doesn't exist.
if ~exist(movie_outdir, 'dir')
        mkdir(movie_outdir)
end


% =================================================================================================
% ================= Smoothing gamma_all ===========================================================
% =================================================================================================
% this is an exponential moving avg == it's causal / only looks at the past.
% gamma_smoothed = gamma_all;
% alpha = 0.3;
% for k = 2:size(gamma_all, 3)
%     gamma_smoothed(:,:,k) = alpha * gamma_all(:,:,k) + (1 - alpha) * gamma_smoothed(:,:,k-1);
% end

% this is a centered moving avg == it looks at the past and present.
% use a moving avg plot to smooth the gammas in gamma_all by removing high-frequency noise.
% adding the 'gaussian' parameter to smooth in time (different from adding gaussian to dbar algorithm because that was smoothing in space).
gamma_smoothed = smoothdata(gamma_all, 3, 'gaussian', 5);






% % =========================================================================
% % DIAGNOSTIC TOOL: Find the "Ghost" Lung
% % =========================================================================
% 
% % 1. Create a map of "Change." The lungs will be the brightest spots.
% activity_map = std(gamma_all, 0, 3); 
% 
% figure('Name', 'CLICK THE STUCK LUNG');
% imagesc(activity_map);
% colorbar;
% title('Click on the center of the Right Lung now...');
% axis square;
% colormap(jet);
% 
% % 2. Let you click on the pixel you are worried about
% [y_click, x_click] = ginput(1); % Returns (x, y) coordinates
% c_idx = round(y_click); % Column index (x-axis click)
% r_idx = round(x_click); % Row index (y-axis click)
% 
% close(gcf); % Close the map
% 
% % 3. Extract and Plot the Time Course for THAT pixel
% pixel_raw = squeeze(gamma_all(r_idx, c_idx, :));
% 
% % Re-calculate smooth with your window of 7
% if exist('smoothdata', 'file')
%     pixel_smooth = smoothdata(pixel_raw, 'gaussian', 7);
% else
%     pixel_smooth = conv(pixel_raw, ones(1,7)/7, 'same');
% end
% 
% % 4. Plot the TRUTH
% figure('Name', 'Lung Truth Graph');
% plot(pixel_raw, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw Data'); hold on;
% plot(pixel_smooth, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Smoothed (Window 7)');
% yline(1.0, 'k--', 'Baseline (1.0)', 'LineWidth', 2); 
% legend;
% grid on;
% title(['Conductivity at Row ' num2str(r_idx) ', Col ' num2str(c_idx)]);
% xlabel('Frame Number');
% ylabel('Conductivity (Gamma)');










% % Pick a pixel index that corresponds to the Right Lung
% % (You might need to adjust these indices based on your grid size)
% % For a 20x20 grid, row 10, col 15 is often right lung.
% r_idx = 10; 
% c_idx = 15; 
% 
% % Extract the time course for this pixel
% pixel_raw = squeeze(gamma_all(r_idx, c_idx, :));
% pixel_smooth = squeeze(gamma_smoothed(r_idx, c_idx, :));
% 
% figure;
% plot(pixel_raw, 'Color', [0.7 0.7 0.7], 'DisplayName', 'Raw Data'); hold on;
% plot(pixel_smooth, 'LineWidth', 2, 'Color', 'r', 'DisplayName', 'Smoothed (Window 20)');
% yline(1.0, 'k--', 'Baseline (1.0)'); % The theoretical "empty" line
% legend;
% title('Right Lung Conductivity Over Time');
% xlabel('Frame Number');
% ylabel('Conductivity');
% 
% % IF the red line never touches the dotted black line, 
% % your data simply does not return to baseline (Reference Frame issue).
% % IF the grey line touches the dotted line but the red line doesn't,
% % your smoothing window (20) is too big.





% =================================================================================================
% ================= Standardize Colorbar for MOVIE plotting ===============================
% =================================================================================================
% truncate the colorbar. enter an integer between 0 and 10. if displaying a small number of images, choose something smaller
% percent_to_truncate_colorbar = 0;

%(option 1) og min/max method.
% max_gamma_all = max(max(max(gamma_all)));
% min_gamma_all = min(min(min(gamma_all)));
% range_gamma = max_gamma_all - min_gamma_all;
% cmax_gamma = max_gamma_all - 0.1*range_gamma;
% cmin_gamma = min_gamma_all + 0.1*range_gamma;

% (option 2) robust percentile method.
all_vals = gamma_smoothed(:);
cmin_gamma = prctile(all_vals, 3);
cmax_gamma = prctile(all_vals, 97);



% ==================================================================================================
% =================================== Set up plotting params and figure for movie ==================
% ==================================================================================================
cmap = 'jet'; % select colormap for plots.

hh = 0.2;
xx = -1:hh:1;

% create video writer object in the output directory.
writerObj = VideoWriter([movie_outdir, '/', movie_avi_fname]);

% set the frame rate to one frame per second
set(writerObj,'FrameRate',5);

% open the writer object.
open(writerObj);

% create a figure window
if display_movie_to_screen == 1
    figure = figure('visible', 'on');
else
    figure = figure('visible', 'off');
end

% initialize the plot with 1 frame
image = imagesc(xx, xx, rot90(gamma_smoothed(:,:,1)));

% SET params for image tiles (doing this outside the loop so it's standardized across tiles).
colormap(cmap)
clim([cmin_gamma, cmax_gamma]);
colorbar
axis square
set(gca, 'Ydir', 'normal')
title_for_tile = title('Frame number: 1');


% ==================================================================================================
% =================================== Generate Movie ==================================================
% ==================================================================================================
all_frames = 1:size(gamma_smoothed, 3);
num_frames = length(all_frames);

% iterate over all frames
for frame_num = all_frames
    
    disp("Frame number: " + frame_num)
    
    curr_data = rot90(gamma_smoothed(:,:,frame_num));

    set(image, 'CData', curr_data);
 
    set(title_for_tile, 'String', sprintf('Frame Number: %d', frame_num));

    drawnow;

    frame_pick = getframe(figure);

    writeVideo(writerObj, frame_pick);

end % END PLOTTING MOVIE


%% save movie to file
% choose [yes/no] to save the movie params to a .mat file.
if save_movie_as_mat == 1
    save([movie_outdir '/' movie_mat_fname, '.mat'], 'gamma_all', 'hh', 'xx');
end

% close the video writer object
close(writerObj);

disp("Movie has been generated.")
