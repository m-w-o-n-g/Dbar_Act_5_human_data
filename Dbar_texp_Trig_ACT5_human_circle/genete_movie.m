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
data_fname = 'dec19_1933_gamma_cond_distributions_57_1';

load([data_dir, data_fname]); % load gamma_all

movie_outdir = 'Dbar_human_recons_movies';                 % directory for the .mat measurements file and .avi recon movie file.
movie_mat_fname = 'dec19_2114_Dbar_movie_sbj001_57_1';                % filename for output .mat measurements file.
movie_avi_fname = [movie_mat_fname, '.avi'];               % filename for output .avi reconstruction movie.

% create the movie outdir if it doesn't exist.
if ~exist(movie_outdir, 'dir')
        mkdir(movie_outdir)
end


% =================================================================================================
% ================= Standardize Colorbar for Individual Frame Plots ===============================
% =================================================================================================
% truncate the colorbar. enter an integer between 0 and 10. if displaying a small number of images, choose something smaller
percent_to_truncate_colorbar = 0;

%(option 1) og min/max method.
max_gamma_all = max(max(max(gamma_all)));
min_gamma_all = min(min(min(gamma_all)));
range_gamma = max_gamma_all - min_gamma_all;
cmax_gamma = max_gamma_all - 0.1*range_gamma;
cmin_gamma = min_gamma_all + 0.1*range_gamma;

% (option 2) robust percentile method.
% all_vals = gamma_all(:);
% cmin_gamma = prctile(all_vals, 1);
% cmax_gamma = prctile(all_vals, 99);



% ==================================================================================================
% =================================== Create a Movie with the Image Reconstructions ================
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

% initialize variable to keep track of the current frame # in the for-loop.
frame_idx = 1; 

all_frames = 1:size(gamma_all, 3);

% iterate over all frames
for frame_num = all_frames
    
    disp("frame number: " + frame_num)

    % choose [yes/no] to display movie to screen.
    if display_movie_to_screen == 1
        figure('visible','on');
    else
        figure('Visible','off');
    end

    colormap(cmap)

    % generate the pretty image reconstruction.
    imagesc(flipud(gamma_all(:,:,frame_idx))) % this one was og
    
    caxis([cmin_gamma,cmax_gamma])
    colorbar
    axis square
    set(gca, 'Ydir', 'normal')

    title(['Frame number = ',num2str(frame_num)]) % add title to figure for reference frame number.
    
    frame_num_double = double(frame_num); % this conversion is somehow needed for title.
    
    % Convert frame_num to string.
    frame_str = ['Frame Number: ' num2str(frame_num_double)];
        
    frame_pick = getframe(gcf);
    
    writeVideo(writerObj, frame_pick);
    
    frame_idx = frame_idx + 1; % to iterate through all frames in gamma all.

end % END PLOTTING MOVIE


%% save movie to file
% choose [yes/no] to save the movie to a .avi file.
if save_movie_as_mat == 1
    save([movie_outdir '/' movie_mat_fname, '.mat'], 'gamma_all', 'hh', 'xx');
end

% close the video writer object
close(writerObj);
