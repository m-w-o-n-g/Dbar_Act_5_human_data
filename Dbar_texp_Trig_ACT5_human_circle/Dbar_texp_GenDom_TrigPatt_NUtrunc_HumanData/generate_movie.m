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
data_fname = 'dec20_1125_gamma_cond_distributions_57_1';

load([data_dir, data_fname]); % load gamma_all

movie_outdir = 'Dbar_human_recons_movies';                 % directory for the .mat measurements file and .avi recon movie file.
movie_mat_fname = 'dec20_1151_%_Dbar_movie_sbj001_57_1_gauss_trunc';    % filename for output .mat measurements file.
movie_avi_fname = [movie_mat_fname, '.avi'];               % filename for output .avi reconstruction movie.

% create the movie outdir if it doesn't exist.
if ~exist(movie_outdir, 'dir')
        mkdir(movie_outdir)
end



% =================================================================================================
% ================= Smooth Gamma_all ==============================================================
% =================================================================================================
gamma_all_smoothed = gamma_all;

alpha = 0.7;  % 0.6–0.8 works well
for k = 2:size(gamma_all,3)
    gamma_all_smoothed(:,:,k) = ...
        alpha*gamma_all(:,:,k) + (1-alpha)*gamma_all_smoothed(:,:,k-1);
end


% =================================================================================================
% ================= Standardize Colorbar for Individual Frame Plots ===============================
% =================================================================================================
% truncate the colorbar. enter an integer between 0 and 10. if displaying a small number of images, choose something smaller
percent_to_truncate_colorbar = 0;

%(option 1) og min/max method.
% max_gamma_all = max(max(max(gamma_all)));
% min_gamma_all = min(min(min(gamma_all)));
% range_gamma = max_gamma_all - min_gamma_all;
% cmax_gamma = max_gamma_all - 0.1*range_gamma;
% cmin_gamma = min_gamma_all + 0.1*range_gamma;

% (option 2) robust percentile method.
all_vals = gamma_all_smoothed(:);
cmin_gamma = prctile(all_vals, 2);
cmax_gamma = prctile(all_vals, 98);

% JUST REMOVED 12/19...
% % Standardizing the colorbar for the image reconstruction.
% max_gamma_all = max(max(max(gamma_all)));
% min_gamma_all = min(min(min(gamma_all)));
% range_gamma = max_gamma_all - min_gamma_all;
% cmax_gamma = max_gamma_all - 0.2*range_gamma; % should this be .2 or 0?
% cmin_gamma = min_gamma_all + 0.2*range_gamma;

% % recompute datamin and datamax from gamma_all for consistent colorbar scaling. 
% datamin = min(gamma_all(:));
% datamax = max(gamma_all(:)); 
% datarange = datamax-datamin;
% colorbartrunc = percent_to_truncate_colorbar * .01;
% datamin = datamin + colorbartrunc * datarange;
% datamax = datamax - colorbartrunc * datarange;




% ==================================================================================================
% =================================== Create a Movie with the Image Reconstructions ================
% ==================================================================================================

% Select colormap for figures
cmap = 'jet';

hh = 0.2; 
xx = -1:hh:1;

% create video writer object in the output directory.
writerObj = VideoWriter([movie_outdir, '/', movie_avi_fname]);

% set the frame rate to one frame per second
set(writerObj,'FrameRate',5);

% open the writer object.
open(writerObj);


% initialize variable to keep track of the current frame # in the for-loop.
% frame_idx = 1; 


%% Plot movie
% iterate over all frames (MINUS the reference frame)
% for frame_num = all_frames
for frame_idx = 1:size(gamma_all_smoothed, 3)

    disp("frame number: " + frame_idx)
    
    % choose [yes/no] to display each frame to screen.
    if display_movie_to_screen == 1 
        figure('visible','on');
    else
        figure('Visible','off');
    end

    colormap(cmap)

    % generate the pretty image reconstruction (this is from new code)
    % imagesc(xx,xx,squeeze(gam_real(frame_num,:,:)),[datamin, datamax]);
    % imagesc(xx,xx,squeeze(gam_real(frame_num,:,:))); % added

    % imagesc(xx, xx, rot90(squeeze(gamma_all(:,:, frame_idx))), [datamin, datamax]);
    imagesc(rot90(gamma_all_smoothed(:, :, frame_idx)));

    caxis([cmin_gamma,cmax_gamma]);
    colorbar
    axis square
    set(gca, 'Ydir', 'normal');
    % axis([-1 1 -1 1 ]);
    % colorbar; % added
    % axis off
    % axis square;

    % add title to each frame slide.
    title(['Frame number = ', num2str(frame_idx)])


    % generate the pretty image reconstruction (this is from og movie script)
    % imagesc(flipud(rot90(gamma_all(:,:,frame_idx))))
    % set(gca, 'Ydir', 'normal')
    % caxis([cmin_gamma,cmax_gamma])
    % colorbar
    % axis square

    % add title to movie figure
    % title(['Frame number = ',num2str(frame_num)]); 
    % title(['Frame number = ',num2str(frame_idx)]); 
    
    % frame_num_double = double(frame_num); % this conversion is somehow needed for title.
    frame_num_double = double(frame_idx); % this conversion is somehow needed for title.

    % Convert frame_num to string.
    frame_str = ['Frame Number: ' num2str(frame_num_double)];
        
    frame_pick = getframe(gcf);
    
    writeVideo(writerObj, frame_pick);
    
    % frame_idx = frame_idx + 1; % to iterate through all frames in gamma all.

end % END PLOTTING MOVIE

% close the video writer object
close(writerObj);


% choose [yes/no] to save the .mat file with measurement data.
if save_movie_as_mat == 1
    % save([movie_outstr, '.mat'], 'gam_real', 'init_trunc', 'max_trunc', 'Mk', 'hz', 'xx', 'numz',  'refframe', 'texpmat' );
    % save([movie_outstr, '.mat']);
    % save([movie_outstr '.mat'], [output_directory, dataset_name, '_GenDom_Skip2Trig', num2str(refimg), '_texp_R', '_',  num2str(max_trunc), '_M', num2str(M),'_zstep', num2str(hh), '.mat'], 'gamma_all', 'max_trunc', 'M', 'hh', 'xx', 'numz' )
    save([movie_outdir '/' movie_mat_fname, '.mat'], 'gamma_all', 'hh', 'xx');
end
