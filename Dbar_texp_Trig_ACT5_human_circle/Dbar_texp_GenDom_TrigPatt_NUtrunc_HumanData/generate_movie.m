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
data_fname = 'gamma_cond_distributions_57_1';

load([data_dir, data_fname]); % load gamma_all

movie_outdir = 'Dbar_human_recons_movies';                 % directory for the .mat measurements file and .avi recon movie file.
movie_mat_fname = 'TESTING_Dbar_movie_sbj001_57_1_gauss_trunc';    % filename for output .mat measurements file.
movie_avi_fname = [movie_mat_fname, '.avi'];               % filename for output .avi reconstruction movie.

% create the movie outdir if it doesn't exist.
if ~exist(movie_outdir, 'dir')
        mkdir(movie_outdir)
end


% ==================================================================================================
% =================================== Create a Movie with the Image Reconstructions ================
% ==================================================================================================
% all_frames = 1:size(gamma_all, 3);

% Truncates colorbar for display purposes. Enter an integer from 0 to 10.
% If displaying a small number of images, choose something smaller.
percent_to_truncate_colorbar = 0;

% Select colormap for figures
cmap = 'jet';

hh = 0.015; 
xx = -1:hh:1;

% create video writer object in the output directory.
writerObj = VideoWriter([movie_outdir, '/', movie_avi_fname]);

% set the frame rate to one frame per second
set(writerObj,'FrameRate',5);

% open the writer object.
open(writerObj);

% Standardizing the colorbar for the image reconstruction.
max_gamma_all = max(max(max(gamma_all)));
min_gamma_all = min(min(min(gamma_all)));
range_gamma = max_gamma_all - min_gamma_all;
cmax_gamma = max_gamma_all - 0.2*range_gamma; % should this be .2 or 0?
cmin_gamma = min_gamma_all + 0.2*range_gamma;

% recompute datamin and datamax from gamma_all for consistent colorbar scaling. 
datamin = min(gamma_all(:));
datamax = max(gamma_all(:)); 
datarange = datamax-datamin;
colorbartrunc = percent_to_truncate_colorbar * .01;
datamin = datamin + colorbartrunc * datarange;
datamax = datamax - colorbartrunc * datarange;

% initialize variable to keep track of the current frame # in the for-loop.
% frame_idx = 1; 


%% Plot movie
% iterate over all frames (MINUS the reference frame)
% for frame_num = all_frames
for frame_idx = 1:size(gamma_all, 3)

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
    imagesc(xx, xx, rot90(squeeze(gamma_all(:,:, frame_idx))), [datamin, datamax]);
    set(gca, 'Ydir', 'normal');
    % axis([-1 1 -1 1 ]);
    % caxis([cmin_gamma,cmax_gamma]); % added
    colorbar; % added
    axis off
    axis square;

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
    % frame_num_double = double(frame_idx); % this conversion is somehow needed for title.

    % Convert frame_num to string.
    % frame_str = ['Frame Number: ' num2str(frame_num_double)];
        
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
    save([movie_outdir '/' movie_mat_fname, '.mat'], 'gamma_all', 'max_trunc', 'M', 'hh', 'xx', 'numz');
end
