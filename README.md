# About This Repository
This repository has scripts to perform Dbar reconstruction on human data using the Act 5 machine under the assumption of a circular domain and manual truncation. 



# Repository Contents:
- ACT5_humanData folder:
    - contains all the human data collected with the ACT5 machine.
- Dbar_human_recons_movies folder:
    - Contains the reconstruction movie and corresponding .mat file (with dataset info) for each dataset.
- siiri_boundary_lower_ring_bdryCoords.txt
    - The human boundary .txt file to be used for generating reconstructions with a physiologically-accurate boundary shape.
- trim_dataset_to_16x16.m:
    - Reshapes sbj02 data matrices from 32x32 --> 16x16
    - Necessary step for sbj02 reconstructions, as data for sbj02 was collected on 16 electrodes only. 
- find_highest_avg_cond_frame.m:
    - Finds the **single** frame with the highest avg conductivity.
    - This frame will be used as the reference frame for ventilation datasets and the sbj001 perfusion dataset. 
- fill_gamma_all.m
    - Reconstructs each frames' gamma (aka, conductivity distribution) within a dataset and saves the real part of it to a .mat file (gamma_all). 
    - Reference frame used for the Dbar algorithm = the **single** best reference frame found using the 'Dbar_texp_Trig_ACT5_human_circle_Lhat.m' script (our 'find-best-reference-frame' script)
    - Includes functionality to plot individual frames of a datset for analysis or debugging purposes. 
- gammas folder
    - Contains the gamma_all .mat file for each dataset.
    - The gamma_all is used for movie plotting in the 'generate_movie.m' script.
- generate_movie.m
    - Generates a movie by plotting each gamma saved in a dataset's gamma_all .mat file.
- make_separate_sbj2_perf_chunk.m
    - Isolate the perfusion chunk frames from the total 2843 frames in the Sbj02_2D_16e_24_10_16_12_38_03_93750 dataset.
- find_start_systoles
    - Find all start-systole in the ECG signal 
    - Uses findpeaks() MATLAB function.
- start_systole_as_refframe.m
    - Reconstructs each frames' gamma (aka, conductivity distribution) within a dataset and PLOTS the movie.
    - Reference frames used for the Dbar algorithm = the **multiple** start-systole points within the ECG signal.
- plot_single_frame_dynamic_refframe.m
    - Reconstruct and plot **one** human data frame as a difference image by selecting a dynamic reference frame (e.g., using a frame at which a mid-systole occurs as the reference frame).
- plot_single_frame_sbj002_perf_chunk.m
    - Reconstruct and plot **one** human data frame as a difference image by selecting a single reference frame: the frame at which a start-systole occurs.
- plot_voltages_ACT5_Subj2.m
    - This script plots the power waveform to reveal ventilation vs perfusion chunks within the sjb002 dataset (because the ventilation and perfusion data was merged into one file)
- movie_dbar_ctboundary.m
    - Functionally the same as fill_gamma_all.m, but uses a CT realistic human boundary instead of a circle.
- movie_Dbar_texp_Trip_ACT5_human_circle_Lhat.m
    - Reconstructs each frames' gamma (aka, conductivity distribution) within a dataset and PLOTS the movie 
    - Reference frame used for the Dbar algorithm = the single best reference frame found using the 'Dbar_texp_Trig_ACT5_human_circle_Lhat.m' script (our "find-best-reference-frame' script)
    - Note: we do not use this script anymore since I made separate "fill gamma_all" and "plot movie" scripts.
    - Includes functionality to plot individual frames of a datset for analysis or debugging purposes. 




# Datasets Chosen from subject001 and subject002

### Sbj001_35kHz_vent_24_10_15_10_45_29_1
- Ventilation set
- Reconstruction method: no gaussian truncation
- Reference frame: 361
- Chosen truncation radius: 4.0-4.4
- Frames to use for reconstruction movie: 55-220
- Frame with the best image of heart and lungs: 70

### Sbj001_93kHz_perf_24_10_15_11_05_09_1
- Perfusion set
- Reconstruction method: no gaussian truncation
- Reference frame: 399
- Chosen truncation radius: 3.6-4.2
- Frames to use for reconstruction movie: 230-300
- Frame with best image of heart and lungs: 250

### Sbj001_93kHz_vent_24_10_15_10_51_57_1
- Ventilation set
- Reconstruction method: gaussian truncation
- Reference frame: 37
- Chosen truncation radius: 4.1-4.6
- Frames to use for reconstruction movie: 130-350 (2 breaths) OR 136-236 (1 breath)
- Frame with best image of heart and lungs: 96

### Sbj02_2D_16e_24_10_16_12_38_03_93750 
- Perfusion set (frames 2440-2700)
- Reconstruction method: no gaussian truncation
- Reference frames: start-systoles
- Chosen truncation radius: 3.6-3.8
- Frame with best image of heart and lungs: 92

### Sbj02_2D_16e_24_10_16_12_39_39_93750
- Ventilation Set
- Reconstruction method: no gaussian truncation
- Reference frame: 513
- Chosen truncation radius: 3.8-4.2
- Frames to use for reconstruction movie: 1200-1500
- Frame with best image of heart and lungs: 1335
