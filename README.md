# About This Repository
This repository has scripts to perform Dbar reconstruction on human data using the Act 5 machine under the assumption of a circular domain. 



# Repository Contents:
- ACT5_humanData folder:
    - contains all the human data collected with the Act5 machine.
- Dbar_human_recons_movies folder:
    - Contains all the reconstruction movies and corresponding .mat files (with dataset info) for each dataset.
- Subject_2_Dbar_texp_Trig_ACT5_human_circle_Lhat.m:
    - Reshapes data matrices from 32x32 --> 16x16 
    - Use when 16 electrodes were used for data collection (not 32), but the data was entered into a 32x32 matrix ==> NaN entries that cause issues in the reconstruction scripts if not resolved.
- Dbar_texp_Trig_ACT5_human_circle_Lhat.m:
    - Finds the best reference frame to use for each dataset.
- movie_Dbar_texp_Trip_ACT5_human_circle_Lhat.m:
    - Generates a reconstruction movie for a datset.
    - Utilizes the original method of truncation (using init_trunc, and max_trunc)
    - Includes functionality to plot individual frames of a datset for analysis or debugging purposes. 
- gauss_trunc_movie_Dbar_texp_Trip_ACT5_human_circle_Lhat.m:
    - Generates a reconstruction movie for the dataset.
    - Utilizes the Gaussian truncation method.
    - Includes functionality to plot individual frames of a datset for analysis or debugging purposes.



# Datasets Chosen from subject001 and subject002

### Sbj001_35kHz_vent_24_10_15_10_45_29_1
- Reference frame: 361
- Chosen truncation radisu: 4.0-4.4
- Frames to use for reconstruction movie: 55-220
- Frame with the best image of heart and lungs: 70

### Sbj001_93kHz_perf_24_10_15_11_05_09_1
- Reference frame: 399
- Chosen truncation radius: 3.6-4.2
- Frames to use for reconstruction movie: 230-300
- Frame with best image of heart and lungs: 250

### Sbj001_93kHz_vent_24_10_15_10_51_57_1
- Reference frame: 37
- Chosen truncation radius: 4.1-4.6
- Frames to use for reconstruction movie: 130-350 (2 breaths) OR 136-236 (1 breath)
- Frame with best image of heart and lungs: 96

### Sbj02_2D_16e_24_10_16_12_39_39_93750
- Reference frame: 513
- Chosen truncation radius: 3.8-4.2
- Frames to use for reconstruction movie: 1200-1500
- Frame with best image of heart and lungs: 1335
