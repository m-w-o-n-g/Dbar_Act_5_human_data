%===================================================================================================
% This script runs the D-bar algorithm, calling the necessary functions to compute the approximate 
% scattering transform texp and solve the Dbar equation
%
% This code is set up to reconstruct human data as difference images by selecting frames at 
% which diastoles occur from a multiframe dataset.
%
% This is for…
% - ACT5 human data
% - circular domain
% - trig patterns
% - Gaussian truncation
% - reference frames used: frames at which diastoles occur in ECG signal.
%
% Plotting Reconstructions:
% - option to plot ONE frame (one gamma) 
% - option to plot movie (iterate over all gammas in gamma_all)
%
% Note: this code is a little messy and could be cleaned up some, but it does seem to work fine
%
% Authors:              Lydia Lonzarich
% Date Modified:        January 16, 2025
%
%===================================================================================================

% clear 
% close all
timeStampstr = strrep(strrep(datestr(now,0),' ', '-'), ':', '-');  % Create timestamp string
timestart = tic; 


% ==================================================================================================
% ================================= Choose What to Plot and Save ===================================
% ==================================================================================================
save_dbar_output_as_mat_file = 0;
display_images_to_screen = 0;
save_images_as_jpg_files = 0;
plot_movie = 0;
saved = 0;
save_movie = 1;

%===================================================================================================
%======================================== Specify External Data ====================================
%===================================================================================================
% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name of .mat file containing EIT data.
% datafname = 'modified_16x16_Sbj02_2D_16e_24_10_16_12_38_03_93750';
datafname = 'Sbj001_93kHz_perf_24_10_15_11_05_09_1';

% Directory for the program output to be saved to. If it doesn't exist, we'll create it later.
outdir = 'HumanRecons';


% ==================================================================================================
% ================================ Specify Mesh Size Parameters  ===================================
% ==================================================================================================
Mk = 8;                % Size of k-grid is Mk x Mk
hz = 0.06;              % z-grid step size used to create the z grid (xx=-1:hz:1). Smaller value => finer mesh.


%===================================================================================================
%======================================== Specify Reconstruction Parameters ========================
%===================================================================================================
startframe = 1;    
endframe = 201;

% determine the total # of frames to reconstruct (we must ignore the reference frame when it's in the range (startframe,endframe))
% if refframe >= startframe & refframe <= endframe
%     total_reconstruct_frames = endframe - startframe;
% else
%     total_reconstruct_frames = endframe - startframe + 1;
% end
% all_frames = endframe - startframe + 1;
% 
% % "grab" all selected frames in the dataset.
% all_frames = startframe:endframe;
% 
% % initialize gamma_all for storing all frame reconstructions. 
% xx = -1:hz:1;
% N = numel(xx);
% gamma_all = zeros(N, N, all_frames);
% 
% frame_idx = 1; % initialize the main for-loop indexing variable to 1.
% 
% % "grab" all selected frames in the dataset.
% % all_frames = startframe:endframe; % 261

gamma_best = 300;

cmap = 'jet';             % Select colormap for figures

init_trunc = 3.6;         % Initial trunc. radius. Used to trunc scattering transform. Choose something smallish
max_trunc = 4.2;          % Final max trunc. radius. Used to trunc scattering transform. Choose something bigger

global_frame_idx = 1;     % to count the reconstructed frames?

active_elecs = 1:32; % because current is only applied to electrodes 1-16, according to Jennifer.

%===================================================================================================
%======================== Load and Extract External Data & Physical Parameters =====================
%===================================================================================================
% Load measured data. We will pull various physical parameters from this.  
load([datadir, datafname])

% trim voltage measurement .mat files.
Vmulti = real(frame_voltage);                     % 32 x 32 x 2843 (raw voltage matrix)
Vmulti = Vmulti(active_elecs, active_elecs, :);   % 16 x 16 x 2843 (remove 0 voltage row)
Vmulti_perf = Vmulti(:,:,startframe:endframe);    % 16 x 16 x num_frames=261

% extract diastoles. These will be used as reference frames.
diastoles = find_diastole(Vmulti_perf);
disp(diastoles) % there are 26.

% "grab" all selected frames in the dataset.
all_frames = startframe:endframe; % 261

% remove the diastole frames from the list of frames to reconstruct.
frames_to_reconstruct = setdiff(all_frames, diastoles);
disp("number of frames to reconstruct: " + num2str(length(frames_to_reconstruct)))

% initialize gamma_all for storing all frame reconstructions. 
xx = -1:hz:1;
N = numel(xx);
total_frames = endframe - startframe + 1; 
gamma_all = NaN(N, N, total_frames); % used to be zeros

% define the number of reference frames to be used.
num_refframes = length(diastoles);


%===================================================================================================
%======================== Generate Reconstructions With the Dbar Algorithm =========================
%===================================================================================================
for cycle_idx = 1:num_refframes
    curr_mid_systole = diastoles(cycle_idx);
    disp("Current diastole: " + curr_mid_systole)
    disp("Reference frame #: " + cycle_idx)

    refframe = diastoles(cycle_idx);  % use the i-th diastole as the current iteration's reference frame.
    % disp("refframe idx: " + refframe)

    Vref = Vmulti_perf(:,:,refframe);    % grab the reference frame's voltages.
    
    % define the frames to be used against the current iteration's reference frame.
    curr_frames = frames_to_reconstruct( ...
        frames_to_reconstruct > diastoles(max(cycle_idx-1,1)) & ...
        frames_to_reconstruct < diastoles(min(cycle_idx+1, end)) );

    num_frames = length(curr_frames);
    disp("num frames: " + num_frames)

    % % reconstruct frames around the current diastole: 5 before diastole --> 5 after diastole
    % % disp("all frames: " + all_frames)
    % % disp("startframe: " + startframe)
    % % disp("curr mid systole: " + curr_mid_systole)
    % frames_to_reconstruct = find(all_frames >= startframe & all_frames <= curr_mid_systole + 5);
    % disp("number of frames to reconstruct: " + length(frames_to_reconstruct))
    % 
    % % increment the starting frame for next iteration
    % startframe = curr_mid_systole + 5;


    % START OF DBAR ALGORITHM
    % Iterate over all frames in dataset subset and fill gamma_all with frame reconstructions. 
    % for frame = frames_to_reconstruct
    for frame = curr_frames
    
        V = Vmulti_perf(:,:,frame); % measured voltages for the current frame. 
    
        % % Voltages (use these to derive DN map)
        % Vmulti = real(frame_voltage);   % Voltages for all frames in .mat file
        % V = Vmulti(:,:,frame);          % Target frame voltage matrix. Selecting the measured voltage at 'frame' index
        % Vref = Vmulti(:,:,refframe);    % Reference frame voltage matrix
        
        % Current pattern matrix (unnormalized and including all columns) (use these to derive DN map)
        J0 = cur_pattern;
        J0 = J0(active_elecs, active_elecs);
        
        L = length(J0);    % Number of electrodes
        numCP = L-1;       % Number of linearly independent current patterns
        
        % L = 32;
        % dtheta = 2*pi/L;
        % theta = (dtheta:dtheta:2*pi)';
        % currents = zeros(L,L-1);
        % for j = 1:(L/2)
        %    currents(:,j) = cos(j*theta);
        % end
        % for j = 1:(L/2)
        %    currents(:,L/2+j) = sin(j*theta);
        % end
        % J0 = currents * current_amp;
        
        % Format of J0: Each column corresponds to a current pattern, as follows
        % curr_amp * [ cos(theta), cos(2*theta), ... ,cos(16*theta), sin(theta), sin(2*theta), ..., sin(16*theta) ]; 
        
        % extract electrode geometry params from the loaded .mat file to create circular domain.
        eheight = elec_height; ewidth = elec_width;     % Electrode height, width (in meters)
        perim = circumference * 0.0254;                 % Domain perimeter in meters (the variable circumference is loaded in inches)
        dtheta = 2*pi/L;                                % We assume equal electrode spacing
        etheta = dtheta:dtheta:2*pi;                    % Angular positions of electrode centers
        eArea = eheight * ewidth;                       % Simple area of electrode (meters^2)
        
        x_bdry = cos(etheta)'; y_bdry = sin(etheta)'; 
        
        
        %========================Set up numerical parameters=======================
        
        % Grid and computational parameters
        s = max_trunc;          % trunc radius (to be used in k-grid)
        h = 2*s/(Mk-1);         % k-grid step size
        
        Mdiv2 = Mk/2;
        Mtimes2 = 2*Mk;
        
        %================ Set up Boundary Data and Arclength function==============
        
        % store position of electrodes as Lx2 matrix.
        coords = zeros(L,2);        
        coords(:,1) = x_bdry;  coords(:,2) = y_bdry;
        
        % Get polygon data: geom = [ area   X_cen  Y_cen  perimeter ]
        [geom,~,~] = polygeom(x_bdry,y_bdry);
        
        % Move origin to centroid of the boundary
        x_bdry = x_bdry - geom(2); y_bdry = y_bdry - geom(3);
        
        % Scale to match physical parameters
        domScaleFactor = perim / geom(4);
        
        %x_bdry = x_bdry * domScaleFactor; y_bdry = y_bdry * domScaleFactor;
        
        %figure
        %plot(x_bdry,y_bdry,'*')  % Plots the boundary
        %axis square
        %title('Scaled boundary shape')
        
        [~,bdry_r] = cart2pol(x_bdry,y_bdry);
        bdry_rmax = max(bdry_r);
        
        % Rotate the boundary so that e1 is at the angular position 0 + dtheta.
        % rot_angle = 0;  
        % Rmat = [cos(rot_angle), -sin(rot_angle); sin(rot_angle), cos(rot_angle)];
        % rot_coords = Rmat*[x_bdry'; y_bdry'];
        % x_bdry = rot_coords(1,:)'; y_bdry = rot_coords(2,:)';
        %figure
        %plot(x_bdry,y_bdry,'o')  % Plots the boundary
        %axis square
        %title('Scaled and rotated boundary shape')
        
        Ldiv2 = floor(L/2);  % half hte number of electrodes. Note if L is odd, this means more sines than cosines
        
        [A1,A2,B1,B2,C1,C2,D1,D2,theta,r_th,a,M_pts,rmax]=Fourier_coefficients_gen_MODIFIED(coords,perim,L,Ldiv2);
        
        
        %======================Set up computational grids==========================
        
        %..........................................................................
        % Construct mesh of z-values representing physical domain. We can throw out
        % z-values outside the domain; these will not be needed in the computation.
        %..........................................................................
        [z1,z2] = meshgrid(xx,xx);
        
        z = z1 + 1i*z2;
        z = reshape(z,N*N,1);                       % Set of z-vals is now a vector
        [IN, ON]=inpolygon(z1,z2,x_bdry,y_bdry);    % Find indices of z-vals in domain
        zidx = find(IN + ON);                       % Indices of z-vals in domain
        z = z(zidx);                                % Get rid of z-vals outside domain
        numz = numel(zidx);                         % Number of in-domain z-vals
        conjz = conj(z);                            % Complex conj of in-domain z-vals
        
        %..........................................................................
        % Construct computational grid with M x M elements & complex variable k
        % We will need to eliminate the k-values outside the truncation radius and
        % within a small radius around k=0, but we will need to keep the original
        % k-grid for some comutations.
        %..........................................................................
        x = -s:h:s;
        [K1, K2] = meshgrid(x,x);
        k = K1 + 1i*K2;                                     % The set of all k-vals (matrix)
        numk = Mk*Mk;                                       % Total number of k-vals
        
        kidx_init = find(abs(k)<init_trunc & abs(k)>0.1);
        kidx_max = find(abs(k)<max_trunc & abs(k)>0.1);     % Indices of k-vals in trunc area
        ktrunc_max = k(kidx_max);
        numktrunc_max = numel(ktrunc_max);
        conjktrunc_max = conj(ktrunc_max);
        conjk = conj(k);                                    % conj of all k-vals (matrix)
        
        % The k-grid for the Green's function beta needs to be larger to accomodate
        % the convolution.
        xBig            = [-(s+((Mdiv2):-1:1)*h),x,s+(1:(Mdiv2))*h];
        [K1Big, K2Big]  = meshgrid(xBig,xBig);
        k_Big           = K1Big + 1i*K2Big;
        
        
        %======================Define Green's function beta========================
        
        beta = h*h ./(pi * k_Big); % Mult by h^2 for when we compute the convolution
        beta(Mk+1,Mk+1)=0;           % Set beta(0,0) = 0 to avoid singularity.
        
        % Take the fast fourier transform of the Green's function.
        % This is an averaging (Andreas's form)
        p1 = shiftr(1:(Mtimes2),0,Mk+1,1);
        p  = shiftr(1:(Mtimes2),0,Mk,1);
        fft_beta  = fftn(beta(p,p)+beta(p1,p1)+beta(p1,p)+beta(p,p1))/4;
        
        %======================= Construct Current Matrix J =======================
        J = J0(:,1:numCP); 
        for kk = 1:numCP
         % J = J/norm(J(:,kk),2);
            J(:,kk) = J(:,kk) / norm(J(:,kk), 2);
        end
        
        %============= Precompute some values necessary for Dbar eqn ==============
        
        %..........................................................................
        % EXP encodes the factor exp(-i(kz+conj(kz))) / (4pi*conj(k)) used in Dbar
        % eqn. This will be multiplied by the scattering transform later to form
        % the pointwise multiplication operator TR.
        %..........................................................................
        EXP = zeros(Mk,Mk, numz);
        for ii = 1:numz
            EXP(:,:,ii) = exp(-1i*(k*z(ii) + conjk*conjz(ii)))./ conjk;
        end
        EXP = EXP / (4*pi);
        
        
        %..................Values necessary for linear solver......................
        % Construct rhs of eqn DBO*m = 1. We also use rhs for the init guess.
        rhs = [ones(numk,1); zeros(numk,1)];
        bnrm2 = Mk;  				 % Norm of rhs.
        tol = 1e-5;                         % Error tolerance
        maxit = 10;                         % Max number of iterations
        restrt = 5;                         % Max iterations before GMRES restart
        e1 = [1; zeros(2*numk-1,1)];        % First basis vector for R^n
        
        
        % num_frames = 1;      % Number of frames to reconstruct
        % num_frames = length(curr_frames);
        
        %================ Construct DN matrix for reference data ==================
        
        % % Normalize the entries so that the voltages sum to zero in each col.
        % % Vref = Vref(1:numCP, :)';      % 16 x 15 ==> 
        % Vref = Vref(:, 1:numCP);         % Vref was 16x15, but we want to keep only the 15 linearly independent current patterns while keeping all electrodes, so drop one col ==> 16 x 15.
        % adj = sum(Vref)/L;               % 1 x 15
        % % Vref = Vref - adj(ones(L,1),:);
        % Vref = Vref - adj;               % 16 x 15 x num_frames (implicit expansion).
        % Vref = Vref';                    % transpose Vref: 16 x 15 ==> 15 x 16
        


        % trying new stuff
        Vref = Vref(1:numCP, :); % 16x16 ==> 15x16
        Vref = Vref';            % 15x16 ==> 16x15
        adj = sum(Vref)/L;       
        Vref = Vref - adj;
        Vref = Vref';            % 16x15 ==> 15x16
        % end of trying new stuff



        refLambda = inv(Vref * J);       % DN map for the reference frame, size numCP x numCP
        
        refLhat1 = (A1-1i*B1)*refLambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
        refLhat2 = (A1-1i*B1)*refLambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
        refLhat3 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
        refLhat4 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');
        
        refLhat = [refLhat1, refLhat2; refLhat3, refLhat4];
        
        %======Loop through all data sets, reconstruct conductivity for each ======
        
        % gamma is the conductivity we will reconstruct. Outside domain will be NaN
        gamma = ones(num_frames,N*N) * NaN;
        
        
        for jj = 1:num_frames
            % disp("iteration: " + jj)
            % V = Vmulti_perf(:,:,jj); 

            actual_frame = curr_frames(jj);
            V = Vmulti_perf(:, :, actual_frame);
            
            gammatemp = zeros(1,numz);
            %================= Construct DN matrix for measured data ==============
            
            
            % % Normalize the entries so that the voltages sum to zero in each col.
            % % V = V(1:numCP,:)'; 
            % V = V(:, 1:numCP);          % 16 x 16 ==> 16 x 15
            % adj = sum(V)/L;             % 1 x 15
            % % V = V - adj(ones(L,1),:);
            % V = V - adj;                % 16 x 15
            % V = V';                     % 16 x 15 ==> 15 x 16



            % new - trying stuff
            V = V(1:numCP, :); % 16x16 ==> 15x16
            V = V';            % 15x16 ==> 16x15
            adj = sum(V)/L;       
            V = V - adj;
            V = V';            % 16x15 ==> 15x16
            % end of trying new stuff


            Lambda = inv(V * J);   % DN map for the measured frame
            
            Lhat1 = (A1-1i*B1)*Lambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
            Lhat2 = (A1-1i*B1)*Lambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
            Lhat3 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
            Lhat4 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');
            Lhat = [Lhat1, Lhat2; Lhat3, Lhat4];
            
            dLambda = Lhat - refLhat; % transformed DN map, size numCP x numCP 
            %dLambda = -(Lhat - refLhat); % transformed DN map, size numCP x numCP 
            %dLambda = Lambda - refLambda;
            
            %==================Compute approx. scattering transform================
            
            texp = zeros(numk,1);
            
            L2=Ldiv2;
            ak_L2=((1i*ktrunc_max).^L2)/factorial(L2);
            akbar_L2=((1i*conjktrunc_max).^L2)/factorial(L2);
            sumjk=zeros(size((1i*conjktrunc_max)));
            sumk=zeros(size((1i*conjktrunc_max)));
            sumj=zeros(size((1i*conjktrunc_max)));
            
            % Compute sums from Jutta's paper and Ethan's work
            % looping over indices to compute sums from dlambda
            for j=1:L2-1
                akbar_j=((1i*conjktrunc_max).^j)/factorial(j);
                sumj=sumj+ak_L2.*akbar_j.*(dLambda(j,L2)+dLambda(L2+j,L2));
                for k=1:L2-1
                    ak_k=((1i*ktrunc_max).^k)/factorial(k);
                    sumk=sumk+akbar_L2.*ak_k.*(dLambda(L2,k)+dLambda(L2,L2+k));
                    sumjk=sumjk+akbar_j.*ak_k.*(dLambda(j,k)+dLambda(L2+j,L2+k)+dLambda(j,L2+k)+dLambda(L2+j,k));
                end
            end
            texp(kidx_max)=sumjk+sqrt(2)*(sumk+sumj)+2*akbar_L2.*ak_L2.*dLambda(L2,L2);
            
            % Implement nonuniform truncation of scattering data
            max_real_texp = max(real(texp(kidx_init)));
            max_imag_texp = max(imag(texp(kidx_init)));
            min_real_texp = min(real(texp(kidx_init)));
            min_imag_texp = min(imag(texp(kidx_init)));
        
            imagtexp = imag(texp); 
            realtexp = real(texp); 
            realtexp( realtexp>max_real_texp | realtexp<min_real_texp ) = 0; 
            imagtexp( imagtexp>max_imag_texp | imagtexp<min_imag_texp ) = 0; 
        
            texp= realtexp + 1i * imagtexp; 
        
            scaling_factor = rmax * dtheta / (eArea * gamma_best);
            texp = texp * scaling_factor; 
            
            texpmat = reshape(texp,Mk,Mk);
            
            % figure
            % subplot(1,2,1)
            % imagesc(real(texpmat))
            % colormap jet
            % title('Real texp');
            % colorbar;
            % axis square;
            
            % subplot(1,2,2)
            % imagesc(imag(texpmat))
            % colormap jet
            % title('Imag texp');
            % colorbar;
            % axis square;
            
            
            % This is the pointwise multiplication operator used in the Dbar eqn.
            TR = repmat(reshape(texp,Mk,Mk),[1,1,numz]) .* EXP;
            
            %==========================Solve Dbar Equation=========================
            
            % Loop through all z values in domain
            for ii = 1:numz
                T = TR(:,:,ii);
                m = rhs;   % Computation result, init guess is rhs
                done = 0;
                
                % Vee is orthog. projector onto Krylov subspace K, cols are ONB for K
                Vee = zeros(2*numk,restrt+1);
                
                % Upper Hessenberg matrix, H = Vee'A Vee
                H = zeros(restrt+1,restrt);
                
                % Needed for Givens rotation
                cs = zeros(restrt,1);
                sn = zeros(restrt,1);
                
                %=================== Inline Code for GMRES =======================
                % This block replaces the function call
                % [m,~] = GMRES(@DBop, rhs, 10, 1e-5, 10, [], [], rhs);
                %==================================================================
                f = m;
                %--------------------- Inline code for DBO ------------------------
                % This block replaces fxn call r = rhs - DBO(m,M,numk,T,fft_beta);
                %------------------------------------------------------------------
                f = f(1:numk) + 1i * f(numk+1:2*numk);
                
                % Construct conj(matf) .* T with zero padding to accommodate convolution
                temp_fT = conj(reshape(f,Mk,Mk)) .* T;
                temp_fT_Big  = zeros(Mtimes2);
                temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)) = temp_fT;  % zero padding
                
                % Compute action of operator on f. Note, the h^2 is already included in
                % beta, but it could go here instead
                temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
                tmp = temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)); %Remove zero padding
                f = f - tmp(:);
                
                %Stack real and imaginary parts
                r = rhs - [real(f); imag(f)];
                %------------------------------------------------------------------
                
                normr = norm(r);
                error = normr / bnrm2;
                
                % Begin first iteration
                Vee(:,1) = r / normr;
                ess = normr*e1;
                
                for kk = 1:restrt   % construct ONB using Gram-Schmidt
                    
                    f = Vee(:,kk);
                    %--------------------- Inline code for DBO ------------------------
                    % This block replaces fxn call w = DBO(Vee(:,kk),M,numk,T,fft_beta);
                    %------------------------------------------------------------------
                    f = f(1:numk) + 1i * f(numk+1:2*numk);
                    
                    % Construct conj(matf) .* T with zero padding to accommodate convolution
                    temp_fT = conj(reshape(f,Mk,Mk)) .* T;
                    temp_fT_Big  = zeros(Mtimes2);
                    temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)) = temp_fT;  % zero padding
                    
                    % Compute action of operator on f. Note, the h^2 is already included in
                    % beta, but it could go here instead
                    temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
                    tmp = temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)); %Remove zero padding
                    f = f - tmp(:);
                    
                    %Stack real and imaginary parts
                    w = [real(f); imag(f)];
                    %------------------------------------------------------------------
                    
                    if kk == 1  % This will be the only kk value in many cases.
                        H(1,1) = w'*Vee(:,1);
                        w = w - H(1,1)*Vee(:,1);
                        H(2,1) = norm(w);
                        
                        % form 1st Givens rotation matrix.
                        % We assume H(2,1) <= H(1,1).
                        temp = H(2,1) / H(1,1);
                        cs(1) = 1.0 / sqrt( 1.0 + temp*temp );
                        sn(1) = temp * cs(1);
                        
                        temp   = cs(1)*ess(1);     	% approximate residual norm
                        ess(2) = -sn(1)*ess(1);
                        ess(1) = temp;
                        H(1,1) = cs(1)*H(1,1) + sn(1)*H(2,1);
                        
                        error  = abs(ess(2)) / bnrm2;
                        
                        if ( error <= tol )
                            % update approximation and exit loop
                            m = m + Vee(:,1)* (ess(1) / H(1,1));
                            done = 1;
                            break
                        end
                        
                        Vee(:,2) = w / H(2,1);
                        H(2,1) = 0.0;
                        
                    else % kk > 1
                        for ll = 1:kk
                            H(ll,kk)= w'*Vee(:,ll);
                            w = w - H(ll,kk)*Vee(:,ll);
                        end
                        
                        H(kk+1,kk) = norm(w);
                        Vee(:,kk+1) = w / H(kk+1,kk);
                        
                        for ll = 1:kk-1            % apply Givens rotation
                            temp     =  cs(ll)*H(ll,kk) + sn(ll)*H(ll+1,kk);
                            H(ll+1,kk) = -sn(ll)*H(ll,kk) + cs(ll)*H(ll+1,kk);
                            H(ll,kk)   = temp;
                        end
                        
                        % form kk-th Givens rotation matrix
                        AA = H(kk,kk); BB = H(kk+1,kk);
                        if( BB ~= 0.0 && abs(BB) <= abs(AA) )
                            temp = BB / AA;
                            cs(kk) = 1.0 / sqrt( 1.0 + temp*temp );
                            sn(kk) = temp * cs(kk);
                        elseif( BB == 0.0 )
                            cs(kk) = 1.0;
                            sn(kk) = 0.0;
                        else
                            temp = AA / BB;
                            sn(kk) = 1.0 / sqrt( 1.0 + temp*temp );
                            cs(kk) = temp * cs(kk);
                        end
                        
                        temp   = cs(kk)*ess(kk);     	% approximate residual norm
                        ess(kk+1) = -sn(kk)*ess(kk);
                        ess(kk)   = temp;
                        H(kk,kk) = cs(kk)*H(kk,kk) + sn(kk)*H(kk+1,kk);
                        H(kk+1,kk) = 0.0;
                        error  = abs(ess(kk+1)) / bnrm2;
                        if ( error <= tol )                        % update approximation
                            Y = H(1:kk,1:kk) \ ess(1:kk);           % and exit
                            m = m + Vee(:,1:kk)*Y;
                            break
                        end
                    end
                end
                
                if ( done == 1 )
                    % Do nothing, we're done after one iteration.
                    
                else
                    % Probably won't get this far... usually converges after 1
                    % iteration unless error tolerance is very small
                    for iter = 2:maxit                 % begin iteration
                        f = m;
                        %--------------------- Inline code for DBO ------------------------
                        % This block, along with the lines immediately above and below,
                        % replace function call r = rhs - DBO(m,M,numk,T,fft_beta);
                        %------------------------------------------------------------------
                        f = f(1:numk) + 1i * f(numk+1:2*numk);
                        
                        % Construct conj(matf) .* T with zero padding to accommodate convolution
                        temp_fT = conj(reshape(f,Mk,Mk)) .* T;
                        temp_fT_Big  = zeros(Mtimes2);
                        temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)) = temp_fT;  % zero padding
                        
                        % Compute action of operator on f. Note, the h^2 is already included in
                        % beta, but it could go here instead
                        temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
                        tmp = temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)); %Remove zero padding
                        f = f - tmp(:);
                        
                        %Stack real and imaginary parts
                        result = [real(f); imag(f)];
                        %------------------------------------------------------------------
                        r = rhs - result;
                        Vee(:,1) = r / norm( r );
                        ess = norm(r)*e1;
                        
                        for kk = 1:restrt  % construct ONB using Gram-Schmidt
                            
                            f = Vee(:,kk);
                            %--------------------- Inline code for DBO ------------------------
                            % This block, along with the lines immediately above and below,
                            % replace function call w = DBO(Vee(:,kk),M,numk,T,fft_beta);
                            %------------------------------------------------------------------
                            f = f(1:numk) + 1i * f(numk+1:2*numk);
                            
                            % Construct conj(matf) .* T with zero padding to accommodate convolution
                            temp_fT = conj(reshape(f,Mk,Mk)) .* T;
                            temp_fT_Big  = zeros(Mtimes2);
                            temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)) = temp_fT;  % zero padding
                            
                            % Compute action of operator on f. Note, the h^2 is already included in
                            % beta, but it could go here instead
                            temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
                            tmp = temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)); %Remove zero padding
                            f = f - tmp(:);
                            
                            %Stack real and imaginary parts
                            result = [real(f); imag(f)];
                            %------------------------------------------------------------------
                            w = result;
                            
                            for ll = 1:kk
                                H(ll,kk)= w'*Vee(:,ll);
                                w = w - H(ll,kk)*Vee(:,ll);
                            end
                            
                            H(kk+1,kk) = norm( w );
                            Vee(:,kk+1) = w / H(kk+1,kk);
                            
                            for ll = 1:kk-1            % apply Givens rotation
                                temp =  cs(ll)*H(ll,kk) + sn(ll)*H(ll+1,kk);
                                H(ll+1,kk) = -sn(ll)*H(ll,kk) + cs(ll)*H(ll+1,kk);
                                H(ll,kk) = temp;
                            end
                            
                            % form kk-th Givens rotation matrix
                            AA = H(kk,kk); BB = H(kk+1,kk);
                            if( BB ~= 0.0 && abs(BB) <= abs(AA) )
                                temp = BB / AA;
                                cs(kk) = 1.0 / sqrt( 1.0 + temp*temp );
                                sn(kk) = temp * cs(kk);
                            elseif( BB == 0.0 )
                                cs(kk) = 1.0;
                                sn(kk) = 0.0;
                            else
                                temp = AA / BB;
                                sn(kk) = 1.0 / sqrt( 1.0 + temp*temp );
                                cs(kk) = temp * cs(kk);
                            end
                            
                            temp   = cs(kk)*ess(kk);     	% approximate residual norm
                            ess(kk+1) = -sn(kk)*ess(kk);
                            ess(kk)   = temp;
                            H(kk,kk) = cs(kk)*H(kk,kk) + sn(kk)*H(kk+1,kk);
                            H(kk+1,kk) = 0.0;
                            error  = abs(ess(kk+1)) / bnrm2;
                            if ( error <= tol )                        % update approximation
                                Y = H(1:kk,1:kk) \ ess(1:kk);           % and exit
                                m = m + Vee(:,1:kk)*Y;
                                break
                            end
                        end
                        
                        if ( error <= tol )
                            break
                        end
                        
                        Y = H(1:restrt,1:restrt) \ ess(1:restrt);
                        m = m + Vee(:,1:restrt)*Y;              % update approximation
                        f = m;
                        %--------------------- Inline code for DBO ------------------------
                        % This block, along with the lines immediately above and below,
                        % replace function call r = rhs - DBO(m,M,numk,T,fft_beta);
                        %------------------------------------------------------------------
                        f = f(1:numk) + 1i * f(numk+1:2*numk);
                        
                        % Construct conj(matf) .* T with zero padding to accommodate convolution
                        temp_fT = conj(reshape(f,Mk,Mk)) .* T;
                        temp_fT_Big  = zeros(Mtimes2);
                        temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)) = temp_fT;  % zero padding
                        
                        % Compute action of operator on f. Note, the h^2 is already included in
                        % beta, but it could go here instead
                        temp_fT_Big = ifftn(fft_beta.*fftn(temp_fT_Big));
                        tmp = temp_fT_Big((Mdiv2+1):(3*Mdiv2), (Mdiv2+1):(3*Mdiv2)); %Remove zero padding
                        f = f - tmp(:);
                        
                        %Stack real and imaginary parts
                        result = [real(f); imag(f)];
                        %------------------------------------------------------------------
                        r = rhs - result;                        % compute residual
                        ess(kk+1) = norm(r);
                        error = ess(kk+1) / bnrm2;        % check convergence
                        if ( error <= tol )
                            break
                        end
                    end
                end        %=================== END Inline Code for GMRES ====================
                
                sqrtgamma = m((numk+Mk)/2 +1) + 1i * m( (3*numk + Mk)/2 + 1);
                gammatemp(ii) = sqrtgamma * sqrtgamma;
            end
            gamma(jj,zidx) = gammatemp;
            
        end % All images have now been processed.
        %======================================================================
        
        % Make each conductivity distribution into a matrix, one for each image.
        gamma = reshape(gamma,num_frames,N,N);
        % total_runtime = toc(timestart)
        
        gamma_real = real(gamma);
        
        % gamma_all(:,:,global_frame_idx) = squeeze(gamma_real(jj,:,:));
        % global_frame_idx = global_frame_idx + 1;

        for jj = 1:num_frames
            frame_for_gamma_all = curr_frames(jj);
            gamma_all(:, :, frame_for_gamma_all) = squeeze(gamma_real(jj, :, :));
            % gamma_all(:,:,global_frame_idx) = squeeze(gamma_real(jj,:,:));
            % global_frame_idx = global_frame_idx + 1;
        end
        % gamma_all(:, :, global_frame_idx) = gamma_real;
        % global_frame_idx = global_frame_idx + 1;
        
        % Switch to DICOM orientation
        % for jj = 1:num_frames
        %     gamma_real(jj,:,:) = fliplr(squeeze(gamma_real(jj,:,:)));
        % end
        
        frames_to_plot = 1:num_frames;
        datamin = min(min(min(gamma_real(frames_to_plot,:,:))));
        datamax = max(max(max(gamma_real(frames_to_plot,:,:))));
        datarange = datamax-datamin;
        
        
        % ==================================================================================================
        % ============================== Set Up the Output Directory and Filename for Each Frame ===========
        % ==================================================================================================
        if ~exist(outdir, 'dir')
               mkdir(outdir)        
        end
        
        outstr = [outdir, '/', datafname, '_R', num2str(init_trunc),'_',  num2str(max_trunc), '_Mk', num2str(Mk), '_recontime_', timeStampstr]; 
        
        
        % ==================================================================================================
        % ================================= Plot and Save Individual Image Reconstruction ==================
        % ==================================================================================================
        if(display_images_to_screen == 1 || save_images_as_jpg_files == 1 )
            for jj = frames_to_plot
                
                % choose [yes/no] to display individual image reconstruction images to screen.
                if( display_images_to_screen == 1 )
                    h = figure;                    % create a blank figure window.
                else
                    h = figure('visible', 'off');  % Suppress display to screen.
                end
                
                colormap(cmap);
                
                % generate the pretty image reconstruction.
                imagesc(xx,xx,flipud(squeeze(gamma_all(:,:,1))),[datamin, datamax]);
                
                set(gca, 'Ydir', 'normal');
                
                colorbar;
                axis([-1 1 -1 1 ]);
                axis square;
                
                title(['Frame number = ',num2str(frame), ', init trunc = ', num2str(init_trunc), ', max trunc = ', num2str(max_trunc), ', refframe = ', num2str(refframe)]); % add title to figure for reference frame number.
                
                % choose [yes/no] to save image individual reconstruction image as a .jpg file.
                if( save_images_as_jpg_files == 1)
                    print(h,'-djpeg', [outstr '.jpg']);
                end
            end
        end
        
        % choose [yes/no] to save a .mat file with the variables used for generating the frame reconstruction.
        if( save_dbar_output_as_mat_file == 1)
            save([outstr, '.mat'],'gamma_real', 'init_trunc', 'max_trunc', 'Mk', 'hz', 'xx', 'numz',  'refframe', 'texpmat' );
        end
        
        fclose('all');
    
    end % END MAIN FOR-LOOP ==> gamma_all has been completely filled with 'total_frames'-# of reconstructions.
end


%%
if save_movie == 1
% ==================================================================================================
% =================================== Create a Movie with the Image Reconstructions ================
% ==================================================================================================
    % set up movie output directories.
    movie_outdir = 'Dbar_human_recons_movies';                 % directory for the .avi recon movie file.
    movie_outFname = ['TESTING_short_Dbar_movie_refframe=diastole_frames1500_1700_colorflip', datafname];         % output filename for recon movie.
    movie_outstr = [movie_outdir, '/', movie_outFname];        % directory for the recon movie's corresponding .mat file. 
    
    % create the movie outdir if it doesn't exist.
    if ~exist(movie_outdir, 'dir')
            mkdir(movie_outdir)
    end
    
    % create video writer object in the output directory.
    writerObj = VideoWriter([movie_outdir, '/', movie_outFname]);
    
    % set the frame rate to one frame per second
    set(writerObj,'FrameRate',5);
    
    % open the writer object.
    open(writerObj);
    
    % Standardizing the colorbar for the image reconstruction.
    max_gamma_all = max(max(max(gamma_all)));
    min_gamma_all = min(min(min(gamma_all)));
    % max_gamma_all = max(gamma_all(:), [], 'omitnan');
    % min_gamma_all = min(gamma_all(:), [], 'omitnan');
    range_gamma = max_gamma_all - min_gamma_all;
    cmax_gamma = max_gamma_all - 0.2*range_gamma;
    cmin_gamma = min_gamma_all + 0.2*range_gamma;
    
    % if cmin_gamma >= cmax_gamma
    %     cmin_gamma = min_gamma_all;
    %     cmax_gamma = max_gamma_all;
    % end
    
    % initialize variable to keep track of the current frame # in the for-loop.
    % frame_idx = 1; 
    
    % num_movie_frames = size(gamma_all, 3);
    movie_frames = frames_to_reconstruct;
    
    % Plot movie
    % iterate over all frames (MINUS the reference frame)
    for ii = 1:length(movie_frames)
    % for frame_num = all_frames
        % choose [yes/no] to display movie to screen.
        if plot_movie == 1
            figure('visible','on');
        else
            figure('Visible','off');
        end
    
        frame_idx = movie_frames(ii);
    
        colormap(cmap)
    
        % generate the pretty image reconstruction.
        imagesc(flipud(gamma_all(:,:,frame_idx)))
        
        caxis([cmin_gamma,cmax_gamma])
        colorbar
        axis square
        set(gca, 'Ydir', 'normal')
    
        % caxis([min(gamma_all(:)), max(gamma_all(:))])
        % caxis([0.9951, 1.0117]);
    
        title(['Frame number = ',num2str(frame_idx), ...
               ', init trunc = ', num2str(init_trunc), ...
               ', max trunc = ', num2str(max_trunc)]); % add title to figure for reference frame number.
        
        % frame_num_double = double(frame_num); % this conversion is somehow needed for title.
        frame_num_double = double(frame_idx); % this conversion is somehow needed for title.
    
        % Convert frame_num to string.
        frame_str = ['Frame Number: ' num2str(frame_num_double)];
            
        frame_pick = getframe(gcf);
        
        writeVideo(writerObj, frame_pick);
        
        % frame_idx = frame_idx + 1; % to iterate through all frames in gamma all.
    
    end % END PLOTTING MOVIE


    % save movie to file
    % choose [yes/no] to save the movie to a .avi file.
    if saved==1
        save([movie_outstr, '.mat'], 'gamma_real', 'init_trunc', 'max_trunc', 'Mk', 'hz', 'xx', 'numz', 'texpmat' );
    end
    
    % close the video writer object
    close(writerObj);

end
