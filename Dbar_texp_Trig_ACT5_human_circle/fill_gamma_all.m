%===================================================================================================
% This script runs the D-bar algorithm, calling the necessary functions to compute the approximate 
% scattering transform texp and solve the Dbar equation
%
% This code is set up to reconstruct human data as difference images by selecting one reference
% frame from a multiframe dataset. 
%
% This is for…
% - ACT5 human data
% - circular domain
% - trig patterns
% - manual truncation (no gaussian truncation)
% - reference frame used: frame found in our "find-reference-frame" script.
%
% This code JUST meant to compute the conductivity distribution (gamma) for each frame in a 
% multiframe dataset, and add the real component of each to a matrix called gamma_all. 
%
% Plotting Reconstructions:
% - option to plot INDIVIDUAL frames (individual gammas) by setting
%   'display_images_to_screen' = 1
%
% Note: this uses the transformed DN map, so general domain functionality could
% be implemented by importing electrode positions from a boundary file, but this is not currently
% done in this code. 
%
% Note: this code is a little messy and could be cleaned up some, but it does seem to work fine
%
% Authors:          Melody Alsaker, Jennifer Mueller, Peter Muller
%                   Modifications made by Drew Fitzpatrick, Lydia Lonzarich, and Matthew Wong
% Date Modified:    December 2025
%
% Edits made in 2025...
% Editor:           Lydia Lonzarich
% Edits made:       separated 'generate gamma for each frame + option to
%                   plot individual frames' and 'generate movie reconstruction' logic into
%                   two separate scripts
% Date Modified:    December 2025
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
display_images_to_screen = 1;
save_images_as_jpg_files = 0;
save_gam_real_as_mat_file = 0;
% plot_movie = 0;
% saved = 1;


%===================================================================================================
%======================================== Specify External Data ====================================
%===================================================================================================
% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name of .mat file containing EIT data.
% datafname = 'modified_16x16_Sbj02_2D_16e_24_10_16_12_38_03_93750';
datafname = 'Sbj001_93kHz_perf_24_10_15_11_05_09_1';

% Directory for the program output to be saved to. If it doesn't exist, we'll create it later.
output_directory = [];
outdir = 'gammas/';
out_fname = 'jan14_1930_gamma_cond_distributions_09_1.mat';


% ==================================================================================================
% ================================ Specify Mesh Size Parameters  ===================================
% ==================================================================================================
Mk = 32;                % Size of k-grid is Mk x Mk
hz = 0.02;              % z-grid step size used to create the z grid (xx=-1:hz:1). Smaller value => finer mesh.


%===================================================================================================
%======================================== Specify Reconstruction Parameters ========================
%===================================================================================================
refframe = 37;    % Reference frame (e.g. at max expiration)
startframe = 200;    
endframe = 200;  

% determine the total # of frames to reconstruct (we must ignore the reference frame when it's in the range (startframe,endframe))
if refframe >= startframe & refframe <= endframe
    total_reconstruct_frames = endframe - startframe;
else
    total_reconstruct_frames = endframe - startframe + 1;
end

% initialize gamma_all for storing all frame reconstructions. 
xx = -1:hz:1;
N = numel(xx);
gamma_all = zeros(N, N, total_reconstruct_frames);

frame_idx = 1; % initialize the main for-loop indexing variable to 1.
  
% "grab" all frames in the dataset (MINUS the reference frame).
all_frames = startframe:endframe;
all_frames(all_frames == refframe) = []; % remove refframe index from the list of frames (that we'll iterate over).

gamma_best = 300;

cmap = 'jet';           % Select colormap for figures

init_trunc = 4.1;       % Initial trunc. radius. Used to trunc scattering transform. Choose something smallish
max_trunc = 4.6;          % Final max trunc. radius. Used to trunc scattering transform. Choose something bigger


%===================================================================================================
%======================== Load and Extract External Data & Physical Parameters =====================
%===================================================================================================
% Load measured data. We will pull various physical parameters from this.  
load([datadir, datafname])


%===================================================================================================
%======================== Generate Reconstructions With the Dbar Algorithm =========================
%===================================================================================================
% START OF MAIN FOR-LOOP (the dbar algorithm)
% Iterate over all frames in dataset (MINUS reference frames) and fill gamma_all with frame reconstructions. 
for frame = all_frames
    disp("frame: " + frame)


    % Voltages (use these to derive DN map)
    Vmulti = real(frame_voltage);   % Voltages for all frames in .mat file
    V = Vmulti(:,:,frame);          % Target frame voltage matrix. Selecting the measured voltage at 'frame' index
    Vref = Vmulti(:,:,refframe);    % Reference frame voltage matrix
    
    % Current pattern matrix (unnormalized and including all columns) (use these to derive DN map)
    J0 = cur_pattern; 
    
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
     J = J/norm(J(:,kk),2);
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
    
    
    num_frames = 1;      % Number of frames to reconstruct
    
    %================ Construct DN matrix for reference data ==================
    
    % Normalize the entries so that the voltages sum to zero in each col.
    Vref = Vref(1:numCP,:)'; 
    adj = sum(Vref)/L;
    Vref = Vref - adj(ones(L,1),:);
    
    
    refLambda = inv(Vref' * J);           % DN map for the reference frame, size numCP x numCP
    
    refLhat1 = (A1-1i*B1)*refLambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
    refLhat2 = (A1-1i*B1)*refLambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
    refLhat3 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
    refLhat4 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');
    
    refLhat = [refLhat1, refLhat2; refLhat3, refLhat4];
    
    %======Loop through all data sets, reconstruct conductivity for each ======
    
    % gamma is the conductivity we will reconstruct. Outside domain will be NaN
    gamma = ones(num_frames,N*N) * NaN;
    
    
    for jj = 1:num_frames
        gammatemp = zeros(1,numz);
        %================= Construct DN matrix for measured data ==============
        
        
        % Normalize the entries so that the voltages sum to zero in each col.
        V = V(1:numCP,:)'; 
        adj = sum(V)/L;
        V = V - adj(ones(L,1),:);
        Lambda = inv(V' * J);   % DN map for the measured frame
        
        Lhat1 = (A1-1i*B1)*Lambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
        Lhat2 = (A1-1i*B1)*Lambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
        Lhat3 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
        Lhat4 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');
        Lhat = [Lhat1, Lhat2; Lhat3, Lhat4];
        
        dLambda = Lhat - refLhat; % transformed DN map, size numCP x numCP 
        % dLambda = -(Lhat - refLhat); % transformed DN map, size numCP x numCP 
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
    
    gamma_real = real(gamma);
    
    gamma_all(:,:,frame_idx) = gamma_real;
    
    frame_idx = frame_idx + 1;
    
    % Switch to DICOM orientation
    % for jj = 1:num_frames
    %     gamma_real(jj,:,:) = fliplr(squeeze(gamma_real(jj,:,:)));
    % end

    if ~exist(outdir, 'dir')
        mkdir(outdir);
    end

    % % choose [yes/no] to save gamma conductivity reconstructions to a .mat file to make a movie with (in a separate script)
    % if save_gam_real_as_mat_file == 1
    %     save([outdir out_fname], 'gamma_all');
    % end 
    



    % ==============================================================
    % ==== Standardize Colorbar for Individual Frame Plots =========
    % ==============================================================
    % (option 1) og min/max method.
    frames_to_plot = 1:num_frames;
    datamin = min(min(min(gamma_real(frames_to_plot,:,:))));
    datamax = max(max(max(gamma_real(frames_to_plot,:,:))));
    datarange = datamax-datamin;
    
    % (option 2) robust percentile method.
    % all_vals = gamma_all(:);
    % datamin = prctile(all_vals, 2);
    % datamax = prctile(all_vals, 98);
    
    
    
    % ==================================================================================================
    % ============================== Set Up the Output Directory and Filename for Each Frame ===========
    % ==================================================================================================
    % if ~exist(outdir, 'dir')
    %        mkdir(outdir)        
    % end
    % 
    % outstr = [outdir, '/', datafname, '_R', num2str(init_trunc),'_',  num2str(max_trunc), '_Mk', num2str(Mk), '_recontime_', timeStampstr]; 
    
    


    % ==================================================================================================
    % =================== Plot and Save Individual Image Reconstruction ================================
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
            % imagesc(xx,xx,flipud(squeeze(gamma_all(:,:,1))),[datamin, datamax]);
            imagesc(xx,xx,flipud(squeeze(gamma_all(1,:,:))),[datamin, datamax]);


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
disp("All frames have now been reconstructed")



% ==================================================================================================
% =================== Save Gamma As A .mat File ====================================================
% ==================================================================================================
% choose [yes/no] to save gamma conductivity reconstructions to a .mat file to make a movie with (in a separate script)
if save_gam_real_as_mat_file == 1
    save([outdir out_fname], 'gamma_all');
end 

disp("gamma_all has been saved into a .mat file")