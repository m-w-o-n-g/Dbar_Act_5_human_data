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
% - Gaussian truncation
% - reference frame used: ONE frame found using our “find-best-refframe” script
%
% This code JUST meant to compute the conductivity distribution (gamma) for each frame in a multiframe dataset, and adds the real component
% of each to a matrix called gamma_all 
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
% original...
% Authors:          Melody Alsaker, Jennifer Mueller, Peter Muller
% Date Modified:    
%
% Edits made in 2025...
% Author: Sean Cowan
% Edits: - Added Gaussian truncation and removed intial truncation
%        - Made code ready for ACT5 Human Data
%
% Edits made in 2025...
% Authors: Drew Fitzpatrick, Lydia Lonzarich, and Matthew Wong
% Edits: Made code actually work with ACT5 human data.
%
%===================================================================================================

% clear all
% close all

timestart = tic;
total_runtime = 0;

%**************************************************************************
%===================== Set up Dataset-Specific Info =======================
% ATTN USER: Make necessary changes to all the parameters in this section.
% Only values within this section need to be modified between datasets.
%==========================================================================

% ==================================================================================================
% ================================= Choose What to Plot and Save ===================================
% ==================================================================================================
save_dbar_output_as_mat_file = 0;
display_images_to_screen = 0;
save_gam_real_as_mat_file = 0;


%===================================================================================================
%======================================== Load External Data ====================================
%===================================================================================================
% directory where data is stored.
data_dir = 'ACT5_humanData/';

% .mat file containing EIT data.
data_fname = 'perf_chunk_Sbj02_2D_16e_24_10_16_12_39_39_93750';

% File containing list of bdry pts and the directory where it is stored
bdry_file = 'siiri_boundary_lower_ring_bdryCoords.txt';
bdry_directory = [];



% ==================================================================================================
% ======================= Specify Mesh Size Parameters =============================================
% ==================================================================================================
% Perim_inches = 40;              % Perimeter of boundary (inches).
% Perim = Perim_inches * 0.0254;  % Perimeter of boundary (meters).
Perim = 0.76;                     % Perimeter of Siiri's boundary (cm).

M = 32;                         % Size of k-grid for Fourier domain is MxM. Enter a power of 2.
                                % M=16 is nice and fast. M=32 is more accurate. M=64 is super great but slowish

hh = 0.03;                      % Spatial z-grid step size. This changes the number of pixels in your reconstruction.  Smaller value => finer mesh.
                                % Choose 0.01 <= hh <= 0.065 for best results.

ee = 0.15;                       % Used to compute width of Gaussian window in FT. Smaller = more truncation

max_trunc = 4.6;                  % Specify a circular truncation region for the low-pass Fourier domain filter.
                                % This determines how much high-freq content we allow.
                                % Final max trunc radius. Choose something bigger.
                        
% Truncates colorbar for display purposes. Enter an integer from 0 to 10.
% If displaying a small number of images, choose something smaller.
percent_to_truncate_colorbar = 0;

% Select colormap for figures
cmap = 'jet';

% active_elecs = 1:32; 
% L = 32;
active_elecs = 1:16;
L = 16;



%===================================================================================================
%======================================== Specify Reconstruction Parameters ========================
%===================================================================================================
ref_frame = 10;
startframe = 1; 
endframe = 297;
highest_cond = 0; 
best_frame = 0;

% determine the total # of frames to reconstruct (we must ignore the reference frame when it's in the range (startframe,endframe))
if ref_frame >= startframe && ref_frame <= endframe
    total_reconstruct_frames = endframe - startframe;
else
    total_reconstruct_frames = endframe - startframe + 1;
end

% intialize a .mat file to store gammas (these are NOT reconstructed yet)
xx = -1:hh:1;
N = numel(xx);
gamma_all = zeros(N, N, total_reconstruct_frames);

frame_idx = 1; % initialize the main for-loop indexing variable to 1.
  
all_frames = startframe:endframe; % "grab" all frames in the dataset (MINUS the reference frame).
all_frames(all_frames == ref_frame) = []; % remove ref_frame index from the list of frames (that we'll iterate over).

gamma_best = 300;



%===================================================================================================
%============================ Generate Recons with Dbar Algorithm ==================================
%===================================================================================================

% begin main for-loop (the dbar algorithm begins)
% iterate over all frames in dataset (MINUS the reference frame)
for frame = all_frames 
    disp("Frame: " + frame)

    % ============================ set up voltage matrices ==============================================
    load([data_dir, data_fname]);
    
    % Load the voltage data for current frame -- i.e. the thing we want to reconstruct
    V_total = frame_voltage;            % grab voltages for all frames in .mat file
    V_total = V_total(active_elecs, active_elecs, :);
    V = V_total(:,:,frame);             % grab the current frame voltages.
    V(:,L) = [];                        % drop the 32nd column. This fixes the singular matrix problem.
    V = V.*1000;                        % scale voltages by 1000
    
    J = cur_pattern;
    J = J(active_elecs, active_elecs);
    J(:,L) = [];                        % drop the 32nd column.
    % total_num_frames = 1;
    
    % Load the voltage data for the reference frame -- i.e. the homogenous voltage data.
    Vref = V_total(:,:,ref_frame);      % grab the ref_frame voltages.
    Vref(:,L) = [];                     % remove the 32nd column.
    Vref = Vref.*1000;                  % scale voltages by 1000
    
    
    coords = load([bdry_directory, bdry_file], '-ascii'); % Load bdry points. Odd-indixed pts are electrode ctrs
    
    numCP = L - 1;                  % Number of linearly independent current patterns
    
    
    %======================== set up numerical parameters =======================
    % Grid and computational parameters
    s = max_trunc;         % truncation radius (to be used in k-grid)
    h = 2*s/(M-1);         % k-grid step size
    
    Mdiv2 = M/2;
    Mtimes2 = 2*M;
    
    
    %======================== set up electrode geometry parameters =======================
    % extract electrode geometry params from the loaded .mat file to create circular domain.
    eheight = 0.0254; ewidth = 0.0254;      % Electrode height, width (meters)
    eArea = ewidth * eheight;               % Area of electrode (meters^2)
    dtheta = 2*pi/L;                        % We assume equal electrode spacing.
    
    
    %================ set up Boundary Data and Arclength function ========================
    x_bdry = coords(:,2); y_bdry = coords(:,1);
    
    % Get polygon data: geom = [ area   X_cen  Y_cen  perimeter ]
    [geom,~,~] = polygeom(x_bdry,y_bdry);
    
    % Move origin to centroid of the boundary
    x_bdry = x_bdry - geom(2); y_bdry = y_bdry - geom(3);
    
    % Scale to match physical parameters
    ScaleFactor = Perim / geom(4);
    
    x_bdry = x_bdry * ScaleFactor; y_bdry = y_bdry * ScaleFactor;
    %figure
    %plot(x_bdry,y_bdry,'*')  % Plots the boundary
    %axis square
    %title('Scaled boundary shape')
    
    [bdry_theta,bdry_r] = cart2pol(x_bdry,y_bdry);
    bdry_rmax = max(bdry_r);
    
    % Rotate the boundary so that e1 is at the angular position 0 + dtheta.
    rot_angle = 0;  
    Rmat = [cos(rot_angle), -sin(rot_angle); sin(rot_angle), cos(rot_angle)];
    rot_coords = Rmat*[x_bdry'; y_bdry'];
    x_bdry = rot_coords(1,:)'; y_bdry = rot_coords(2,:)';
    %figure
    %plot(x_bdry,y_bdry,'o')  % Plots the boundary
    %axis square
    %title('Scaled and rotated boundary shape')
    
    % Angles corresponding to electrode centers
    etheta = dtheta:dtheta:2*pi;  % assume equal angular spacing
    
    Ldiv2 = floor(L/2);  % Note if L is odd, this means more sines than cosines
    
    [A1,A2,B1,B2,C1,C2,D1,D2,theta,r_th,a,M_pts,rmax]=Fourier_coefficients_gen(bdry_directory,bdry_file,Perim,L,Ldiv2);
    
    
    %======================Set up computational grids==========================
    %..........................................................................
    % Construct mesh of z-values representing physical domain. We can throw out
    % z-values outside the domain; these will not be needed in the computation.
    %..........................................................................
    xx = -1:hh:1;
    N = numel(xx);
    [z1,z2] = meshgrid(xx,xx);
    
    
    x_bdry = x_bdry / bdry_rmax;
    y_bdry = y_bdry / bdry_rmax;
    z = z1 + 1i*z2;
    z = reshape(z,N*N,1);                   % Set of z-vals is now a vector
    [IN, ON]=inpolygon(z1,z2,x_bdry,y_bdry); % Find indices of z-vals in domain
    zidx = find(IN + ON);                   % Indices of z-vals in domain
    z = z(zidx);                            % Get rid of z-vals outside domain
    numz = numel(zidx);                     % Number of in-domain z-vals
    conjz = conj(z);                        % Complex conj of in-domain z-vals
    
    %..........................................................................
    % Construct computational grid with M x M elements & complex variable k
    % We will need to eliminate the k-values outside the truncation radius and
    % within a small radius around k=0, but we will need to keep the original
    % k-grid for some comutations.
    %..........................................................................
    x = -s:h:s;
    [K1, K2] = meshgrid(x,x);
    k = K1 + 1i*K2;                         % The set of all k-vals (matrix)
    numk = M*M;                             % Total number of k-vals
    
    % kidx_max = find(abs(k)<max_trunc & abs(k)>0.1); % Indices of k-vals in trunc area
    kidx_max = find(abs(k)<max_trunc & abs(k)>1e-6); % Indices of k-vals in trunc area
    ktrunc_max = k(kidx_max);
    numktrunc_max = numel(ktrunc_max);
    conjktrunc_max = conj(ktrunc_max);
    conjk = conj(k);                        % conj of all k-vals (matrix)
    
    % The k-grid for the Green's function beta needs to be larger to accomodate
    % the convolution.
    xBig            = [-(s+((Mdiv2):-1:1)*h),x,s+(1:(Mdiv2))*h];
    [K1Big, K2Big]  = meshgrid(xBig,xBig);
    k_Big           = K1Big + 1i*K2Big;
    
    
    %======================Define Green's function beta========================
    beta = h*h ./(pi * k_Big); % Mult by h^2 for when we compute the convolution
    beta(M+1,M+1)=0;           % Set beta(0,0) = 0 to avoid singularity.
    
    % Take the fast fourier transform of the Green's function.
    % This is an averaging (Andreas's form)
    p1 = shiftr(1:(Mtimes2),0,M+1,1);
    p  = shiftr(1:(Mtimes2),0,M,1);
    fft_beta  = fftn(beta(p,p)+beta(p1,p1)+beta(p1,p)+beta(p,p1))/4;
    
    
    %======================= Construct Current Matrix J =======================
    CurrAmp = max(max(J));
    
    % Now normalize columns of J with respect to the L2 norm.
    J = J * sqrt(2/L)/CurrAmp;
    J(:,L/2) = J(:,L/2) * sqrt(1/2); % The L/2 col. gets different treatment
    
    
    %============= Precompute some values necessary for Dbar eqn ==============
    %..........................................................................
    % EXP encodes the factor exp(-i(kz+conj(kz))) / (4pi*conj(k)) used in Dbar
    % eqn. This will be multiplied by the scattering transform later to form
    % the pointwise multiplication operator TR.
    %..........................................................................
    EXP = zeros(M,M, numz);
    for ii = 1:numz
        EXP(:,:,ii) = exp(-1i*(k*z(ii) + conjk*conjz(ii)))./ conjk;
    end
    EXP = EXP / (4*pi);
    
    
    %..................Values necessary for linear solver......................
    % Construct rhs of eqn DBO*m = 1. We also use rhs for the init guess.
    rhs = [ones(numk,1); zeros(numk,1)];
    bnrm2 = M;  				 % Norm of rhs.
    tol = 1e-5;                         % Error tolerance
    maxit = 10;                         % Max number of iterations
    restrt = 5;                         % Max iterations before GMRES restart
    e1 = [1; zeros(2*numk-1,1)];        % First basis vector for R^n
    
    
    num_frames = 1;      % Number of frames to reconstruct
    
    
    %================ Construct DN matrix for reference data ==================
    % Normalize the entries so that the voltages sum to zero in each col.
    %  The voltages you get after synthesizing trig are different if you do not do this step first
    adj = sum(Vref)/L;
    Vref = Vref - adj(ones(L,1),:);
    
    
    refLambda = inv(Vref' * J);           % DN map, size L-1 x L-1
    
    refLhat1 = (A1-1i*B1)*refLambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
    refLhat2 = (A1-1i*B1)*refLambda(1:Ldiv2,Ldiv2+1:L-1)*(C2.'+1i*D2.');
    refLhat3 = (A2-1i*B2)*refLambda(Ldiv2+1:L-1,1:Ldiv2)*(C1.'+1i*D1.');
    refLhat4 = (A2-1i*B2)*refLambda(Ldiv2+1:L-1,Ldiv2+1:L-1)*(C2.'+1i*D2.');
    
    refLhat = [refLhat1, refLhat2; refLhat3, refLhat4];
    
    
    %======Loop through all data sets, reconstruct conductivity for each ======
    % gamma is the conductivity we will reconstruct. Outside domain will be NaN
    gamma = ones(num_frames,N*N) * NaN;
    
    
    for jj = 1:num_frames
        gammatemp = zeros(1,numz);
    
    
        %================= Construct DN matrix for measured data ==============
        % Normalize the entries so that the voltages sum to zero in each col.
        adj = sum(V)/L;
        V = V - adj(ones(L,1),:);
        Lambda = inv(V' * J);
        
        Lhat1 = (A1-1i*B1)*Lambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
        Lhat2 = (A1-1i*B1)*Lambda(1:Ldiv2,Ldiv2+1:L-1)*(C2.'+1i*D2.');
        Lhat3 = (A2-1i*B2)*Lambda(Ldiv2+1:L-1,1:Ldiv2)*(C1.'+1i*D1.');
        Lhat4 = (A2-1i*B2)*Lambda(Ldiv2+1:L-1,Ldiv2+1:L-1)*(C2.'+1i*D2.');
        Lhat = [Lhat1, Lhat2; Lhat3, Lhat4];
        
        dLambda = (Lhat - refLhat); % transformed DN map, size L-1 x L-1. use for datasets 1,3, . 
        % dLambda = -(Lhat - refLhat); % transformed DN map, size L-1 x L-1. use for datasets 2, 5.
    
        %==================Compute approx. scattering transform================
        texp = zeros(numk,1);
        
        L2=Ldiv2;
        ak_L2=((1i*ktrunc_max).^L2)/factorial(L2);
        akbar_L2=((1i*conjktrunc_max).^L2)/factorial(L2);
        sumjk=zeros(size((1i*conjktrunc_max)));
        sumk=zeros(size((1i*conjktrunc_max)));
        sumj=zeros(size((1i*conjktrunc_max)));
        
        % Compute sums from Jutta's paper and Ethan's work
        for j=1:L2-1
            akbar_j=((1i*conjktrunc_max).^j)/factorial(j);
            sumj=sumj+ak_L2.*akbar_j.*(dLambda(j,L2)+dLambda(L2+j,L2));
            for r=1:L2-1
                ak_k=((1i*ktrunc_max).^r)/factorial(r);
                sumk=sumk+akbar_L2.*ak_k.*(dLambda(L2,r)+dLambda(L2,L2+r));
                sumjk=sumjk+akbar_j.*ak_k.*(dLambda(j,r)+dLambda(L2+j,L2+r)+dLambda(j,L2+r)+dLambda(L2+j,r));
            end
        end
        texp(kidx_max)=sumjk+sqrt(2)*(sumk+sumj)+2*akbar_L2.*ak_L2.*dLambda(L2,L2);
        
        texp = texp * rmax* dtheta / (eArea); % JM ADDED THE FACTOR RMAX HERE - SEE LINE (5.18) OF ETHAN'S THESIS.  I think a FACTOR CurrAmp^2 does not need to go IN THE DENOMINATOR because our trig patterns do not include the current amplitude
        % for plotting only
        
        % Construct Gaussian window function adjusted to max_trunk
        a = -log(ee)/max_trunc^2;
        imaginary_k = imag(k);
        real_k = real(k);
        g_window = exp(-a*(real_k.^2 + imaginary_k.^2));
        
    
        tmat = reshape(texp,M,M);
    
        tmat_trunc = tmat.*g_window;
     
        % figure(2)
        % surf(real(tmat_trunc))
        % title('Real texp');
        % colorbar;
        % axis square;
        % 
        % figure(3)
        % surf(imag(tmat_trunc))
        % title('Imag texp');
        % colorbar;
        % axis square;
    
        %                 break
        
        % This is the pointwise multiplication operator used in the Dbar eqn.
        TR = repmat(tmat_trunc,[1,1,numz]) .* EXP;
        
    
        %========================== solve Dbar Equation =========================
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
            temp_fT = conj(reshape(f,M,M)) .* T;
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
                temp_fT = conj(reshape(f,M,M)) .* T;
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
                    temp_fT = conj(reshape(f,M,M)) .* T;
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
                        temp_fT = conj(reshape(f,M,M)) .* T;
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
                    temp_fT = conj(reshape(f,M,M)) .* T;
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
            
            sqrtgamma = m((numk+M)/2 +1) + 1i * m( (3*numk + M)/2 + 1);
            gammatemp(ii) = sqrtgamma * sqrtgamma;
        end
        gamma(jj,zidx) = gammatemp;
        
    end % All images have now been processed.
    disp("All images have now been processed")
    %======================================================================
    
    % Make each conductivity distribution into a matrix, one for each image.
    % gamma is 1 x 134 x 134 (num_frames=1, N=134) (1 is to specify that this gamma is just for 1 frame).
    gamma = reshape(gamma, num_frames, N, N); 
    
    % total_runtime = toc(timestart)
    
    % save only the real part of gamma. gam_real is same size as gamma.
    gam_real = real(gamma);  

    % squeeze(gam_real) ==> 134 x 134 matrix (the conductivity distribution matrix for a single frame).
    % this is where we add the conductivity distribution (gamma) for the current frame to a matrix, in position 'frame_idx', to a matrix with the rest of the gammas. 
    % gamma_all size is 134 x 134 x num_frames
    gamma_all(:, :, frame_idx) = squeeze(gam_real); % JUST REMOVED 12/19/2025
    % gamma_all(:, :, frame_idx) = gam_real; 

    frame_idx = frame_idx + 1;

    % logic for finding best refframe iteratively.
    gam_real(isnan(gam_real))=0;
    avg = mean(mean(gam_real));
    if avg > highest_cond
        highest_cond = avg;
        best_frame = frame;
    end


end % END MAIN FOR-LOOP ==> gamma_all has been completely filled with 'total_frames'-# of reconstructions.
disp("All frames have now been reconstructed")

disp(['Frame number ',num2str(best_frame),' to be chosen for reference frame.'])