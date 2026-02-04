%===================================================================================================
% This script runs the D-bar algorithm, calling the necessary functions to compute the approximate 
% scattering transform texp and solve the Dbar equation
%
% This code is set up to reconstruct and plot 1 human data frame as a difference
% image by selecting a frame at which a mid-systole occurs as the reference
% frame.
%
% Note: fill_gamma_all.m can also reconstruct and plot single frames, but
% the reference frame used in that script is the frame with the highest
% average conductivity. In contrast, this script chooses the mid-systole
% frame closest to the frame chosen to reconstruct (e.g., say the
% mid-systoles in a dataset occur at frames 5, 11, 29, ... and I want to reconstruct frame 7. Then, 
% I would choose frame 5 as the reference frame).
%
% This is for…
% - ACT5 human data
% - circular domain
% - trig patterns
% - Gaussian truncation
% - reference frame used: frame at which a mid-systole occurs in ECG signal.
%
% Note: this code is a little messy and could be cleaned up some, but it does seem to work fine
%
% Authors:              Lydia Lonzarich
% Date Modified:        January 13, 2025
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
plot_movie = 1;
saved = 0;

%===================================================================================================
%======================================== Specify External Data ====================================
%===================================================================================================
% Directory where data is stored:
datadir = 'ACT5_humanData/';

% File name of .mat file containing EIT data.
% datafname = 'modified_16x16_Sbj02_2D_16e_24_10_16_12_38_03_93750';
datafname = 'perf_chunk_2_Sbj02_2D_16e_24_10_16_12_38_03_93750';

% File containing list of bdry pts and the directory where it is stored
bdry_file = 'EllipseBdry_1_952p6mm.txt';
bdry_directory = [];


% ==================================================================================================
% ================================ Specify Mesh Size Parameters  ===================================
% ==================================================================================================
perim_inches = 40;              % Perimeter of boundary (inches).
perim = perim_inches * 0.0254;  % Perimeter of boundary (meters).

Mk = 32;                        % Size of k-grid is Mk x Mk
hz = 0.02;                      % z-grid step size used to create the z grid (xx=-1:hz:1). Smaller value => finer mesh.


%===================================================================================================
%======================================== Specify Reconstruction Parameters ========================
%===================================================================================================
% must include ALL frames in startframe-endframe range to ensure mid-systoles are found correctly.   
startframe = 1;  
endframe = 261;

ee = 0.3;

max_trunc = 4.6; 

cmap = 'jet';

% L = 16;                   % FOR SBJ002 DATASETS...because current is only applied to electrodes 1-16, accordig to Jennifer.
active_elecs = 1:16;      % FOR SBJ002 DATASETS...because current is only applied to electrodes 1-16, according to Jennifer.

gamma_best = 300;

global_frame_idx = 1;     % to count the reconstructed frames?

%===================================================================================================
%======================== Load and Extract External Data & Physical Parameters =====================
%===================================================================================================
% Load measured data. We will pull various physical parameters from this.  
load([datadir, datafname])

% trim voltage measurement .mat files.
Vmulti = real(frame_voltage);                     % 32 x 32 x 2843 (raw voltage matrix)
% Vmulti = frame_voltage;
Vmulti = Vmulti(active_elecs, active_elecs, :);   % 16 x 16 x 2843 (remove 0 voltage row)

figure;
subplot(1,2,1); imagesc(Vmulti(:,:,1)); title('Raw Data'); 
xlabel('Columns'); ylabel('Rows');
colorbar



Vmulti_perf = Vmulti(:,:,startframe:endframe);    % 16 x 16 x num_frames

% extract mid-systoles.
mid_systoles = find_mid_systoles(Vmulti_perf);
disp(mid_systoles)

% "grab" all selected frames in the dataset.
all_frames = startframe:endframe;

% remove the mid-systole frames from the list of frames to reconstruct.
frames_to_reconstruct = setdiff(all_frames, mid_systoles);
disp("Number of frames available to reconstruct: " + num2str(length(frames_to_reconstruct)))


%===================================================================================================
%===================== Generate SINGLE Reconstructions With the Dbar Algorithm =====================
%===================================================================================================
refframe = mid_systoles(9);          % choose reference frame (one of the mid-systoles).
disp("Mid-systole chosen as reference frame: " + refframe)
Vref = Vmulti_perf(:,:,refframe);    % grab the reference frame's voltages.
Vref = Vref.*1000;

target_frame = refframe + 3;         % choose the frame to reconstruct relative to the chosen reference frame.
disp("Frame to reconstruct " + target_frame)

num_frames = 1;

% V = Vmulti_perf(:,:,frame);         % measured voltages for the target frame. 

% Current pattern matrix (unnormalized and including all columns) (use these to derive DN map)
J = cur_pattern;
J = J(active_elecs, active_elecs); % keep only 16 electrodes and the 16 current patterns.

L = length(J);    % Number of electrodes
numCP = L-1;       % Number of linearly independent current patterns

coords = load([bdry_directory, bdry_file], '-ascii'); % Load bdry points. Odd-indixed pts are electrode ctrs



%========================Set up numerical parameters=======================

% Grid and computational parameters
s = max_trunc;          % trunc radius (to be used in k-grid)
h = 2*s/(Mk-1);         % k-grid step size

Mdiv2 = Mk/2;
Mtimes2 = 2*Mk;


%======================== Set up electrode geometry parameters =======================

% extract electrode geometry params from the loaded .mat file to create circular domain.
% eheight = elec_height; ewidth = elec_width;     % Electrode height, width (in meters)
% perim = circumference * 0.0254;                 % Domain perimeter in meters (the variable circumference is loaded in inches)
% dtheta = 2*pi/L;                                % We assume equal electrode spacing
% etheta = dtheta:dtheta:2*pi;                    % Angular positions of electrode centers
% eArea = eheight * ewidth;                       % Simple area of electrode (meters^2)

% x_bdry = cos(etheta)'; y_bdry = sin(etheta)'; 

eheight = 0.0254; ewidth = 0.0254;
eArea = ewidth * eheight;
dtheta = 2*pi/L;


%================ Set up Boundary Data and Arclength function==============

% store position of electrodes as Lx2 matrix.
x_bdry = coords(:,1); y_bdry = coords(:,2);

% Get polygon data: geom = [ area   X_cen  Y_cen  perimeter ]
[geom,~,~] = polygeom(x_bdry,y_bdry);

% Move origin to centroid of the boundary
x_bdry = x_bdry - geom(2); y_bdry = y_bdry - geom(3);

% Scale to match physical parameters
domScaleFactor = perim / geom(4);

x_bdry = x_bdry * domScaleFactor; y_bdry = y_bdry * domScaleFactor;
%figure
%plot(x_bdry,y_bdry,'*')  % Plots the boundary
%axis square
%title('Scaled boundary shape')

[~,bdry_r] = cart2pol(x_bdry,y_bdry);
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

% angles corresponding to electrode centers
% etheta = dtheta:dtheta:2*pi;

Ldiv2 = floor(L/2);  % half hte number of electrodes. Note if L is odd, this means more sines than cosines

% [A1,A2,B1,B2,C1,C2,D1,D2,theta,r_th,a,M_pts,rmax]=Fourier_coefficients_gen_MODIFIED(coords,perim,L,Ldiv2);
[A1,A2,B1,B2,C1,C2,D1,D2,theta,r_th,a,M_pts,rmax]=Fourier_coefficients_gen(bdry_directory,bdry_file,perim,L,Ldiv2);



%======================Set up computational grids==========================
%..........................................................................
% Construct mesh of z-values representing physical domain. We can throw out
% z-values outside the domain; these will not be needed in the computation.
%..........................................................................
xx = -1:hz:1;
N = numel(xx);
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

% kidx_init = find(abs(k)<init_trunc & abs(k)>0.1);
% kidx_max = find(abs(k)<max_trunc & abs(k)>0.1);     % Indices of k-vals in trunc area
kidx_max = find(abs(k)<max_trunc & abs(k)>1e-6); % Indices of k-vals in trunc area
ktrunc_max = k(kidx_max);
numktrunc_max = numel(ktrunc_max);
conjktrunc_max = conj(ktrunc_max);
conjk = conj(k);                                    % conj of all k-vals (matrix)

% The k-grid for the Green's function beta needs to be larger to accomodate the convolution.
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
CurrAmp = max(max(J));

J = J(:,1:numCP); 
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
% num_frames = length(target_frame);


%================ Construct DN matrix for reference data ==================

Vref = Vref(1:numCP, :); % 16x16 ==> 15x16 (currentpatterns x electrodes). keep only the 15 linearly ind. current patterns.
Vref = Vref';            % 16x16 (electrodes x currentpatterns)
adj = sum(Vref)/L;       
Vref = Vref - adj;
Vref = Vref';            % 16x15 ==> 15x16 (currentpatterns x electrodes)

refLambda = inv(Vref * J);       % DN map for the reference frame, size numCP x numCP

refLhat1 = (A1-1i*B1)*refLambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
refLhat2 = (A1-1i*B1)*refLambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
refLhat3 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
refLhat4 = (A2-1i*B2)*refLambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');

refLhat = [refLhat1, refLhat2; refLhat3, refLhat4];


% ====== Reconstruct conductivity distribution for each frame (gamma) ======
gamma = ones(num_frames,N*N) * NaN;

for jj = 1:num_frames
    actual_frame = target_frame(jj);     % grab frame to reconstruct.
    V = Vmulti_perf(:, :, actual_frame);
    V = V.*1000;
    V = V(1:numCP, :);    % 16x16 ==> 15x16 (currentpatterns x electrodes). Keep only the 15 linearly ind. current patterns
    V = V';               % 15x16 ==> 16x15 (electrodes x currentpatterns).
    adj = sum(V)/L;       
    V = V - adj;
    V = V';               % 16x15 ==> 15x16 (currentpatterns x electrodes).

    gammatemp = zeros(1,numz);

    Lambda = inv(V * J);   % DN map for the measured frame
    
    Lhat1 = (A1-1i*B1)*Lambda(1:Ldiv2,1:Ldiv2)*(C1.'+1i*D1.');
    Lhat2 = (A1-1i*B1)*Lambda(1:Ldiv2,Ldiv2+1:numCP)*(C2.'+1i*D2.');
    Lhat3 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,1:Ldiv2)*(C1.'+1i*D1.');
    Lhat4 = (A2-1i*B2)*Lambda(Ldiv2+1:numCP,Ldiv2+1:numCP)*(C2.'+1i*D2.');
    Lhat = [Lhat1, Lhat2; Lhat3, Lhat4];
    
    % dLambda = Lhat - refLhat;  % transformed DN map, size numCP x numCP 
    dLambda = -(Lhat - refLhat); % transformed DN map, size numCP x numCP 
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
    
end % target frame's conductivity distribution has now been reconstructed.
%======================================================================

% Make each conductivity distribution into a matrix, one for each image.
gamma = reshape(gamma,num_frames,N,N);

gamma_real = real(gamma);

gamma_single = squeeze(gamma_real(1,:,:));

% Switch to DICOM orientation
% for jj = 1:num_frames
%     gamma_real(jj,:,:) = fliplr(squeeze(gamma_real(jj,:,:)));
% end

frames_to_plot = 1:num_frames;
datamin = min(min(min(gamma_real(frames_to_plot,:,:))));
datamax = max(max(max(gamma_real(frames_to_plot,:,:))));
datarange = datamax-datamin;

figure;
colormap(cmap);
imagesc(xx, xx, flipud(gamma_single)); 
axis square;
set(gca, 'YDir', 'normal')
colorbar;

title(['Target frame = ', num2str(target_frame), ...
    ', reference (mid-systole) = ', num2str(refframe), ...
    ', init_trunc = ', num2str(init_trunc), ...
    ', max trunc = ', num2str(max_trunc)]);

