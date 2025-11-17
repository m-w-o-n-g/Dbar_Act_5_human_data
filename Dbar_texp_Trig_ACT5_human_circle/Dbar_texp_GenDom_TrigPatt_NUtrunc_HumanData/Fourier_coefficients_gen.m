function [A1,A2,B1,B2,C1,C2,D1,D2,theta,r_th,a,M_pts,rmax]=Fourier_coefficients_gen(bdry_directory, bdry_file, Perim,L,Ldiv2)
% This function computes the Fourier coefficients needed to transform the
% Dirichlet-to-Neumann map when the domain is closed and non-symmetric
%----------Inputs-------------------
% bdry_directory gives the local directory for the boundary info file
% bdry_file is the subject-specific .txt file containing the points on the
%   boundary
% Perim is the measured perimeter of the subject (in meters)

%---------Outputs-------------------
% A1,A2,B1,B2,C1,C2,D1,D2 are the resulting matrices of the Fourier
%   coefficients
% The following is from Melody Alsaker's code for retrieving boundary shape
% from a photo or CT scan

%global num_elec % number of electrodes
%L=num_elec;
theta=linspace(0,2*pi,300);
%if ~isempty(bdry_directory)
    % Load bdry pts. Odd-indexed pts are electrode ctrs
    coords = load([bdry_directory, bdry_file], '-ascii');
    [L_on_bdry,~] = size(coords); % Num. elecs. applied, according to bdry file. 
    L_on_bdry = L_on_bdry/2;
    dtheta=2*pi/L; % Uniform angle between each electrode (approx. true)
    %================ Set up Boundary Data and Arclength function==============
    x_bdry = coords(:,1); y_bdry = coords(:,2);
    % Get polygon data: geom = [ area   X_cen  Y_cen  perimeter ]
    [geom,~,~] = polygeom(x_bdry,y_bdry);
    % Move origin to centroid of the boundary
    x_bdry = x_bdry - geom(2); y_bdry = y_bdry - geom(3);
    % Scale to match physical parameters
    ScaleFactor = Perim / geom(4);
    x_bdry = x_bdry * ScaleFactor; y_bdry = y_bdry * ScaleFactor;
    [bdry_theta,bdry_r] = cart2pol(x_bdry,y_bdry);
    bdry_rmax = max(bdry_r);
    npts = 200; % Number of points to use
    coords_npts = interparc(npts, x_bdry, y_bdry); % Interpolate extra points along boundary
    % These now go from 0 to 2*pi
    M_pts = 15; % Approx. r(theta) as a series of cos,sin(m*theta) for m=1:M_pts
    a = Param_Bdry_Func_new(coords_npts(:,1),coords_npts(:,2),M_pts);
    % function handles for s(theta)*r^ii(theta) and r^ii(theta)
    srii='gen_int_sr_th';
    rii='gen_int_r_th';
    % Compute r(theta) to to get max(r) later
    r_th=a(1)*ones(size(theta));
    for m=1:M_pts
        r_th=r_th+a(m+1)*cos(m*theta)+a(m+M_pts+1)*sin(m*theta);
    end
    % Normalize boundary to be r=1 at largest
    rmax=max(r_th);
    r_th=r_th./rmax; 
    a=a./rmax;
%end
% The following comes from my work following Ethan Murphy and my work on
% non-circular domains
% Compute Fourier coefficients
L2=Ldiv2;  % Need to be careful for odd number of electrodes
A1=zeros(L2);A2=A1;B1=A1;B2=A1;C1=A1;C2=C1;D1=C1;D2=C1;

% Number of subintervals for integration
N_th=1000;
for ii=1:L2
    for jj=1:L2
        funA1=@(th) feval(srii,th,a,M_pts,ii).*cos(ii*th).*cos(jj*th);
        funA2=@(th) feval(srii,th,a,M_pts,ii).*cos(ii*th).*sin(jj*th);
        funB1=@(th) feval(srii,th,a,M_pts,ii).*sin(ii*th).*cos(jj*th);
        funB2=@(th) feval(srii,th,a,M_pts,ii).*sin(ii*th).*sin(jj*th);
        funC1=@(th) feval(rii,th,a,M_pts,ii).*cos(ii*th).*cos(jj*th);
        funC2=@(th) feval(rii,th,a,M_pts,ii).*cos(ii*th).*sin(jj*th);
        funD1=@(th) feval(rii,th,a,M_pts,ii).*sin(ii*th).*cos(jj*th);
        funD2=@(th) feval(rii,th,a,M_pts,ii).*sin(ii*th).*sin(jj*th);
%---------------To test on ellipse data------------------------------------
%         r='r_ellipse';s='s_ellipse';
%         funA1=@(th) feval(s,th).*feval(r,th).^ii.*cos(ii*th).*cos(jj*th);
%         funA2=@(th) feval(s,th).*feval(r,th).^ii.*cos(ii*th).*sin(jj*th);
%         funB1=@(th) feval(s,th).*feval(r,th).^ii.*sin(ii*th).*cos(jj*th);
%         funB2=@(th) feval(s,th).*feval(r,th).^ii.*sin(ii*th).*sin(jj*th);
%         funC1=@(th) feval(r,th).^ii.*cos(ii*th).*cos(jj*th);
%         funC2=@(th) feval(r,th).^ii.*cos(ii*th).*sin(jj*th);
%         funD1=@(th) feval(r,th).^ii.*sin(ii*th).*cos(jj*th);
%         funD2=@(th) feval(r,th).^ii.*sin(ii*th).*sin(jj*th);
%--------------------------------------------------------------------------        
        A1(ii,jj)=(1/pi)*simp1D(funA1,-pi,pi,N_th);
        A2(ii,jj)=(1/pi)*simp1D(funA2,-pi,pi,N_th);
        B1(ii,jj)=(1/pi)*simp1D(funB1,-pi,pi,N_th);
        B2(ii,jj)=(1/pi)*simp1D(funB2,-pi,pi,N_th);
        C1(ii,jj)=(1/pi)*simp1D(funC1,-pi,pi,N_th);
        C2(ii,jj)=(1/pi)*simp1D(funC2,-pi,pi,N_th);
        D1(ii,jj)=(1/pi)*simp1D(funD1,-pi,pi,N_th);
        D2(ii,jj)=(1/pi)*simp1D(funD2,-pi,pi,N_th);
    end
end
% 2 matrices aren't as large as 1's because they correspond to the sine
% patterns which are one fewer than cosine patterns
if (mod(L,2) == 0)
A2=A2(1:L2-1,1:L2-1);
C2=C2(1:L2-1,1:L2-1);
B2=B2(1:L2-1,1:L2-1);
D2=D2(1:L2-1,1:L2-1);
end
% If L is odd, no need to drop a pattern
    
%----- Uncomment following for ellipse test else comment-------------------
% r_th=r_ellipse(theta);a=[];M_pts=[];