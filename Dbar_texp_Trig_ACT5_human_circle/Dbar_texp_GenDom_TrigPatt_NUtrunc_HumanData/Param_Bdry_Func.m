% x,y are coords of points on boundary
% M is number of points for Fourier expansion

function [r_theta, s_theta, theta, rmax] = Param_Bdry_Func(x,y,M)

[theta, r] = cart2pol(x,y);

% theta(theta<0) = theta(theta<0) + 2* pi; % Make theta values in (0,2pi).

% This makes theta values within (0,2pi) and in order!!
a = find(theta<0);
b = a(1);
theta(b:end) = theta(b:end) + 2* pi; 
[N,~] = size(theta);

% Compute the matrix Q and the derivative matrix Qprime
Q=ones((2*M+1),N);
Qprime=zeros((2*M+1),N);
[T,scalar] = meshgrid(theta,1:M);
angle = T.* scalar;
Q(2:M+1,:) = cos(angle);
Q(M+2:2*M+1,:) = sin(angle);
Qprime(2:M+1,:) = -scalar .* sin(angle);
Qprime(M+2:2*M+1,:) = scalar .* cos(angle);

% Compute the vector a of unknowns
a=(Q*Q')\ Q*r;                             

% Compute the parameterized boundary function r(theta)
r_theta = Q'*a;
rmax = max(r_theta);
r_theta = r_theta / rmax;  % Scaling to maximum radius of 1.00

% Compute the arclength function s(theta)
r_theta_prime = (Qprime'*a) / rmax;
s_theta = sqrt( r_theta.^2 + r_theta_prime.^2 );





