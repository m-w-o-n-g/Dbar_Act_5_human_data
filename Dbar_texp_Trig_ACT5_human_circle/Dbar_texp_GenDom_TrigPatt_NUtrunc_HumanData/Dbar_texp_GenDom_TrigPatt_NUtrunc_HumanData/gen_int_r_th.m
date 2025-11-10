function rii=gen_int_r_th(theta,a,M,ii)
% This functions builds r(theta) and s(theta) based on the vector a found using Ethan
% Murphy's ideas for parametrizing a boundary. This function will be called
% by a quadrature rule for highly oscillatory integrands
%-----------Inputs-----------------
% a is a vector of coefficient for the Fourier expansion of r(theta)
% M is the number of terms taken in the expansion
% ii index number for the expansion of s(theta)r^ii(theta)cos(ii*theta)
%---------Output--------------------
% rii computes r^ii(theta) for the given boudary information
r_th=a(1)*ones(size(theta));
for m=1:M
    r_th=r_th+a(m+1)*cos(m*theta)+a(m+M+1)*sin(m*theta);
end
rii=r_th.^ii;